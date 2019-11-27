import matplotlib
matplotlib.use('TkAgg')
import specfile4 as sp
import netCDF4 as nc
import makedatalist as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime as dt
import astropy.stats

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, LSR, AltAz, ICRS
import astropy.units as u

import importlib
import scipy.optimize

from copy import copy
from math import pi
#from regex import re
import re

#from IPython.html.widgets import FloatProgress
#from IPython.display import display

importlib.reload(md)
importlib.reload(sp)

plt.rcParams['font.size'] = 12

def getEnvironment():
    try:
        env = get_ipython().__class__.__name__
        if env == 'ZMQInteractiveShell':
            return 'Jupyter'
        elif env == 'TerminalInteractiveShell':
            return 'IPython'
        else:
            return 'OtherShell'
    except NameError:
        return 'Interpreter'



isnoteBook = getEnvironment() == 'Jupyter'
if isnoteBook:
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm



globBaseDir = '.'

def setDir(basedir) :
    if basedir != None : globBaseDir = basedir

def gaussianFunc(x, pk, pos, sigm) :
    return pk * np.exp( - ( ( ( x - pos)/sigm)**2/2.))


def binnedArray(arr, nbin) :
    if nbin <=1 : return arr
    else        : return np.array( [np.average(arr[(i)*nbin:(i+1)*nbin],axis=0) for i in range(int(len(arr)/nbin)) ])

def binnedArray2(arr, nbin) :
    if nbin <=1 : return arr
    else        : return np.array( [np.average(arr[:,(i)*nbin:(i+1)*nbin],axis=1) for i in range(int(arr.shape[1]/nbin)) ]).T


def modelBeam(xymesh, posx, posy, widthx, widthy, peak, sky) : #modelPar) :
    gsbeam = lambda x, w :  np.exp( -(x/w)**2/2.)
    axis_x = xymesh[0][0]
    axis_y = xymesh[1][:,0]
    modelx = gsbeam( axis_x - posx, widthx)[:,np.newaxis]
    modely = gsbeam( axis_y - posy, widthy)[:,np.newaxis]
    return (np.dot(modelx, modely.T) * peak + sky).flatten()


def showFileList(date='', source='', mode='', freq='', fmlo='') :
    utc = dt.timezone(dt.timedelta(hours=+0), 'utc')
    for x in md.getFileList(globBaseDir, ffilter=date) :
        p = md.getNCprops(x)
        if p == None : continue
        if re.compile(source).search(p['srcname']) == None : continue
        if re.compile(mode)  .search(p['obsmode']) == None : continue
        if re.compile(freq)  .search(str(p['obsfreq'])) == None : continue
        if re.compile(fmlo)  .search(str(p['fmlofile'])) == None : continue

        print(p['obsnum'], dt.datetime.fromtimestamp(p['tm'][0], tz=utc), '{:3.1f}'.format(p['tm'][1]-p['tm'][0]), p['srcname'], p['obsfreq'], p['obsmode'], p['obsgoal'])

def returnFileList(date='', source='', mode='', freq='', obsgoal='') :
    utc = dt.timezone(dt.timedelta(hours=+0), 'utc')
    info_list = []

    for x in md.getFileList(globBaseDir, ffilter=date) :
        p = md.getNCprops(x)
        if p == None : continue
        if re.compile(source).search(p['srcname']) == None : continue
        if re.compile(mode)  .search(p['obsmode']) == None : continue
        if re.compile(freq)  .search(str(p['obsfreq'])) == None : continue
        if re.compile(obsgoal)  .search(str(p['obsgoal'])) == None : continue

        info = [[str(p['obsnum']).zfill(6), str(dt.datetime.fromtimestamp(p['tm'][0], tz=utc)), '{:3.1f}'.format(p['tm'][1]-p['tm'][0]), p['srcname'], str(p['obsfreq']), p['obsmode'], p['obsgoal']]]
        info_list = info_list + info

    return info_list

class LinePointing :

    def __init__(self, obsnum, calnum=-1) :
#        def.beInf  = {}
#        def.ncInf  = {}

#        def.bedat  = None
        #        def.obsnum = =1
        self.obsnum = obsnum
        self.calnum = calnum

        self.ncfname = md.getLMTtpmByObsNum(self.obsnum, globBaseDir)
        self.befname = md.getXFFTSbyObsNum(obsnum, globBaseDir)
        self.calbename= None if calnum == -1 else md.getXFFTSbyObsNum(calnum, globBaseDir)


        print(f'   nc file = {self.ncfname}')
        print(f'   BE file = {self.befname}')


    def loadlmttpm(self) :

        ncfile= nc.Dataset(self.ncfname)
        ncprop = md.getNCprops(self.ncfname)

        utc = dt.timezone(dt.timedelta(hours=+0), 'utc')

        nctime = ncfile.variables['Data.TelescopeBackend.TelTime'][:]
        bufpos = ncfile.variables['Data.TelescopeBackend.BufPos'][:]

        self.az_map = ncfile['Data.TelescopeBackend.TelAzMap'][:]
        self.el_map = ncfile['Data.TelescopeBackend.TelElMap'][:]

        self.fif    = [ncfile['Header.B4r.If1Freq'][:][0], ncfile['Header.B4r.If2Freq'][:][0]]
        self.fline  = ncfile['Header.B4r.LineFreq'][0]
        self.srcename = ncprop['srcname']

        self.dAZuser = ncfile['Header.PointModel.AzUserOff'][0]
        self.dELuser = ncfile['Header.PointModel.ElUserOff'][0]

        self.antonpos = bufpos == 0
        self.antoffpos= bufpos == 1

        self.t_ant = np.array([dt.datetime.fromtimestamp(t, tz=utc) for t in nctime], 'datetime64[us]')

        ncfile.close()


    def loadxffts(self,
                  nIF = 4,
                  sbdef  = [ 'lsb', 'usb', 'lsb', 'usb' ],
                  poldef = [ 'pol1', 'pol1', 'pol2', 'pol2' ],
                  usePol= ['pol1','pol1'],
                  useFreq=  129.36324 ,  #GHz
                  useVel =  0 ,        #km/s
                  useVelRange = 500,     #km/s
                  timediffTorr= 2) :

        ##############
        # read xffts

        fif    = self.fif
        fline  = self.fline
        print( f"   obsFreq = {fline} GHz")
        print( f"   IF freq = {fif} GHz")

        df_org = 2.5/32768        #GHz
        iffreq_ch_org = fif[1] + np.arange(0,-2.5, -df_org)
        fusb = fline + iffreq_ch_org
        flsb = fline - iffreq_ch_org

        print( "   freq range (usb) : {f} GHz".format(f=[fusb[-1], fusb[0]]))
        print( "              (lsb) : {f} GHz".format(f=[flsb[ 0], flsb[-1]]))

        self.LineSkipFlag = False

        if  fusb[-1] < useFreq < fusb[0]:
            useSB ='usb'
            fIF_org = fusb
        elif flsb[0] < useFreq < flsb[-1]:
            useSB = 'lsb'
            fIF_org = flsb
        else :
            self.LineSkipFlag = True
            #raise Exception('requested frequency {useFreq} is out of observation range')
#            return False

        if not self.LineSkipFlag:
            usespw = [ i for i in range(nIF) if (sbdef[i] == useSB ) and ( poldef[i] in usePol) ]
            print(f'   {useSB} {usePol} (spw={usespw}) will be used for pointing meas.')


            vrange = np.array([-useVelRange/2.+useVel, useVelRange/2.+useVel])
            frange = useFreq * ( 1 - vrange/3.e5)

            ch_st, ch_ed = sorted([np.abs( frange[0] - fIF_org).argmin(), np.abs( frange[-1] - fIF_org).argmin()])
            print("   BE channel range to use = {d1} : {d2}".format( d1=ch_st, d2=ch_ed))

            xfftsf = sp.SpecFile(self.befname, chrange=[ch_st,ch_ed], spw=usespw)
            fIF     = fIF_org[ch_st : ch_ed]
            xfftsf.analyScanPattern()

            if self.calbename != None :
                xffts_cal = sp.SpecFile(self.calbename, chrange=[ch_st,ch_ed], spw=usespw)
                xffts_cal.analyScanPattern()


            ################################
            # load header inf
            ################################

            ndat = xfftsf.getSize()
            betime     = np.array([ xfftsf.binHead(i)['date'][:-4] for i in range(ndat) ],'datetime64[us]')   #+  np.timedelta64(9, 'h')
            integtime  = np.array([ xfftsf.binHead(i)['integTime'] for i in range(ndat) ],'timedelta64[us]')
            bebuf  = np.array([ xfftsf.binHead(i)['bufPos'] for i in range(ndat) ])
            beonpos  = bebuf == 'ON'
            beoffpos = bebuf == 'REF'

            ntpsync = np.array( [ xfftsf.binHead(i)['date'][-4:-1] == b'GPS' for i in range(ndat) ])

            tm_PCsync = copy(betime)
            tm_PCsync[ ntpsync == False ] -=  np.timedelta64(9, 'h')

            blnktime = np.timedelta64(1000,'us')

            ntpsync_work        = copy(ntpsync)

            print('### BE time corretion')
            i_needscorrect, = np.where(ntpsync_work==False)
            print('  step1')
            for i in tqdm(i_needscorrect) :
                if i<1 : continue    #just in case
                tmincr = tm_PCsync[i-1] + integtime[i-1] + blnktime
                tmdiff = abs((tmincr-tm_PCsync[i]).astype(float)) * 1.e-6    #microsec->sec
    #            print('   {s}({d:.1f})'.format(s=tm_PCsync[i], d=tmdiff), end="")    #, tmdiff, ntpsync_work[i-1])
                if tmdiff < timediffTorr and ntpsync_work[i-1] == True :
                    tm_PCsync[i] = tmincr
                    ntpsync_work[i] = True
    #                print( '=> {s}'.format(s=tmincr))
    #            else :
    #                print( ' : correction failed')

            i_needscorrect, = np.where(ntpsync_work==False)
            #try:
            print('  step2')
            for i in tqdm(i_needscorrect[::-1]) :
                if i> ndat -1 : continue
                tmincr = tm_PCsync[i+1] - integtime[i] - blnktime
                tmdiff = abs((tmincr-tm_PCsync[i]).astype(float))*1.e-6
    #            print('   {s}({d:.1f})'.format(s=tm_PCsync[i], d=tmdiff), end="")    #, tmdiff, ntpsync_work[i-1])
    #            print('   {s}'.format(s=tm_PCsync[i]), end="")    #, tmdiff, ntpsync_work[i-1])
                if tmdiff< timediffTorr and ntpsync_work[i+1] == True :
                    tm_PCsync[i]= tmincr
                    ntpsync_work[i] = True
    #                print( '=> {s}'.format(s=tmincr))
    #            else :
    #                print( ' : correction failed')

    #                print( 'corrected : {s}'.format(s=tmincr))

    #        print('### stats ###')
            #except:
            #    print('step2 fail')


            ncorrupted = ntpsync.tolist().count(False)
            nrecovered = ((ntpsync == False) & (ntpsync_work == True)).tolist().count(True)
            print(f'   num of corrupted timestamps = {ncorrupted}/{ndat}')
            print(f'   recovered timestamps        = {nrecovered}/{ncorrupted}')
            print( '   flagged time                = {f}%'.format( f = np.sum( integtime[ntpsync_work==False])/np.sum(integtime)))



            self.nbedat = ndat
            self.betime = tm_PCsync
            self.integtime = integtime
            self.beonpos   = beonpos  & ntpsync_work
            self.beoffpos  = beoffpos & ntpsync_work
            self.corrupted = ntpsync_work==False
            self.xffts      = xfftsf
            self.fIF        = fIF
            self.f_use      = useFreq
            self.p_use      = usePol
            if self.calbename != None :
                self.xffts_cal  = xffts_cal

    def showAntInf(self):
        fig = plt.figure(figsize=(8,4))
        plt.subplot(1,2,1);  plt.scatter( self.azon_raw, self.elon_raw );     plt.title('dAz-dEl (org)')
        plt.subplot(1,2,2);  plt.scatter( self.azon_intrp, self.elon_intrp ); plt.title('dAz-dEl (intrp)')


    def getOnOffSpecs(self, binning=8, interp='linear', calib=True, yfactor=2., Tr=300.):
        tant_on = self.t_ant.astype(float)[self.antonpos]
        tbe_on  = self.betime.astype(float)[self.beonpos]
        azon_raw = self.az_map[self.antonpos] * 3600 * 180/pi     #rad -> sec
        elon_raw = self.el_map[self.antonpos] * 3600 * 180/pi

        azon_intrp  = interp1d(tant_on - tant_on[0], azon_raw, fill_value='extrapolate')(tbe_on - tant_on[0])
        elon_intrp  = interp1d(tant_on - tant_on[0], elon_raw, fill_value='extrapolate')(tbe_on - tant_on[0])





        #### create off-pos spectra #####


        # off pos

        print('### integrating off position spectra ... ')
        offscans =  np.unique(np.array(self.xffts.scanIds)[self.beoffpos] )
        offtime = [None] * len(offscans)
        offspecs= [None] * len(offscans)
        for ix in tqdm(range(len(offscans))) :
            scanid = offscans[ix]
            #time
            avgtime = np.average(self.betime[ np.array(self.xffts.scanIds==scanid) ].astype(float))
            offtime[ix] = avgtime

            #spec
            timesum = 0
            specsum = 0
            for idat in np.arange(self.nbedat)[self.xffts.scanIds==scanid] :
                theData = self.xffts.readData(idat)
                tint    = self.integtime[idat].astype(float)
                timesum += tint
                specsum += binnedArray2(theData,binning)  # * tint
#                theSpec = self.xffts.readSpec0(idat)
#                tint = theSpec['integTime']
#                timesum += tint
#                specsum += binnedArray2(theSpec['data'],binning)  # * tint
            offspecs[ix] = specsum/timesum


        # calib

        print('### R-Sky spectra ###')
        if calib==False or self.calbename == None :
            print(f'   using given Yfactor {yfactor} and Tr={Tr}K')
            sky_spec = np.average(offspecs,axis=0)
            r_spec   = sky_spec * 10**(yfactor/1.)
            rskySpec     = r_spec - sky_spec
        else :
            print(f'   using calib BE data {self.calbename}(self.calnum) and given Tr={Tr}K')
            ndat_cal = self.xffts_cal.getSize()
            calRpos = np.array([x == ('CAL', 'R') for x in self.xffts_cal.scanPattern])
            calSpos = np.array([x == ('CAL', 'S') for x in self.xffts_cal.scanPattern])
            calRspec = []
            calRint  = 0
            calSspec = []
            calSint  =0
            print('   collecting R-pos data ...')
            for ical in tqdm(np.arange(ndat_cal)[calRpos]) :
                calTint = self.xffts_cal.binHead(ical)['integTime']
                calSpec = binnedArray2(self.xffts_cal.readData(ical), binning)
                calRspec.append(calSpec)
                calRint += calTint
            print('   collecting S-pos data ...')
            for ical in tqdm(np.arange(ndat_cal)[calSpos]) :
                calTint = self.xffts_cal.binHead(ical)['integTime']
                calSpec = binnedArray2(self.xffts_cal.readData(ical), binning)
                calSspec.append(calSpec)
                calSint += calTint

            rspec = np.sum( calRspec, axis=0)/calRint
            sspec = np.sum( calSspec, axis=0)/calSint
            print('   measured Y-factor is {yf:1.3}'.format(yf = np.log10(np.median(rspec/sspec))*10))

            rskySpec = rspec - sspec

        # on pos

        onspecIx = np.arange(len(self.beonpos))[self.beonpos]
        ontime = self.betime.astype(float)[self.beonpos]

        print('### collecting on-position spectra ... ')
        onspecs = [None] * len(onspecIx)
        for i in tqdm(range(len(onspecIx))) :
#            dat = self.xffts.readSpec0(onspecIx[i])
#            onspecs[i]= binnedArray2(dat['data'],binning)/dat['integTime']
            dat = self.xffts.readData(onspecIx[i])
            onspecs[i]= binnedArray2(dat,binning)/self.integtime[onspecIx[i]].astype(float)

        print('### interpolating off-position spectra ... ')
        if interp=='linear' :
            offspec_interp = interp1d( offtime - offtime[0], offspecs, axis=0) (ontime - offtime[0])
        elif interp=='nearest' :
            offspec_interp = np.array([  offspecs[np.abs( offtime - ontime[i]).argmin()] for i in range(len(ontime))])
        else :
            raise Exception(f'unsupported interpolation {interp}')

        onoffspecs_raw = (onspecs - offspec_interp)/rskySpec * Tr
        onoffspecs = np.array([ [ispec[j] - astropy.stats.biweight_location(ispec[j]) for j in range(self.xffts.narray)] for ispec in onoffspecs_raw ])



        self.onoffspecs = onoffspecs
        self.freq       = binnedArray(self.fIF, binning)
        self.vel        = ( 1 - self.freq/self.f_use) * 2.99792458e5
        self.azon_intrp = azon_intrp
        self.elon_intrp = elon_intrp
        self.azon_raw   = azon_raw
        self.elon_raw   = elon_raw






    def showSpectrum(self, nTimeAvg=16, vrange =[ None, None], sig_thres=3.) :

        ##

        if vrange ==[ None, None] :
            fitch = [ 0, len(self.freq)-1 ]
        else :
#            frange = self.f_use * ( 1 - np.array(vrange)/2.99792458e5 )
            fitch = sorted([ np.abs( self.vel - vrange[0]).argmin(), np.abs( self.vel - vrange[1]).argmin() ])
#            print(frange, fitch)
        ## mk time average

        specTimeAvg = binnedArray(self.onoffspecs, nTimeAvg) #np.array( [ np.average( self.onoffspecs[:,i*nTimeAvg:(i+1)*nTimeAvg,:], axis=1) for i in range(int(self.onoffspecs.shape[1]/nTimeAvg))])
        specMax = np.max(specTimeAvg, axis=0)

        ##
        #fig = plt.figure(figsize=(8,4))
        plt.subplot(4,1,1)
        for i in range(specMax.shape[0]):
            plt.step( self.vel, specMax[i], where='mid')
        plt.xlabel('V_topo (km/s)')
        plt.legend(self.p_use)

        ##
        # fitting

        specAvg = np.average( specMax[:, fitch[0]:fitch[1]+1], axis=0)
        bl = astropy.stats.biweight_location(specAvg)
        specx, specy = np.arange(fitch[0],fitch[1]+1), specAvg - bl
        noise = astropy.stats.biweight_scale(specAvg)

        valid_ix = specy > noise* sig_thres
        pkpos_est  = sum( (specx*specy)[valid_ix] ) / sum(specy[valid_ix])
        sigm_est = np.sqrt(np.sum( ((specx - pkpos_est)**2 * specy)[valid_ix]) / np.sum(specy[valid_ix]))
        pkval_est  = np.max(specy[valid_ix])

        #print(pkpos_est, sigm_est, pkval_est)
        popt, trash = scipy.optimize.curve_fit( gaussianFunc, specx, specy, p0=(pkval_est, pkpos_est, sigm_est))

        plt.subplot(4,1,2)
        plt.step( self.vel[fitch[0]:fitch[1]+1], specy, where='mid')
        yrange = (min(specy), max(specy))
        plt.ylim( (yrange[0]*1.1 - yrange[1]*0.1, yrange[1]*1.1 - yrange[0]*0.1 ))
#        plt.plot ( self.vel[fitch[0]:fitch[1]+1], gaussianFunc(specx, *popt), c='r')
#        plt.plot ( self.vel[fitch[0]:fitch[1]+1], gaussianFunc(specx, *popt), c='r')
        plt.xlabel('V_topo (km/s)')
        plt.grid()

        ch_fm = int(popt[1] - 2*abs(popt[2]))
        ch_to = int(popt[1] + 2*abs(popt[2])) + 1

        plt.fill_between( self.vel[ch_fm:ch_to], np.average(specMax[:,ch_fm:ch_to],axis=0)-bl, color='r', step='mid')

        dvCH = self.vel[1]-self.vel[0]
        dfCH = self.freq[1]-self.freq[0]

        #
        print('### spec fitting result :')
        print('   peak/rms   = {f:1.2}'.format( f= popt[0]/noise ))
        print('   location   = {f:3.5}GHz' .format( f= self.freq[0]+ popt[1]*dfCH))
        print('              = {f:3.5}km/s'.format( f= self.vel[0]+ popt[1]*dvCH))
        print('   FWHM       = {f:3.2}MHz' .format( f= abs(popt[2]*1.e3*dfCH)*2.35))
        print('              = {f:3.2}km/s'.format( f= abs(popt[2]*dvCH)*2.35))

        self.chrangeMap = [ ch_fm, ch_to ]

    def mkMap(self, size=[None, None], grid=[1.,1.], smooth=[5.,5.], velRange=[None, None], searchRadius=2. , mode='int'):


        #################
        ### making map

        v0   = self.vel[0]
        dvCH = self.vel[1]-self.vel[0]
        if velRange !=[None,None] :
            chrangeMap = sorted([ int((v-v0)/dvCH) for v in velRange])
            chrangeMap[1] +=1
        else :
            chrangeMap = self.chrangeMap

        print('### mapping channel range = {f}'.format(f=self.chrangeMap))


        self.smoothFWHM = np.array(smooth)
        smoothSTD = self.smoothFWHM/2.35

        ##########
        # determine mapping range

        mpx = self.azon_intrp
        mpy = self.elon_intrp

        if size[0] == None :
            xmax, xmin = np.floor(mpx.max()/grid[0]) * grid[0],  np.floor(mpx.min()/grid[0])*grid[0]
        else :
            xmax, xmin = size[0]/2., -size[0]/2.

        if size[1] == None :
            ymax, ymin = np.floor(mpy.max()/grid[1]) * grid[1],  np.floor(mpy.min()/grid[1])*grid[1]
        else :
            ymax, ymin = size[1]/2., -size[1]/2.

        nx = int(abs(xmax-xmin)/grid[0]) + 1
        ny = int(abs(ymax-ymin)/grid[1]) + 1

        nz = self.chrangeMap[1] - self.chrangeMap[0]


        ###########
        # create axis

        xaxis = np.linspace( xmin, xmax, num=nx, endpoint=True)
        yaxis = np.linspace( ymin, ymax, num=ny, endpoint=True)

        ###########
        #

        imagemap  = np.zeros((nx,ny,nz))
        weightmap = np.zeros((nx,ny))

        ###########
        # gaussian

        gsbeam = lambda x, w :  np.exp( -(x/w)**2/2.)
        obs_sp = np.sum(self.onoffspecs[:,:, self.chrangeMap[0]:self.chrangeMap[1]], axis=1)
        kpix    = float(self.smoothFWHM[0])/grid[0], float(self.smoothFWHM[1])/grid[1]
        search  = int(kpix[0]*searchRadius), int(kpix[1]*searchRadius)
        print(kpix, search)

        ##########

        for i in tqdm(range(obs_sp.shape[0])):
            xix = np.abs( mpx[i] - xaxis).argmin()
            yix = np.abs( mpy[i] - yaxis).argmin()

            xr = max([xix-search[0],0]), min([xix+search[0],nx-1])
            yr = max([yix-search[1],0]), min([yix+search[1],ny-1])

            wtx = gsbeam( xaxis[xr[0]:xr[1]] - mpx[i], smoothSTD[0])[:,np.newaxis]
            wty = gsbeam( yaxis[yr[0]:yr[1]] - mpy[i], smoothSTD[1])[:,np.newaxis]

            wtmap_i = np.dot(wtx, wty.T)
            weightmap[xr[0]:xr[1],yr[0]:yr[1]] += wtmap_i

            imagemap_i = wtmap_i[:,:,np.newaxis] * obs_sp[i]
            imagemap[xr[0]:xr[1],yr[0]:yr[1]] += imagemap_i

            #print(xix, yix, xr, yr)
 #           plt.subplot(1,2,1); plt.imshow( np.rot90(wtmap_i), cmap='jet')
 #           plt.subplot(1,2,2); plt.imshow( np.rot90(weightmap), cmap='jet')

        masked = (weightmap==0.)

        weightmap[masked] = 0.1

        self.imCube  = imagemap/weightmap[:,:,np.newaxis]
        if mode == 'int' :   #integrated intensity
            self.imObj   = np.sum(self.imCube,axis=2) * abs(dvCH)
        elif mode == 'peak' :  #peak intensity
            self.imObj   = np.max    (self.imCube,axis=2)
        else :
            raise Exception(f'unsupported mode {mode}')

#        self.imObj[masked] = None

        self.obspk = np.max(self.imObj)
        obspix= np.where( self.imObj == self.obspk )
        self.obspos= xaxis[obspix[0][0]], yaxis[obspix[1][0]]
        self.obsrms= astropy.stats.biweight_scale( self.imObj[masked==False] )
        self.obssky= astropy.stats.biweight_location( self.imObj[masked==False] )
        self.xaxis=xaxis
        self.yaxis=yaxis
        self.mask = masked

        print("   pos        = {f:}".format(f=self.obspos))
        print("   pk         = {f:.2f}".format(f=self.obspk))
        print("   source/rms = {f:.2f}".format(f=self.obspk/self.obsrms))
        print("   source/sky = {f:.2f}".format(f=self.obspk/self.obssky))

        #fit = plt.figure(figsize=(12,5))
        plt.subplot(4,2,5)
        plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.title('{name}'.format(name={'int':'Integrated Intensity', 'peak':'peak intensity'}[mode]))
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')
        plt.subplot(4,2,6)
        plt.imshow(np.rot90(weightmap)[:,::-1], extent=(xmax, xmin, ymin, ymax), cmap='jet')
        plt.title('weight')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

    def mkFit(self) :
        xmax,xmin = np.max(self.xaxis), np.min(self.xaxis)
        ymax,ymin = np.max(self.yaxis), np.min(self.yaxis)

        p0=(self.obspos[0], self.obspos[1],
            10/2.35, 10/2.35,
            self.obspk , self.obssky )

        popt, trash = scipy.optimize.curve_fit( modelBeam, np.meshgrid(self.xaxis,self.yaxis), self.imObj.flatten() , p0=p0, maxfev=100000)
        modelim = modelBeam(np.meshgrid(self.xaxis, self.yaxis), *popt).reshape(self.imObj.shape)

        #fit = plt.figure(figsize=(12,5))

#        plt.subplot(1,3,1)
#        plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
#        plt.title('obs')
#        plt.xlabel('dAz (arcsec)')
#        plt.ylabel('dEl (arcsec)')

        plt.subplot(4,2,7)
        plt.title('model')
        plt.imshow(np.rot90(modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        plt.subplot(4,2,8)
        plt.title('residual')
        plt.imshow(np.rot90(self.imObj - modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        print('### beam fitting result :')
        self.fitloc  = np.array([popt[0], popt[1]])
        self.fitbeam = np.array([popt[2], popt[3]])*2.35
        self.fitbeam_corr = np.sqrt( self.fitbeam**2 - self.smoothFWHM**2)
        self.fitpeak = (popt[4] - popt[5])
        self.fitpeak_err = np.sqrt(np.sqrt(np.diag(trash))[4]**2 + np.sqrt(np.diag(trash))[5]**2)
        print('   location       : {x:1.3f}",  {y:1.3f}"'.format(x=self.fitloc[0], y=self.fitloc[1]))
        print('   FWHM(raw)      : {x:1.3f}" x {y:1.3f}"'.format(x=self.fitbeam[0], y=self.fitbeam[1]))
        print('   FWHM(corrected): {x:1.3f}" x {y:1.3f}"'.format(x=self.fitbeam_corr[0], y=self.fitbeam_corr[1]))

    def mkContFit(self) :
        xmax,xmin = np.max(self.xaxis), np.min(self.xaxis)
        ymax,ymin = np.max(self.yaxis), np.min(self.yaxis)

        p0=(0.0, 0.0,
            10/2.35, 10/2.35,
            self.obspk , self.obssky )

        popt, trash = scipy.optimize.curve_fit( modelBeam, np.meshgrid(self.xaxis,self.yaxis), self.imObj.flatten() , p0=p0,maxfev=100000)
        modelim = modelBeam(np.meshgrid(self.xaxis, self.yaxis), *popt).reshape(self.imObj.shape)

        #fit = plt.figure(figsize=(12,5))

#        plt.subplot(1,3,1)
#        plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
#        plt.title('obs')
#        plt.xlabel('dAz (arcsec)')
#        plt.ylabel('dEl (arcsec)')

        plt.subplot(2,2,3)
        plt.title('model')
        plt.imshow(np.rot90(modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        plt.subplot(2,2,4)
        plt.title('residual')
        plt.imshow(np.rot90(self.imObj - modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        print('### beam fitting result :')
        self.fitloc  = np.array([popt[0], popt[1]])
        self.fitbeam = np.array([popt[2], popt[3]])*2.35
        self.fitbeam_corr = np.sqrt( self.fitbeam**2 - self.smoothFWHM**2)
        self.fitpeak = (popt[4] - popt[5])
        self.fitpeak_err = np.sqrt(np.sqrt(np.diag(trash))[4]**2 + np.sqrt(np.diag(trash))[5]**2)
        print('   location       : {x:1.3f}",  {y:1.3f}"'.format(x=self.fitloc[0], y=self.fitloc[1]))
        print('   FWHM(raw)      : {x:1.3f}" x {y:1.3f}"'.format(x=self.fitbeam[0], y=self.fitbeam[1]))
        print('   FWHM(corrected): {x:1.3f}" x {y:1.3f}"'.format(x=self.fitbeam_corr[0], y=self.fitbeam_corr[1]))

    def PolSB2freq(self,Pol='B',SB='LSB'):
        if SB == 'LSB':
            self.useFreq = self.fline - self.fif[1] + 2.5/2.
        elif SB == 'USB':
            self.useFreq = self.fline + self.fif[1] - 2.5/2.

        if Pol == 'A':
            self.usePol = 'pol2'
        elif Pol == 'B':
            self.usePol = 'pol1'


    def LinePointingQlook(self, useFreq=129.36324):
        self.loadlmttpm()
        self.loadxffts(useVel=0, useVelRange=500, useFreq=useFreq)
        if not self.LineSkipFlag:
            self.getOnOffSpecs(binning=8 )
            self.showSpectrum(nTimeAvg=4,sig_thres=1,vrange=[None,None])
            self.mkMap()
            self.mkFit()
            self.outputfig = plt.gcf()


    def ContPointingQlook(self, useFreq=129.36324, Pol='B', SB='LSB'):
        self.loadlmttpm()

        self.PolSB2freq(Pol=Pol,SB=SB)
        self.loadxffts(useVel=0, useVelRange=12000, useFreq=self.useFreq,usePol= [self.usePol,])

        self.getOnOffSpecs(binning=256)
        #self.showContSpectrum(nTimeAvg=4,sig_thres=1,vrange=[None,None])
        self.mkContMap()
        self.mkContFit()
        self.outputfig = plt.gcf()
















######
