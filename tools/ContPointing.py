import xarray as xr
import numpy as np
from astropy.stats import sigma_clipped_stats
#import fmflow as fm
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
import netCDF4 as nc
import datetime as dt
import Lib_b4r as Lib
import makedatalist as md
from tqdm import tqdm
import os
import sys
import glob
from scipy.interpolate import interp1d
import astropy.stats
import fmflow as fm

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, LSR, AltAz, ICRS
import astropy.units as u
import scipy

globBaseDir = '.'

def gaussianFunc(x, pk, pos, sigm) :
    return pk * np.exp( - ( ( ( x - pos)/sigm)**2/2.))

def modelBeam(xymesh, posx, posy, widthx, widthy, peak, sky) : #modelPar) :
    gsbeam = lambda x, w :  np.exp( -(x/w)**2/2.)
    axis_x = xymesh[0][0]
    axis_y = xymesh[1][:,0]
    modelx = gsbeam( axis_x - posx, widthx)[:,np.newaxis]
    modely = gsbeam( axis_y - posy, widthy)[:,np.newaxis]
    return (np.dot(modelx, modely.T) * peak + sky).flatten()

def convrt2nc(path):
    if os.path.exists(path+'.nc'):
        return path+'.nc'
    else:
        Lib.xffts2netcdf(path,path+'.nc')
        return path+'.nc'


class ContPointing:

    def __init__(self, obsnum, calnum, XFFTSnum) :
        self.obsnum = obsnum
        self.calnum = calnum

        xfftsnum = str(XFFTSnum).zfill(2)

        print(globBaseDir+ 'xffts_links/*'+str(calnum)+'*.xfftsx*.'+xfftsnum)
        self.data_sci = xr.open_dataset(convrt2nc(glob.glob(globBaseDir + 'xffts_links/*'+str(obsnum)+'*.xfftsx*.'+xfftsnum)[0]))
        self.data_cal = xr.open_dataset(convrt2nc(glob.glob(globBaseDir + 'xffts_links/*'+str(calnum)+'*.xfftsx*.'+xfftsnum)[0]))
        self.data_ant = glob.glob(globBaseDir + 'lmttpm/lmttpm_*'+str(obsnum)+'_*.nc')[0]

    def findSkyLevel(self):
        Psky = fm.array(self.data_sci['array'].copy().values/(self.data_sci['integtime'].copy().values)[:,None]).mean('ch')
        date = self.data_sci['date'].copy().values
        t_sci = np.array([s[:-4] for s in date], 'datetime64[us]')
        flag1 = t_sci < t_sci[-1]
        #flag2 = data_sci['bufpos'] == 'ON'

        flag = flag1 #& flag2

        model_sky = models.Polynomial1D(2)
        model_sky.c0 = 0.0
        model_sky.c1 = -1.0
        model_sky.c2 = 1.0
        pfit = fitting.LinearLSQFitter()
        opfit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip,niter=15, sigma=2.0)
        #fitted_sky = pfit(model_sky, t_sci[flag], Psky[flag])
        fitted_sky,filtered_data = opfit(model_sky, t_sci[flag]-t_sci[0], Psky[flag])
        #opfit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip,niter=3, sigma=3.0)
        #fitted_sky,filtered_data = opfit(model_sky,t_sci[flag][~filtered_data],Psky[flag][~filtered_data])
        #print(filtered_data)
        #print(fitted_sky.c1,fitted_sky.c0)
        self.t_sci = t_sci
        self.fitted_sky = fitted_sky
        self.Psky = Psky
        self.flag = flag
        self.filtered_data = filtered_data

    def loadlmttpm(self) :
        print(self.data_ant)
        ncfile= nc.Dataset(self.data_ant)
        ncprop = md.getNCprops(self.data_ant)

        utc = dt.timezone(dt.timedelta(hours=+0), 'utc')

        nctime = ncfile.variables['Data.TelescopeBackend.TelTime'][:]
        bufpos = ncfile.variables['Data.TelescopeBackend.BufPos'][:]

        self.az_map = ncfile['Data.TelescopeBackend.TelAzMap'][:]
        self.el_map = ncfile['Data.TelescopeBackend.TelElMap'][:]

        self.fif    = [ncfile['Header.B4r.If1Freq'][:][0], ncfile['Header.B4r.If2Freq'][:][0]]
        self.fline  = ncfile['Header.B4r.LineFreq'][0]
        self.Tamb = 273. + ncfile['Header.Weather.Temperature'][0]
        self.srcename = ncprop['srcname']

        self.dAZ = np.rad2deg(ncfile['Header.PointModel.AzUserOff'][0])*3600
        self.dEL = np.rad2deg(ncfile['Header.PointModel.ElUserOff'][0])*3600.

        self.antonpos = bufpos == 0
        self.antoffpos= bufpos == 1

        self.t_ant = np.array([dt.datetime.fromtimestamp(t, tz=utc) for t in nctime], 'datetime64[us]')

        ncfile.close()

    def RSkyCal(self):

        self.findSkyLevel()
        self.P_cal = fm.array(self.data_cal['array'].copy().values/(self.data_cal['integtime'].copy().values)[:,None]).mean('ch')


        #Psky = findSkyLevel(data_sci)

        # TSYS

        self.Tsys = self.Tamb * (self.P_cal[self.data_cal['chopperpos']=='S']).mean('t')/((self.P_cal[self.data_cal['chopperpos']=='R']).mean('t')-(self.P_cal[self.data_cal['chopperpos']=='S']).mean('t'))

        # R-sky
        self.Ta = self.Tsys * (self.Psky - self.fitted_sky(self.t_sci-self.t_sci[0]))/self.fitted_sky(self.t_sci-self.t_sci[0])

    def mkMap(self, size=[None, None], grid=[3.,3.], smooth=[3.,3.], velRange=[None, None], searchRadius=2. , mode='int'):
        azon_raw = self.az_map[self.antonpos] * 3600 * 180/np.pi     #rad -> sec
        elon_raw = self.el_map[self.antonpos] * 3600 * 180/np.pi
        tant_on = self.t_ant.astype(float)[self.antonpos]
        self.smoothFWHM = np.array(smooth)
        smoothSTD = self.smoothFWHM/2.35
        #self.chrangeMap = [0,2**15-1]

        mpx  = interp1d(tant_on - tant_on[0], azon_raw, fill_value='extrapolate')(self.t_sci - self.t_sci[0])
        mpy  = interp1d(tant_on - tant_on[0], elon_raw, fill_value='extrapolate')(self.t_sci - self.t_sci[0])

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

        #nz = self.chrangeMap[1] - self.chrangeMap[0]
        nz = 1


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
        #print(self.Ta.copy().values.shape)
        gsbeam = lambda x, w :  np.exp( -(x/w)**2/2.)
        obs_sp = np.array([self.Ta.copy().values]).T
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
            self.imObj   = np.sum(self.imCube,axis=2)
        elif mode == 'peak' :  #peak intensity
            self.imObj   = np.max    (self.imCube,axis=2)
        else :
            raise Exception(f'unsupported mode {mode}')

#        self.imObj[masked] = None
        print('done')
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
        plt.subplot(3,2,3)
        plt.imshow(np.rot90(self.imObj)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.title('{name}'.format(name={'int':'Integrated Intensity', 'peak':'peak intensity'}[mode]))
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')
        plt.subplot(3,2,4)
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

        plt.subplot(3,2,5)
        plt.title('model')
        plt.imshow(np.rot90(modelim)[:,::-1], extent=(xmax, xmin, ymin, ymax), vmin=self.obssky, vmax=self.obspk, cmap='jet')
        plt.xlabel('dAz (arcsec)')
        plt.ylabel('dEl (arcsec)')

        plt.subplot(3,2,6)
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

    def efficiency(self,D = 50.,SB='LSB',pol='Bpol'):
        if SB == 'LSB':
            self.useFreq = self.fline - self.fif[1] + 2.5/2.
        elif SB == 'USB':
            self.useFreq = self.fline + self.fif[1] - 2.5/2.

        HPBW_x = self.fitbeam_corr[0]
        HPBW_y = self.fitbeam_corr[1]
        lobs = 3.0e8/self.useFreq/1000.
        #Tpeak = 3.84641
        #amp = 4.13259

        # Brightness temperature of Uranus (Griffin & Orton 1993)
        a0 = -795.694
        a1 = 845.179
        a2 = -288.946
        a3 = 35.200
        Turanus = a0 + a1*np.log10(lobs) + a2*(np.log10(lobs))**2 + a3*(np.log10(lobs))**3 # [K]
        Turanus = Turanus*0.931 # Sayers+12

        # main beam efficiency
        theta_eq = 3.72 # [arcsec]
        theta_pol = 3.64 # [arcsec]

        #bff = 1 - np.exp(-np.log(2) * theta_eq * theta_pol / fwhm_x / fwhm_y)
        #eta = Tpeak/Turanus/bff
        bff = 1 - np.exp(-np.log(2) * theta_eq * theta_pol / HPBW_x/HPBW_y)
        eta = self.fitpeak/Turanus/bff
        eta_err = self.fitpeak_err/Turanus/bff

        # aperture efficiency
        etaa = eta/1.2 # assuming -12 dB edge taper
        etaa_err = eta_err/1.2

        # Jy/K
        JpK = HPBW_x*HPBW_y/13.6/(lobs*1e-3)**2

        self.eta = eta
        self.etaa = etaa
        self.JpK = JpK
        self.eta_err = eta_err
        self.etaa_err = etaa_err
        self.Uranus_K = Turanus
        self.Uranus_Jy = Turanus*JpK*bff
        self.bff = bff


    def ContPointingQlook(self):
        self.findSkyLevel()
        self.loadlmttpm()
        self.RSkyCal()

        self.mkMap()
        self.mkFit()
        self.efficiency()
        plt.subplot(3,1,1)
        POS = '{x:1.3f}",{y:1.3f}"'.format(x=self.fitloc[0], y=self.fitloc[1])
        FWHM = '{x:1.3f}"x{y:1.3f}"'.format(x=self.fitbeam[0], y=self.fitbeam[1])
        FWHM_corr = '{x:1.3f}"x{y:1.3f}"'.format(x=self.fitbeam_corr[0], y=self.fitbeam_corr[1])
        dAZEL = '{x:1.3f}",{y:1.3f}"'.format(x=self.dAZ+self.fitloc[0], y=self.dEL+self.fitloc[1])
        etaa = '{0:1.3f}'.format(self.etaa)+r'$\pm$'+'{0:1.3f}'.format(self.etaa_err)
        intens = '(Uranus: '+'{0:1.3f}'.format(self.Uranus_K)+'K/'+'{0:1.3f}'.format(self.Uranus_Jy)+'Jy) bff='+'{0:1.3f}'.format(self.bff)

        titleinfo = 'ID'+str(self.obsnum)+' Pos:'+POS+' FWHM:'+FWHM+' FWHM(corr):'+FWHM_corr+'\n'+r'$\eta_a=$'+etaa+' '+intens+' NewdAZEL: '+dAZEL
        plt.plot(self.t_sci[self.flag]-self.t_sci[self.flag][0],self.Psky[self.flag]/1.0e7,'o')
        plt.plot(self.t_sci[self.flag][self.filtered_data]-self.t_sci[0],self.Psky[self.flag][self.filtered_data]/1.0e7,'o')
        plt.plot(self.t_sci[self.flag]-self.t_sci[0],self.fitted_sky(self.t_sci[self.flag]-self.t_sci[0])/1.0e7)
        plt.xlabel('time [sec]')
        plt.title(titleinfo)
        print('done')
        self.outputfig = plt.gcf()
