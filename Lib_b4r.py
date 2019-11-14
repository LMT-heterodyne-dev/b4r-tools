#!/usr/bin/env python3

### common modules ###
import numpy as np
import os

### functions ###

################################################################################
# date2mjs

def date2mjs(t):
	# t = '2018-10-05T05:37:50.994900'
	y  = int(str(t)[0:4])
	m  = int(str(t)[5:7])
	d  = int(str(t)[8:10])
	hh = int(str(t)[11:13])
	mm = int(str(t)[14:16])
	ss = float(str(t)[17:])

	return (int(365.25*y)+y//400-y//100.+int(30.59*(m-2))+d-678912)*24.*3600.+hh*3600+mm*60+ss

# end date2mjs
################################################################################
# xffts2netcdf

def xffts2netcdf(infile,outfile):
	#os.system('python xffts2netcdf.py '+infile+' '+outfile)

	"""Convert XFFTS binary file to netCDF.
	Usage: $ python xffts2netcdf.py <xffts binary file> <netcdf>
	"""

	# standard library
	import os
	import re
	import sys
	from collections import deque, OrderedDict
	from pathlib import Path
	from struct import Struct


	# dependent packages
	import numpy as np
	from netCDF4 import Dataset
	from tqdm import tqdm


	# module constants
	CONFIG = (
	    ('date', '28s'),
	    ('junk1', '4s'),
	    ('obsnum', 'q'),
	    ('bufpos', '8s'),
	    ('scanmode', '8s'),
	    ('chopperpos', '8s'),
	    ('scancount', 'q'),
	    ('speccount', 'q'),
	    ('integtime', 'l'),
	    ('junk2', '172s'),
	    ('array', 'f', 32768)
	)
	''' Old config
	CONFIG = [
	    ('date', '28s'),
	    ('junk1', '4s'),
	    ('obsnum', 'l'),
	    ('bufpos', '8s'),
	    ('scanmode', '8s'),
	    ('chopperpos', '8s'),
	    ('scancount', 'l'),
	    ('speccount', 'l'),
	    ('integtime', 'l'),
	    ('junk2', '184s'),
	    ('array', 'f', 32768)
	]
	'''


	# functions and classes
	class Struct2NetCDF:
	    def __init__(self, path, config, overwrite=False,
	                 byteorder='<', unlimited_dim='t'):
	        """Initialize structure-to-netCDF converter.
	        Args:
	            path (str or path object): Path of netCDF to be created.
	            config (tuple of items): Tuple that contains tuple items
	                (name, format character, shape) for member in a structure.
	            overwrite (bool, optional): If True, the converter overwrites
	                an existing netCDF of the same name. Default is False.
	            byteorder (str, optional): Format character that indicates
	                the byte order of a binary file. Default is for little
	                endian ('<'). Use '>' for big endian instead.
	        """
	        # instance variables
	        self.path = Path(path).expanduser()
	        self.byteorder = byteorder
	        self.unlimited_dim = unlimited_dim

	        # check netCDF existance
	        if self.path.exists() and not overwrite:
	            raise FileExistsError(f'{self.path} already exists')

	        # initialization
	        self.config = self.parse_config(config)
	        self.struct = self.create_struct()
	        self.dataset = self.create_empty_dataset()

	        # initialize writing counter
	        self.n_write = 0

	        # aliases
	        self.close = self.dataset.close
	        self.readsize = self.struct.size

	    def write(self, binary):
	        """Convert binary data compatible to netCDF format and write it."""
	        if len(binary) == 0:
	            raise EOFError('Reached the end of file')

	        assert len(binary) == self.readsize
	        data = deque(self.struct.unpack(binary))

	        for name, (fmt, shape) in self.config.items():
	            variable = self.dataset[name]

	            # in the case of no additional dimensions
	            if shape == (1,):
	                variable[self.n_write] = data.popleft()
	                continue

	            # otherwise
	            flat = [data.popleft() for i in range(np.prod(shape))]
	            variable[self.n_write] = np.reshape(flat, shape)

	        self.n_write += 1

	    def create_struct(self):
	        """Create structure object for unpacking binary string."""
	        joined = ''.join(c*np.prod(s) for c, s in self.config.values())
	        return Struct(self.byteorder + joined)

	    def create_empty_dataset(self):
	        """Create empty netCDF dataset according to structure config."""
	        # add unlimited dimension
	        empty = Dataset(self.path, 'w')
	        empty.createDimension(self.unlimited_dim)

	        # add variables and additional dimensions
	        for name, (fmt, shape) in self.config.items():
	            dtype = self.convert_fmt_to_dtype(fmt)

	            # in the case of no additional dimensions
	            if shape == (1,):
	                dims = (self.unlimited_dim,)
	                empty.createVariable(name, dtype, dims)
	                continue

	            # otherwise
	            n_dims = len(shape)
	            dims = [f'{name}_dim{i}' for i in range(n_dims)]

	            for i in range(n_dims):
	                dim, size = dims[i], shape[i]
	                empty.createDimension(dim, size)

	            dims = (self.unlimited_dim,) + tuple(dims)
	            empty.createVariable(name, dtype, dims)

	        return empty

	    @staticmethod
	    def parse_config(config):
	        """Parse structure config to ordered dictionary."""
	        parsed = OrderedDict()

	        for item in config:
	            if not 2 <= len(item) <= 3:
	                raise ValueError(item)

	            name, fmt = item[:2]
	            shape = (item[2:] or (1,))[0]
	            if isinstance(shape, int):
	                shape = (shape,)

	            assert isinstance(name, str)
	            assert isinstance(fmt, str)
	            assert isinstance(shape, tuple)

	            parsed[name] = fmt, shape

	        return parsed

	    @staticmethod
	    def convert_fmt_to_dtype(fmt):
	        """Convert format character to NumPy dtype object."""
	        if re.search('[bhil]', fmt, re.I):
	            return np.int32
	        elif re.search('[q]', fmt, re.I):
	            return np.int64
	        elif re.search('[ef]', fmt):
	            return np.float32
	        elif re.search('[d]', fmt):
	            return np.float64
	        elif re.search('[csp]', fmt):
	            return np.str
	        elif re.search('[?]', fmt):
	            return np.bool
	        else:
	            raise ValueError(fmt)

	    def __enter__(self):
	        """Special method for with statement."""
	        return self

	    def __exit__(self, exc_type, exc_value, traceback):
	        """Special method for with statement."""
	        self.close()


	def main(infile,outfile):
		"""Main function for command line tool."""
		#args = sys.argv[1:]
		args = [infile,outfile]

		xffts  = Path(args[0]).expanduser()
		if len(args) == 1:
			netcdf = Path(f"{xffts}.nc")
		elif len(args) == 2:
			netcdf = Path(args[1]).expanduser()
		else:
			print(__doc__)
			sys.exit(0)

		with xffts.open('rb') as f, Struct2NetCDF(netcdf, CONFIG) as g:
			filesize = xffts.stat().st_size
			readsize = g.readsize

			assert not filesize % readsize
			n_struct = int(filesize / readsize)

			for i in tqdm(range(n_struct)):
			    binary = f.read(readsize)
			    g.write(binary)

	main(infile,outfile)

# end xffts2netcdf
################################################################################
# pathOFnc

def pathOFnc(path_cal_raw,path_sci_raw):

    if path_cal_raw[-3:] == '.nc':
        path_cal = path_cal_raw
    else:
        path_cal  = path_cal_raw+'.nc'
        #os.system('rm -rf '+path_cal)
        xffts2netcdf(path_cal_raw,path_cal)

    if path_sci_raw[-3:] == '.nc':
        path_sci = path_sci_raw
    else:
        path_sci  = path_sci_raw+'.nc'
        #os.system('rm -rf '+path_sci)
        xffts2netcdf(path_sci_raw,path_sci)

    return path_cal,path_sci

# end pathOFnc
################################################################################
# create_freq

def create_freq(path_ant,sideband,tBW=2.5,nchan=2**15):

    # modules
    import xarray as xr

    # get frequency array
    ant = xr.open_dataset(path_ant)
    LO1 = ant['Header.B4r.LineFreq'].values[0]
    LO2 = ant['Header.B4r.If2Freq'].values[0]
    if (sideband=='USB'):
        return np.linspace(LO1+LO2, LO1+LO2-tBW,nchan)
    elif (sideband=='LSB'):
        return np.linspace(LO1-LO2, LO1-LO2+tBW,nchan)

# end create_freq
################################################################################
# read_rawdata_science

def read_rawdata_science(path_sci,freq,timestampFlag=True,cal=True,Tsys=None,threshold_despiking=100.):

    # modules
    import xarray as xr
    import fmflow as fm

    with xr.open_dataset(path_sci) as ds:
        array     = ds['array'].copy().values
        integtime = ds['integtime'].copy().values
        scantype  = ds['bufpos'].copy().values
        scanno    = ds['scancount'].copy().values
        date      = ds['date'].copy().values

	# flagging invalid datetime
    if timestampFlag:
        t_sci = np.array([s[:-4] for s in date], 'datetime64[us]')
        flag = t_sci < t_sci[-1]
        array     = array[flag]
        integtime = integtime[flag]
        scantype  = scantype[flag]
        scanno    = scanno[flag]

    state_id = (scantype == 'ON').astype('int32')
    tcoords  = {'scanno': scanno, 'scantype': scantype, 'state_id': state_id}
    chcoords = {'fsig': freq}
    state_id = {'state_id': state_id}

    P = fm.array(array/integtime[:,None],
    			tcoords=tcoords, chcoords=chcoords)

    #print(P)
    if cal:
        Pon  = P[scantype=='ON']
        Poff = P[scantype=='REF']

        Tcal = fm.zeros_like(Pon)
        for no in np.unique(Pon.scanno):
            # ON array of single scan
            Pon_ = Pon[Pon.scanno==no]

            # OFF array for single scan
            Poff_l = Poff[Poff.scanno==no-1].mean('t')
            Poff_r = Poff[Poff.scanno==no+1].mean('t')

            if np.all(np.isnan(Poff_r)):
            	Poff_ = Poff_l
            else:
            	Poff_ = (Poff_l+Poff_r) / 2

            #print(Tsys)
            #print(Pon_)
            #print(Poff_)
            # calibrated array
            Tcal_ = Tsys * (Pon_-Poff_) / (Poff_)

            Tcal[Tcal.scanno==no] = Tcal_.T

        Tcal.values[np.abs(Tcal)>threshold_despiking] = 0
        return Tcal
    else:
        return P

# end read_rawdata_science
################################################################################
# cal_Tsys_from_rawdata_R

def cal_Tsys_from_rawdata_R(path_cal,path_ant,freq,CASAoffset=3506716797.):

    # modules
    import xarray as xr
    import fmflow as fm

    with xr.open_dataset(path_cal) as ds:
        array      = ds['array'].copy().values
        integtime  = ds['integtime'].copy().values
        chopperpos = ds['chopperpos'].copy().values
        date       = ds['date'].copy().values

    ant = xr.open_dataset(path_ant)
    Tamb = 273.0 + ant['Header.Weather.Temperature'].copy().values
    #Tamb = 273.0 + ant['Header.B4r.LoadAmbientTemp'].copy().values

    t_tsys = np.array([s[:-4] for s in date], 'datetime64[us]')[0]
    Tsys_time = CASAoffset + (t_tsys - np.array('1970-01-01T00:00:00.000000',dtype='datetime64[us]'))/np.timedelta64(1,'s')


    Pr = fm.array((array/np.array([integtime]).T)[chopperpos=='R']).mean('t')
    Psky = fm.array((array/np.array([integtime]).T)[chopperpos=='S']).mean('t')
    chcoords = {'fsig': freq}

    Tsys = Tamb*Psky/(Pr-Psky)

    return Tsys, Tsys_time

# end cal_Tsys_from_rawdata_R
################################################################################
# getOBSinfo

def getOBSinfo(path_ant):

    # modules
    import xarray as xr

    ant = xr.open_dataset(path_ant)
    sourcename_ = ant['Header.Source.SourceName'].copy().values.astype('|U')
    sourcename = np.array([sourcename_,sourcename_])[0].replace(' ','')
    ra = ant['Header.Source.Ra'].copy().values[0]/np.pi*180.
    dec = ant['Header.Source.Dec'].copy().values[0]/np.pi*180.
    sysvel = ant['Header.Source.Velocity'].copy().values + 0.
    project_ = ant['Header.Dcs.ProjectId'].copy().values.astype('|U')
    project = np.array([project_,project_])[0].replace(' ','')
    observer_ = ant['Header.Telescope.Operator'].copy().values.astype('|U')
    observer = np.array([observer_,observer_])[0].replace(' ','')

    return ra,dec,sysvel,sourcename,project,observer

# end getOBSinfo
################################################################################
# data2coord

def data2coord(path_sci,path_ant,nrows=False,cal=True):

	# modules
	import xarray as xr
	from scipy.interpolate import interp1d
	from datetime import datetime as dt
	from datetime import timezone as tz

	with xr.open_dataset(path_sci) as ds:
	    array     = ds['array'].copy().values
	    integtime = ds['integtime'].copy().values
	    scantype  = ds['bufpos'].copy().values
	    scanno    = ds['scancount'].copy().values
	    date      = ds['date'].copy().values

	# flagging invalid datetime
	t_sci = np.array([s[:-4] for s in date], 'datetime64[us]')
	flag = t_sci < t_sci[-1]

	array     = array[flag]
	integtime = integtime[flag]
	scantype  = scantype[flag]
	scanno    = scanno[flag]
	date      = date[flag]
	t_sci     = t_sci[flag]

	#ON
	#t_sci = t_sci[scantype == 'ON']

	with xr.open_dataset(path_ant) as ds:
	    x_ant = np.rad2deg(ds['Data.TelescopeBackend.TelAzSky'].values)
	    y_ant = np.rad2deg(ds['Data.TelescopeBackend.TelElSky'].values)
	    unix  = ds['Data.TelescopeBackend.TelTime'].values # unix time
	    SourceName = bytes(ds['Header.Source.SourceName'].values).decode()

	t_ant = np.array([dt.fromtimestamp(t,tz.utc) for t in unix], 'datetime64[us]')

	dt_sci = (t_sci - t_ant[0]).astype(np.float64)
	dt_ant = (t_ant - t_ant[0]).astype(np.float64)

	x_sci = interp1d(dt_ant, x_ant, fill_value='extrapolate')(dt_sci)
	y_sci = interp1d(dt_ant, y_ant, fill_value='extrapolate')(dt_sci)

	if cal:
		x_sci = x_sci[scantype=='ON']
		y_sci = y_sci[scantype=='ON']
		t_sci = t_sci[scantype=='ON']

	if nrows:
	    return x_sci, y_sci, t_sci, x_sci.shape[0]
	else:
	    return x_sci, y_sci, t_sci

# end data2coord
################################################################################
# AZEL2RADEC

def AZEL2RADEC(az,
               el,
               utc,
               lon = -97.31461167440136,
               lat = 18.98607157170612,
               height = 4640.):

	# modules
	from astropy import coordinates
	from astropy.time import Time
	import astropy.units as u

	site = coordinates.EarthLocation.from_geodetic(lon=lon,lat=lat,height=height)
	c = str(az) + ' ' + str(el)

	radec = coordinates.SkyCoord(c, unit=(u.deg,u.deg),frame='altaz',obstime=Time(utc.astype('|U')),location=site).transform_to('icrs')
	return radec.ra.value, radec.dec.value

################################################################################
# PointingINFO

def PointingINFO(path_sci,path_ant,CASAoffset=3506716797.,cal=True):

	# modules

	x,y,t,nrows = data2coord(path_sci,path_ant,nrows=True,cal=cal)

	direction = np.zeros((nrows,2),dtype='float64')
	time = np.full(nrows,5.04545e+09,dtype='float64')

	for i in range(nrows):
		ra,dec = AZEL2RADEC(x[i],y[i],t[i])
		direction[i][0] = ra
		direction[i][1] = dec
		#time[i] = CASAoffset + (t[i] - np.array('1970-01-01T00:00:00.000000',dtype='datetime64[us]'))/np.timedelta64(1,'s')
		time[i] = date2mjs(t[i])

	return direction, time

# end PointingINFO
################################################################################
# loadB4Rdata

def loadB4Rdata(path_cal_raw,path_sci_raw,path_ant,sideband,cal=True,blsub=False):

	# modules
	import xarray as xr

	# load
	path_cal, path_sci = pathOFnc(path_cal_raw,path_sci_raw)
	freq = create_freq(path_ant,sideband)

	Tsys, Tsys_time = cal_Tsys_from_rawdata_R(path_cal,path_ant,freq)
	P = read_rawdata_science(path_sci,freq,timestampFlag=True,Tsys=Tsys,cal=cal)

	if blsub:
		for no in np.unique(P.scanno):
		    P_ = P[P.scanno==no]
		    N = P_.shape[0]
		    t = xr.DataArray(np.arange(N), dims='t')
		    Pbl_ = ((N-1-t)/(N-1)*P_[0]).values + (t/(N-1)*P_[-1]).values

		    P[P.scanno==no] -= Pbl_

	return P, Tsys, Tsys_time

# end loadB4Rdata
################################################################################
# np.array(['Hz'],dtype='|S3'data

def loadB4Rfulldata(path_cal_raw_head,path_sci_raw_head,path_antlog,chbin,cal=True,blsub=False):

    # modules
	import fmflow as fm

    # Params
	num = ['.01','.03','.02','.04']
	sideband = ['LSB','LSB','USB','USB']

	for i in [0,1]:
		if os.path.exists(path_cal_raw_head+num[i]+'.nc'):
			path_cal_raw = path_cal_raw_head+num[i]+'.nc'
		else:
			path_cal_raw = path_cal_raw_head+num[i]

		if os.path.exists(path_sci_raw_head+num[i]+'.nc'):
			path_sci_raw = path_sci_raw_head+num[i]+'.nc'
		else:
			path_sci_raw = path_sci_raw_head+num[i]

		P_, Tsys_, Tsys_time = loadB4Rdata(path_cal_raw,path_sci_raw,path_antlog,sideband[i],cal=cal,blsub=blsub)
		if i==0:
			P_LSB_X = fm.chbinning(P_, chbin).copy().values
			Tsys_LSB_X = fm.chbinning(Tsys_, chbin).copy().values

			state_id = P_['state_id'].copy().values
			freq_LSB = fm.chbinning(P_, chbin)['fsig'].copy().values
			nrows = P_LSB_X.shape[0]
			nchan = P_LSB_X.shape[1]

		else:
			P_LSB_Y = fm.chbinning(P_, chbin).copy().values
			Tsys_LSB_Y = fm.chbinning(Tsys_, chbin).copy().values

	for i in [2,3]:
		if os.path.exists(path_cal_raw_head+num[i]+'.nc'):
			path_cal_raw = path_cal_raw_head+num[i]+'.nc'
		else:
			path_cal_raw = path_cal_raw_head+num[i]

		if os.path.exists(path_sci_raw_head+num[i]+'.nc'):
			path_sci_raw = path_sci_raw_head+num[i]+'.nc'
		else:
			path_sci_raw = path_sci_raw_head+num[i]

		P_, Tsys_, Tsys_time = loadB4Rdata(path_cal_raw,path_sci_raw,path_antlog,sideband[i],cal=cal,blsub=blsub)
		if i==2:
			P_USB_X = fm.chbinning(P_, chbin).copy().values
			Tsys_USB_X = fm.chbinning(Tsys_, chbin).copy().values
			freq_USB = fm.chbinning(P_, chbin)['fsig'].copy().values
		else:
			P_USB_Y = fm.chbinning(P_, chbin).copy().values
			Tsys_USB_Y = fm.chbinning(Tsys_, chbin).copy().values

	P = np.zeros((nrows,2,2,nchan),dtype='float64')
	P[:,0,0] = P_LSB_X
	P[:,0,1] = P_LSB_Y
	P[:,1,0] = P_USB_X
	P[:,1,1] = P_USB_Y


	Tsys = np.zeros((2,2,nchan),dtype='float64')
	Tsys[0,0] = Tsys_LSB_X
	Tsys[0,1] = Tsys_LSB_Y
	Tsys[1,0] = Tsys_USB_X
	Tsys[1,1] = Tsys_USB_Y

	freq = np.zeros((2,nchan),dtype='float64')
	freq[0] = freq_LSB
	freq[1] = freq_USB

	return P, Tsys, freq, state_id, Tsys_time

# end loadB4Rfulldata
################################################################################
# addKeywords
def addKeywords(colname,keyword_list,type_list,value_list,f):

    if not isinstance(keyword_list, list):
        keyword_list = [keyword_list]

    if not isinstance(type_list, list):
        type_list = [type_list]

    if not isinstance(value_list, list):
        value_list = [value_list]

    f.write('.keywords '+colname+'\n')

    for keyword,type,value in zip(keyword_list,type_list,value_list):
        f.write(keyword+' '+type+' '+value+'\n')

    f.write('.endkeywords'+'\n')


# end addKeywords
################################################################################
# makeMAIN

def makeMAIN(tablename,outputfilename,specdata,time,state_id,texp=0.2,tBW=2.5e9):

	'''
	specdata:    (nrow,nspw,npol,nchan) array
	time:        (nrow,)
	'''

	# modules
	import casacore.tables as tb

    # params
	nrow  = specdata.shape[0]
	nspw  = specdata.shape[1]
	npol = specdata.shape[2]
	nchan  = specdata.shape[3]

	weight = tBW/float(nchan) * texp
	sigma = (tBW/float(nchan) * texp)**-0.5

	ind_spw = (np.linspace(0,2*nrow-1,2*nrow,dtype='int32') % 2)

    #header
	f = open(outputfilename+'.header','w')

	header1 = 'UVW;WEIGHT;SIGMA;ANTENNA1;ANTENNA2;ARRAY_ID;DATA_DESC_ID;EXPOSURE;FEED1;FEED2;FIELD_ID;FLAG_ROW;INTERVAL;OBSERVATION_ID;PROCESSOR_ID;SCAN_NUMBER;STATE_ID;TIME;TIME_CENTROID;FLAG_CATEGORY;FLAG;FLOAT_DATA'
	header2 = 'D3;R2;R2;I;I;I;I;D;I;I;I;B;D;I;I;I;I;D;D;B3;B2,'+str(nchan)+';R2,'+str(nchan)

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	f = open(outputfilename+'.dat','w')

    # table
	for i in range(nrow):
		for i_spw in range(nspw):
			UVW = '0;0;0'
			WEIGHT = str(weight)+';'+str(weight)
			SIGMA = str(sigma)+';'+str(sigma)
			ANTENNA1 = '0'
			ANTENNA2 = '0'
			ARRAY_ID = '0'
			DATA_DESC_ID = str(i_spw)
			EXPOSURE = str(texp)
			FEED1 = '0'
			FEED2 = '0'
			FIELD_ID = '0'
			FLAG_ROW = '0'
			INTERVAL = str(texp)
			OBSERVATION_ID = '0'
			PROCESSOR_ID = '0'
			SCAN_NUMBER = str(i)
			STATE_ID = str(state_id[i])
			TIME = str(time[i])
			TIME_CENTROID = str(time[i])
			FLAG_CATEGORY = ''
			FLAG = ''
			FLOAT_DATA = ''

			f.write(UVW+';')
			f.write(WEIGHT+';')
			f.write(SIGMA+';')
			f.write(ANTENNA1+';')
			f.write(ANTENNA2+';')
			f.write(ARRAY_ID+';')
			f.write(DATA_DESC_ID+';')
			f.write(EXPOSURE+';')
			f.write(FEED1+';')
			f.write(FEED2+';')
			f.write(FIELD_ID+';')
			f.write(FLAG_ROW+';')
			f.write(INTERVAL+';')
			f.write(OBSERVATION_ID+';')
			f.write(PROCESSOR_ID+';')
			f.write(SCAN_NUMBER+';')
			f.write(STATE_ID+';')
			f.write(TIME+';')
			f.write(TIME_CENTROID+';')
			f.write(FLAG_CATEGORY)
			f.write(FLAG)
			f.write(FLOAT_DATA+'\n')

	f.close()

    # make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = np.zeros_like(np.concatenate([specdata[:,0],specdata[:,1]],axis=0),dtype='bool').transpose(0,2,1)
	returned_table.putcol('FLAG',value)

	value  = np.zeros_like(np.concatenate([specdata[:,0],specdata[:,1]],axis=0),dtype='float64')

	for i in range(2):
		value[ind_spw==i] = specdata[:,i].copy()

	value = value.transpose(0,2,1)
	returned_table.putcol('FLOAT_DATA',value)

	#value = np.zeros([2*nrow,3],dtype='bool')
	#returned_table.putcol('FLAG_CATEGORY',value)

	value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
	         'MEASINFO': {'type': 'uvw', 'Ref': 'ITRF'},
	        }
	returned_table.putcolkeywords('UVW',value)

	value = {'CATEGORY': np.array([],dtype='|S1')}
	returned_table.putcolkeywords('FLAG_CATEGORY',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2')}
	returned_table.putcolkeywords('EXPOSURE',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2')}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
	         'MEASINFO': {'type': 'epoch', 'Ref': 'UTC'}
	        }
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
	         'MEASINFO': {'type': 'epoch', 'Ref': 'UTC'},
	        }
	returned_table.putcolkeywords('TIME_CENTROID',value)

	value = {'UNIT': 'K'}
	returned_table.putcolkeywords('FLOAT_DATA',value)

	returned_table.flush()
	returned_table.close()

# end makeMAIN
################################################################################
# makeMAIN2

def makeMAIN2(tablename,specdata,time,state_id,texp=0.2,tBW=2.5e9):

	# modules
	import casacore.tables as tb

	# params
	nrow  = specdata.shape[0]
	nspw  = specdata.shape[1]
	npol = specdata.shape[2]
	nchan  = specdata.shape[3]

	weight = tBW/float(nchan) * texp
	sigma = (tBW/float(nchan) * texp)**-0.5

	ind_spw = (np.linspace(0,2*nrow-1,2*nrow,dtype='int32') % 2)

	# tables
	colnames = ['UVW',
				'FLAG',
				'FLAG_CATEGORY',
				'WEIGHT',
				'SIGMA',
				'ANTENNA1',
				'ANTENNA2',
				'ARRAY_ID',
				'DATA_DESC_ID',
				'EXPOSURE',
				'FEED1',
				'FEED2',
				'FIELD_ID',
				'FLAG_ROW',
				'INTERVAL',
				'OBSERVATION_ID',
				'PROCESSOR_ID',
				'SCAN_NUMBER',
				'STATE_ID',
				'TIME',
				'TIME_CENTROID',
				'FLOAT_DATA'
				]

	colkeywords = [
		{'MEASINFO': {'Ref': 'ITRF', 'type': 'uvw'},'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2')},
		{},
		{'CATEGORY': np.array([],dtype='|S1')},
		{},{},{},{},{},{},
		{'QuantumUnits': np.array(['s'],dtype='|S2')},
		{},{},{},{},
		{'QuantumUnits': np.array(['s'],dtype='|S2')},
		{},{},{},{},
		{'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
		{'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
		{'UNIT': 'K'}
				  ]

	ndims = [1,2,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2]
	isarrays = [True,True,True,True,True,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,True]
	valuetypes = ['double','bool','bool','float','float','int','int','int','int','double','int','int','int','bool','double','int','int','int','int','double','double','float']

	descs = []
	for colname,colkeyword,ndim,valuetype,isarray in zip(colnames,colkeywords,ndims,valuetypes,isarrays):
		if colname=='UVW':
			descs.append(tb.makearrcoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',ndim=ndim,keywords=colkeyword,valuetype=valuetype,options=5,shape=np.array([3], dtype='int32')))
		elif isarray:
			descs.append(tb.makearrcoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',ndim=ndim,keywords=colkeyword,valuetype=valuetype))
		else:
			descs.append(tb.makescacoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',keywords=colkeyword,valuetype=valuetype))

	td = tb.maketabdesc(descs=descs)

	returned_table = tb.table(tablename,tabledesc=td,nrow=2*nrow,readonly=False)

	# put values
	value = np.zeros([2*nrow,3],dtype='float64')
	returned_table.putcol('UVW',value)

	value = np.full([2*nrow,2],sigma)
	returned_table.putcol('SIGMA',value)

	value = np.full([2*nrow,2],weight)
	returned_table.putcol('WEIGHT',value)

	value = np.zeros([2*nrow],dtype='int32')
	returned_table.putcol('ANTENNA1',value)
	returned_table.putcol('ANTENNA2',value)
	returned_table.putcol('ARRAY_ID',value)
	returned_table.putcol('FEED1',value)
	returned_table.putcol('FEED2',value)
	returned_table.putcol('FIELD_ID',value)
	returned_table.putcol('OBSERVATION_ID',value)
	returned_table.putcol('PROCESSOR_ID',value)

	value = np.zeros([2*nrow],dtype='bool')
	returned_table.putcol('FLAG_ROW',value)

	value = np.full([2*nrow],texp,dtype='float64')
	returned_table.putcol('EXPOSURE',value)
	returned_table.putcol('INTERVAL',value)

	value = np.zeros_like(np.concatenate([specdata[:,0],specdata[:,1]],axis=0),dtype='bool')
	value = value.transpose(0,2,1)
	returned_table.putcol('FLAG',value)

	value  = np.zeros_like(np.concatenate([specdata[:,0],specdata[:,1]],axis=0),dtype='float64')

	for i in range(2):
		value[ind_spw==i] = specdata[:,i].copy()

	value = value.transpose(0,2,1)
	returned_table.putcol('FLOAT_DATA',value)

	value = np.zeros(2*nrow,dtype='int32')
	for i in range(2):
		value[ind_spw==i] = state_id

	value = np.zeros(2*nrow,dtype='int32')
	for i in range(2):
		value[ind_spw==i] = i

	returned_table.putcol('DATA_DESC_ID',value)

	value = np.zeros(2*nrow,dtype='int32')
	for i in range(2):
		value[ind_spw==i] = state_id
	returned_table.putcol('STATE_ID',value)

	value = np.zeros(2*nrow,dtype='int32')
	for i in range(2):
		value[ind_spw==i] = np.linspace(0,nrow-1,nrow,dtype='int32')

	returned_table.putcol('SCAN_NUMBER',value)

	value = np.zeros(2*nrow,dtype='float64')
	for i in range(2):
		value[ind_spw==i] = time.copy()
	returned_table.putcol('TIME',value)
	returned_table.putcol('TIME_CENTROID',value)

	returned_table.flush()
	returned_table.close()

# end makeMAIN2
################################################################################
# makeANTENNA

def makeANTENNA(tablename,
                outputfilename,
                dish_diameter = 50.,
                nbeam = 1,
                beamname = 'LMT-b4r-beam',
                lon = -97. + 18./60. + 53./3600.,
                lat = 18. + 59./60. + 6./3600.,
                height = 4600.,
                type = 'GROUND-BASED',
                mounit='ALT-AZ'):

	# modules
	from astropy import coordinates
	import casacore.tables as tb

	# params
	c = coordinates.EarthLocation.from_geodetic(lon=lon,lat=lat,height=height)
	x = c.value[0]
	y = c.value[1]
	z = c.value[2]
	beamname_list = []
	for i in range(nbeam):
		beamname_list.append(beamname+str(i))

	# header
	header1 = 'OFFSET;POSITION;TYPE;DISH_DIAMETER;FLAG_ROW;MOUNT;NAME;STATION'
	header2 = 'D3;D3;A;D;B;A;A;A'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nbeam):
	    OFFSET = '0;0;0'
	    POSITION = str(x)+';'+str(y)+';'+str(z)
	    TYPE = '"'+type+'"'
	    DISH_DIAMETER = str(dish_diameter)
	    FLAG_ROW = '0'
	    MOUNT = '"'+mounit+'"'
	    NAME = '"'+';'.join(beamname_list)+'"'
	    STATION = '""'

	    f.write(OFFSET+';')
	    f.write(POSITION+';')
	    f.write(TYPE+';')
	    f.write(DISH_DIAMETER+';')
	    f.write(FLAG_ROW+';')
	    f.write(MOUNT+';')
	    f.write(NAME+';')
	    f.write(STATION+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
	         'MEASINFO': {'type': 'position','Ref': 'ITRF'},
	        }
	returned_table.putcolkeywords('OFFSET',value)

	value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
	         'MEASINFO': {'type': 'position','Ref': 'ITRF'},
	        }
	returned_table.putcolkeywords('POSITION',value)

	value = {'QuantumUnits': np.array(['m'],dtype='|S2')}
	returned_table.putcolkeywords('DISH_DIAMETER',value)

# end makeANTENNA
################################################################################
# makeANTENNA

def makeDATA_DESCRIPTION(tablename,outputfilename,nspw=2):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'FLAG_ROW;POLARIZATION_ID;SPECTRAL_WINDOW_ID'
    header2 = 'B;I;I'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    for i in range(nspw):
        FLAG_ROW = '0'
        POLARIZATION_ID = '0'
        SPECTRAL_WINDOW_ID = str(i)
        f.write(FLAG_ROW+';')
        f.write(POLARIZATION_ID+';')
        f.write(SPECTRAL_WINDOW_ID+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makeDATA_DESCRIPTION
################################################################################
# makeDOPPLER

def makeDOPPLER(tablename,outputfilename):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'DOPPLER_ID;SOURCE_ID;TRANSITION_ID;VELDEF'
    header2 = 'I;I;I;D'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    f.write('')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['m/s'],dtype='|S4'),
             'MEASINFO': {'Ref': 'RADIO', 'type': 'doppler'},
            }
    returned_table.putcolkeywords('VELDEF',value)

# end makeDOPPLER
################################################################################
# makeFEED

def makeFEED(tablename,outputfilename,nspw=2):

    #modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'POSITION;BEAM_OFFSET;POLARIZATION_TYPE;POL_RESPONSE;RECEPTOR_ANGLE;ANTENNA_ID;BEAM_ID;FEED_ID;INTERVAL;NUM_RECEPTORS;SPECTRAL_WINDOW_ID;TIME'
    header2 = 'D3;D2,2;A2;X2,2;D2;I;I;I;D;I;I;D'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    for i in range(nspw):
        POSITION = '0;0;0'
        BEAM_OFFSET = '0;0;0;0'
        POLARIZATION_TYPE = '"X";"Y"'
        POL_RESPONSE = '0;0;0;0;0;0;0;0'
        RECEPTOR_ANGLE = '0;0'
        ANTENNA_ID = '0'
        BEAM_ID = '0'
        FEED_ID = '0'
        INTERVAL = '0'
        NUM_RECEPTORS = '0'
        SPECTRAL_WINDOW_ID = str(i)
        TIME = '0'

        f.write(POSITION+';')
        f.write(BEAM_OFFSET+';')
        f.write(POLARIZATION_TYPE+';')
        f.write(POL_RESPONSE+';')
        f.write(RECEPTOR_ANGLE+';')
        f.write(ANTENNA_ID+';')
        f.write(BEAM_ID+';')
        f.write(FEED_ID+';')
        f.write(INTERVAL+';')
        f.write(NUM_RECEPTORS+';')
        f.write(SPECTRAL_WINDOW_ID+';')
        f.write(TIME+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['m', 'm', 'm'],dtype='|S2'),
             'MEASINFO': {'type': 'position','Ref': 'ITRF'},
            }
    returned_table.putcolkeywords('POSITION',value)

    value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
             'MEASINFO': {'type': 'direction','Ref': 'J2000'},
            }
    returned_table.putcolkeywords('BEAM_OFFSET',value)

    value = {'QuantumUnits': np.array(['rad'],dtype='|S4'),}
    returned_table.putcolkeywords('RECEPTOR_ANGLE',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
    returned_table.putcolkeywords('INTERVAL',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME',value)

# end makeFEED
################################################################################
# makeFIELD

def makeFIELD(tablename,outputfilename,ra,dec,fieldname,time):

    #modules
    import casacore.tables as tb

    # params
    ra_rad  = ra/180.*np.pi
    dec_rad = dec/180.*np.pi
    time_mean = np.mean(time)

    # header
    header1 = 'DELAY_DIR;PHASE_DIR;REFERENCE_DIR;CODE;FLAG_ROW;NAME;NUM_POLY;SOURCE_ID;TIME'
    header2 = 'D2,1;D2,1;D2,1;A;B;A;I;I;D'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    DELAY_DIR = str(ra_rad)+';'+str(dec_rad)
    PHASE_DIR = str(ra_rad)+';'+str(dec_rad)
    REFERENCE_DIR = str(ra_rad)+';'+str(dec_rad)
    CODE = '""'
    FLAG_ROW = '0'
    NAME = '"'+fieldname+'"'
    NUM_POLY = '0'
    SOURCE_ID = '0'
    TIME = str(time_mean)

    f.write(DELAY_DIR+';')
    f.write(PHASE_DIR+';')
    f.write(REFERENCE_DIR+';')
    f.write(CODE+';')
    f.write(FLAG_ROW+';')
    f.write(NAME+';')
    f.write(NUM_POLY+';')
    f.write(SOURCE_ID+';')
    f.write(TIME+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
             'MEASINFO': {'type': 'direction','Ref': 'J2000'},
            }
    returned_table.putcolkeywords('DELAY_DIR',value)

    value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
             'MEASINFO': {'type': 'direction','Ref': 'J2000'},
            }
    returned_table.putcolkeywords('PHASE_DIR',value)

    value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
             'MEASINFO': {'type': 'direction','Ref': 'J2000'},
            }
    returned_table.putcolkeywords('REFERENCE_DIR',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME',value)

# end makeFIELD
################################################################################
# makeFLAG_CMD

def makeFLAG_CMD(tablename,outputfilename):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'APPLIED;COMMAND;INTERVAL;LEVEL;REASON;SEVERITY;TIME;TYPE'
    header2 = 'B;A;D;I;A;I;D;A'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    f.write('')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
    returned_table.putcolkeywords('INTERVAL',value)

# end makeFLAG_CMD
################################################################################
# makeFREQ_OFFSET

def makeFREQ_OFFSET(tablename,outputfilename):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'ANTENNA1;ANTENNA2;FEED_ID;INTERVAL;OFFSET;SPECTRAL_WINDOW_ID;TIME'
    header2 = 'I;I;I;D;D;I;D'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    f.write('')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
    returned_table.putcolkeywords('INTERVAL',value)

    value = {'QuantumUnits': np.array(['Hz'],dtype='|S3'),
             'MEASINFO': {'type': 'frequency','Ref': 'LSRK'},
            }
    returned_table.putcolkeywords('OFFSET',value)

# end makeFREQ_OFFSET
################################################################################
# makeHISTORY

def makeHISTORY(tablename,outputfilename):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'APP_PARAMS;CLI_COMMAND;APPLICATION;MESSAGE;OBJECT_ID;OBSERVATION_ID;ORIGIN;PRIORITY;TIME'
    header2 = 'A1;A1;A;A;I;I;A;A;D'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    f.write('')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME',value)

# end makeHISTORY
################################################################################
# makeLMT_ARRAY

def makeLMT_ARRAY(tablename,outputfilename,nbeam=1,npol=2,nspw=2):

    # modules
    import casacore.tables as tb

    # params
    n_array = nbeam*npol*nspw

    # header
    header1 = 'ARRAY;BEAM;POLARIZATION;SPECTRAL_WINDOW'
    header2 = 'I;I;I;I'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    for i in range(n_array):
        ARRAY = str(i)
        BEAM = str((i // npol) % nbeam)
        POLARIZATION = str((i % npol) * 3 + 9)
        SPECTRAL_WINDOW = str((i // (npol*nbeam)) % nspw)

        f.write(ARRAY+';')
        f.write(BEAM+';')
        f.write(POLARIZATION+';')
        f.write(SPECTRAL_WINDOW+'\n')

    f.close()
    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makeLMT_ARRAY
################################################################################
# makeOBSERVATION

def makeOBSERVATION(tablename,outputfilename,time,project,observer,telescope_name='NRO'):

    # modules
    import casacore.tables as tb

    # params
    time_start = time.min()
    time_end   = time.max()

    # header
    header1 = 'TIME_RANGE;LOG;SCHEDULE;FLAG_ROW;OBSERVER;PROJECT;RELEASE_DATE;SCHEDULE_TYPE;TELESCOPE_NAME'
    header2 = 'D2;A1;A1;B;A;A;D;A;A'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    TIME_RANGE = str(time_start)+';'+str(time_end)
    LOG = ''
    SCHEDULE = ''
    FLAG_ROW = '0'
    OBSERVER = '"'+observer+'"'
    PROJECT = '"'+project+'"'
    RELEASE_DATE = '0'
    SCHEDULE_TYPE = ''
    TELESCOPE_NAME = '"'+telescope_name+'"'

    f.write(TIME_RANGE+';')
    f.write(LOG+';')
    f.write(SCHEDULE+';')
    f.write(FLAG_ROW+';')
    f.write(OBSERVER+';')
    f.write(PROJECT+';')
    f.write(RELEASE_DATE+';')
    f.write(SCHEDULE_TYPE+';')
    f.write(TELESCOPE_NAME+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME_RANGE',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('RELEASE_DATE',value)

# end makeOBSERVATION
################################################################################
# makePOINTING

def makePOINTING(tablename,outputfilename,direction,time,nbeam=1,texp=0.2):

	'''
	direction: (nscan*nbeam) array
	'''

	# modules
	import casacore.tables as tb

	# params
	direction_rad = direction/180.*np.pi
	n_direction = direction.shape[0]
	time_end   = time.max()

	# header
	header1 = 'DIRECTION;ANTENNA_ID;INTERVAL;NAME;NUM_POLY;TARGET;TIME;TIME_ORIGIN;TRACKING'
	header2 = 'D2,1;I;D;A;I;D2,1;D;D;B'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(n_direction):
		#DIRECTION = str(direction_rad[i][0])+';'+str(direction_rad[i][1])
		DIRECTION = ';'
		ANTENNA_ID = str((n_direction // (n_direction/nbeam))-1)
		INTERVAL = str(texp)
		NAME = ''
		NUM_POLY= '0'
		TARGET = ';'
		TIME = time[i].astype('|U')
		TIME_ORIGIN = '0'
		TRACKING = '0'

		f.write(DIRECTION+';')
		f.write(ANTENNA_ID+';')
		f.write(INTERVAL+';')
		f.write(NAME+';')
		f.write(NUM_POLY+';')
		f.write(TARGET+';')
		f.write(TIME+';')
		f.write(TIME_ORIGIN+';')
		f.write(TRACKING+'\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value  = np.zeros([n_direction,1,2],dtype='float64')

	value[:,0] = direction_rad.copy()

	returned_table.putcol('DIRECTION',value)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
	         'MEASINFO': {'type': 'direction','Ref': 'J2000'},
	        }
	returned_table.putcolkeywords('DIRECTION',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
	         'MEASINFO': {'type': 'direction','Ref': 'J2000'},
	        }
	returned_table.putcolkeywords('TARGET',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
	         'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
	        }
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
	         'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
	        }
	returned_table.putcolkeywords('TIME_ORIGIN',value)

	returned_table.flush()
	returned_table.close()

# end makePOINTING
################################################################################
# makePOINTING2

def makePOINTING2(tablename,direction,time,nbeam=1,texp=0.2):

	# modules
	import casacore.tables as tb

	# params
	direction_rad = direction/180.*np.pi
	n_direction = direction.shape[0]
	time_end   = time.max()

	# make table
	colnames = ['DIRECTION',
				'ANTENNA_ID',
				'INTERVAL',
				'NAME',
				'NUM_POLY',
				'TARGET',
				'TIME',
				'TIME_ORIGIN',
				'TRACKING']
	colkeywords = [{'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4')},
				   {},
				   {'QuantumUnits': np.array(['s'],dtype='|S2')},
				   {},{},
				   {'MEASINFO': {'Ref': 'J2000', 'type': 'direction'},'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4')},
				   {'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
				   {'MEASINFO': {'Ref': 'UTC', 'type': 'epoch'}, 'QuantumUnits': np.array(['s'],dtype='|S2')},
				   {}
				  ]

	#ndims = [2,1,1,1,1,2,1,1,1]
	#isarrays = [True,False,False,False,False,True,False,False,False]
	valuetypes = ['double','int','double','string','int','double','double','double','bool']
	ndims = [2,1,1,1,1,-1,1,1,1]
	descs = []

	for colname,colkeyword,valuetype,ndim in zip(colnames,colkeywords,valuetypes,ndims):
		if (colname=='DIRECTION' or colname=='TARGET'):
			descs.append(tb.makearrcoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',ndim=ndim,keywords=colkeyword,valuetype=valuetype,options=0))
		else:
			descs.append(tb.makescacoldesc(colname,0.0,datamanagertype='StandardStMan',datamanagergroup='StandardStMan',keywords=colkeyword,valuetype=valuetype))

	td = tb.maketabdesc(descs=descs)

	returned_table = tb.table(tablename,tabledesc=td,nrow=n_direction,readonly=False)

	value  = np.zeros([n_direction,1,2],dtype='float64')
	value[:,0] = direction_rad.copy()
	returned_table.putcol('DIRECTION',value)

	value = np.zeros([n_direction],dtype='int32')
	returned_table.putcol('ANTENNA_ID',value)
	returned_table.putcol('NUM_POLY',value)

	value = np.zeros([n_direction],dtype='float64')
	returned_table.putcol('TIME_ORIGIN',value)

	value = np.zeros([n_direction],dtype='bool')
	returned_table.putcol('TRACKING',value)

	value = np.full([n_direction],texp,dtype='float64')
	returned_table.putcol('INTERVAL',value)

	value = time
	returned_table.putcol('TIME',value)

	returned_table.flush()
	returned_table.close()

# end makePOINTING2
################################################################################
# makePOINTING

def makePOLARIZAION(tablename,outputfilename,npol=2):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'CORR_TYPE;CORR_PRODUCT;FLAG_ROW;NUM_CORR'
    header2 = 'I'+str(npol)+';I'+str(npol)+','+str(npol)+';B;I'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    CORR_TYPE = ''
    for i in range(npol):
        CORR_TYPE = CORR_TYPE + str(3*i+9) + ';'
    CORR_TYPE = CORR_TYPE[0:-1]
    CORR_PRODUCT = '0;'*(npol**2-1)+'0'
    FLAG_ROW = '0'
    NUM_CORR = str(npol)

    f.write(CORR_TYPE+';')
    f.write(CORR_PRODUCT+';')
    f.write(FLAG_ROW+';')
    f.write(NUM_CORR+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makePOLARIZAION
################################################################################
# makePROCESSOR

def makePROCESSOR(tablename,outputfilename):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'FLAG_ROW;MODE_ID;TYPE;TYPE_ID;SUB_TYPE'
    header2 = 'B;I;A;I;A'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    f.write('0;-1;"";-1;""')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

# end makePROCESSOR
################################################################################
# makeSOURCE

def makeSOURCE(tablename,outputfilename,ra,dec,sysvel,sourcename,time,freq):

    # modules
    import casacore.tables as tb

    # params
    ra_rad  = ra/180.*np.pi
    dec_rad = dec/180.*np.pi
    time_interval = time.max() - time.min()
    time_mean = np.mean(time)
    nspw = freq.shape[0]

    # header
    header1 = 'DIRECTION;PROPER_MOTION;CALIBRATION_GROUP;CODE;INTERVAL;NAME;NUM_LINES;SOURCE_ID;SPECTRAL_WINDOW_ID;TIME;TRANSITION;REST_FREQUENCY;SYSVEL'
    header2 = 'D2;D2;I;A;D;A;I;I;I;D;A1;D1;D1'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    for i in range(nspw):
        DIRECTION = str(ra_rad)+';'+str(dec_rad)
        PROPER_MOTION = '0;0'
        CALIBRATION_GROUP ='-1'
        CODE = '""'
        INTERVAL = str(time_interval)
        NAME = '"'+sourcename+'"'
        NUM_LINES = '1'
        SOURCE_ID = '0'
        SPECTRAL_WINDOW_ID = str(i)
        TIME = str(time_mean)
        TRANSITION = ''
        REST_FREQUENCY = str(np.mean(freq[i]))
        SYSVEL = str(sysvel)

        f.write(DIRECTION+';')
        f.write(PROPER_MOTION+';')
        f.write(CALIBRATION_GROUP+';')
        f.write(CODE+';')
        f.write(INTERVAL+';')
        f.write(NAME+';')
        f.write(NUM_LINES+';')
        f.write(SOURCE_ID+';')
        f.write(SPECTRAL_WINDOW_ID+';')
        f.write(TIME+';')
        f.write(TRANSITION+';')
        f.write(REST_FREQUENCY+';')
        f.write(SYSVEL+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['rad', 'rad'],dtype='|S4'),
             'MEASINFO': {'type': 'direction','Ref': 'J2000'},
            }
    returned_table.putcolkeywords('DIRECTION',value)

    value = {'QuantumUnits': np.array(['rad/s'],dtype='|S6')}
    returned_table.putcolkeywords('PROPER_MOTION',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2')}
    returned_table.putcolkeywords('INTERVAL',value)

    value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
             'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
            }
    returned_table.putcolkeywords('TIME',value)


    value = {'QuantumUnits': np.array(['Hz'],dtype='|S3'),
             'MEASINFO': {'type': 'frequency','Ref': 'LSRK'},
            }
    returned_table.putcolkeywords('REST_FREQUENCY',value)

    value = {'QuantumUnits': np.array(['m/s'],dtype='|S4'),
             'MEASINFO': {'Ref': 'LSRK', 'type': 'radialvelocity'},
            }
    returned_table.putcolkeywords('SYSVEL',value)


# end makeSOURCE
################################################################################
# makeSPECTRAL_WINDOW

def makeSPECTRAL_WINDOW(tablename,outputfilename,freq,nspw=2):

	# modules
	import casacore.tables as tb

	# params
	nchan = freq.shape[1]
	tBW = abs(freq[0].max() - freq[0].min())

	# header
	header1 = 'MEAS_FREQ_REF;REF_FREQUENCY;FLAG_ROW;FREQ_GROUP;FREQ_GROUP_NAME;IF_CONV_CHAIN;NAME;NET_SIDEBAND;NUM_CHAN;TOTAL_BANDWIDTH;CHAN_FREQ;CHAN_WIDTH;EFFECTIVE_BW;RESOLUTION'
	header2 = 'I;D;B;I;A;I;A;I;I;D;D'+str(nchan)+';D'+str(nchan)+';D'+str(nchan)+';D'+str(nchan)

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(nspw):
		MEAS_FREQ_REF = '1'
		#CHAN_FREQ = ';'.join(map(str,freq[i]))
		#CHAN_FREQ = ''
		REF_FREQUENCY = str(freq[i][0])
		#CHAN_WIDTH = (str(np.diff(freq[i])[0])+';')*(nchan-1)+str(np.diff(freq[i])[0])
		#CHAN_WIDTH = ''
		#EFFECTIVE_BW = (str(np.abs(np.diff(freq[i]))[0])+';')*(nchan-1)+str(np.abs(np.diff(freq[i]))[0])
		#EFFECTIVE_BW = ''
		#RESOLUTION = (str(np.abs(np.diff(freq[i]))[0])+';')*(nchan-1)+str(np.abs(np.diff(freq[i]))[0])
		#RESOLUTION = ''
		FLAG_ROW = '0'
		FREQ_GROUP = '0'
		FREQ_GROUP_NAME = ''
		IF_CONV_CHAIN = '0'
		NAME = ''
		NET_SIDEBAND = str(i)
		NUM_CHAN = str(nchan)
		TOTAL_BANDWIDTH = str(tBW)

		f.write(MEAS_FREQ_REF+';')
		#f.write(CHAN_FREQ+';')
		f.write(REF_FREQUENCY+';')
		#f.write(CHAN_WIDTH+';')
		#f.write(EFFECTIVE_BW+';')
		#f.write(RESOLUTION+';')
		f.write(FLAG_ROW+';')
		f.write(FREQ_GROUP+';')
		f.write(FREQ_GROUP_NAME+';')
		f.write(IF_CONV_CHAIN+';')
		f.write(NAME+';')
		f.write(NET_SIDEBAND+';')
		f.write(NUM_CHAN+';')
		f.write(TOTAL_BANDWIDTH+';\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = freq.copy()
	returned_table.putcol('CHAN_FREQ',value)

	value = np.zeros_like(freq)
	for i in range(2):
		value[i] = np.full(nchan,np.diff(freq[i])[0])

	returned_table.putcol('CHAN_WIDTH',value)

	value = np.zeros_like(freq)
	for i in range(2):
		value[i] = np.full(nchan,abs(np.diff(freq[i])[0]))

	returned_table.putcol('EFFECTIVE_BW',value)
	returned_table.putcol('RESOLUTION',value)

	value = {'MEASINFO': {'TabRefCodes': np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8, 64], dtype='uint32'),
	                      'TabRefTypes': np.array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP','CMB', 'Undefined'],dtype='|S10'),
	                      'VarRefCol': 'MEAS_FREQ_REF',
	                      'type': 'frequency'},
	         'QuantumUnits': np.array(['Hz'],dtype='|S3')}
	returned_table.putcolkeywords('CHAN_FREQ',value)

	value = {'MEASINFO': {'TabRefCodes': np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8, 64], dtype='uint32'),
	                      'TabRefTypes': np.array(['REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP','CMB', 'Undefined'],dtype='|S10'),
	                      'VarRefCol': 'MEAS_FREQ_REF',
	                      'type': 'frequency'},
	         'QuantumUnits': np.array(['Hz'],dtype='|S3')}
	returned_table.putcolkeywords('REF_FREQUENCY',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('CHAN_WIDTH',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('EFFECTIVE_BW',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('RESOLUTION',value)

	value = {'QuantumUnits': np.array(['Hz'],dtype='S3')}
	returned_table.putcolkeywords('TOTAL_BANDWIDTH',value)

	returned_table.flush()
	returned_table.close()

# end makeSPECTRAL_WINDOW
################################################################################
# makeSTATE

def makeSTATE(tablename,outputfilename):

    # modules
    import casacore.tables as tb

    # params

    # header
    header1 = 'CAL;FLAG_ROW;LOAD;OBS_MODE;REF;SIG;SUB_SCAN'
    header2 = 'D;B;D;A;B;B;I'

    f = open(outputfilename+'.header','w')

    f.write(header1+'\n')
    f.write(header2+'\n')
    f.close()

    # table
    f = open(outputfilename+'.dat','w')

    f.write('0;0;0;"OBSERVE_TARGET#OFF_SOURCE,POSITION_SWITCH";1;0;2'+'\n')
    f.write('0;0;0;"OBSERVE_TARGET#ON_SOURCE,POSITION_SWITCH";1;0;1'+'\n')

    f.close()

    # make table
    returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

    value = {'QuantumUnits': np.array(['K'],dtype='S2')}
    returned_table.putcolkeywords('CAL',value)

    value = {'QuantumUnits': np.array(['K'],dtype='S2')}
    returned_table.putcolkeywords('LOAD',value)

# end makeSTATE
################################################################################
# makeSYSCAL

def makeSYSCAL(tablename,outputfilename,Tsys,Tsys_time,time,nbeam=1,npol=2,nspw=2):

	'''
	Tsys: (nspw,npol,nchan) array
	'''

	# modules
	import casacore.tables as tb

	# params
	nspw = Tsys.shape[0]
	npol = Tsys.shape[1]
	nchan = Tsys.shape[2]
	n_array = nbeam*nspw
	time_interval = time.max() - Tsys_time

	# header
	header1 = 'ANTENNA_ID;FEED_ID;INTERVAL;SPECTRAL_WINDOW_ID;TIME;TSYS' #;TSYS_SPECTRUM'
	header2 = 'I;I;D;I;D;R'+str(npol) #+';R'+str(npol)+','+str(nchan)

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	for i in range(n_array):
		ANTENNA_ID = str(i % nbeam)
		FEED_ID = '0'
		INTERVAL = str(time_interval)
		SPECTRAL_WINDOW_ID = str(i // nbeam)
		TIME = str(Tsys_time)
		#TSYS_SPECTRUM = ';'.join(map(str,Tsys[i // nbeam].T.reshape(npol*nchan)))
		TSYS_SPECTRUM = ''
		TSYS = ';'.join(map(str,np.median(Tsys[(i // nbeam)],axis=(1,))))

		f.write(ANTENNA_ID+';')
		f.write(FEED_ID+';')
		f.write(INTERVAL+';')
		f.write(SPECTRAL_WINDOW_ID+';')
		f.write(TIME+';')
		f.write(TSYS+'\n')
		#f.write(TSYS_SPECTRUM+';\n')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	#value = Tsys.copy().transpose(0,2,1)
	#returned_table.putcol('TSYS_SPECTRUM',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
	         'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
	        }
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	#value = {'QuantumUnits': np.array(['K'],dtype='S2')}
	#returned_table.putcolkeywords('TSYS_SPECTRUM',value)

	value = {'QuantumUnits': np.array(['K'],dtype='S2')}
	returned_table.putcolkeywords('TSYS',value)

# end makeSYSCAL
################################################################################
# makeWEATHER

def makeWEATHER(tablename,outputfilename):

	# modules
	import casacore.tables as tb

	# params

	# header
	header1 = 'ANTENNA_ID;INTERVAL;TIME;TEMPERATURE;PRESSURE;REL_HUMIDITY;WIND_SPEED;WIND_DIRECTION'
	header2 = 'I;D;D;R;R;R;R;R'

	f = open(outputfilename+'.header','w')

	f.write(header1+'\n')
	f.write(header2+'\n')
	f.close()

	# table
	f = open(outputfilename+'.dat','w')

	f.write('')

	f.close()

	# make table
	returned_table = tb.tablefromascii(tablename,outputfilename+'.dat',headerfile=outputfilename+'.header',sep=';',readonly=False)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),
	         'MEASINFO': {'type': 'epoch','Ref': 'UTC'},
	        }
	returned_table.putcolkeywords('TIME',value)

	value = {'QuantumUnits': np.array(['s'],dtype='|S2'),}
	returned_table.putcolkeywords('INTERVAL',value)

	value = {'QuantumUnits': np.array(['K'],dtype='|S2'),}
	returned_table.putcolkeywords('TEMPERATURE',value)

	value = {'QuantumUnits': np.array(['hPa'],dtype='|S4'),}
	returned_table.putcolkeywords('PRESSURE',value)

	value = {'QuantumUnits': np.array(['%'],dtype='|S2'),}
	returned_table.putcolkeywords('REL_HUMIDITY',value)

	value = {'QuantumUnits': np.array(['m/s'],dtype='|S4'),}
	returned_table.putcolkeywords('WIND_SPEED',value)

	value = {'QuantumUnits': np.array(['rad'],dtype='|S4'),}
	returned_table.putcolkeywords('WIND_DIRECTION',value)

	returned_table.flush()
	returned_table.close()

# end makeWEATHER
################################################################################
# createMS2

def createMS2(MS2name,specdata,time,ra,dec,sysvel,sourcename,project,observer,direction,freq,Tsys,Tsys_time,state_id,path_temp='./temp/',removetemp=True):

	# modules
	import casacore.tables as tb

	# params
	os.system('rm -rf '+path_temp)
	os.system('rm -rf '+MS2name)
	os.system('mkdir -p '+path_temp)

	# make tables
	#makeMAIN(MS2name,path_temp+'MAIN',specdata,time,state_id)
	makeMAIN2(MS2name,specdata,time,state_id)
	makeANTENNA(MS2name+'/ANTENNA',path_temp+'ANTENNA')
	makeDATA_DESCRIPTION(MS2name+'/DATA_DESCRIPTION',path_temp+'DATA_DESCRIPTION')
	makeDOPPLER(MS2name+'/DOPPLER',path_temp+'DOPPLER')
	makeFEED(MS2name+'/FEED',path_temp+'FEED')
	makeFIELD(MS2name+'/FIELD',path_temp+'FIELD',ra,dec,sourcename,time)
	makeFLAG_CMD(MS2name+'/FLAG_CMD',path_temp+'FLAG_CMD')
	makeFREQ_OFFSET(MS2name+'/FREQ_OFFSET',path_temp+'FREQ_OFFSET')
	makeHISTORY(MS2name+'/HISTORY',path_temp+'HISTORY')
	makeOBSERVATION(MS2name+'/OBSERVATION',path_temp+'OBSERVATION',time,project,observer)
	#makePOINTING(MS2name+'/POINTING',path_temp+'POINTING',direction,time)
	makePOINTING2(MS2name+'/POINTING',direction,time)
	makePOLARIZAION(MS2name+'/POLARIZATION',path_temp+'POLARIZATION')
	makePROCESSOR(MS2name+'/PROCESSOR',path_temp+'PROCESSOR')
	makeSOURCE(MS2name+'/SOURCE',path_temp+'SOURCE',ra,dec,sysvel,sourcename,time,freq)
	makeSPECTRAL_WINDOW(MS2name+'/SPECTRAL_WINDOW',path_temp+'SPECTRAL_WINDOW',freq)
	makeSTATE(MS2name+'/STATE',path_temp+'STATE')
	makeSYSCAL(MS2name+'/SYSCAL',path_temp+'SYSCAL',Tsys,Tsys_time,time)
	makeWEATHER(MS2name+'/WEATHER',path_temp+'WEATHER')
	makeLMT_ARRAY(MS2name+'/LMT_ARRAY',path_temp+'LMT_ARRAY')

	# put keywords
	abs_path = os.path.abspath(MS2name)

	keywords = {'MS_VERSION': 2.0,
				'ANTENNA': 'Table: '+abs_path+'/ANTENNA',
				'DATA_DESCRIPTION': 'Table: '+abs_path+'/DATA_DESCRIPTION',
				'DOPPLER': 'Table: '+abs_path+'/DOPPLER',
				'FEED': 'Table: '+abs_path+'/FEED',
				'FIELD': 'Table: '+abs_path+'/FIELD',
				'FLAG_CMD': 'Table: '+abs_path+'/FLAG_CMD',
				'FREQ_OFFSET': 'Table: '+abs_path+'/FREQ_OFFSET',
				'HISTORY': 'Table: '+abs_path+'/HISTORY',
				'OBSERVATION': 'Table: '+abs_path+'/OBSERVATION',
				'POINTING': 'Table: '+abs_path+'/POINTING',
				'POLARIZATION': 'Table: '+abs_path+'/POLARIZATION',
				'PROCESSOR': 'Table: '+abs_path+'/PROCESSOR',
				'SOURCE': 'Table: '+abs_path+'/SOURCE',
				'SPECTRAL_WINDOW': 'Table: '+abs_path+'/SPECTRAL_WINDOW',
				'STATE': 'Table: '+abs_path+'/STATE',
				'SYSCAL': 'Table: '+abs_path+'/SYSCAL',
				'WEATHER': 'Table: '+abs_path+'/WEATHER',
				'LMT_ARRAY': 'Table: '+abs_path+'/LMT_ARRAY',
				}
	returnedTable = tb.table(MS2name,readonly=False)
	returnedTable.putkeywords(keywords)

	returnedTable.flush(recursive=True)
	returnedTable.close()

	if removetemp:
	    os.system('rm -rf '+path_temp)

	#return returnedTable

# end createMS2
################################################################################





###
