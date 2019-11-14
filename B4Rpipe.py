#!/usr/bin/env python3

# modules
import Lib_b4r as Lib
from optparse import OptionParser

import glob
import os

# option
parser = OptionParser()
parser.add_option("-s", "--sci", dest="sci", help="write Science FILENAME", metavar="SCIFILE")
parser.add_option("-c", "--cal", dest="cal", help="write Clibration FILENAME", metavar="CALFILE")
parser.add_option("-n", "--nchbin", dest="chbin", help="numer of channel for binning", metavar="CALFILE")
parser.add_option("-a", "--antlog", dest="antlog", help="write Antenna log FILENAME", metavar="CALFILE")
parser.add_option("-o", "--output", dest="output", help="write output MS2 FILENAME", metavar="CALFILE")
parser.add_option("-t", "--obstype", dest="obstype", help="write obstype (cont or spec)", metavar="CALFILE")

(options, args) = parser.parse_args()

path_cal_raw_head = options.cal
path_sci_raw_head = options.sci
chbin = int(options.chbin)
path_antlog = options.antlog
MS2name = options.output

if options.obstype == 'cont':
    blsub=True
elif options.obstype == 'spec':
    blsub=False
else:
    raise ValueError("[B4Rpipe error] Please select cont or spec in obstype.")

# main
print('################')
print('# step 1/4     #')
print('# load data    #')
print('################')
specdata, Tsys, freq_, state_id, Tsys_time = Lib.loadB4Rfulldata(path_cal_raw_head,path_sci_raw_head,path_antlog,chbin,cal=True,blsub=blsub)
print('')
print('################')
print('# step 2/4     #')
print('# get obs info #')
print('################')
ra,dec,sysvel,sourcename,project,observer = Lib.getOBSinfo(path_antlog)
print('')
print('################')
print('# step 3/4     #')
print('# get map info #')
print('################')
direction, time = Lib.PointingINFO(path_sci_raw_head+'.01.nc',path_antlog,cal=True)

freq = freq_*1.0e9

print('')
print('################')
print('# step 4/4     #')
print('# create MS2   #')
print('################')
returned_table = Lib.createMS2(MS2name,specdata,time,ra,dec,sysvel,sourcename,project,observer,direction,freq,Tsys,Tsys_time,state_id)

print('')
print('################')
print('# finished!    #')
print('################')
