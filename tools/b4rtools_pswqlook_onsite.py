import matplotlib
matplotlib.use('TkAgg')
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np
import glob
import os
import PySimpleGUI as sg
from astropy.stats import sigma_clipped_stats

if os.getenv('B4RTOOLS_PATH') == None or os.getenv('B4RTOOLS_PATH') == '':
	b4rtools_path = '.'
else:
	b4rtools_path = os.getenv('B4RTOOLS_PATH')

def makelogfile(inputlogfile,obsid,calid,chbin_in):
	inputfiles = np.array([obsid,calid,chbin_in])
	np.save(inputlogfile,inputfiles)

def updateSymlinks(path_xffts_teslx6,path_xffts_links_teslx6):
    import os
    import re
    import struct

    XFFTS_DIR = os.path.expanduser(path_xffts_teslx6)
    LINKS_DIR = os.path.expanduser(path_xffts_links_teslx6)
    os.system('mkdir -p '+LINKS_DIR)
    XFFTS_PATTERN = r"^xffts20[0-9]{12}\.xfftsx\.0[1-4]$"

    def is_exists(path):
        """"""
        path = os.path.expanduser(path)

        if not os.path.exists(path):
            raise FileNotFoundError("{path}: not found".format(path=path))

        return path


    def create_symlink(path):
        """
        """
        path = is_exists(path)

        with open(path, "rb") as f:
            f.seek(32)
            obsnum_bin = f.read(8)

        try:
            obsnum = struct.Struct("q").unpack(obsnum_bin)[0]
        except struct.error:
            obsnum = "no_lmt"

        fname = os.path.basename(path)
        #new_fname = "{obsnum}_{fname}".format(obsnum=obsnum, fname=fname)
        new_fname = str(obsnum).zfill(6) + '_' + fname

        new_path = os.path.expanduser(
            "{links_dir}/{new_fname}".format(
                links_dir=LINKS_DIR,
                new_fname=new_fname
            )
        )

        try:
            os.symlink(path, new_path)
            #os.system('ln -sf '+path+' '+new_path)
        except OSError:
            # rprint("{new_path}: already exists".format(new_path=new_path))
            pass

        return

    pattern = re.compile(XFFTS_PATTERN)
    for xffts in os.listdir(XFFTS_DIR):
        if pattern.match(xffts):
            create_symlink(XFFTS_DIR + "/" + xffts)



def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

def dataCal(path_sci,path_cal,path_ant,chbin,xffts_num,nocaldata=True):
	import Lib_b4r as Lib
	import fmflow as fm

	sideband_list = ['LSB','USB','LSB','USB']
	Pol_list = ['B','B','A','A']
	sideband = sideband_list[xffts_num]
	Pol = Pol_list[xffts_num]

	if os.path.exists(path_sci+'.nc'):
	    path_sci_in = path_sci+'.nc'
	else:
	    path_sci_in = path_sci

	if not nocaldata:
		if os.path.exists(path_cal+'.nc'):
		    path_cal_in = path_cal+'.nc'
		else:
		    path_cal_in = path_cal

	layout = [[sg.Text('')],]

	winp = sg.Window('processing...',layout,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)

	ra,dec,sysvel,sourcename,project,observer = Lib.getOBSinfo(path_ant)
	if nocaldata:
		P, Tsys, Tsys_time = Lib.loadB4Rdata('',path_sci_in,path_ant,sideband,cal=True,blsub=False,qlook=nocaldata)
	else:
		P, Tsys, Tsys_time = Lib.loadB4Rdata(path_cal_in,path_sci_in,path_ant,sideband,cal=True,blsub=False,qlook=nocaldata)

	integtime = P['scantype'].shape[0]

	winp.close()

	Pmean = fm.chbinning(P,chbin).copy().mean('t')
	Tsys_value = '{0:4.2f}'.format(Tsys.median('ch').values)

	return Pmean, Tsys_value, sourcename, integtime, Pol, sideband

def plotspec(path_sci_head,path_cal_head,path_ant,chbin,obsid):
	subplot_order = [221,222,223,224]
	RMS_list = np.zeros(4)

	for xffts_num in range(4):
		path_sci = path_sci_head + '.' + str(xffts_num+1).zfill(2)
		if path_cal_head == None:
			Pmean, Tsys_value, sourcename, integtime, Pol, sideband = dataCal(path_sci,None,path_ant,chbin,xffts_num,nocaldata=True)
			dammy_tsys = '(dammy)'
		else:
			path_cal = path_cal_head + '.' + str(xffts_num+1).zfill(2)
			Pmean, Tsys_value, sourcename, integtime, Pol, sideband = dataCal(path_sci,path_cal,path_ant,chbin,xffts_num,nocaldata=False)
			dammy_tsys = ''

		mean,median,rms = sigma_clipped_stats(Pmean)
		RMS_list[xffts_num] = rms

		titlename = 'ID'+str(obsid)+' '+sourcename+' Pol:'+Pol+' Sideband:'+sideband+'\n'

		titlename = titlename +'Tsys:'+'{0:4.2f}'.format(float(Tsys_value))+dammy_tsys+' ton:'+str(integtime)+'s RMS:'+'{0:2.4f}'.format(rms)+'K '

		plt.subplot(2,2,xffts_num+1)
		plt.title(titlename)
		plt.step(Pmean['fsig'],Pmean-median) #blsub
		plt.xlabel('Frequency [GHz]')
		plt.ylabel(r'Ta$^*$ [K]')
		plt.ylim(-5*rms,Pmean.max()+rms)

	fig = plt.gcf()
	return fig, integtime, Tsys_value

def plotspecall(fig,integtime,Tsys_value,obsid,chbin):

    os.system('mkdir -p '+os.path.abspath(b4rtools_path)+'/log_qlookpsw/')
    path_Out = os.path.abspath(b4rtools_path)+'/log_qlookpsw/'+obsid+'.qlook.chbin'+str(chbin)
    figure_x, figure_y, figure_w, figure_h = fig.bbox.bounds

    stayhere = True
    while stayhere:
    	layoutp = [[sg.Text('')],]
    	winp = sg.Window('Wait a moment...',layoutp,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)

        # create the form and show it without the plot
    	layout = [[sg.OK(),sg.Cancel(r'Save & Quit')],
    			  [sg.Text('Path to output file'), sg.In(path_Out, key='path_Out'), sg.Text('.png')],
				  [sg.Canvas(size=(figure_w, figure_h), key='canvas')],
    			  ]
    	window = sg.Window('Spectrum - PSW Qlook - B4R tools',
    	                   layout, finalize=True)

    	# add the plot to the window
    	winp.close()
    	fig_canvas_agg = draw_figure(window['canvas'].TKCanvas, fig)

    	event, values = window.read()
    	path_out = values['path_Out']+'.png'

    	if event == r'Save & Quit' and path_out == '.png':
    		sg.popup('[WARN] Please enter path to output.')
    		window.Close()
    		continue

    	elif event is None or event == r'Save & Quit':
    		if os.path.exists(path_out):
    			overwrite = sg.PopupYesNo(path_out+'\n'+'aleady exits.'+'\n'+'Anyway overwrite?')
    			if overwrite == 'Yes':
    				plt.savefig(path_out)
    				break
    			else:
    				window.Close()
    				continue
    		else:
    			plt.savefig(path_out)
    			break

    	stayhere = False

    window.close()

def main(path_xffts_teslx6,path_lmttpm_teslx6,xffts_links_name='xffts_links'):
	chbin_in = '1'
	obsid = '77777'
	calid = '77778'
	inputlogfile = b4rtools_path+'/.b4rtools.pswqlook.onsite.input.npy'

	path_xffts_links_teslx6 = os.path.dirname(path_xffts_teslx6) + '/' + xffts_links_name

	#layoutp = [[sg.Text('')],]
	#winp = sg.Window('updating symlinks...',layoutp,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)
	#updateSymlinks(path_xffts_teslx6,path_xffts_links_teslx6)
	#winp.close()

	if os.path.exists(inputlogfile):
		inputfiles = np.load(inputlogfile)
		obsid = inputfiles[0]
		calid = inputfiles[1]
		chbin_in = inputfiles[2]


	ContinueReduc = 'Yes'

	while ContinueReduc == 'Yes':
		stayhere = True
		while stayhere:
			layout = [
				[sg.Text('Please enter Obs ID and # of ch binning.', size=(30,1), font=('Any',12),text_color='#1c86ee' ,justification='left')],
				[sg.Text('Obs ID'), sg.In(obsid,size=(7,1), key='Obsid')],
				[sg.Text('Cal ID (optional)'), sg.In(calid,size=(7,1), key='Calid')],
				[sg.Text('# of ch binning'), sg.In(chbin_in,size=(7,1), key='Chbin'),],
				[sg.OK(), sg.Button('Return to select'), sg.Button('Update symlinks'), sg.Cancel('Quit')],]
			win = sg.Window('PSW Qlook - B4R tools',default_element_size=(30,1),text_justification='right',auto_size_text=False).Layout(layout)
			event, values = win.Read()
			obsid = values['Obsid']
			calid = values['Calid']
			chbin_in = values['Chbin']
			makelogfile(inputlogfile,obsid,calid,chbin_in)

			if event is None or event == 'Quit':
			    stayhere = False
			    raise SystemExit()

			elif event == 'Return to select':
			    stayhere = False
			    nextstep = 'select'
			    win.close()
			    break

			elif event == 'Update symlinks':
			    win.close()
			    layoutp = [[sg.Text('')],]
			    winp = sg.Window('updating symlinks...',layoutp,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)
			    updateSymlinks(path_xffts_teslx6,path_xffts_links_teslx6)
			    winp.close()
			    nextstep = 'skip'
			    break

			else:
			    nextstep = 'quit'

			xffts_name_01  = glob.glob(path_xffts_links_teslx6+'/' + obsid.zfill(6) + '_xffts*sx.01')
			xffts_cal_name_01  = glob.glob(path_xffts_links_teslx6+'/' + calid.zfill(6) + '_xffts*sx.01')
			lmttpm_name_list = glob.glob(path_lmttpm_teslx6+'/lmttpm_*' + obsid.zfill(6) + '*.nc')

			if len(xffts_name_01) < 1:
			    sg.popup('[WARN] XFFTS data of Obs ID '+obsid.zfill(6)+' dose not exist.' )
			    win.Close()
			    continue

			if len(lmttpm_name_list) < 1:
			    sg.popup('[WARN] Antenna log of Obs ID '+obsid.zfill(6)+' dose not exist.' )
			    win.Close()
			    continue

			if len(xffts_cal_name_01) < 1:
				xffts_cal_name_01 = None


			stayhere = False
			win.Close()

		if nextstep == 'select':
			break

		if not nextstep == 'skip':
			path_sci_head = xffts_name_01[0][:-3]
			if xffts_cal_name_01 == None:
				path_cal_head = None
			else:
				path_cal_head = xffts_cal_name_01[0][:-3]

			path_ant = lmttpm_name_list[0]
			chbin = int(chbin_in)

			### main ###

			plt.close()
			plt.rcParams["font.size"] = 7
			plt.figure(figsize=[8,5])
			plt.subplots_adjust(wspace=0.4, hspace=0.4,top=0.9,bottom=0.1)
			fig, integtime, Tsys_value = plotspec(path_sci_head,path_cal_head,path_ant,chbin,obsid)
			plotspecall(fig,integtime,Tsys_value,obsid,chbin)
			plt.close()

			############

			ContinueReduc = sg.PopupYesNo('Do you want to reduce another data?'+'\n'+'(otherwise quit)')

	return nextstep



































###
