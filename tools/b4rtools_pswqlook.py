import matplotlib
matplotlib.use('TkAgg')
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np
import os
import PySimpleGUI as sg
import Lib_b4r as Lib
from astropy.stats import sigma_clipped_stats
import fmflow as fm


def makelogfile(inputlogfile,path_cal,path_sci,path_ant,chbin_in):
	inputfiles = np.array([path_cal,path_sci,path_ant,chbin_in])
	np.save(inputlogfile,inputfiles)

def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

def plotspec(Pmean,Tsys_value,sourcename,path_cal,Pol,sideband,chbin,integtime):

	layoutp = [[sg.Text('Input files')],]
	winp = sg.Window('Wait a moment...',layoutp,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)

	mean,median,rms = sigma_clipped_stats(Pmean)

	plt.close()
	plt.title(os.path.basename(path_cal))
	plt.step(Pmean['fsig'],Pmean-median) #blsub
	plt.xlabel('Frequency [GHz]')
	plt.ylabel(r'Ta$^*$ [K]')
	fig = plt.gcf()      # if using Pyplot then get the figure from the plot
	figure_x, figure_y, figure_w, figure_h = fig.bbox.bounds

	winp.close()

	stayhere = True
	while stayhere:
		layoutp = [[sg.Text('Input files')],]
		winp = sg.Window('Wait a moment...',layoutp,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)
		# create the form and show it without the plot
		layout = [[sg.Text(' Source name: '+sourcename, font='Any 14',text_color='black')],
				  [sg.Text(' Pol: '+Pol+'; Sideband: '+sideband+'; # of ch binning: '+str(chbin), font='Any 14',text_color='black')],
				  [sg.Text(' Tsys: '+Tsys_value+' K; '+'Integration time: '+str(integtime)+' sec; '+'RMS: '+'{0:4.4f}'.format(rms)+' K', font='Any 14',text_color='black')],
				  [sg.Canvas(size=(figure_w, figure_h), key='canvas')],
				  [sg.Text('Path to output file'), sg.In('', key='path_Out'), sg.Text('.png')],
		          [sg.OK(),sg.Cancel(r'Save & Quit')],
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

	plt.close()
	window.close()

def main():
	path_cal = ''
	path_sci = ''
	path_ant = ''
	chbin_in = '1'
	inputlogfile = '.b4rtools.pswqlook.input.npy'

	if os.path.exists(inputlogfile):
		inputfiles = np.load(inputlogfile)
		path_cal = inputfiles[0]
		path_sci = inputfiles[1]
		path_ant = inputfiles[2]
		chbin_in = inputfiles[3]


	ContinueReduc = 'Yes'

	while ContinueReduc == 'Yes':
		stayhere = True
		while stayhere:
			layout = [
				[sg.Text('Please enter input files', size=(18,1), font=('Any',18),text_color='#1c86ee' ,justification='left')],
				[sg.Text('Path to calibration data (binary)'), sg.In(path_cal,size=(100,1), key='path_Cal'), sg.FileBrowse()],
				[sg.Text('Path to science data (binary)'), sg.In(path_sci,size=(100,1), key='path_Sci'), sg.FileBrowse()],
				[sg.Text('Path to antenna log (netcdf)'), sg.In(path_ant,size=(100,1), key='path_Ant'), sg.FileBrowse()],
				[sg.Text('# of ch binning'), sg.In(chbin_in,size=(5,1), key='Chbin'),],
				[sg.OK(), sg.Button('Return to select'), sg.Cancel('Quit')],
					]

			win = sg.Window('PSW Qlook - B4R tools',
							default_element_size=(30,1),
							text_justification='right',
							auto_size_text=False).Layout(layout)


			event, values = win.Read()

			makelogfile(inputlogfile,path_cal,path_sci,path_ant,chbin_in)

			if event is None or event == 'Quit':
				stayhere = False
				raise SystemExit()

			elif event == 'Return to select':
				stayhere = False
				nextstep = 'select'
				win.close()
				break
			else:
				nextstep = 'quit'

			path_cal = values['path_Cal']
			path_sci = values['path_Sci']
			path_ant = values['path_Ant']
			chbin_in = values['Chbin']
			chbin = int(chbin_in)

			try:
				xffts_num_cal = int(path_cal[-1])-1
			except:
				sg.popup('[WARN] Path to calibration seems to be incorrect.'+'\n'+'(Please enter a binary file (*.0[1-4]), not netcdf file (*.nc))')
				win.Close()
				continue

			try:
				xffts_num_sci = int(path_sci[-1])-1
			except:
				sg.popup('[WARN] Path to science seems to be incorrect.'+'\n'+'(Please enter a binary file (*.0[1-4]), not netcdf file (*.nc))')
				win.Close()
				continue

			if not os.path.exists(path_cal):
				sg.popup('[WARN] Path to calibration does not exist.')
				win.Close()
				continue

			if not os.path.exists(path_sci):
				sg.popup('[WARN] Path to science does not exist.')
				win.Close()
				continue

			if not os.path.exists(path_ant):
				sg.popup('[WARN] Path to antenna log does not exist.')
				win.Close()
				continue

			if not xffts_num_cal == xffts_num_sci:
				sg.popup('[WARN] xffts No. is different between CAL and SCI.')
				win.Close()
				continue

			stayhere = False
			win.Close()

		if nextstep == 'select':
			break

		sideband_list = ['LSB','USB','LSB','USB']
		Pol_list = ['B','B','A','A']
		sideband = sideband_list[xffts_num_cal]
		Pol = Pol_list[xffts_num_cal]

		if os.path.exists(path_cal+'.nc'):
			path_cal_in = path_cal+'.nc'
		else:
			path_cal_in = path_cal
		if os.path.exists(path_sci+'.nc'):
			path_sci_in = path_sci+'.nc'
		else:
			path_sci_in = path_sci

		layout = [
			[sg.Text('Input files')],
				]

		winp = sg.Window('processing...',layout,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)

		ra,dec,sysvel,sourcename,project,observer = Lib.getOBSinfo(path_ant)
		P, Tsys, Tsys_time = Lib.loadB4Rdata(path_cal_in,path_sci_in,path_ant,sideband,cal=True,blsub=False)

		integtime = P['scantype'].shape[0]

		print(Tsys)
		print(Tsys_time)

		winp.close()

		Pmean = fm.chbinning(P,chbin).copy().mean('t')
		Tsys_value = '{0:4.2f}'.format(Tsys.median('ch').values)
		plotspec(Pmean,Tsys_value,sourcename,path_cal,Pol,sideband,chbin,integtime)

		ContinueReduc = sg.PopupYesNo('Do you want to reduce another data?'+'\n'+'(otherwise quit)')

	return nextstep



































###
