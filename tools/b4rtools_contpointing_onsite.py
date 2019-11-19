import matplotlib
matplotlib.use('TkAgg')
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import os
import numpy as np
import glob
import PySimpleGUI as sg
from astropy.stats import sigma_clipped_stats
import importlib
import LinePointing as lp
importlib.reload(lp)

if os.getenv('B4RTOOLS_PATH') == None or os.getenv('B4RTOOLS_PATH') == '':
	b4rtools_path = '.'
else:
	b4rtools_path = os.getenv('B4RTOOLS_PATH')

def makelogfile(inputlogfile,obsnum_in,calnum_in):
	inputfiles = np.array([obsnum_in,calnum_in])
	np.save(inputlogfile,inputfiles)

def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

def plotspecall(obj):
	plt.close()
	#plt.figure(figsize=[10,10])

	os.system('mkdir -p '+os.path.abspath(b4rtools_path)+'/log_contpointing/')
	path_Out = os.path.abspath(b4rtools_path)+'/log_linepointing/'+str(obj.obsnum)+'.linpointing'

	plt.rcParams["font.size"] = 7
	plt.figure(figsize=[6,6.5])

	#obj.LinePointingQlook()
	obj.ContPointingQlook()

	if not obj.LineSkipFlag:

		#plt.subplot(4,1,1)
		#POS = '{x:1.3f}",{y:1.3f}"'.format(x=obj.fitloc[0], y=obj.fitloc[1])
		#FWHM = '{x:1.3f}"x{y:1.3f}"'.format(x=obj.fitbeam[0], y=obj.fitbeam[1])
		#FWHM_corr = '{x:1.3f}"x{y:1.3f}"'.format(x=obj.fitbeam_corr[0], y=obj.fitbeam_corr[1])
		#titleinfo = 'ID'+str(obj.obsnum)+' Pos:'+POS+' FWHM:'+FWHM+' FWHM(corr):'+FWHM_corr
		#plt.title(titleinfo)

		#plt.subplot(4,1,2)
		#plt.title('Sourcename: '+obj.srcename)

		plt.subplots_adjust(wspace=0.0, hspace=0.62,top=0.95,bottom=0.1)
		fig = obj.outputfig
		figure_x, figure_y, figure_w, figure_h = fig.bbox.bounds

		stayhere = True
		while stayhere:
			layoutp = [[sg.Text('')],]
			winp = sg.Window('Wait a moment...',layoutp,size=(300,10),disable_close=True,default_element_size=(20,1),auto_size_text=False,finalize=True)

		    # create the form and show it without the plot
			layout = [#[sg.Text('Tsys is dummy (200K).'),sg.Text('On source time: '+str(integtime)+' sec')],
					  [sg.Text('Path to output file'), sg.In(path_Out, key='path_Out'), sg.Text('.png')],
			          [sg.OK(),sg.Cancel(r'Save & Quit')],
		              [sg.Canvas(size=(figure_w, figure_h), key='canvas')],
					  ]
			window = sg.Window('Outputs - Line Pointing - B4R tools',
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
		plt.close()


	else:
		sg.popup('SiO v=1 may not be included in the data.')

def main(path_data):
    inputlogfile = b4rtools_path +'/.b4rtools.ContPointing.onsite.input.npy'

    lp.globBaseDir = path_data

    ContinueReduc = 'Yes'

    while ContinueReduc == 'Yes':
        stayhere = True

        obsnum_in = '77777'
        calnum_in = ''

        if os.path.exists(inputlogfile):
        	inputfiles = np.load(inputlogfile)
        	obsnum_in = inputfiles[0]
        	calnum_in = inputfiles[1]

        while stayhere:
            layout = [
            [sg.Text('Please enter Obs ID and # of ch binning.', size=(30,1), font=('Any',12),text_color='#1c86ee' ,justification='left')],
            [sg.Text('Obs ID'), sg.In(obsnum_in,size=(7,1), key='Obsid'),],
            [sg.Text('Cal ID (optional)'), sg.In(calnum_in,size=(7,1), key='Calid'),],
            [sg.OK(), sg.Button('Return to select'),sg.Cancel('Quit')],
            ]

            win = sg.Window('Results - Line Pointing - B4R tools',default_element_size=(30,1),text_justification='right',auto_size_text=False).Layout(layout)


            event, values = win.Read()
            obsnum_in = values['Obsid']
            obsnum = int(obsnum_in)
            if not values['Calid']=='':
                calnum_in = int(values['Calid'])
                calnum = int(calnum_in)
            else:
                calnum_in = ''
                calnum = 0

            makelogfile(inputlogfile,obsnum_in,calnum_in)

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

            stayhere = False
            win.Close()

        if nextstep == 'select':
            break

        if not nextstep == 'skip':

            ### main ###

            if calnum == 0:
                obj = lp.LinePointing(obsnum=obsnum)
            else:
                obj = lp.LinePointing(obsnum=obsnum,calnum=calnum)
            plotspecall(obj)

            ############

            ContinueReduc = sg.PopupYesNo('Do you want to reduce another data?'+'\n'+'(otherwise quit)')

    return nextstep
