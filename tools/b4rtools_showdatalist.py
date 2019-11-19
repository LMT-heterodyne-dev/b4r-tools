import LinePointing as lp
import importlib
importlib.reload(lp)
import os
import sys
import numpy as np
import datetime as dt

import PySimpleGUI as sg

if os.getenv('B4RTOOLS_PATH') == None or os.getenv('B4RTOOLS_PATH') == '':
	b4rtools_path = '.'
else:
	b4rtools_path = os.getenv('B4RTOOLS_PATH')

def makelogfile(inputlogfile,scrname,obsgoal):
	inputfiles = np.array([scrname,obsgoal])
	np.save(inputlogfile,inputfiles)

def showFileListGUI(date='',mode='',obsgoal='',path_rsync=''):
    showlist = True
    while showlist:
        FileList = lp.returnFileList(date=date,mode=mode,obsgoal=obsgoal)
        listGUI = []
        for File in FileList:
            listGUI = listGUI + [' '.join(File)]

        layout = [[sg.Text('#ObsID, Date, UTC, T_tot[s], Source name, Freq[GHz], Mode, Obs Goal')],
                  [sg.Listbox(values=listGUI,size=(100,50),key='File')],
                  [sg.OK(),sg.Button('rsync & update'), sg.Button('update')],]

        win = sg.Window('FileList - Show Filelist - B4R tools').Layout(layout)

        event, values = win.read()
        if event == 'rsync & update':
            rsync = os.system('bash '+path_rsync)
            print('rsync...')
            if not rsync == 0:
                print('rsync fail')
                sg.popup('rsync fail')
            else:
                print('rsync done')

        if event == 'rsync & update' or event == 'update':
            win.close()
        else:
            win.close()
            showlist = False
            break


def main(path_data,path_rsync):

	inputlogfile = b4rtools_path+'/.b4rtools.showfilelist.input.npy'

	layoutp = [[sg.Text('')],]

	scrname = ''
	obsgoal = ''
	if os.path.exists(inputlogfile):
		inputfiles = np.load(inputlogfile)
		scrname = inputfiles[0]
		obsgoal = inputfiles[1]

	ContinueReduc = 'Yes'

	while ContinueReduc == 'Yes':
	    stayhere = True

	    while stayhere:
	        layout = [
	        [sg.Text('Search', size=(30,1), font=('Any',12),text_color='#1c86ee' ,justification='left')],
	        [sg.Text('Date (yyyy-mm-dd)'), sg.In(dt.date.today().isoformat(),size=(14,1), key='Date'),],
	        [sg.Text('Source name'), sg.In(scrname,size=(14,1), key='Scrname'),],
	        [sg.Text('Obs goal'), sg.In(obsgoal,size=(14,1), key='Obsgoal'),],
	        [sg.Text('Mode'), sg.Checkbox('OBS',default=True)],
	        [sg.Text(''),sg.Checkbox('CALIB',default=False)],
	        [sg.OK(), sg.Button('Return to select'), sg.Cancel('Quit')],
	        ]

	        win = sg.Window('Show Filelist - B4R tools',default_element_size=(30,1),text_justification='right',auto_size_text=False).Layout(layout)


	        event, values = win.Read()
	        scrname = str(values['Scrname'])
	        date = str(values['Date'])

	        dOBS = values[0]
	        dCALIB = values[1]
	        obsgoal = values['Obsgoal']

	        makelogfile(inputlogfile,scrname,obsgoal)

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
	        lp.globBaseDir = path_data

	        if dOBS and dCALIB:
	            mode = ''
	        elif dOBS:
	            mode = 'OBS'
	        elif dCALIB:
	            mode = 'CALIB'
	        else:
	            mode = 'Other'

	        showFileListGUI(date=date,mode=mode,obsgoal=obsgoal,path_rsync=path_rsync)

	        ############

	        ContinueReduc = sg.PopupYesNo('Do you want to search again?'+'\n'+'(otherwise quit)')

	return nextstep
