import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import os
import sys
import PySimpleGUI as sg

if os.getenv('B4RDATA_PATH') == None or os.getenv('B4RDATA_PATH') == '':
	path_data = '.'
else:
	path_data = os.getenv('B4RDATA_PATH')

if os.getenv('B4RDATA_RSYNC_LOC') == None or os.getenv('B4RDATA_RSYNC_LOC') == '':
	path_rsync = './rsync_from_b4r1.sh'
else:
	path_rsync = os.getenv('B4RDATA_RSYNC_LOC')

if os.getenv('B4RTOOLS_PATH') == None or os.getenv('B4RTOOLS_PATH') == '':
	b4rtools_path = '.'
else:
	b4rtools_path = os.getenv('B4RTOOLS_PATH')
sys.path.append(b4rtools_path+'/tools')

refuctype_list = ['Show Filelist', 'PSW Qlook', 'Line Pointing', 'Cont Pointing', 'MS2 Pipeline']

# confing
#path_xffts_teslx6 = '/home/tescam/b4r/log/xffts_links/'
#path_lmttpm_teslx6 = '/home/tescam/b4r/log/lmttpm/'

# main
stayhere = True
while stayhere:
	layout = [
		[sg.Text('Please select tool', size=(12,1), font=('Any',18),text_color='#1c86ee' ,justification='left')],
		[sg.Text('Path to data:'),],
		[sg.In(path_data,size=(40,1),key='Path_data')],
		[sg.Text(''),],
		[sg.Text('Select:'),],
		[sg.Listbox(default_values='Show Filelist',values=refuctype_list,size=(40,6), key='tool')],
		[sg.OK(), sg.Cancel('Quit')],
			]

	win = sg.Window('B4R tools',default_element_size=(30,1),text_justification='left',auto_size_text=False).Layout(layout)

	event, values = win.Read()
	if event is None or event =='Quit':
		raise SystemExit()

	try:
		reductype = values['tool'][0]
	except:
		sg.popup('[WARN] please select tool.')
		win.close()
		continue

	path_data = values['Path_data']
	if not path_data[-1] == '/':
		path_data = path_data + '/'
	print(path_data)
	path_xffts_teslx6  = path_data + 'xffts/'
	path_lmttpm_teslx6 = path_data + 'lmttpm/'

	if reductype == 'PSW Qlook':
		if not os.path.exists(path_xffts_teslx6):
			sg.popup('[WARN] '+'"'+path_xffts_teslx6+'"'+' dose not exist.')
			nextstep = 'skip'
			win.close()
			continue

		elif not os.path.exists(path_lmttpm_teslx6):
			sg.popup('[WARN] '+'"'+path_lmttpm_teslx6+'"'+' dose not exist.')
			nextstep = 'skip'
			win.close()
			continue

		else:
			import b4rtools_pswqlook_onsite as LibPswQlookOnsite
			win.close()
			nextstep = LibPswQlookOnsite.main(path_xffts_teslx6,path_lmttpm_teslx6)

	elif reductype == 'Line Pointing':
		import b4rtools_linepointing_onsite as LibLinepointingOnsite
		#sg.popup('[WARN] "Line Pointing" mode is under development.')
		win.close()
		nextstep = LibLinepointingOnsite.main(path_data)

	elif reductype == 'Cont Pointing':
		import b4rtools_contpointing_onsite as LibContpointingOnsite
		#sg.popup('[WARN] "Line Pointing" mode is under development.')
		win.close()
		nextstep = LibContpointingOnsite.main(path_data)
		#sg.popup('[WARN] "Cont Pointing" mode is under development.')
		#win.Close()
		#continue
		#import b4rtools_contpointing_onsite as LibContPointingOnsite
		#win.Close()
		#nextstep = LibContPointingOnsite.main(path_data)

	elif reductype == 'Show Filelist':
		import b4rtools_showdatalist as LibLShowDataList
		win.Close()
		nextstep = LibLShowDataList.main(path_data,path_rsync)

	elif reductype == 'MS2 Pipeline':
		sg.popup('[WARN] "MS2 pipeline" mode is under development.')
		win.Close()
		continue

	if nextstep == 'select':
		staythere = True
	else:
		stayhere=False
#########
