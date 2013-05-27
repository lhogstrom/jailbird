#! /usr/bin/env python
'''
loop through HOG plates submit dose_plate_jobs
'''
import os
import glob as glob
import subprocess
import os
import time

baseDir = '/xchip/obelix/pod/brew'
outBase = '/xchip/cogs/projects/HOG/dose_plate_output-by_pert_id_pert_dose'
cellLst = ['A549','MCF7']
tpLst = ['6H','24H']
plateLst = ['HOG001','HOG002']
cntrlT = 'pc'

#specifications for subprocess
processes = set()
max_processes = 5

#make system call - run dose_plate_tool
for cellL in cellLst:
	for tp in tpLst:
		for platePrefix in plateLst:
			plateName = '_'.join([platePrefix,cellL,tp])
			print plateName
			inPut = glob.glob('/'.join([baseDir,
										cntrlT,
										plateName,
										'by_pert_id_pert_dose',
										plateName+'_COMPZ.MODZ_SCORE_LM_n*x978.gctx']))[0]
			# brewFolder=''
			outdir = '/'.join([outBase,plateName,cntrlT])
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			#system call to run dose tool
			cmd = ' '.join(['python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py',
					 inPut,
					 '-o ' + outdir,
					 '--no_summary_self_connect'])
			# os.system(cmd)
			# processes.add(subprocess.call([cmd]))
			processes.add(subprocess.Popen(cmd,shell=True))
			if len(processes) >= max_processes:
				os.wait()
				processes.difference_update(
					p for p in processes if p.poll() is not None)

#create hyperlink 
#make system call - run dose_plate_tool
for cellL in cellLst:
	for tp in tpLst:
		for platePrefix in plateLst:
			plateName = '_'.join([platePrefix,cellL,tp])
			outdir = '/'.join([outBase,plateName,cntrlT])
			outpath = glob.glob('/'.join([outdir,'may24/*']))[0]
			hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/HOG_dose_plates/' + plateName
			cmd = ' '.join(['ln -s',
					 outpath,
					 hyperLnkPath])
			os.system(cmd)


			