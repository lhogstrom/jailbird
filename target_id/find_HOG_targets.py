import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo

# work_dir = '/xchip/cogs/projects/target_id/HOG_25June2013'
# if not os.path.exists(work_dir):
#     os.mkdir(work_dir)

### make target_dict
# targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_short.txt'
# '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
targetSheetF = '/xchip/cogs/projects/target_id/7June2014/A2_DrugBank_targets_tab.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			splt2 = line.split('\t')
			pID = splt2[0] #the pert_id listed the line
			pDesc = splt2[1]
			targets = splt2[2:]
			targets = [x for x in targets if x != '']
			# targets = targets.split(';')
			if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
				continue
			else:
				targetDict[pID] = targets
				pDescDict[pID] = pDesc

targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			splt2 = line.split('\t')
			pID = splt2[0] #the pert_id listed the line
			pDesc = splt2[1]
			targets = splt2[2]
			targets = targets.split(';')
			targets = [x for x in targets if x != '']
			if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
				continue
			else:
				targetDict[pID] = targets
				pDescDict[pID] = pDesc