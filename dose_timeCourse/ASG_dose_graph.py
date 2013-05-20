#!/bin/py
#run file like this
#execfile('/xchip/cogs/hogstrom/scripts/dose_timeCourse/ASG_dose_graph.py')

import cmap.analytics.sc as sc
import os.path
import cmap.analytics.signature_strength as ss
import numpy
import scipy

cellLine = 'MCF7'
timeP = '24H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 
outdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/%s_sigs/ssCCplots_%s' % (cellLine,timeP)

#creaete outdir if it does not exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

sco = sc.SC()
sco.add_sc_from_gctx_meta(gctfile)
sco.set_thresh_by_specificity(0.8)
#sco.plot(include='trich') #basic plot
dose = [float(x.split('::')[0].split(':')[2]) for x in sco.pid]
cmpds = [x.split(':::')[0].split('::')[1] for x in sco.pid]
#sco.plot(include=cmpds[0],size=dose, title='trichostatin A', out='/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/tmp.png')

### calculate alternative ss ###
#SS = ss.SigStrength()
#SS.sig_strength_from_gct_file(gctfile,do_zthresh=True)
##SS.sig_strength_from_gct_file(gctfile,do_zthresh=False)
#ssList = SS.ss
#ssArray = scipy.array(ssList) #convert to array
#sco.s = ssArray #overwrite sco sig strength with new metric


uCmpds = set(cmpds)
uCmpds = list(uCmpds)
uCmpds.remove('DMSO') #remove dmso because it doesn't have multiple doses

##loop though only unique cmpds - make sc plot

for ic, c in enumerate(uCmpds):
    #fname = '%s.png' % cmpds[ic]
    fname = '%s.png' % c
    fout = os.path.join(outdir, fname)
    #sco.plot(include=cmpds[ic],size=dose, title=cmpds[ic], out=fout)
    graphTitle = '%s - %s %s' % (c,cellLine,timeP)
    sco.plot(include=c,size=dose, title=graphTitle, out=fout)

## build a dose response clasifier
#how to get column info from gctObj?

#SS = ss.SigStrength()
#SS.sig_strength_from_gct_file(gctfname)
#mtrx = SS.GCTObj.matrix
#n_col = mtrx.shape[1]
#for i in range(n_col):
	#col = mtrx[:,i]
	