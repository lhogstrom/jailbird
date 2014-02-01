#! /usr/bin/env python

'''
Example script to run clique_significane module

LH 02/2014
'''
import os
import cmap.analytics.clique_significance as cls

wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/meta_compound_classes_vs_dmso_null'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

reload(cls)
cDir = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_groups_n147/summly'
cr = cls.CliqueResult(cDir)
cr.load_clique_res(group_min=2)
nullMtrx = '/xchip/cogs/projects/connectivity/null/dmso/lass_n1000x7147.gctx'
cr.load_null_matrix(nullMtrx)
