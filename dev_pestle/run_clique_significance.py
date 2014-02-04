#! /usr/bin/env python

'''
Example script to run clique_significane module

LH 02/2014
'''
import os
import cmap.analytics.clique_significance as cls

# wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/meta_compound_classes_vs_dmso_null_uQ'
# wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/significance/cpd_groups_n147_vs_dmso_uQ'
wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/significance/cpd_targets_n368_vs_random_median'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

reload(cls)
# cDir = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_groups_n147/summly'
cDir = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_targets_n368/summly'
cr = cls.CliqueResult(cDir)
cr.set_output_dir(wkdir)
cr.load_clique_res(group_min=2)
# nullMtrx = '/xchip/cogs/projects/connectivity/null/dmso/lass_n1000x7147.gctx'
nullMtrx = '/xchip/cogs/projects/connectivity/null/random/lass_n1000x7147.gctx'
cr.load_null_matrix(nullMtrx)
cr.group_p_values(test_statistic='median',n_perm=1000,graph=True) # upper_quartile_median
cr.fdr_correction()
cr.store_parameters_rpt()
cr.make_summary_table()

# to-do:
# make a seperate figures directory
# fix quartile connection
# what about sym vs. non-sym

#run matlab html page generator
# tbl = parse_tbl('/xchip/cogs/projects/connectivity/null/clique_analysis/significance/cpd_groups_n147_vs_dmso_uQ/clique_group_qval_summary.txt','outfmt','record')
# fig = parse_tbl('/xchip/cogs/projects/connectivity/null/clique_analysis/significance/cpd_groups_n147_vs_dmso_uQ/figure_table.txt','outfmt','record')
# mk_html_table('/xchip/cogs/projects/connectivity/null/clique_analysis/significance/cpd_groups_n147_vs_dmso_uQ/index.html',tbl,'figures',fig)

# ln -s /xchip/cogs/projects/connectivity/null/clique_analysis/significance/cpd_groups_n147_vs_dmso_uQ /xchip/cogs/web/icmap/hogstrom/clique_significance/cpd_groups_n147_vs_dmso_uQ




