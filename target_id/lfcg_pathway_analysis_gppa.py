#! /usr/bin/env python
'GPPA (genomic perturbation pathay analyzer) - cluster genomic perturbation data acording to pathways'

import numpy as np
import os
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.io.gmt as gmt
from matplotlib import cm
from cmap.analytics.queryer import Queryer
import cmap.analytics.gppa as gppa
from cmap.io import queryresult

wkdir = '/xchip/cogs/projects/target_id/pathway_clustering/lfcg_gppa_KD_spearman'
if not os.path.exists(wkdir):
    os.makedirs(wkdir)

#pathway annotations from reactome 
pathGMT = '/xchip/cogs/projects/target_id/KD_pathway_clustering/ReactomePathways.gmt'
gmtDict = gmt.read(pathGMT)
pathwayDict = {}
for dict1 in gmtDict:
    pathwayDict[dict1['id']] = dict1['sig']

# aprioriPathways = ['Cholesterol biosynthesis', 
#   'p53-Dependent G1 DNA Damage Response', 
#   'p53-Dependent G1/S DNA damage checkpoint', 
#   'Antigen processing: Ubiquitination & Proteasome degradation', 
#   'Regulation of activated PAK-2p34 by proteasome mediated degradation', 
#   'Signaling by TGF-beta Receptor Complex',
#   'mTOR signalling']

# Lessons from the cancer genome pathways:
lfcgPathways = ['p38MAPK events',
  'ERK/MAPK targets',
  'MAPK targets/ Nuclear events mediated by MAP kinases',
  'Negative regulation of the PI3K/AKT network',
  'PI3K Cascade',
  'PI3K/AKT Signaling in Cancer',
  'Constitutive PI3K/AKT Signaling in Cancer',
  'Signaling by NOTCH1 in Cancer',
  'FBXW7 Mutants and NOTCH1 in Cancer',
  'Signaling by NOTCH',
  'Signaling by Wnt',
  'WNT ligand biogenesis and trafficking',
  'Signaling by TGF-beta Receptor Complex',
  'Downregulation of TGF-beta receptor signaling',
  'TAK1 activates NFkB by phosphorylation and activation of IKKs complex',
  'Methylation',
  'mRNA Splicing',
  'mRNA Splicing - Major Pathway',
  'Antigen processing: Ubiquitination & Proteasome degradation',
  'Extension of Telomeres',
  'Intrinsic Pathway for Apoptosis',
  'Regulation of Apoptosis',
  'Apoptosis']

# mystring.replace (" ", "_")
# 'Cholesterol biosynthesis' # vs statin signature?
for pathway1 in lfcgPathways:
    pathGenes = pathwayDict[pathway1]  
    pathName = pathway1.replace(" ","_")
    pathName = pathName.replace("/","_")
    pathName = pathName.replace("&","_")
    pathName = pathName.replace(":","_")
    ### run analysis with gppa object
    gp_type = 'KD'
    metric = 'spearman'
    out = wkdir
    gppaObj = gppa.GPPA(pathGenes, 
                        pathName,
                        gp_type, 
                        metric, 
                        out, 
                        row_space = 'lm')
    gppaObj.find_cell_lines()
    for cell in gppaObj.cell_lines:
        if not os.path.exists(out + '/null_queries/' + cell):
            #run new null query
            gppaObj.null_queryer(cell,n_rand=200)
        else: 
            #load exisiting null query
            queryer = Queryer()
            queryer.set_params(out=out + '/null_queries/' + cell,
                                metric = metric)
            fout = queryer.get_result_file()
            qres =  queryresult.QueryResult()
            qres.read(fout)
            gppaObj.randScores = qres.score
        gppaObj.pathway_queryer(cell)
        gppaObj.connection_dist(cell)
    gppaObj.make_html_report()

#make general index file
indexfile = os.path.join(out,'index.html')
with open(indexfile,'w') as f:
    lineWrite = '<h2>Pathways: </h2>'
    f.write(lineWrite + '\n')
    for pathway1 in lfcgPathways:
        pathName = pathway1.replace(" ","_")
        pathName = pathName.replace("/","_")
        pathName = pathName.replace("&","_")
        pathName = pathName.replace(":","_")
        lineWrite =  '<a href="' + pathName + '/index.html">' + pathName + '</a> <BR>'
        f.write(lineWrite + '\n')

# to-do:
# heatmap threshold
# test spearman
# OE - display gene names on axis
# test cancer genome pathways (+ all genes together)

