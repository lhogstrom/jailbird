#! /usr/bin/env python
'''
Acording to STITCH, what are the genes that relate to the to the TP53 drugs?

how do the summly results relate to these?

'''
import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd
import matplotlib

work_dir = '/xchip/cogs/hogstrom/analysis/summly/TP53'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

stitchF = '/xchip/cogs/hogstrom/notes/stitch_db/STITCH_drug_target_800%2B_activators.txt'

