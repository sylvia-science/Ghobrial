#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:01:44 2020

@author: sujwary
"""
import anndata
from pylab import rcParams
import matplotlib as mpl

import rcParams
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

base_input = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/NK/Cluster/PCA30/res3/Data/'

sample = anndata.read_loom(base_input + "data.loom")
 

sample_obs = pd.read_csv(base_input + "cellID_obs.csv")
umap_cord = pd.read_csv(base_input + "cell_embeddings.csv")
cell_clusters = pd.read_csv(base_input + "clusters.csv")
