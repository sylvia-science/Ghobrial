#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:21:34 2020

@author: sujwary
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 18:43:56 2020

@author: sujwary
"""
import pandas as pd
import numpy as np
from scipy import io
import scipy.sparse
from numpy import savetxt
import bbknn
from scipy.sparse import csr_matrix
import os
import scanpy as sc
import anndata
#import igraph as ig
import scanpy.external as sce
import louvain
import igraph as ig
import pdb

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]

#sample_list = c('GL3404BM')
i = 0


patient1 = 20
patient2 = 30
folder_input = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/baseline/' + 'Patients' + str(patient1)+ '_'+str(patient2) + '/SNN_Umap/'
path = folder_input +'baseline_PCA.csv'
pca_matrix = pd.read_csv (path,index_col=0)
pca_matrix  = pca_matrix.to_numpy()
batch_list = np.genfromtxt((folder_input  +'baseline_sample.txt'), dtype=str)

path = folder_input +'MergeCounts.csv'
counts = pd.read_csv(path,index_col=0)


batch_list_df = pd.DataFrame(batch_list, index = counts.T.index).rename(columns={0: "batch"})
batch_list_df.index = batch_list_df.index.str.replace('.','-')
#adata = anndata.AnnData(X = counts)
adata = anndata.AnnData(X=pca_matrix, obs=batch_list_df)
#sc.tl.pca(adata)

adata.obsm['X_pca'] = pca_matrix

adata_bbknn = bbknn.bbknn(adata, neighbors_within_batch=5, n_pcs=30, trim=0, copy=True, batch_key = 'batch')

#sce.pp.bbknn(adata, batch_key='batch',n_pcs = 30)

sc.pp.neighbors(adata_bbknn)
sc.tl.louvain(adata_bbknn,resolution = 1)
sc.tl.umap(adata_bbknn)



sc.pl.umap(adata_bbknn, color=['batch','louvain'])


####
distances, connectivities = bbknn.bbknn_pca_matrix(pca_matrix, batch_list)


 
    
savetxt(folder_input + 'distances.csv', distances.todense(), delimiter=',')
savetxt(folder_input + 'connectivities.csv',  connectivities.todense(), delimiter=',')


## 13 samples

folder_input = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/Samples13/'

path = folder_input  +'PCA.csv'
pca_matrix = pd.read_csv (path,index_col=0)
pca_matrix  = pca_matrix.to_numpy()
sample_list = np.genfromtxt((folder_input +'sample.txt'), dtype=str)
sample_list = sample_list.reshape(-1,1)
kit_list = np.genfromtxt((folder_input +'kit.txt'), dtype=str)
kit_list = kit_list.reshape(-1,1)

batch_list = np.genfromtxt((folder_input +'batch.txt'), dtype=str)
batch_list = batch_list.reshape(-1,1)

#batch_list = np.concatenate((sample_list,kit_list),axis=1)

#batch_list= kit_list
batch_list = batch_list
distances, connectivities = bbknn.bbknn_pca_matrix(pca_matrix, batch_list)

#folder_output = folder_input + 'Batch_Batch/'+'SNN_Umap/'
folder_output = folder_input + 'Batch_Batch/'+'SNN_Umap/'
    
    
savetxt(folder_output + 'distances.csv', distances.todense(), delimiter=',')
savetxt(folder_output + 'connectivities.csv',  connectivities.todense(), delimiter=',')


####################
## All samples

#folder = 'Intra-v3_1'
folder = 'AllSamples'
folder = 'AllSamplesDownsample'
folder_input = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/' + folder +'/'

downsample = pd.read_csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv',index_col=0)
downsample = downsample['x'].to_list()

path = folder_input  +'PCA.csv'
pca_matrix = pd.read_csv (path,index_col=0)
pca_matrix  = pca_matrix.to_numpy()
sample_list = np.genfromtxt((folder_input +'sample.txt'), dtype=str)
sample_list = sample_list.reshape(-1,1)
kit_list = np.genfromtxt((folder_input +'kit.txt'), dtype=str)
kit_list = kit_list.reshape(-1,1)

path = "/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/AllSamples/Batch_Sample/SNN_Umap/GeneralCellType.csv"
GeneralCellType = pd.read_csv (path,index_col=0)

colnames = np.genfromtxt((folder_input +'rownames.txt'), dtype=str)


mask = colnames['x'].isin(downsample).to_list()
pca_matrix = pca_matrix[mask,:]
sample_list = sample_list[mask]
kit_list = kit_list[mask]
GeneralCellType = GeneralCellType[mask]
colnames = colnames[mask]
    
#batch_list = np.genfromtxt((folder_input +'batch.txt'), dtype=str)
#batch_list = batch_list.reshape(-1,1)

#batch_list = np.concatenate((sample_list,kit_list),axis=1)

#batch_list= kit_list
batch_list = sample_list

batch_list_df = pd.DataFrame(batch_list, index = colnames).rename(columns={0: "batch"})
batch_list_df.index = batch_list_df.index.str.replace('.','-')
#adata = anndata.AnnData(X = counts)
adata = anndata.AnnData(X=pca_matrix, obs=batch_list_df)
#sc.tl.pca(adata)

adata.obsm['X_pca'] = pca_matrix


adata_bbknn = bbknn.bbknn(adata, neighbors_within_batch=5, n_pcs=30, trim=0, copy=True, batch_key = 'batch')

cell_type_vector = adata_bbknn.obs['cell_type']


cell_type_vector.index = range(0,len(cell_type_vector))

N_cell_types = len(adata_bbknn.obs['cell_type'].astype('category').cat.categories)    
    
adata_bbknn.obs['cell_type'] = GeneralCellType['x'].to_list()
adata_bbknn.obs['kit'] = kit_list.tolist()

distribute_datasets(adata_bbknn)

distances = adata_bbknn.obsp['distances']
#sce.pp.bbknn(adata, batch_key='batch',n_pcs = 30)

#sc.pp.neighbors(adata_bbknn)
sc.tl.louvain(adata_bbknn,resolution = 3.4)
sc.tl.umap(adata_bbknn)


sc.pl.umap(adata_bbknn, color=['batch','louvain'])


Umap_list = adata_bbknn.obsm['X_umap']
cluster_labels = adata_bbknn.obs['louvain']
####

#distances, connectivities = bbknn.bbknn_pca_matrix(pca_matrix, batch_list)

#folder_output = folder_input + 'Batch_Batch/'+'SNN_Umap/'
folder_output = folder_input + 'Batch_Sample/'+'SNN_Umap/'
    
try:
    os.makedirs(folder_output)
except OSError:
    print ("Creation of the directory %s failed" % folder_output)
else:
    print ("Successfully created the directory %s" % folder_output)
   
    
savetxt(folder_output + 'Umap_coord.csv', Umap_list, delimiter=',')
cluster_labels.to_csv(folder_output + 'cluster_label.csv')
savetxt(folder_output + 'connectivities.csv',  connectivities.todense(), delimiter=',')
#savetxt(folder_output + 'distances.csv', distances.todense(), delimiter=',')
    
    