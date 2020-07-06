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
#folder = 'AllSamples'
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


GeneralCellType = pd.read_csv (folder_input +'GeneralCellType.csv',index_col=0)
GeneralCellType  = GeneralCellType.to_numpy()

#path = "/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/AllSamples/Batch_Sample/SNN_Umap/GeneralCellType.csv"
#GeneralCellType = pd.read_csv (path,index_col=0)

colnames = np.genfromtxt((folder_input +'rownames.txt'), dtype=str)




batch_list = kit_list

batch_list_df = pd.DataFrame(batch_list, index = colnames).rename(columns={0: "batch"})
batch_list_df.index = batch_list_df.index.str.replace('.','-')
#adata = anndata.AnnData(X = counts)
adata = anndata.AnnData(X=pca_matrix, obs=batch_list_df)
adata.obsm['X_pca'] = pca_matrix
adata.obs['cell_type'] = GeneralCellType
adata.obs['kit'] = kit_list.tolist()
adata.obs['sample'] = sample_list
adata_bbknn_kit = bbknn.bbknn(adata, neighbors_within_batch=5, n_pcs=30, trim=0, copy=True, batch_key = 'batch')


batch_list = sample_list

batch_list_df = pd.DataFrame(batch_list, index = colnames).rename(columns={0: "batch"})
batch_list_df.index = batch_list_df.index.str.replace('.','-')
#adata = anndata.AnnData(X = counts)
adata = anndata.AnnData(X=pca_matrix, obs=batch_list_df)
adata.obsm['X_pca'] = pca_matrix
adata.obs['cell_type'] = GeneralCellType
adata.obs['kit'] = kit_list.tolist()
adata.obs['sample'] = sample_list
adata_bbknn_sample = bbknn.bbknn(adata, neighbors_within_batch=5, n_pcs=30, trim=0, copy=True, batch_key = 'batch')
neighbor = adata_bbknn_sample.uns['neighbors'] 
connectivities = adata_bbknn_sample.obsp['connectivities'] 

folder_output_kit = folder_input + 'Batch_Kit/'+'SNN_Umap/'
folder_output_sample = folder_input + 'Batch_Sample/'+'SNN_Umap/'
distribute_datasets(adata_bbknn_kit, folder_output_kit)
distribute_datasets(adata_bbknn_sample, folder_output_sample)

#############3
#distances = adata_bbknn.obsp['distances']
connectivities_kit = adata_bbknn_kit.obsp['connectivities']
connectivities_sample = adata_bbknn_sample.obsp['connectivities']
(connectivities_kit!=connectivities_sample).nnz==0 
#sce.pp.bbknn(adata, batch_key='batch',n_pcs = 30)

#sc.pp.neighbors(adata_bbknn)
sc.tl.louvain(adata_bbknn_kit,resolution = 3.4)
sc.tl.umap(adata_bbknn_kit)
sc.pl.umap(adata_bbknn_kit, color=['batch','louvain'])
Umap_list_kit = adata_bbknn_kit.obsm['X_umap']
cluster_labels_kit = adata_bbknn_kit.obs['louvain']

sc.tl.louvain(connectivities_sample,resolution = 3.4)
sc.tl.umap(connectivities_sample)
sc.pl.umap(connectivities_sample, color=['batch','louvain'])
Umap_list_sample = connectivities_sample.obsm['X_umap']
cluster_labels_sample = connectivities_sample.obs['louvain']
####

#distances, connectivities = bbknn.bbknn_pca_matrix(pca_matrix, batch_list)


#folder_output = folder_input + 'Batch_Batch/'+'SNN_Umap/'

    
try:
    os.makedirs(folder_output)
except OSError:
    print ("Creation of the directory %s failed" % folder_output)
else:
    print ("Successfully created the directory %s" % folder_output)
   
    
savetxt(folder_output_kit + 'Umap_coord.csv', Umap_list_kit, delimiter=',')
cluster_labels_kit.to_csv(folder_output_kit + 'cluster_label.csv')
savetxt(folder_output_kit + 'connectivities.csv',  connectivities_kit.todense(), delimiter=',')
#savetxt(folder_output + 'distances.csv', distances.todense(), delimiter=',')

    
savetxt(folder_output_sample + 'Umap_coord.csv', Umap_list_sample, delimiter=',')
cluster_labels_sample.to_csv(folder_output_sample + 'cluster_label.csv')
savetxt(folder_output_sample + 'connectivities.csv',  connectivities_sample.todense(), delimiter=',')
#savetxt(folder_output + 'distances.csv', distances.todense(), delimiter=',')
    
connectivities_kit.todense()[10:20,10:20]
connectivities_sample.todense()[1:10,10:20]
    
(connectivities_kit == connectivities_sample).all()


connectivities_kit = connectivities_kit.todense()
connectivities_sample = connectivities_sample.todense()
all.equal(connectivities_kit.todense() ,connectivities_sample.todense() )
