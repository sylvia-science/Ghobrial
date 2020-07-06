#Compute entropy over AnnData objects

import argparse
import numpy as np
import pandas as pd
import scanpy as sc 
import anndata
import scipy
from math import log


def shannon_entropy (x, b_vec, N_b):
    #pdb.set_trace()
    
    tabled_values = b_vec[x > 0].value_counts()/ len(b_vec[x >0]) #class 'pandas.core.series.Series'

    tabled_val = tabled_values.tolist() 
    
    entropy = 0.0
    for element in tabled_val:
        if element != 0:
            entropy += element * log(element)
            
    entropy /= log(N_b)

    return(-entropy) #the entropy formula is the -sum, this is why we include the minus sign

def save_file_to_csv(results,folder_output):
    results.to_csv( (folder_output + 'BBKNN_k5_entropy.csv'), header = True, index = False)

def compute_entropy(df, **kwargs):
    #apply function
    
    print(kwargs)
    sample_entropy = df.apply(shannon_entropy, axis=0, args=(kwargs['sample_vector'],kwargs['N_samples']))
    kit_entropy = df.apply(shannon_entropy, axis=0, args=(kwargs['kit_vector'],kwargs['N_kits']))
    celltype_entropy = df.apply(shannon_entropy, axis=0, args=(kwargs['cell_type_vector'] ,kwargs['N_cell_types']))
    print("Entropy calculated!")
    
    #pdb.set_trace()
    results = {'sample_entropy': sample_entropy,'kit_entropy': kit_entropy, "celltype_entropy":celltype_entropy}
    results = pd.concat(results, axis = 1, keys = ['sample_entropy','kit_entropy', 'celltype_entropy'])

    save_file_to_csv(results,kwargs['folder_output'] )

def distribute_datasets(dataset, folder_output):
    kwargs = {}
    
    kwargs['folder_output'] = folder_output
    #pdb.set_trace()
    batch_key = 'kit'
    celltype_key = 'cell_type'
    print(batch_key)
    
    
    kwargs['sample_vector'] = dataset.obs['sample']
    print(celltype_key)
    #modify index of batch vector so it coincides with matrix's index
    kwargs['sample_vector'].index = range(0,len(kwargs['sample_vector']))
    #number of batches
    kwargs['N_samples'] = len(dataset.obs['sample'].astype('category').cat.categories)


    #batch vector(batch id of each cell)
    kwargs['kit_vector'] = dataset.obs[batch_key]
    print(celltype_key)
    #modify index of batch vector so it coincides with matrix's index
    kwargs['kit_vector'].index = range(0,len(kwargs['kit_vector']))
    #number of batches
    kwargs['N_kits'] = len(dataset.obs['batch'].astype('category').cat.categories)

    print(celltype_key)
    
    #cell_type vector( betch id of each cell)
    kwargs['cell_type_vector'] = dataset.obs[celltype_key]
    #modify index of cell_type vector so it coincides with matrix's index
    kwargs['cell_type_vector'].index = range(0,len(kwargs['cell_type_vector']))
    #number of cell_types
    kwargs['N_cell_types'] = len(dataset.obs['cell_type'].astype('category').cat.categories)    
    

    try:
        knn_graph = dataset.uns['neighbors']
        print('BBKNN corrected object!') 
    
    except KeyError:
        #Both: pre corrected logcounts and Scanorama counts enter through this way.
        #compute neighbors
        sc.tl.pca(dataset, n_comps = args.n_pcs)
        sc.pp.neighbors(dataset, n_neighbors = args.n_neighbors, knn = True)
                        
    #knn graph
    knn_graph = dataset.uns['neighbors']['connectivities']
    #transforming csr_matrix to dataframe
    df = pd.DataFrame(knn_graph.toarray())
    
    compute_entropy(df, **kwargs)

def read_h5ad(dataset):
    
    #read input h5ad
    return sc.read(dataset)
    print("File read!")
    
