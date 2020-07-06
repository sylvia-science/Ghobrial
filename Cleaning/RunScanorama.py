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
import scanorama
from scipy.sparse import csr_matrix
import os


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = pd.read_excel(filename_testIntRun)


#folder = 'Inter-version'
#sample_list = (Samples_runs['Samples'][Samples_runs['Folder']== folder]).iloc[0]
#sample_list = [x.strip() for x in sample_list.split(',')]

folder = 'AllSamples'
folder = 'AllSamplesDownsample'
sample_list = metaData['Sample'].tolist()


#sample_list = c('GL3404BM')
i = 0

patient_list = [5, 6, 10]
patient_list = [12, 16, 20]
for i in range(0,len(patient_list)):
    patient = patient_list[i]
    
    datasets = [None]*2 #List of data sets (matrices of cells-by-genes):
    genes_list = [None]*2 # List of gene lists:
    sample_name_baseline = metaData['Sample'][(metaData['Patient Number'] == patient) &  (metaData['Treatment'] == 'baseline')]
    sample_name_baseline=  sample_name_baseline.iloc[0]
    print(sample_name_baseline)
    folder_input = ('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/' + sample_name_baseline + '/')
    sparsematrix = io.mmread(folder_input +sample_name_baseline +'_matrix.txt')
    row_names = np.genfromtxt((folder_input + sample_name_baseline +'_rownames.txt'), dtype=str)
    col_names = np.genfromtxt((folder_input + sample_name_baseline +'_colnames.txt'), dtype=str)
    sparsematrix_T = csr_matrix(np.transpose(sparsematrix))
    datasets[0] = sparsematrix_T
    genes1 =  pd.read_csv(folder_input + sample_name_baseline +'_var_feature.csv' ,index_col=0)
    genes1['x'] = genes1['x'].apply(str)
    genes1 = list(genes1['x'].to_numpy())
    genes_list[0] = row_names
    
    sample_name_C9D1 = metaData['Sample'][(metaData['Patient Number'] == patient) &  (metaData['Treatment'] == 'C9D1')]
    sample_name_C9D1=  sample_name_C9D1.iloc[0]
    print(sample_name_C9D1)
    folder_input = '/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/' + sample_name_C9D1 + '/'
    sparsematrix = io.mmread(folder_input + sample_name_C9D1 +'_matrix.txt')
    row_names = np.genfromtxt((folder_input +sample_name_C9D1 +'_rownames.txt'), dtype=str)
    col_names = np.genfromtxt((folder_input +sample_name_C9D1 +'_colnames.txt'), dtype=str)
    sparsematrix_T = csr_matrix(np.transpose(sparsematrix))
    datasets[1] = sparsematrix_T
    genes2 =  pd.read_csv(folder_input + sample_name_C9D1 +'_var_feature.csv' ,index_col=0)
    genes2['x'] = genes2['x'].apply(str)
    genes2 = list(genes2['x'].to_numpy())
    
    genes_list[1] = row_names
    
    integrated, corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=True)
    
    path_base = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/Patient'+ str(patient) + '/'
    try:
        os.makedirs(path_base)
    except OSError:
        print ("Creation of the directory %s failed" % path_base)
    else:
        print ("Successfully created the directory %s" % path_base)
       
    corrected[0] = corrected[0].todense()
    corrected[1] = corrected[1].todense()
    file0 = path_base + sample_name_baseline + '.csv' 
    file1= path_base + sample_name_C9D1 + '.csv' 
    file_gene= path_base + 'genes' + '.csv'
    
    savetxt(file0,  corrected[0], delimiter=',')
    savetxt(file1,  corrected[1], delimiter=',')
    
    file0 = path_base + sample_name_baseline + '_integrated.csv' 
    file1= path_base + sample_name_C9D1 + '_integrated.csv' 

    savetxt(file0,  integrated[0], delimiter=',')
    savetxt(file1,  integrated[1], delimiter=',')



    #savetxt(file_gene ,  genes, delimiter=',')
    pd.DataFrame(genes).to_csv(file_gene,header=None)
    
    
patient1 = 20
patient2 = 30

datasets = [None]*2 #List of data sets (matrices of cells-by-genes):
genes_list = [None]*2 # List of gene lists:
sample_name_baseline1 = metaData['Sample'][(metaData['Patient Number'] == patient1) &  (metaData['Treatment'] == 'baseline')]
sample_name_baseline1=  sample_name_baseline1.iloc[0]
print(sample_name_baseline1)
folder_input = ('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/' + sample_name_baseline1 + '/')
sparsematrix = io.mmread(folder_input +sample_name_baseline1 +'_matrix.txt')
row_names = np.genfromtxt((folder_input + sample_name_baseline1 +'_rownames.txt'), dtype=str)
col_names = np.genfromtxt((folder_input + sample_name_baseline1 +'_colnames.txt'), dtype=str)
sparsematrix_T = csr_matrix(np.transpose(sparsematrix))
datasets[0] = sparsematrix_T
genes_list[0] = list(row_names)

sample_name_baseline2 = metaData['Sample'][(metaData['Patient Number'] == patient2) &  (metaData['Treatment'] == 'baseline')]
sample_name_baseline2=  sample_name_baseline2.iloc[0]
print(sample_name_baseline2)
folder_input = '/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/' + (sample_name_baseline2) + '/'
sparsematrix = io.mmread(folder_input + sample_name_baseline2 +'_matrix.txt')
row_names = np.genfromtxt((folder_input +sample_name_baseline2 +'_rownames.txt'), dtype=str)
col_names = np.genfromtxt((folder_input +sample_name_baseline2 +'_colnames.txt'), dtype=str)
sparsematrix_T = csr_matrix(np.transpose(sparsematrix))
datasets[1] = sparsematrix_T
genes_list[1] = list(row_names)

integrated, corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=True)
path_base = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/baseline/' + 'Patients' + str(patient1)+ '_'+str(patient2) + '/'

#path_base = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/Baseline' + '/'
try:
    os.makedirs(path_base)
except OSError:
    print ("Creation of the directory %s failed" % path_base)
else:
    print ("Successfully created the directory %s" % path_base)
   

#integrated[0] = integrated[0].todense()
#integrated[1] = integrated[1].todense()
corrected[0] = corrected[0].todense()
corrected[1] = corrected[1].todense()

file0 = path_base + sample_name_baseline1 + '_integrated.csv' 
file1= path_base + sample_name_baseline2 + '_integrated.csv' 

savetxt(file0,  integrated[0], delimiter=',')
savetxt(file1,  integrated[1], delimiter=',')

file0 = path_base + sample_name_baseline1 + '_corrected.csv' 
file1= path_base + sample_name_baseline2 + '_corrected.csv' 

savetxt(file0,  corrected[0], delimiter=',')
savetxt(file1,  corrected[1], delimiter=',')

#savetxt(file_gene ,  genes, delimiter=',')
file_gene= path_base + 'genes' + '.csv'
pd.DataFrame(genes).to_csv(file_gene,header=None)



###############
## Many samples
###############


#folder = 'Samples13   
 
downsample = pd.read_csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv',index_col=0)
downsample = downsample['x'].to_list()

downsample = [s.split("_")[0] for s in downsample]

datasets = [None]*len(sample_list) #List of data sets (matrices of cells-by-genes):
genes_list = [None]*len(sample_list) # List of gene lists:

for i in range(0,len(sample_list)):
    sample_name = sample_list[i]
    print(sample_name)
    folder_input = ('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/' + sample_name + '/')

    sparsematrix = pd.read_csv (folder_input +sample_name +'_NormData.csv',index_col=0)
    sparsematrix  = sparsematrix.to_numpy()

    row_names = pd.read_csv((folder_input + sample_name +'_rownames.csv'), index_col=0)['x'].to_list()
    col_names = pd.read_csv((folder_input + sample_name +'_colnames.csv'), index_col=0)
    var_genes= pd.read_csv((folder_input + sample_name +'_varFeatures.csv'), index_col=0)

    mask = col_names['x'].isin(downsample).to_list()
    sparsematrix = sparsematrix[:,mask]
    col_names = col_names[mask]
    
    sparsematrix_T = (np.transpose(sparsematrix)) # Turn into cell by genes
    datasets[i] = sparsematrix_T
    genes_list[i] = row_names
    #genes_list[i] = list(var_genes['x'])

integrated, corrected, genes = scanorama.correct(datasets, genes_list, return_dimred=True)

path_base = '/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/'+folder+'/SNN_Umap/'

try:
    os.makedirs(path_base)
except OSError:
    print ("Creation of the directory %s failed" % path_base)
else:
    print ("Successfully created the directory %s" % path_base)
   

#integrated[0] = integrated[0].todense()
#integrated[1] = integrated[1].todense()
    
for i in range(0,len(sample_list)):
    sample_name = sample_list[i]
    print(sample_name)
    corrected[i] = corrected[i].todense()
    
    
    file0 = path_base + sample_name + '_integrated.csv' 
    
    
    savetxt(file0,  integrated[i], delimiter=',')
    
    
    file0 = path_base + sample_name + '_corrected.csv' 
    
    
    savetxt(file0,  corrected[i], delimiter=',')


#savetxt(file_gene ,  genes, delimiter=',')
file_gene= path_base + 'genes' + '.csv'
pd.DataFrame(genes).to_csv(file_gene,header=None)
