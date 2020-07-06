#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:18:11 2020

@author: sujwary
"""
import pandas as pd
import numpy as np
import scrublet as scr
import matplotlib.pyplot as plt
from scipy import io
import os
filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]

#metaData =metaData['Run'] == 1
filename_sample_parameters = '/home/sujwary/Desktop/scRNA/Data/sample_parameters.xlsx'
sample_parameters = pd.read_excel(filename_sample_parameters)
i = 4
# i = 18 doesn't have auto estimate
for i in range(0,(metaData.shape[0])):
    
    #i = 3
    sample_name = metaData['Sample'][i]
    threshold = sample_parameters['Scrublet_threshold'][sample_parameters['Sample'] == sample_name].iloc[0]

    print(i)
    print(sample_name)
    print(threshold)
    folder_input = '/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/' + sample_name + '/'
    sparsematrix = io.mmread(folder_input +sample_name +'_matrix.txt')
    col_names = np.genfromtxt((folder_input +sample_name +'_colnames.txt'), dtype=str)

    sparsematrix_T = np.transpose(sparsematrix)
    
    scrub = scr.Scrublet(sparsematrix_T)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    path_base = '/home/sujwary/Desktop/scRNA/Output/Soup_MT_Scrublet/Umap Plots/'+ sample_name + '/' + str(threshold) + '/'
    try:
        os.makedirs(path_base)
    except OSError:
        print ("Creation of the directory %s failed" % path_base)
    else:
        print ("Successfully created the directory %s" % path_base)
    
    if i != 18:
        histo = scrub.plot_histogram()
        path= path_base +sample_name + 'histoDefault.png'
        histo[0].savefig(path)
    else:
        predicted_doublets = [0]*len(doublet_scores)
        
        
    output = scrub.call_doublets(threshold=threshold)
    histo = scrub.plot_histogram()
    histo[0].savefig(path_base + sample_name + 'histo' + str(threshold) + '.png')
    
    #np.savetxt(('/home/sujwary/Desktop/scRNA/Output/Doublets/' + sample_name + '_predicted_doublets' + threshold,'.csv'), predicted_doublets, delimiter=',')
    
    data = np.array([col_names, output,doublet_scores,predicted_doublets])
    df = pd.DataFrame(data=data).T
    
    #df.rename(columns =  {0:'Cell',1:'Scrublet_Boolean',2:'doublet_scores',3:'predicted_doublet'})
    
    path = path_base + sample_name + '_scrublet' + str(threshold) + '.csv'
    df.rename(columns =  {0:'Cell',1:'Scrublet_Boolean',2:'doublet_scores',3:'predicted_doublet'}).to_csv(path, index = False)
