import pandas as pd
import numpy as np
import scrublet as scr


folder_base_input = '/home/sujwary/Desktop/scRNA/Output/'
folder_base = '/home/sujwary/Desktop/scRNA/Output/'

filename_sampleParam = '/home/sujwary/Desktop/scRNA/Data/sample_parameters.xlsx'
sampleParam = pd.read_excel(filename_sampleParam)

filename_sampleParam_integrate ='/home/sujwary/Desktop/scRNA/Data/sample_Combine_parameters.xlsx'
sampleParam_integrate = pd.read_excel(filename_sampleParam_integrate)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1,]
sample_type = 'PrePostEOTNBM_MT15'
integrate_merge = 'Integrate'
folder_base_output = ('/home/sujwary/Desktop/scRNA/Output/' + integrate_merge + ' All/' + sample_type +'/')

#data_i_filtered = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
#data_i_filtered = CreateSeuratObject(counts = data_i_filtered, project = "BM", min.cells = 3,min.features = 100)
  

for i in range(1,(metaData.shape[1])):
    sample_name = metaData['Sample'][i]
  
    print(sample_name)

    # Run on filtered with 100 min counts
    sparsematrix = io.mmread('/home/sujwary/Desktop/scRNA/Data/RawMatrix/'+sample_name +'_matrix.txt')

    sparsematrix_T = np.transpose(sparsematrix)
   # m_dense = sparsematrix.toarray()
    
    #var_names = np.genfromtxt(('/home/sujwary/Desktop/scRNA/Data/RawMatrix/'+sample_name +'_rownames.txt'), dtype=str)
    #col_names = np.genfromtxt(('/home/sujwary/Desktop/scRNA/Data/RawMatrix/'+sample_name +'_colnames.txt'), dtype=str)
    
    # Export to txt:
    #df = pd.DataFrame(m_dense, columns=col_names, index=var_names)
    
    scrub = scr.Scrublet(sparsematrix_T)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    np.savetxt(('/home/sujwary/Desktop/scRNA/Output/Doublets/' + sample_name + '_doublet_scores.csv'), doublet_scores, delimiter=',')
    np.savetxt(('/home/sujwary/Desktop/scRNA/Output/Doublets/' + sample_name + '_predicted_doublets.csv'), predicted_doublets, delimiter=',')
