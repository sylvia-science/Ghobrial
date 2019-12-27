
# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Subset/MainSubset.R')
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(stringr)
library(data.table)


#library(biomaRt)
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_Func.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Subset/runPipelineSubset.R')


folder_base_output <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Merge All/PrePostNBM/'

filename_sampleParam_subset <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Combine_parameters_Tcell.xlsx'
sampleParam_Combine_subset <- read_excel(filename_sampleParam_subset)

filename_sample_Combine_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Combine_pairs.xlsx'
sample_Combine_pairs <- read_excel(filename_sample_Combine_pairs)


cell_type = 'T Cell'
integrate_merge = 'Merge'
sample_type = 'PrePostNBM'

if (integrate_merge == 'Integrate' || integrate_merge == 'Merge'){
  print(paste0('integrate_merge:', integrate_merge))
  filename_sampleParam <- paste0('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_','Combine','_parameters.xlsx')
  sampleParam <- read_excel(filename_sampleParam)
}else{
  filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
  sampleParam <- read_excel(filename_sampleParam)
}

print('Start')
patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all
#patient_list = c(16, 51, 6, 40) # v3
#patient_list = c(34, 28, 21, 31, 16, 51, 6, 40) # all

patient_list_dexaF = c(5,12, 16)
patient_list_dexaT = c(10,20,34,28,21,31,51,6,40)

patient_list = c('MergeAll_PrePostNBM')
run = TRUE

for(patient in patient_list){ # Patient numbers 
  pair_list =  sample_Combine_pairs[ sample_Combine_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre BM')]]
  folder_output_main = folder_base_output #makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,TRUE)
  folder_output = paste0(folder_output_main,'Subset/',cell_type,'/')
  
  dir.create( paste0(folder_output,'QC'), recursive = TRUE)
  dir.create( paste0(folder_output,'PCA'), recursive = TRUE)
  dir.create( paste0(folder_output,'Cluster'), recursive = TRUE)
  dir.create( paste0(folder_output,'Cell Type'), recursive = TRUE)
  dir.create( paste0(folder_output,'DE'), recursive = TRUE)
  dir.create( paste0(folder_output,'Stats'), recursive = TRUE)

    if (run){
      
      path = paste0(folder_output_main,'data_run_',integrate_merge,'_PCAdim',20,'_',sample_type,'.Robj')
      print(path)
      data = loadRData(path)
  
      print('Starting Run')
      print(paste0('Sample: ', sample_name_pre))
        
      print(paste0('folder: ', folder_output))
        
      cluster_IDs = sampleParam[['Cluster_IDs']][sampleParam['Sample'] == sample_name_pre]
      data = label_cells(data,cluster_IDs)
        
      ##
      ## Cluster T cells seperately and add back to original data
      ##
      #browser()
      data_subset_orig = subset(data, idents = cell_type)
      data_subset = runPipelineSubset(data_subset_orig,folder_output,sample_name_pre,sampleParam,scale_TF = FALSE)
      save(data_subset,file=paste0(folder_output,'data.Robj'))

    }else{

      print('Starting Plotting')
      
      print(paste0('Sample: ', sample_name_pre))
      
      print(paste0('folder: ', folder_output))
      
      data = loadRData(paste0(folder_output,'data.Robj'))
      
      plotAll(data,folder_output,sample_name_pre,sampleParam_Combine_subset,label_TF = FALSE)
      #plotAll(data,folder_output,sample_name_pre,sampleParam_Combine_subset, label_TF = TRUE)
      
      # Plot Markers
      get_cellType(data,data,folder_output,sample_name_pre) 
      
      cluster_IDs <- sampleParam[['Cluster_IDs']][sampleParam['Sample'] == sample_name_pre]
      data = label_cells(data,cluster_IDs)
      
      
    }
  
}

#data = loadRData('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/GL1080BM/Filtered/Regress/Subset/T Cell/data.Robj')
