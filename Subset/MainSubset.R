
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
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Subset/runPipelineSubset.R')


folder_base_output <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_sampleParam_integrate <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters.xlsx'
sampleParam_integrate <- read_excel(filename_sampleParam_integrate)


filename_sampleParam_subset <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters_Tcell.xlsx'
sampleParam_integrate_subset <- read_excel(filename_sampleParam_subset)

filename_sample_Integrate_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Integrate_pairs.xlsx'
sample_Integrate_pairs <- read_excel(filename_sample_Integrate_pairs)


sample_type = 'BM'
#sample_type = 'PB'
cell_type = 'T Cell'

print('Start')
patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all
#patient_list = c(16, 51, 6, 40) # v3
#patient_list = c(34, 28, 21, 31, 16, 51, 6, 40) # all

patient_list_dexaF = c(5,12, 16)
patient_list_dexaT = c(10,20,34,28,21,31,51,6,40)

run = TRUE

for(patient in patient_list){ # Patient numbers 
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$`Patient Number` == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_output_main = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,TRUE)
  folder_output = paste0(folder_output_main,'/Subset/',cell_type,'/')
  
  dir.create( paste0(folder_output,'QC'), recursive = TRUE)
  dir.create( paste0(folder_output,'PCA'), recursive = TRUE)
  dir.create( paste0(folder_output,'Cluster'), recursive = TRUE)
  dir.create( paste0(folder_output,'Cell Type'), recursive = TRUE)
  dir.create( paste0(folder_output,'DE'), recursive = TRUE)
  dir.create( paste0(folder_output,'Stats'), recursive = TRUE)

    if (run){
    
      data <- loadRData(paste0(folder_output_main,'data.Robj'))
        
      
      
      print('Starting Run')
      print(paste0('Sample: ', sample_name_pre))
        
      print(paste0('folder: ', folder_output))
        
      cluster_IDs = sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Sample'] == sample_name_pre]
      data = label_cells(data,cluster_IDs)
        
      ##
      ## Cluster T cells seperately and add back to original data
      ##
      #browser()
      data_subset_orig = subset(data, idents = cell_type)
      data_subset = runPipelineSubset(data_subset_orig,folder_output,sample_name_pre,sampleParam_integrate_subset,scale_TF = FALSE)
      save(data_subset,file=paste0(folder_output,'data.Robj'))

    }else{

      print('Starting Plotting')
      
      print(paste0('Sample: ', sample_name_pre))
      
      print(paste0('folder: ', folder_output))
      
      data = loadRData(paste0(folder_output,'data.Robj'))
      
      plotAll(data,folder_output,sample_name_pre,sampleParam_integrate_subset,label_TF = FALSE)
      #plotAll(data,folder_output,sample_name_pre,sampleParam_integrate_subset, label_TF = TRUE)
      
      # Plot Markers
      get_cellType(data,data,folder_output,sample_name_pre) 
      
      cluster_IDs <- sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Sample'] == sample_name_pre]
      data = label_cells(data,cluster_IDs)
      
      
    }
  
}

#data = loadRData('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/GL1080BM/Filtered/Regress/Subset/T Cell/data.Robj')
