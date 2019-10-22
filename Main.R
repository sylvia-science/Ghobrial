
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

run = TRUE

sample_name <- 'GL1080BM' # 10 Pre treatment
sample_name <- 'GL1374BM' # 10 Post treatment

#sample_name <- 'GL1024BM' # 5 Pre treatment
#sample_name <- 'GL1290BM' # 5 Post treatment



folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

regress_TF = TRUE

for(i in 1:6){
  sample_name <- metaData$Sample[i]
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  if (run == TRUE){
    print('running')
    
    filter <- FALSE
    folder = makeFolders(folder_base,sample_name,filter)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
                        
    filter <- TRUE
    folder = makeFolders(folder_base,sample_name,filter)
    data = run_pipeline(filename,folder,sample_name,sampleParam,filter,regress_TF)
    
    if (regress_TF){
      subfolder = 'Regress/' # Save files in regression folder
    }else{
      subfolder = 'No_Regress/' # Save files in No regression folder
    }
    print('done!')
    print(paste0(folder,subfolder,'data.Robj'))
    
    param_data = sampleParam[ sampleParam$Sample ==sample_name , ] 
    
    save(data,file=paste0(folder,subfolder,'data.Robj'))
    write.csv(param_data, file = paste0(folder,subfolder,'parameters.csv'),row.names=FALSE)
    
    data_orig = load_data(filename)
    get_cellType(data,data_orig,folder,sample_name,filter)
    #label_cells(data,folder,sample_name,sampleParam,filter)
  }else{
    print('here')
    # Load Data
    filter <- TRUE
    
    folder = makeFolders(folder_base,sample_name,filter)
      
    data <- loadRData(paste0(folder,'data.Robj'))
    plotAll(data,folder,sample_name)
    get_cellType(data,folder,sample_name,filter)
  }
}

