
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

run = TRUE

sample_name <- 'GL1080BM' # 10 Pre treatment
#sample_name <- 'GL1374BM' # 10 Post treatment

#sample_name <- 'GL1024BM' # 5 Pre treatment
#sample_name <- 'GL1290BM' # 5 Post treatment

#sample_name <- metaData$Sample[6]

filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'

if (run == TRUE){
  print('running')
  
  filter <- FALSE
  folder = makeFolders(folder_base,sample_name,filter)
  data = run_pipeline(filename,folder,sample_name,sampleParam,filter)
                      
  filter <- TRUE
  folder = makeFolders(folder_base,sample_name,filter)
  data = run_pipeline(filename,folder,sample_name,sampleParam,filter)
  save(data,file=paste0(folder,'data.Robj'))
  #get_cellType(data,folder,sample_name,filter)
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


