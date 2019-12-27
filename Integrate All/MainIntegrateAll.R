# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/MainIntegrateAll.R')
# Libraries
gc()

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)

require(gridExtra)

#library(biomaRt)
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
#source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
#source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

folder_base_input <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'
folder_base = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'


filename_sampleParam_integrate <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters.xlsx'
sampleParam_integrate <- read_excel(filename_sampleParam_integrate)



filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_sampleParam_integrate <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters.xlsx'
sampleParam_integrate <- read_excel(filename_sampleParam_integrate)

filename_sample_Integrate_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Integrate_pairs.xlsx'
sample_Integrate_pairs <- read_excel(filename_sample_Integrate_pairs)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

sample_type = 'NBM'
sample_name = paste0('IntegrateAll_',sample_type)
integrate_merge = 'Merge'
folder_base_output = paste0('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/',integrate_merge ,' All/',sample_type,'/')
dir.create( folder_base_output, recursive = TRUE)

print('Start')


data_list = vector(mode = "list", length = 0)
data_list = vector(length = 0)
sample_name_list = vector(mode = "list", length = 0)
for(i in 1:nrow(metaData)){  #nrow(metaData)
    
  sample_name = metaData$Sample[i]
  filename = paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  Treatment = metaData$Treatment[i]
  dexa = metaData$'Dexa or not'[i]
    
  folder = makeFolders(folder_base,sample_name,filter = TRUE,regress_TF= TRUE,makeFolder_TF = FALSE)
  print(paste0('folder: ', folder))
  data_i = loadRData(paste0(folder,'data.Robj'))
  
  if (sample_type == 'Pre'){
    condition = Treatment == 'pre' && dexa != 'NBM'
  }else if (sample_type == 'Post'){
    condition = Treatment == 'post' && dexa != 'NBM'
  }else if (sample_type == 'PrePost'){
    condition = (Treatment == 'pre' || Treatment == 'post') && dexa != 'NBM'
  }else if (sample_type == 'PreNBM'){
    condition = Treatment == 'pre' || dexa == 'NBM'
  }else if (sample_type == 'PostNBM'){
    condition = Treatment == 'post' || dexa == 'NBM'
  }else if (sample_type == 'PrePostNBM'){
    condition = TRUE
  }else if (sample_type == 'NBM'){
    condition = dexa == 'NBM'
  }
  
  if (condition){
    data_i@meta.data$orig.ident = paste0('data_',Treatment)
    data_i@meta.data$data_pre_renamed = data_i@active.ident
  
    data_i@meta.data$sample_name = sample_name
    data_i@meta.data$dexa = dexa
    data_i@meta.data$'Patient Number' = metaData$'Patient Number'[i]
    data_i@meta.data$Response = metaData$'Response'[i]
    data_i@meta.data$'10X Kit' = metaData$'10X Kit'[i]
    data_i@meta.data$'MM' = TRUE
      
    data_list = c(data_list, data_i)
    sample_name_list= c(sample_name_list,sample_name)
  }
      
}

write.csv(sample_name_list, file =paste0(folder_base_output,'sample_name_list','_',sample_type,'.csv'),row.names = FALSE,col.names = FALSE)

if (integrate_merge == 'Integrate'){
  #anchors = FindIntegrationAnchors(object.list = data_list,k.filter =100)
  #save(anchors,file=paste0(folder_base_output,'anchors_',integrate_merge,'_',sample_type,'.Robj'))
  #data_integrate = IntegrateData(anchorset = anchors)
  #save(data_integrate,file=paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,'.Robj'))
}

if (integrate_merge == 'Merge'){
  data_merge = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
  for (i in 3:length(data_list)){
    print(i)
    data_merge = merge(x =  data_merge,y = data_list[[i]],merge.data = FALSE)
  }
  save(data_merge,file=paste0(folder_base_output,'data_Merge','_',sample_type,'.Robj'))
}


