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

source("/home/sujwary/Desktop/scRNA/Code/Functions.R")
source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')

folder_base = '/home/sujwary/Desktop/scRNA/Output/'

# filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
# sampleParam <- read_excel(filename_sampleParam)
# 
# filename_sampleParam_integrate <-'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Combine_parameters.xlsx'
# sampleParam_integrate = read_excel(filename_sampleParam_integrate)
# 
# filename_sample_Integrate_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Combine_pairs.xlsx'
# sample_Integrate_pairs <- read_excel(filename_sample_Integrate_pairs)



# filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/EloRD Meta.xlsx'
# metaData <- read_excel(filename_metaData)

sample_type = 'PrePostEOTNBM'
integrate_merge = 'Integrate'
folder_base_output = paste0('/home/sujwary/Desktop/scRNA/Output/',integrate_merge ,' All/',sample_type,'/')
dir.create( folder_base_output, recursive = TRUE)

print('Start')


data_list = vector(mode = "list", length = 0)
data_list = vector(length = 0)
sample_name_list = vector(mode = "list", length = 0)

for(i in 1:nrow(metaData)){  #nrow(metaData)
    
  
  sample_name = metaData$Sample[i]
  Treatment = metaData$Treatment[i]
  dexa = metaData$'Dexa or not'[i]
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i = load_data(filename)
  data_i[["percent.mt"]] = PercentageFeatureSet(data_i, pattern = "^MT-")
  
  nFeature_RNA_min = as.numeric(nFeature_RNA_list[1])
  nFeature_RNA_max = as.numeric(nFeature_RNA_list[2])

  # Threshohold
  expr <- FetchData(object = data_i, vars = 'nFeature_RNA')
  data_i = data_i[, which(x = expr > nFeature_RNA_min & expr < nFeature_RNA_max)]
  
  expr <- FetchData(object = data_i, vars = 'percent.mt')
  data_i = data_i[, which(x = expr < percent_mt)]
  data_i = NormalizeData(data_i, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i = FindVariableFeatures(data_i, verbose = FALSE, selection.method = "vst", nfeatures =2000)

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
  }else if (sample_type == 'PrePostEOTNBM'){
    condition = TRUE
  }else if (sample_type == 'NBM'){
    condition = dexa == 'NBM'
  }
  
  if (condition){
    data_i@meta.data$orig.ident = paste0('data_',Treatment)
    data_i@meta.data$data_pre_renamed = data_i@active.ident
  
    data_i@meta.data$sample_name = sample_name
    data_i@meta.data$dexa = dexa
    data_i@meta.data$CHIP =  metaData$'CHIP or not'[i]
    
    data_i@meta.data$'Patient Number' = metaData$'Patient Number'[i]
    data_i@meta.data$Response = metaData$'Response'[i]
    data_i@meta.data$'10X kit' = metaData$'10X kit'[i]
    data_i@meta.data$'MM' = TRUE
    data_i@meta.data$'batch' = metaData$'batch'[i]
    print(data_i)
    data_list = c(data_list, data_i)
    sample_name_list= c(sample_name_list,sample_name)
  }
      
}

write.csv(sample_name_list, file =paste0(folder_base_output,'sample_name_list','_',sample_type,'.csv'),row.names = FALSE,col.names = FALSE)


if (integrate_merge == 'Integrate'){
  

  # anchors = FindIntegrationAnchors(object.list = data_list,k.filter =100,anchor.features = 1000)
  # save(anchors,file=paste0(folder_base_output,'anchors_',integrate_merge,'_',sample_type,'.Robj'))
  # 
  # 
  # path = paste0(folder_base_output,'anchors_',integrate_merge,'_',sample_type,'.Robj')
  # print(path)
  # anchors = loadRData(path)
  # 
  # data_integrate = IntegrateData(anchorset = anchors)
  # save(data_integrate,file=paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,'.Robj'))
  # 

  #browser()
  # run using rpca
  
   features <- SelectIntegrationFeatures(object.list = data_list)
   data_list <- lapply(X = data_list, FUN = function(x) {
     x <- ScaleData(x, features = features, verbose = FALSE)
     x <- RunPCA(x, features = features,npcs = 50, verbose = FALSE)
   })
  
   filter = TRUE
   if (filter == FALSE){
     filterTF = '_filterF'
   }else{
     filterTF = ''
   }
   
  #anchors = FindIntegrationAnchors(object.list = data_list,k.filter =30, k.score = 30,anchor.features = 2000)
  #save(anchors,file=paste0(folder_base_output,'anchors_',integrate_merge,'_rpca_TEST_',sample_type,filterTF,'.Robj'))
  #data_integrate = IntegrateData(anchorset = anchors, dims = 1:5)
  
  anchors_rpca = FindIntegrationAnchors(object.list = data_list,k.filter =100,anchor.features = 2000,dims = 1:30,reduction = "rpca")
  path = paste0(folder_base_output,'anchors_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj')
  save(anchors_rpca,file=path)
  
  anchors_rpca = loadRData(path)
  data_integrate_rpca = IntegrateData(anchorset = anchors_rpca,dims = 1:30)
  
  path = paste0(folder_base_output,'data_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj')
  save(data_integrate_rpca,file=path)
  browser()
  # online fix
  # data_list <- lapply(X = data_list, FUN = function(x) {
  #   x <- SCTransform(x)
  # })
  # 
  # features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
  # data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)
  # data_list <- lapply(X = data_list, FUN = RunPCA, verbose = FALSE, features = features)
  # anchors <- FindIntegrationAnchors(object.list = hca.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
  # data_integrate <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  # path  = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,'.Robj')
  # save(data_integrate,file=path)
  # 
}

browser()
if (integrate_merge == 'Merge'){
  data_merge = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
  for (i in 3:length(data_list)){
    print(i)
    data_merge = merge(x =  data_merge,y = data_list[[i]],merge.data = FALSE)
  }
  save(data_merge,file=paste0(folder_base_output,'data_Merge_filter','_',sample_type,'.Robj'))
}


