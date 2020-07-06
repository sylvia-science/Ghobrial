# source('/home/sujwary/Desktop/scRNA/Code/Integrate All/MainIntegrateAll.R')
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
library(hdf5r)
require(gridExtra)
#library(biomaRt)
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
#source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
#source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')

folder_base_input <- '/home/sujwary/Desktop/scRNA/Output/'
folder_base = '/home/sujwary/Desktop/scRNA/Output/'

filename_sampleParam <- '/home/sujwary/Desktop/scRNA/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_sampleParam_integrate <-'/home/sujwary/Desktop/scRNA/Data/sample_Combine_parameters.xlsx'
sampleParam_integrate = read_excel(filename_sampleParam_integrate)

filename_metaData <- '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData <- read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sample_type = 'Tcell_MergeInt'
integrate_merge = 'Integrate'
folder_base_output = paste0('/home/sujwary/Desktop/scRNA/Output/',integrate_merge ,' All/',sample_type,'/')
dir.create( folder_base_output, recursive = TRUE)

print('Start')

filter = TRUE
if (filter == FALSE){
  filterTF = '_filterF'
}else{
  filterTF = ''
}

data_list = vector(mode = "list", length = 0)
data_list = vector(length = 0)
sample_name_list = vector(mode = "list", length = 0)

for(i in 1:nrow(metaData)){
  
  
  sample_name = metaData$Sample[i]
  Treatment = metaData$Treatment[i]
  dexa = metaData$'Dexa or not'[i]
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  if (metaData$`10X kit`[i] == 'Microwell-seq'){
    filename <- paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_dge.txt",sep = "")
    
    data_i = read.table(file = filename,row.names = 1,header = T)
    data_i = CreateSeuratObject(data_i)
  }else{
    filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
    data_i = load_data(filename)
  }
  data_i[["percent.mt"]] = PercentageFeatureSet(data_i, pattern = "^MT-")
  
  nFeature_RNA_min = as.numeric(nFeature_RNA_list[1])
  nFeature_RNA_max = as.numeric(nFeature_RNA_list[2])
  
  # Threshohold
  expr <- FetchData(object = data_i, vars = 'nFeature_RNA')
  data_i = data_i[, which(x = expr > nFeature_RNA_min & expr < nFeature_RNA_max)]
  
  expr <- FetchData(object = data_i, vars = 'percent.mt')
  data_i = data_i[, which(x = expr < percent_mt)]

  if (length(data_i$orig.ident) > 100){
    
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
    data_i@meta.data$'SampleType' = metaData$'Sample Type'[i]
    
    print(data_i)
    
    
    data_list = c(data_list, data_i)
    sample_name_list= c(sample_name_list,sample_name)
  }
  
}
browser()
if (integrate_merge == 'Merge'){
  data_merge = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
  for (i in 3:length(data_list)){
    print(i)
    data_merge = merge(x =  data_merge,y = data_list[[i]],merge.data = FALSE)
  }
  
  path_merge = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,filterTF,'.Robj')
  
  save(data_merge,file=path_merge)
}


data_merge_1 = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
data_merge_2 = merge(x =  data_list[[3]],y = data_list[[4]], merge.data = FALSE)

data_merge_1 = NormalizeData(data_merge_1, normalization.method = "LogNormalize", scale.factor = 10000)
data_merge_1 = FindVariableFeatures(data_merge_1, verbose = FALSE, selection.method = "vst", nfeatures =2000)

data_merge_2 = NormalizeData(data_merge_2, normalization.method = "LogNormalize", scale.factor = 10000)
data_merge_2 = FindVariableFeatures(data_merge_2, verbose = FALSE, selection.method = "vst", nfeatures =2000)


data_list = vector(mode = "list", length = 2)
data_list[[1]] = data_merge_1
data_list[[2]] = data_merge_2

anchors = FindIntegrationAnchors(object.list = data_list, 
                                 dims = 1:30, k.filter =30, k.score = 30,
                                 anchor.features = 2000)
data_integrate = IntegrateData(anchorset = anchors, 
                               dims = 1:30,k.weight  = 100)
path_integrate = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,filterTF,'.Robj')
save(data_integrate,file=path_integrate)
