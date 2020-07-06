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


sample_type = 'PrePostEOTNBM_MT15'
integrate_merge = 'Integrate'
celltype = 'T Cell'
#############################

folder_base_output = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM_MT15/Regkit/'
rpca = ''
PCA_dim = 30

path = paste0(folder_base_output,
              '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
              '.Robj')
data_integrate_subset = loadRData(path) 
cell_list = c('0', '2', '3', '4', '5', '6', '7', '9', '10', '11', '13', '22', '28', '37', '38')
data_integrate_subset = data_integrate_subset[,Idents(data_integrate_subset) %in% cell_list]


#############################


folder_base_input = '/home/sujwary/Desktop/scRNA/Output/'
folder_base = '/home/sujwary/Desktop/scRNA/Output/'

filename_sampleParam = '/home/sujwary/Desktop/scRNA/Data/sample_parameters.xlsx'
sampleParam = read_excel(filename_sampleParam)

filename_sampleParam_integrate ='/home/sujwary/Desktop/scRNA/Data/sample_Combine_parameters.xlsx'
sampleParam_integrate = read_excel(filename_sampleParam_integrate)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]
sample_type = 'PrePostEOTNBM_MT15'
integrate_merge = 'Integrate'
folder_base_output = paste0('/home/sujwary/Desktop/scRNA/Output/',integrate_merge ,' All/',sample_type,'/')
dir.create( folder_base_output, recursive = TRUE)

print('Start')

data_list = vector(length = 0)
sample_name_list = vector(mode = "list", length = 0)

for(i in 1:nrow(metaData)){  #nrow(metaData)
  
  sample_name = metaData$Sample[i]
  Treatment = metaData$Treatment[i]
  dexa = metaData$'Dexa or not'[i]
  nFeature_RNA_list = list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i = load_data(filename)
  data_i[["percent.mt"]] = PercentageFeatureSet(data_i, pattern = "^MT-")
  
  nFeature_RNA_min = as.numeric(nFeature_RNA_list[1])
  nFeature_RNA_max = as.numeric(nFeature_RNA_list[2])
  
  # Can't use subset function because doesn't work with variables
  # This is a workaround
  expr = FetchData(object = data_i, vars = 'nFeature_RNA')
  data_i = data_i[, which(x = expr > nFeature_RNA_min & expr < nFeature_RNA_max)]
  
  expr = FetchData(object = data_i, vars = 'percent.mt')
  data_i = data_i[, which(x = expr < percent_mt)]
  
  #data_i_normalize = NormalizeData(data_i, normalization.method = "LogNormalize", scale.factor = 10000)
  #data_i_varFeatures = NormalizeData(data_i_normalize, normalization.method = "LogNormalize", scale.factor = 10000)

  
  if (sample_name %in% unique(data_integrate_subset$sample_name)){
    data_integrate_subset_sample = data_integrate_subset[,data_integrate_subset$sample_name == sample_name]
    
    data_integrate_subset_sample$cell_names = data_integrate_subset_sample@assays[["integrated"]]@data@Dimnames[[2]]
    data_integrate_subset_sample$cell_names = gsub("_.*", "", data_integrate_subset_sample$cell_names)
    
    data_i$cell_names = data_i@assays[["RNA"]]@data@Dimnames[[2]]
    
    data_i_subset =  data_i[,data_i$cell_names %in% data_integrate_subset_sample$cell_names]
    
    data_i_subset_normalize = NormalizeData(data_i_subset, normalization.method = "LogNormalize", scale.factor = 10000)
    data_i_subset_var_features= FindVariableFeatures(data_i_subset_normalize, selection.method = "vst", 
                                  nfeatures = 2000, verbose = FALSE)
    
    #data_i_subset_var_features = data_i_subset
    #cell_names = data_integrate_subset_sample@assays[["integrated"]]@data@Dimnames[[2]]
    
    data_i_subset_var_features@meta.data$Treatment = paste0(Treatment)
    data_i_subset_var_features@meta.data$data_pre_renamed = data_i_subset@active.ident
    
    data_i_subset_var_features@meta.data$sample_name = sample_name
    data_i_subset_var_features@meta.data$dexa = dexa
    data_i_subset_var_features@meta.data$'Patient Number' = metaData$'Patient Number'[i]
    data_i_subset_var_features@meta.data$Response = metaData$'Response'[i]
    data_i_subset_var_features@meta.data$'10X kit' = metaData$'10X kit'[i]
    data_i_subset_var_features@meta.data$'MM' = TRUE
    num_cells = length(data_i_subset_var_features@assays[["RNA"]]@data@p )
    if (num_cells> 31){
      #print(data_i_subset_var_features)
      data_list = c(data_list, data_i_subset_var_features)
      sample_name_list= c(sample_name_list,sample_name)
    }
  }
  
}

folder_base_output = paste0('/home/sujwary/Desktop/scRNA/Output/',
                            integrate_merge ,' All/',sample_type,'',clean,'',str,'/')
folder_base_output = paste0(folder_base_output,'SubsetIntegrate/',celltype,'/')
dir.create( folder_base_output, recursive = TRUE)


if (integrate_merge == 'Integrate'){
  
  for (i in 1:length(data_list)){
    print(i)
    data_i = data_list[i][[1]]
    print(length(data_i@assays[["RNA"]]@counts@Dimnames[[2]]))
  }
  
  
  
  # Regular Integration
  # dim:	Which dimensions to use from the CCA to specify the neighbor search space
  anchors = FindIntegrationAnchors(object.list = data_list, 
                                   dims = 1:30, k.filter =30, k.score = 30,
                                   anchor.features = 2000)
  path = paste0(folder_base_output,'anchors_',integrate_merge,'_',sample_type,'.Robj')
  save(anchors,file=path)
  anchors = loadRData(path)
  
  # dim: Number of PCs to use in the weighting procedure
  data_integrate_subset = IntegrateData(anchorset = anchors, 
                                 dims = 1:30,k.weight  = 100)
  path_integrate = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,'.Robj')
  save(data_integrate_subset,file=path_integrate)
  
  # RPCA
  features <- SelectIntegrationFeatures(object.list = data_list)
  data_list_PCA <- lapply(X = data_list, FUN = function(x) {
    print(unique(x$sample_name))
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features,npcs = 30, verbose = T)
  })
  
  anchors_rpca = FindIntegrationAnchors(object.list = data_list_PCA,
                                        dims = 1:30, k.filter =30, k.score = 30,
                                        anchor.features = 2000,
                                        reduction = "rpca",scale = T)
  save(anchors_rpca,file=paste0(folder_base_output,'anchors_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj'))
  anchors_rpca = loadRData(paste0(folder_base_output,'anchors_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj'))
  data_integrate_subset_rpca = IntegrateData(anchorset = anchors_rpca,
                                      dims = 1:30,k.weight  = 100)
  
  
  path = paste0(folder_base_output,'data_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj')
  save(data_integrate_subset_rpca,file=path)
  browser()

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


data_merge_1 = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
data_merge_2 = merge(x =  data_list[[3]],y = data_list[[4]], merge.data = FALSE)

data_merge_1 = NormalizeData(data_merge_1, normalization.method = "LogNormalize", scale.factor = 10000)
#data_merge_1 = FindVariableFeatures(data_merge_1, verbose = FALSE, selection.method = "vst", nfeatures =2000)

data_merge_2 = NormalizeData(data_merge_2, normalization.method = "LogNormalize", scale.factor = 10000)
#data_merge_2 = FindVariableFeatures(data_merge_2, verbose = FALSE, selection.method = "vst", nfeatures =2000)


data_list = vector(mode = "list", length = 2)
data_list[[1]] = data_merge_1
data_list[[2]] = data_merge_2

features_list = read_excel('/home/sujwary/Desktop/scRNA/Data/PaperData/TcellVarGenes.xlsx')
features_list = features_list$`Supplementary Data 2: List of highly variable genes used for T cell clustering`

features_list = features_list[features_list %in% data_merge_1@assays[["RNA"]]@counts@Dimnames[[1]]]
features_list = features_list[features_list %in% data_merge_2@assays[["RNA"]]@counts@Dimnames[[1]]]

anchors = FindIntegrationAnchors(object.list = data_list, 
                                 dims = 1:30, k.filter =30, k.score = 30,
                                 anchor.features = features_list)
data_integrate_subset_output = IntegrateData(anchorset = anchors, 
                               dims = 1:30,k.weight  = 100)
path_integrate = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,'','.Robj')
save(data_integrate_subset_output,file=path_integrate)

(unlist(data_merge_1@raw.data@Dimnames[1]))

