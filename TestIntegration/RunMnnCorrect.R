# merge the samples

# select highly variable genes

# run PCA and then 

# run BKNN instead of findneighbors. 
library(Matrix)
library(readxl)
library(Seurat)
library( batchelor)
library(scran)
library(readxl)
library(ggplot2)
library(sc)
library(scater)
library(dplyr)


source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

patient_list = c(12, 16, 20)

i = 2

for (i in 2:length(patient_list)){
  patient = patient_list[i]
  
  
  sample_name_baseline = metaData$'Sample'[(metaData['Patient Number'] == patient) &  (metaData['Treatment'] == 'baseline')]
  
  print(sample_name_baseline)
  folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name_baseline , '/')
  data1 = loadRData(paste0(folder_input,sample_name_baseline,'_data.Robj'))
  data1$sample = sample_name_baseline
  path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name_baseline,'/cellIdents.csv')
  cellIdents1 = read.csv(path,sep = ',',row.names = 1)
  cellIdents1$x = paste0(cellIdents1$x, ' S1')
  
  sample_name_C9D1 = metaData$'Sample'[(metaData['Patient Number'] == patient) &  (metaData['Treatment'] == 'C9D1')]
  print(sample_name_C9D1)
  folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name_C9D1 , '/')
  data2 = loadRData(paste0(folder_input,sample_name_C9D1,'_data.Robj'))
  data2$sample = sample_name_C9D1
  #data_list = c(data_baseline,data_C9D1)
  path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name_C9D1,'/cellIdents.csv')
  cellIdents2 = read.csv(path,sep = ',',row.names = 1)
  cellIdents2$x = paste0(cellIdents2$x, ' S2')
  
  # The input expression values should generally be log-transformed, e.g., log-counts, see logNormCounts for details. 
  # They should also be normalized
  
  data1$CellType = cellIdents1$x
  data2$CellType = cellIdents2$x
  
  data1 = ScranNorm(data1)
  data2 = ScranNorm(data2)
  #data1 = ScaleData(data1)
  #data2 = ScaleData(data2)
  
  
  # merge the samples
  data_merge = merge(x =  data1,y = data2, merge.data = T)
  data_merge = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
  gene_union = intersect(data1@assays[["RNA"]]@var.features,data2@assays[["RNA"]]@var.features)
  
  gene_union = intersect(gene_union,data_merge@assays[["RNA"]]@var.features)
  mask1 =  rownames( data1@assays[["RNA"]]@scale.data) %in% gene_union
  mask2 =  rownames( data2@assays[["RNA"]]@scale.data) %in% gene_union
  tmp1 = data1@assays[["RNA"]]@scale.data[mask1,]
  tmp2 = data2@assays[["RNA"]]@scale.data[mask2,]
  data_list = list(tmp1, tmp2)
  (rownames(tmp1) %in% rownames(tmp2))
  all(rownames(tmp1) == rownames(tmp2))
  
  # Should set largest batch as first
  correct = mnnCorrect(tmp1,tmp2,
    k = 20,
    prop.k = NULL,
    sigma = 0,
    cos.norm.in = T,
    cos.norm.out = F
  )
  
  tmp =correct@assays@data@listData[["corrected"]]
  tmp = as.matrix(tmp)
  cell_NA= ((is.na(colSums(tmp))))
  
  sum(is.na(correct@assays@data@listData[["corrected"]]))
  sum(!is.na(correct@assays@data@listData[["corrected"]]))
  data_merge@assays[["RNA"]]@scale.data =  correct@assays@data@listData[["corrected"]]
  data_merge@assays[["RNA"]]@var.features = rownames(tmp1)
  
  data_merge = RemoveNan(data_merge, 'scale.data')
  data_merge@assays[["RNA"]]@data = 
    data_merge@assays[["RNA"]]@data[,colnames( data_merge@assays[["RNA"]]@data) %in% colnames(data_merge@assays[["RNA"]]@scale.data )]
  
  resolution_val = 1.4
  data_merge_run = RunPCA(data_merge,npcs = 30)
  data_merge_run = FindNeighbors(data_merge_run, dims = 1:30)
  data_merge_run = FindClusters(data_merge_run, resolution = resolution_val)
  data_merge_run = RunUMAP(data_merge_run,umap.method = 'umap-learn',graph= 'RNA_snn')
  
  # MNN correction places all cells from all batches within the same coordinate system. This means that the corrected values can be freely used to define distances between cells for dimensionality reduction or clustering. However, the correction does not preserve the mean-variance relationship. As such, we do not recommend using the corrected values for studying heterogeneity.
  
  # MNN-corrected values are generally not suitable for differential expression (DE)
}