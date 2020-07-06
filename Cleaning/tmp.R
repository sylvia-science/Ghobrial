library(scRNAseq)
library(scran)
library(Seurat)
library(scater)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readxl)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')

library(Matrix)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

sample_name = metaData$Sample[i]

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')
sample_list = c('GL1497BM', 'GL1160BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')

sample_list = c('GL1497BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N','GL1502BM', 'GL2185BM', 'GL3417BM', 'GL2653BM')


#sample_list = c('GL3404BM')
i = 2
# Soup + MT
run = T

df = as.data.frame(matrix(ncol = 4, nrow = nrow(sampleParam)))
colnames(df) = c('sample','NumCell', 'aveCnt','maxCnt')
df$sample = sampleParam$Sample # Unormalized Mean per cell

for (i in 1:nrow(sampleParam) ){
  sample_name = metaData$Sample[i]
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  RNA_features_min = sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
  RNA_features_max = sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name]
  
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
  
  colSum_list = colSums(data_i_raw ) # Needs to be from Matrix library
  keep = colSum_list >= 100
  data_i_filtered = data_i_raw[,keep]
  data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")
  

  
  data_i_filtered = data_i_filtered[, data_i_filtered$percent.mt < percent_mt]
  data_i_filtered = data_i_filtered[, data_i_filtered$nFeature_RNA > RNA_features_min]
  data_i_filtered = data_i_filtered[, data_i_filtered$nFeature_RNA < RNA_features_max]
  
  
  colSums_list = colSums(data_i_filtered@assays[["RNA"]]@counts)
  print(sample_name)
  print(max(colSums_list))
  
  df$numCell[i] = ncol(data_i_filtered)
  df$maxCnt[i] = max(colSums_list)
  df$meanCnt[i] = mean(colSums_list)
  df = df[order(df$numCell),]
}


write.csv(df, file = paste0('/home/sujwary/Desktop/scRNA/Data/','CellCounts','.csv'))
