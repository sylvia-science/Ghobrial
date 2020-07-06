library(scRNAseq)
library(Seurat)
library(dplyr)
library(ggplot2)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')


sample_name = metaData$Sample[i]

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')

i = 3
# Soup + MT
for (i in 4:length(sample_list) ){
  
  sample_name = sample_list[i]
  print(sample_name)
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_data.Robj')
  data_i_run = loadRData(path)
  
  rowMeans_list = rowMeans(data_i_run)
  
  bin_list = cut(rowMeans_list, 6, include.lowest=TRUE)
    
  # plot the total UMI count per cell on the x axis (aka colSums) and normalized expression on the y axis.
  
}


