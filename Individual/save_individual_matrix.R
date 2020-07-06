library(DropletUtils)
library(Seurat)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

library(Matrix)

library(dplyr)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
library(hdf5r)
require(gridExtra)


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

sc_list = c()
for(i in 1:nrow(metaData)){
  
  sample_name = metaData$Sample[i]
  print(sample_name)

  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3,min.features = 1)
  
  
  data_i_filtered = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_filtered = CreateSeuratObject(counts = data_i_filtered, project = "BM", min.cells = 3,min.features = 100)
  
  
  
  data_matrix = data_i_raw@assays[["RNA"]]@counts

  write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/RawMatrix/',sample_name,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/RawMatrix/',sample_name,'_rownames.txt'))
  writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Data/RawMatrix/',sample_name,'_matrix.txt'))
  
}