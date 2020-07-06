# merge the samples

# select highly variable genes

# run PCA and then 

# run BKNN instead of findneighbors. 
library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
library(sc)
library(scater)
library(dplyr)
library(scran)
library(reshape2)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

cellNames = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
cellNames  =cellNames$x


patient_list = c(12, 16, 20)

i = 1

for (i in 1:length(patient_list)){
  patient = patient_list[i]

  sample_name_baseline = metaData$'Sample'[(metaData['Patient Number'] == patient) &  (metaData['Treatment'] == 'baseline')]

  print(sample_name_baseline)
  folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name_baseline , '/')
  data_baseline = loadRData(paste0(folder_input,sample_name_baseline,'_data.Robj'))
  data_baseline$sample = sample_name_baseline
  
  sample_name_C9D1 = metaData$'Sample'[(metaData['Patient Number'] == patient) &  (metaData['Treatment'] == 'C9D1')]
  
  print(sample_name_C9D1)
  folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name_C9D1 , '/')
  data_C9D1 = loadRData(paste0(folder_input,sample_name_C9D1,'_data.Robj'))
  data_C9D1$sample = sample_name_C9D1
  #data_list = c(data_baseline,data_C9D1)
  
  # merge the samples
  data_merge = merge(x =  data_baseline,y = data_C9D1, merge.data = FALSE)
  
  data_baseline <- NormalizeData(data_baseline, verbose = FALSE)
  data_baseline <- FindVariableFeatures(data_baseline, selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  
  data_C9D1 <- NormalizeData(data_C9D1, verbose = FALSE)
  data_C9D1 <- FindVariableFeatures(data_C9D1, selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
  data_anchors <- FindIntegrationAnchors(object.list = c(data_baseline,data_C9D1), dims = 1:30)
  data_integrate <- IntegrateData(anchorset = data_anchors, dims = 1:30)
  
  # select highly variable genes
  data_merge <- FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
  
  # scale?
  data_merge = ScaleData(data_merge)
  # run PCA 
  data_merge = RunPCA(data_merge,npcs = 30)
  # run BBKNN instead of findneighbors. 
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/BBKNN/patient',patient,'/')
  dir.create(folder, recursive =  T)
  path = paste0(folder,'/','Patient',patient,'_input.Robj')
  save(data_merge,file= path)
  
  data_matrix = as.matrix(data_merge@reductions[["pca"]]@cell.embeddings)
  write(colnames(data_matrix), file = paste0(folder,'Patient',patient,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0(folder,'Patient',patient,'_rownames.txt'))
  write(data_matrix, file = paste0(folder,'Patient',patient,'_matrix.csv',sep = ","))
  
  write.csv(data_matrix, file = paste0(folder,'Patient',patient,'_matrix.csv'))
  
  write(data_merge$sample, file = paste0(folder,'Patient',patient,'_sample.txt'))
  
  
  
}

