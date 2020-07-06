library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)
i = 1
## Save filtered matrices
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  
  
  print(sample_name)
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  data_i_scrublet = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_scrublet = CreateSeuratObject(counts = data_i_scrublet, project = "BM", min.cells = 3,min.features = 1)
  colSum_list = colSums(data_i_scrublet)
  keep = colSum_list >= 100
  
  data_i_scrublet = data_i_scrublet[,keep]
  
  data_matrix = data_i_scrublet@assays[["RNA"]]@counts
  
  write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/FilterMatrix_MinF100/',sample_name,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/FilterMatrix_MinF100/',sample_name,'_rownames.txt'))
  writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Data/FilterMatrix_MinF100/',sample_name,'_matrix.txt'))
  
}

i = 1
# Save filtered + soup + MT
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  data_i_run = data_soup
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Feature100/',sample_name,'/')
  #path = paste0(folder,'/',sample_name,'_Feature100.Robj')
  #data_i_run = loadRData(path)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv'))
  
  data_i_run$emptyProb = br_e_sorted$LogProb
  data_i_run$is_cell = br_e_sorted$is_cell
  
  data_i_run = data_i_run[,data_i_run$is_cell]
  data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
  
  data_i_run = NormalizeData(data_i_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_run = FindVariableFeatures(data_i_run, selection.method = "vst", nfeatures = 2000)
  data_i_run = ScaleData(data_i_run)
  data_i_run = RunPCA(data_i_run,npcs = 30)
  data_i_run = FindNeighbors(data_i_run, dims = 1:30)
  data_i_run = FindClusters(data_i_run)
  data_i_run = RunUMAP(data_i_run, dims = 1:30)
  
  data_matrix = data_i_run@assays[["RNA"]]@counts
  write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/MinF100_Soup_Empty_MT15/',sample_name,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/MinF100_Soup_Empty_MT15/',sample_name,'_rownames.txt'))
  writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Data/MinF100_Soup_Empty_MT15/',sample_name,'_matrix.txt'))
  
  
}

i = 4
# Save + soup + MT
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  data_i_run = data_soup
  
  data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
  
  data_matrix = data_i_run@assays[["RNA"]]@counts
  write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/',sample_name,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/',sample_name,'_rownames.txt'))
  writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/',sample_name,'_matrix.txt'))
  
  
  next
  data_i_run = NormalizeData(data_i_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_run = FindVariableFeatures(data_i_run, selection.method = "vst", nfeatures = 2000)
  data_i_run = ScaleData(data_i_run)
  data_i_run = RunPCA(data_i_run,npcs = 30)
  data_i_run = FindNeighbors(data_i_run, dims = 1:30)
  data_i_run = FindClusters(data_i_run)
  data_i_run = RunUMAP(data_i_run, dims = 1:30)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/')
  dir.create(folder, recursive = T)
  path = paste0(folder,'/',sample_name,'_data.Robj')
  r (i in 1:nrow(metaData) ){
    sample_name = metaData$Sample[i]
    print(sample_name)
    percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
    
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
    path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
    data_soup = loadRData(path)
    data_i_run = data_soup
    
    data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
    
    data_matrix = data_i_run@assays[["RNA"]]@counts
    write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/',sample_name,'_colnames.txt'))
    write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/',sample_name,'_rownames.txt'))
    writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/',sample_name,'_matrix.txt'))
    
    
    next
    data_i_run = NormalizeData(data_i_run, normalization.method = "LogNormalize", scale.factor = 10000)
    data_i_run = FindVariableFeatures(data_i_run, selection.method = "vst", nfeatures = 2000)
    data_i_run = ScaleData(data_i_run)
    data_i_run = RunPCA(data_i_run,npcs = 30)
    data_i_run = FindNeighbors(data_i_run, dims = 1:30)
    data_i_run = FindClusters(data_i_run)
    data_i_run = RunUMAP(data_i_run, dims = 1:30)
    
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/')
    dir.create(folder, recursive = T)
    path = paste0(folder,'/',sample_name,'_data.Robj')
    save(data_i_run,file= path)
    save(data_i_run,file= path)
}


i = 1
# Save soup + MT
for (i in 1:nrow(metaData) ){
  
  
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  data_i_run = data_soup
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv')) # empty cells
  
  data_i_run$is_cell = br_e_sorted$is_cell
  data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
  data_i_run = data_i_run[, data_i_run$is_cell]
  
  
  data_i_run = NormalizeData(data_i_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_run = FindVariableFeatures(data_i_run, selection.method = "vst", nfeatures = 2000)
  data_i_run = ScaleData(data_i_run)
  data_i_run = RunPCA(data_i_run,npcs = 30)
  data_i_run = FindNeighbors(data_i_run, dims = 1:30)
  data_i_run = FindClusters(data_i_run)
  data_i_run = RunUMAP(data_i_run, dims = 1:30)
  
  
  data_matrix = data_i_run@assays[["RNA"]]@counts
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT_empty/',sample_name,'/')
  dir.create(folder, recursive = T)
  
  write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT_empty/',sample_name,'/',sample_name,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT_empty/',sample_name,'/',sample_name,'_rownames.txt'))
  writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT_empty/',sample_name,'/',sample_name,'_matrix.txt'))
  
  
  
  
  path = paste0(folder,'/',sample_name,'_data.Robj')
  save(data_i_run,file= path)
}
