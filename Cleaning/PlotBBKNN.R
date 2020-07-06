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
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

patient_list = c(12, 16, 20)

i = 1

for (i in 1:length(patient_list)){
  patient = patient_list[i]
  
}

################3
folder = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/baseline','/')
data_orig = loadRData(paste0(folder, 'baseline_input.Robj'))
nrow(data_orig)
ncol(data_orig)
connect = read.csv(file = paste0(folder,'connectivities.csv'),header = F)
path = paste0(folder,'baseline_matrix.csv')
pca_matrix = read.csv(file = path,header = T,row.names = 1)
colnames(connect) = rownames(pca_matrix)
rownames(connect) = rownames(pca_matrix)

dis = read.csv(file = paste0(folder,'distances.csv'),header = F)
colnames(dis) = colnames(data_orig)
rownames(dis) = colnames(data_orig)

data_orig_run = data_orig
nrow(connect)
ncol(connect)
#data_orig_run@graphs[["RNA_nn"]] = as.matrix(connect)
data_orig_run@graphs[["RNA_snn"]] = as.matrix(connect)
#data_orig@graphs[["RNA_snn"]] = as.matrix(dis)
names(data_orig_run) = names(data_orig_tmp)

#data_orig_run <- FindVariableFeatures(data_orig, selection.method = "vst", nfeatures = 2000)
#data_orig_run <- ScaleData(data_orig_run)
#data_orig_run <- RunPCA(data_orig_run, features = VariableFeatures(object = data_orig_run),npcs = 30)
data_orig_run <- FindNeighbors(data_orig, dims = 1:30)
data_orig_run = FindClusters(data_orig_run,resolution = 1.4)
data_orig_run = RunUMAP(data_orig_run, dims = 1:30)


data_orig_tmp <- FindNeighbors(data_orig, dims = 1:30)

data_orig_tmp = FindClusters(data_orig_tmp,resolution = 1.4)
data_orig_tmp = RunUMAP(data_orig_tmp, dims = 1:30)




#data_orig_tmp <- FindVariableFeatures(data_orig, selection.method = "vst", nfeatures = 2000)
#data_orig_tmp <- ScaleData(data_orig_tmp)
#data_orig_tmp <- RunPCA(data_orig_tmp, features = VariableFeatures(object = data_orig_tmp),npcs = 30)

nrow(data_orig_tmp@graphs[["RNA_nn"]])
data_orig_tmp@graphs[["RNA_snn"]]
