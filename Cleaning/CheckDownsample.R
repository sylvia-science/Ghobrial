
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
library(stringr)

source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')

source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PipelineIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PlotFunctionIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/IntegrateAll_ClusterUmap.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')

source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/TestIntegration/SaveScranSeuratMergeInt_Helper.R')
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

sample_type = 'Integrate_AllSamples'
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 

folder = 'AllSamples'

folder_output_merge = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Merge/',folder,'/','SNN_Umap','/')

resolution_val = 1.4
filepath_cluster = paste0( folder_output_merge, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
path = paste0(filepath_cluster,'','data','.Robj')
data_merge_run = loadRData(path)
data_merge_run$CellType = ''
data_merge_run = addMetaData(data_merge_run, metaData)
data_merge_run = load_emptyDrops(data_merge_run)
data_merge_run = load_Doublets(data_merge_run)
data_merge_run = load_CellLabel(data_merge_run)
data_merge_run$GeneralCellType = str_match(data_merge_run$CellType, "(^.+)\\s")[, 2]
data_merge_run$kit = data_merge_run$`10X kit`

folder = 'AllSamplesDownsample'

folder_output_merge = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Merge/',folder,'/','SNN_Umap','/')

resolution_val = 1.4
filepath_cluster = paste0( folder_output_merge, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
path = paste0(filepath_cluster,'','data','.Robj')
data_merge_run_DS = loadRData(path)
data_merge_run_DS$CellType = ''
data_merge_run_DS = addMetaData(data_merge_run_DS, metaData)
data_merge_run_DS = load_emptyDrops(data_merge_run_DS)
data_merge_run_DS = load_Doublets(data_merge_run_DS)
data_merge_run_DS = load_CellLabel(data_merge_run_DS)
data_merge_run_DS$kit = data_merge_run_DS$`10X kit`


data = data.frame(as.character(data_merge_run$sample))
colnames(data) = 'sample'
data$kit =as.character(data_merge_run$kit)
data$CellType = as.character(data_merge_run$GeneralCellType)
data$DataType = 'Full Dataset'


data_DS_sample = data.frame(as.character(data_merge_run_DS$sample))
colnames(data_DS_sample) = 'sample'
data_DS_sample$kit =as.character(data_merge_run_DS$kit)
data_DS_sample$CellType =as.character(data_merge_run_DS$GeneralCellType)
data_DS_sample$DataType = 'Downsample'

data = rbind(data,data_DS_sample)


100*table(data_merge_run$sample)/ncol(data_merge_run)
100*table(data_merge_run_DS$sample)/ncol(data_merge_run_DS)

100*table(data_merge_run$kit)/ncol(data_merge_run)
100*table(data_merge_run_DS$kit)/ncol(data_merge_run_DS)

100*table(data_merge_run$GeneralCellType)/ncol(data_merge_run)
100*table(data_merge_run_DS$GeneralCellType)/ncol(data_merge_run_DS)



library(dplyr)
data_sample <- data %>% 
  group_by(DataType,sample) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

data_kit <- data %>% 
  group_by(DataType,kit) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

data_CellType <- data %>% 
  group_by(DataType,CellType) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

ggplot(data_sample, aes(x = factor(DataType), y = perc*100, fill = factor(sample))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "DataType", y = "percent", fill = "sample") +
  theme_minimal(base_size = 14)


ggplot(data_kit, aes(x = factor(DataType), y = perc*100, fill = factor(kit))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "DataType", y = "percent", fill = "kit") +
  theme_minimal(base_size = 14)


ggplot(data_CellType, aes(x = factor(DataType), y = perc*100, fill = factor(CellType))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "DataType", y = "percent", fill = "CellType") +
  theme_minimal(base_size = 14)
