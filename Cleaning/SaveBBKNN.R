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
library(stringr)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')

filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  =downsample$x


patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)


#folder = 'Intra-v3_1'
#folder = 'Inter-version'
folder = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample


data_list = vector(mode = "list",length = length(sample_list))
data_list_norm = vector(mode = "list",length = length(sample_list))
#for (i in 1:nrow(sampleParam)){
for (i in 1:length(sample_list)){
  #sample_name = sampleParam$Sample[i]
  sample_name = sample_list[i]
  print(sample_name)
  folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name , '/')
  data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
  data_i$sample = sample_name
  path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
  cellIdents = read.csv(path,sep = ',',row.names = 1)
  cellIdents$x = paste0(cellIdents$x, ' S', i)
  data_i$CellType = cellIdents
  
  if (!is.na(downsample)){
    downsample = sub("_.*", "", downsample)
    cellnames = colnames(data_i)
    cellnames = sub("_.*", "", cellnames)
    data_i = data_i[,cellnames %in% downsample]
    #browser()
  }
  
  if (ncol(data_i) > 100){
    data_list[[i]] = data_i
    data_list_norm[[i]] = ScranNorm(data_i)
  }
}
data_list_norm =data_list_norm[lengths(data_list_norm) != 0]

data_merge = merge(x = data_list_norm[[1]] ,y = data_list_norm[2:length(data_list_norm)], merge.data = T)
  
resolution_val = 1.4

  
#data_merge_run = ScranNorm(data_merge)
data_merge_run = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
data_merge_run = ScaleData(data_merge_run)
data_merge_run = RunPCA(data_merge_run, features = VariableFeatures(object = data_merge_run),npcs = 30)
  
  
folder = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/',folder,'/')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)

data_merge_run = addMetaData(data_merge_run, metaData)  
dir.create(folder, recursive =  T)
path = paste0(folder,'/','data','.Robj')
save(data_merge_run,file= path)

path = paste0(folder,'/','data','.Robj')
data_merge_run = loadRData(path)


data_merge_run = addMetaData(data_merge_run, metaData)
data_merge_run = load_CellLabel(data_merge_run)
counts_matrix = data_merge@assays[["RNA"]]@counts
write.csv(counts_matrix,file = paste0(folder,'MergeCounts','.csv'))
  
data_matrix = as.matrix(data_merge_run@reductions[["pca"]]@cell.embeddings)
write(colnames(data_matrix), file = paste0(folder,'colnames.txt'))
write(rownames(data_matrix), file = paste0(folder,'rownames.txt'))
  
write.csv(data_matrix, file = paste0(folder,'PCA.csv'),sep = ",")
write(data_merge_run$sample, file = paste0(folder,'sample.txt'))
write(data_merge_run$`10X kit`, file = paste0(folder,'kit.txt'))
write(data_merge_run$Batch, file = paste0(folder,'batch.txt'))
write.csv(data_merge_run$GeneralCellType, file = paste0(folder,'GeneralCellType.csv'))

  
  
  



######################################################
# Compare baseline between 12 and 20

patient1 = 20
patient2 = 30

sample_name_baseline1 = metaData$'Sample'[(metaData['Patient Number'] == patient1) &  (metaData['Treatment'] == 'baseline')]
print(sample_name_baseline1)
folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name_baseline1 , '/')
data1 = loadRData(paste0(folder_input,sample_name_baseline1,'_data.Robj'))
data1$sample = sample_name_baseline1
path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name_baseline1,'/cellIdents.csv')
cellIdents1 = read.csv(path,sep = ',',row.names = 1)
cellIdents1$x = paste0(cellIdents1$x, ' S1')

sample_name_baseline2 = metaData$'Sample'[(metaData['Patient Number'] == patient2) &  (metaData['Treatment'] == 'baseline')]
print(sample_name_baseline2)
folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name_baseline2 , '/')
data2 = loadRData(paste0(folder_input,sample_name_baseline2,'_data.Robj'))
data2$sample = sample_name_baseline2
#data_list = c(data_baseline,data_C9D1)
path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name_baseline2,'/cellIdents.csv')
cellIdents2 = read.csv(path,sep = ',',row.names = 1)
cellIdents2$x = paste0(cellIdents2$x, ' S2')


# merge the samples
data_merge = merge(x =  data1,y = data2, merge.data = FALSE)
data_merge$CellType = ''
cellnames_baseline1 = colnames(data1)
cellnames_baseline2 = colnames(data2)
data_merge$CellType[colnames(data_merge) %in% cellnames_baseline1] = cellIdents1$x
data_merge$CellType[colnames(data_merge) %in% cellnames_baseline2] = cellIdents2$x

data_merge_run = ScranNorm(data_merge)
data_merge_run = FindVariableFeatures(data_merge_run, selection.method = "vst", nfeatures = 2000)
data_merge_run = ScaleData(data_merge_run)
data_merge_run = RunPCA(data_merge_run, features = VariableFeatures(object = data_merge_run),npcs = 30)


folder = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/baseline/',
                'Patients',patient1,'_',patient2,'/SNN_Umap/')

dir.create(folder, recursive =  T)
path = paste0(folder,'/','','baseline','_input.Robj')
save(data_merge_run,file= path)

counts_matrix = data_merge@assays[["RNA"]]@counts
write.csv(counts_matrix,file = paste0(folder,'MergeCounts','.csv'))

data_matrix = as.matrix(data_merge_run@reductions[["pca"]]@cell.embeddings)
write(colnames(data_matrix), file = paste0(folder,'baseline','_colnames.txt'))
write(rownames(data_matrix), file = paste0(folder,'baseline','_rownames.txt'))

write.csv(data_matrix, file = paste0(folder,'baseline','_PCA.csv'),sep = ",")
write(data_merge_run$sample, file = paste0(folder,'baseline','_sample.txt'))

 