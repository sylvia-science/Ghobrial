# merge the samples

# select highly variable genes

# run PCA and then 

# run BKNN instead of findneighbors. 
library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(harmony)
library(ggplot2)
library(SoupX)
library(sc)
library(scater)
library(dplyr)
library(scran)
library(reshape2)
library(stringr)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')

filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  = NA #downsample$x


sample_type = 'Harmony_AllSamples_Sample_Kit'
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 


patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

#folder = 'Intra-v3_1'
#folder = 'Inter-version'
folder_name = 'AllSamples'
#folder_name = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample

folder = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                '/Batch_Sample_Kit/','/')
dir.create(folder,recursive = T)


run = F

if (run){
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
    data_i = data_i[,data_i$nFeature_RNA > 200]
    
    if (!is.na(downsample)){
      downsample = sub("_.*", "", downsample)
      cellnames = colnames(data_i)
      cellnames = sub("_.*", "", cellnames)
      data_i = data_i[,cellnames %in% downsample]
      #browser()
    }
    
    print(ncol(data_i))
    if (ncol(data_i) > 100){
    #if (T){
      
      data_list[[i]] = data_i
      data_list_norm[[i]] = ScranNorm(data_i)
    }
    
  }
  
  data_list_norm =data_list_norm[lengths(data_list_norm) != 0]
  data_merge = merge(x =  data_list_norm[[1]] ,y = data_list_norm[2:length(data_list_norm)], merge.data = T)
  
  #data_merge = merge(x =  data_list_norm[[1]],y = data_list_norm[[2]], merge.data = T)
  #for (i in 3:length(data_list_norm)){
  #  print(i)
  #  data_merge = merge(x =  data_merge,y = data_list_norm[[i]],merge.data = T)
  #}
  
  
  data_merge = addMetaData(data_merge, metaData)
  data_merge = load_emptyDrops(data_merge)
  data_merge = load_Doublets(data_merge)
  data_merge = load_CellLabel(data_merge)
  data_merge$GeneralCellType = str_match(data_merge$CellType, "(^.+)\\s")[, 2]
  data_merge$kit = data_merge$`10X kit`
  data_merge$split_var = ''
  
  data_merge_run = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
  data_merge_run = ScaleData(data_merge_run)
  data_merge_run = RunPCA(data_merge_run,npcs = 40)

  data_merge_run = RunHarmony(data_merge_run,group.by.vars =  c("sample", "10X kit"))
  
  data_merge_run = RunUMAP(reduction = "harmony",data_merge_run, dims = 1:30)
  data_merge_run = FindNeighbors(data_merge_run, reduction = "harmony", dims = 1:30)
  data_merge_run = FindClusters(data_merge_run,resolution = resolution_val)
  
  data_merge_run = addMetaData(data_merge_run, metaData)
  data_merge_run = load_emptyDrops(data_merge_run)
  
  path = paste0(folder,'data_run','.Robj')
  save(data_merge_run,file= path)

}else{
  path = paste0(folder,'data_run','.Robj')
  data_merge_run = loadRData(path)
  
}


#metaData = read_excel(filename_metaData)
data_merge_run = addMetaData(data_merge_run, metaData)
#data_merge_run = load_emptyDrops(data_merge_run)
#data_merge_run = load_Doublets(data_merge_run)
data_merge_run = load_CellLabel(data_merge_run)
print(unique(data_merge_run$GeneralCellType))
#data_merge_run$GeneralCellType = str_match(data_merge_run$CellType, "(^.+)\\s")[, 2]
#data_merge_run$kit = data_merge_run$`10X kit`
#data_merge_run$split_var = ''
#data_merge_run_label = label_cells(data_merge_run,cluster_IDs)


resolution_val = 1.4
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')

groupBy_list = c('sample','Diagnosis','kit',
                 'Treatment','Batch','LowCount',
                 'Doublet','GeneralCellType')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
splitBy_list = NA

plotAll(data_merge_run, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val)


plotAll(data_merge_run_label, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = resolution_val,str = '_label')

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
PlotKnownMarkers(data_merge_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',plotTogether = F)
print(DimPlot(data_merge_run, reduction = "umap", label = TRUE, pt.size = .1))



group = 'CellType'
pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=5000, height=1500)
plot = DimPlot(data_merge_run,pt.size = 1, reduction = "umap",label = FALSE,group.by  = group)

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
)


#plot = plot + scale_color_manual(values=color_list)

print(plot)

dev.off()


filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
k_num = 30
dim_num = 30
folder_output = paste0(filepath_cluster,'Entropy/')
dir.create(folder_output)
folder_output = paste0(folder_output,'Harmony')


#data_merge_run_label_downsample = data_merge_run_label[,cellNames]

compute_entropy_Seurat(data_merge_run, 
                       corrected_assay = 'harmony',folder = folder_output,
                       k_num = k_num, dim_num = dim_num)


file_entropy = paste0(folder_output,'_k',k_num ,'_entropy.csv')
entropy = read.csv(file_entropy,sep = ',')

entropy$Method = 'Harmony'

entropy$batch_entropy <- NULL
plotEntropy(entropy,folder_output)

