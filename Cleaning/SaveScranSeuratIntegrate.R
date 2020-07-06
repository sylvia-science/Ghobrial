
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


filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

cellNames = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
cellNames  =cellNames$x

#folder = 'Intra-v2'
#folder = 'Inter-version'
folder = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")

sample_list = metaData$Sample

folder_output_merge = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Merge/',folder,'/','SNN_Umap','/')
folder_output_integrate = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Integrate/',folder,'/','PCA_Umap','/')

#folder_output_merge = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Merge/PBMC_10X/','SNN_Umap','/')
#folder_output_integrate = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Integrate/PBMC_10X/','SNN_Umap','/')



SaveScranSeuratMergeInt_Helper(sample_list, mode = 'Merge',
                               folder_output_merge,folder_output_integrate, 
                               downsample = cellNames)  




fc_BCell = colorRampPalette(c("orange", "darkorange"))
fc_CD14 = colorRampPalette(c("purple", "purple4"))
fc_CD16 = colorRampPalette(c("lightgreen", "darkgreen"))
fc_CD8 = colorRampPalette(c("goldenrod1", "goldenrod4"))
fc_DC = colorRampPalette(c("yellow", "yellow4"))
fc_Eryth = colorRampPalette(c("red", "red4"))
fc_HSC = colorRampPalette(c("lightcyan", "darkcyan"))
fc_MK = 'darkseagreen'
fc_neutro = colorRampPalette(c("burlywood", "burlywood4"))
fc_NK = colorRampPalette(c("springgreen", "springgreen4"))
fc_pDC = colorRampPalette(c("lightpink", "pink4"))
fc_Plasma = 'plum'
fc_PreB = colorRampPalette(c("seashell", "seashell4"))
fc_T = colorRampPalette(c("lightblue", "darkblue"))




color_list = c(fc_BCell(9),fc_CD14(13),fc_CD16(5),
               fc_CD8(12),fc_DC(6),fc_Eryth(6),
               fc_HSC(10),fc_MK,fc_neutro(3),fc_NK(12),
               fc_pDC(2), fc_Plasma,fc_PreB(2),fc_T(13))


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




resolution_val = 1.4
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
groupBy_list = c('sample','Diagnosis','kit','Treatment','Batch','LowCount','Doublet','GeneralCellType')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
splitBy_list = NA

plotAll(data_merge_run, folder = folder_output_merge,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = resolution_val)


PlotKnownMarkers(data_merge_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '', plotTogether = F)


tmp = sort(unique(data_merge_run$CellType))
tmp = str_replace(tmp, 'HSc', 'HSC')
sum(grepl("B Cell",tmp))


level_list = sort(unique(data_merge_run$CellType))
data_merge_run$CellType = factor( data_merge_run$CellType , levels = level_list)
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

filepath_cluster = paste0( folder_output_merge, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
folder_output = paste0(filepath_cluster,'Entropy/')
dir.create(folder_output)
folder_output = paste0(folder_output,'Merge')
k_num = 30
dim_num = 30
#data_merge_run_downsample= subset(data_merge_run, cells = sample(Cells(data_merge_run), 40000))
#cellNames = colnames(data_merge_run_downsample)
#write.csv(cellNames,'/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')

#data_merge_run_downsample = data_merge_run[,cellNames]
compute_entropy_Seurat(data_merge_run, 
                       corrected_assay = 'pca',folder = folder_output,
                       k_num = k_num, dim_num = dim_num)

file_entropy = paste0(folder_output,'_k',k_num ,'_entropy.csv')
entropy = read.csv(file_entropy,sep = ',')


entropy$Method = 'Merge'

plotEntropy(entropy,folder_output)



###############################################
# Integrate

resolution_val = 1.4
filepath_cluster = paste0( folder_output_integrate, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
path = paste0(filepath_cluster,'','data','.Robj')
data_integrate_run = loadRData(path)

data_integrate_run = addMetaData(data_integrate_run, metaData)
data_integrate_run = load_emptyDrops(data_integrate_run)
data_integrate_run = load_Doublets(data_integrate_run)
data_integrate_run = load_CellLabel(data_integrate_run)

data_integrate_run$split_var = ''
data_integrate_run$kit = data_integrate_run$`10X kit`

data_integrate_run_label = label_cells(data_integrate_run,cluster_IDs)



resolution_val = 1.4
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
groupBy_list = c('sample','Diagnosis','kit','Treatment','Batch','LowCount','Doublet','GeneralCellType')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')

#groupBy_list = c('sample','Diagnosis','Treatment','Batch')
#groupBy_list = c('sample')
splitBy_list = NA
plotAll(data_integrate_run, folder = folder_output_integrate,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list, featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val)


plotAll(data_integrate_run_label, folder = folder_output_integrate,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list, featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val, str = '_label')

PlotKnownMarkers(data_integrate_run, 
                 folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',plotTogether = F)



group = 'CellType'
pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=5000, height=1500)
plot = DimPlot(data_integrate_run,pt.size = 1, reduction = "umap",label = FALSE,group.by  = group)

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

## 
filepath_cluster = paste0( folder_output_integrate, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
folder_output = paste0(filepath_cluster,'Entropy/')
dir.create(folder_output)
folder_output = paste0(folder_output,'Integrate')

k_num = 30
dim_num = 30
#cellNames_int = vapply(strsplit(cellNames, '_'), function(x)paste(x[seq.int(2)], collapse='_'), character(''))

data_integrate_run_label_downsample = data_integrate_run_label[,colnames(data_merge_run) %in% cellNames]

compute_entropy_Seurat(data_integrate_run, 
                       corrected_assay = 'pca',folder = folder_output,
                       k_num = k_num, dim_num = dim_num)

file_entropy = paste0(folder_output,'_k',k_num,'_' ,'entropy.csv')
entropy = read.csv(file_entropy,sep = ',')

entropy$Method = 'Integrate'

entropy$batch_entropy <- NULL
plotEntropy(entropy,folder_output,k_num)
 
 
