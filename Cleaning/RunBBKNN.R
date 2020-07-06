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
library(vroom)
library(Rfast)
library(reshape2)
library(stringr)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

sample_type = 'BBKNN_AllSamples'
sample_type = 'BBKNN_AllSamples_Sample_Kit'
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 

cellNames = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
cellNames  =cellNames$x

patient_list = c(12, 16, 20)

i = 1
# Integrate
patient = patient_list[i]
#folder_name = 'Intra-v3_1'
#folder_name = 'Inter-version'
folder_name = 'AllSamplesDownsample'

BBKNN_batch = 'Sample'
folder = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/',folder_name,
                '/Batch_',BBKNN_batch,'/','SNN_Umap','/')

if (run){

  #file = paste0(folder,'Patient',patient,'_input.Robj')
  file = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/BBKNN/',folder_name,'/','data.Robj')
  data_orig = loadRData(file)
  
  file = paste0(folder,'Umap_coord.csv')
  Umap_coord = read.csv(file = file,header = F)
  file = paste0(folder,'cluster_label.csv')
  cluster_label = read.csv(file = file,header = F)
  
  Umap_coord = data.matrix(Umap_coord)
  rownames(Umap_coord) = colnames(data_orig)
  Umap_coord <- CreateDimReducObject(embeddings = Umap_coord,global = T, key = 'UMAP_',assay  = 'RNA')
  
  data_orig@reductions[["umap"]] = Umap_coord
  Idents(data_orig) = cluster_label$V2
  #file = paste0(folder,'Patient',patient,'_connectivities.csv')
  file = paste0(folder,'connectivities.csv')
  connect = vroom(file, delim = ",",col_names = F)
  
  
  #connect = Rfast:::data.frame.to_matrix(connect)
  #read_csv()
  
  #connect = read.csv(file = file,header = F)
  #connect = as.matrix(connect)
  #colnames(connect) = colnames(data_orig)
  #rownames(connect) = colnames(data_orig)
  
  data_orig_run = data_orig
  
  
  cell_list = Idents(data_orig_run)
  cell_list = as.numeric(levels(cell_list))[cell_list]
  Idents(data_orig_run) = factor(Idents(data_orig_run),levels = (0:max(cell_list)))
  #nrow(connect)
  #ncol(connect)
  
  #data_orig_run@graphs[["RNA_snn"]] <- as.Graph(as.matrix(connect))
  
  #data_orig_run[["RNA_snn"]] <-  as.Graph(as.matrix(connect))
  
  #resolution_val = 1.4

  #data_orig_run = FindClusters(data_orig_run,resolution = resolution_val)
  #data_orig_run = RunUMAP(data_orig_run,umap.method = 'umap-learn',graph= 'RNA_snn')
  
  data_orig_run$split_var = ''
  
  data_orig_run = addMetaData(data_orig_run, metaData)
  data_orig_run = load_emptyDrops(data_orig_run)
  data_orig_run = load_Doublets(data_orig_run)
  data_orig_run = load_CellLabel(data_orig_run)
  
  data_orig_run$split_var = ''
  data_orig_run$kit = data_orig_run$`10X kit`
  
  path = paste0(folder,'data_run','.Robj')
  save(data_orig_run,file= path)

}else{
  
  path = paste0(folder,'data_run','.Robj')
  data_orig_run = loadRData(path)
  
  metaData = read_excel(filename_metaData)

  data_orig_run = addMetaData(data_orig_run, metaData)
  #data_orig_run = load_emptyDrops(data_orig_run)
  #data_orig_run = load_Doublets(data_orig_run)
  data_orig_run = load_CellLabel(data_orig_run)
  data_orig_run$kit = data_orig_run$`10X kit`
  
  #GeneralCellType = data_orig_run$GeneralCellType
  #cellTypeSample = data_orig_run$CellType
  #write.csv(GeneralCellType,paste0(folder,'GeneralCellType.csv'))
  #write.csv(cellTypeSample,paste0(folder,'cellTypeSample.csv'))
  
  
  
  data_orig_run_label = label_cells(data_orig_run,cluster_IDs)
  
  resolution_val = 3.4
  cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')

  groupBy_list = c('sample','Diagnosis','kit',
                   'Treatment','Batch','LowCount',
                   'Doublet','GeneralCellType')
  #groupBy_list = c('sample')
  featurePlot_list = c('percent.mt','nCount_RNA')
  
  #groupBy_list = c('sample')
  splitBy_list = NA
  
  plotAll(data_orig_run, folder = folder,
          sample_name,sampleParam = NA,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = F, keepOldLabels = T, 
          groupBy = groupBy_list, splitBy = splitBy_list, featurePlot_list = featurePlot_list,
          PCA_dim = 30,resolution_val = resolution_val)
  
  
  plotAll(data_orig_run_label, folder = folder,
          sample_name,sampleParam = NA,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = F, keepOldLabels = T, 
          groupBy = groupBy_list, splitBy = splitBy_list, featurePlot_list = featurePlot_list,
          PCA_dim = 30,resolution_val = resolution_val,str = '_label')
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
  PlotKnownMarkers(data_orig_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                   plotType ='FeaturePlotFix' , str = '',plotTogether = F)
  
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
  
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
  cellList = sort(unique(data_orig_run$CellType))
  data_orig_run$CellType = factor(data_orig_run$CellType, levels = cellList)
  group = 'CellType'
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
  png(file=pathName,width=5000, height=1500)
  plot = DimPlot(data_orig_run,pt.size = 1.5, reduction = "umap",label = FALSE,group.by  = group)
  
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
  ################################################
  file = paste0(folder,'connectivities.csv')
  connect = vroom(file, delim = ",",col_names = F)
  
  resolution_val = 3.4
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
  
  k_num=5
  folder_output = paste0(filepath_cluster,'Entropy/','BBKNN','_',BBKNN_batch)
  colnames(connect) = colnames(data_orig_run)
  rownames(connect) = colnames(data_orig_run)
  data_orig_run_label_downsample = data_orig_run_label[,cellNames]
  
  
  #connect_downsample = connect[ colnames(data_orig_run_label)  %in% cellNames,colnames(data_orig_run_label)  %in% cellNames]
  compute_entropy_Seurat(data_orig_run,BBKNN = T, BBKNN_KNN = connect,
                         corrected_assay = 'pca',folder = folder_output,
                         k_num = k_num, dim_num = 30)
  
  file_entropy = paste0(folder_output,'_k',k_num,'_' ,'entropy.csv')
  entropy = read.csv(file_entropy,sep = ',')
  
  entropy$Method = 'BBKNN'
  
  entropy$batch_entropy <- NULL
  
  plotEntropy(entropy,folder_output)
  
  
  
}
