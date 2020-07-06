library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
library(DoubletFinder )
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

#pbmc.data <- Read10X(data.dir = '/home/sujwary/Downloads/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/')
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
i = 4


# Soup Empty cells + MT + Doublet
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  threshold=  sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  path = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_Scrublet/Umap Plots/',sample_name,'/', threshold,'/',sample_name,'_scrublet',threshold,'.csv')
  scrb = read.csv(path,header = T)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv'))
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  data_i_run = data_soup
  
  
  
  data_i_run$emptyProb = br_e_sorted$LogProb
  data_i_run$is_cell = br_e_sorted$is_cell
  
  data_i_run = data_i_run[,data_i_run$is_cell]
  data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
  
  data_i_run$Scrublet_Boolean = as.logical(as.character(toupper(scrb$Scrublet_Boolean)))
  data_i_run$doublet_scores = scrb$doublet_scores
  data_i_run$predicted_doublet = as.logical(as.character(toupper(scrb$predicted_doublet)))
  
  
  data_i_run = NormalizeData(data_i_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_run = FindVariableFeatures(data_i_run, selection.method = "vst", nfeatures = 2000)
  data_i_run = ScaleData(data_i_run)
  data_i_run = RunPCA(data_i_run,npcs = 30)
  data_i_run = FindNeighbors(data_i_run, dims = 1:30)
  data_i_run = FindClusters(data_i_run)
  data_i_run = RunUMAP(data_i_run, dims = 1:30)

  folder_name = 'Soup_Empty_MT_Scrublet'
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/',folder_name,'/',sample_name,'/')
  dir.create(folder, recursive = T)
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','Doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'predicted_doublet'))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  folder_feature = paste0(folder,'Featureplots/Pre_Doublet/' )
  dir.create(folder_feature,recursive = T)
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    plot = FeaturePlotFix(data_i_run, feature = gene, folder =folder,
                          str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=24),
      legend.title=element_text(size=24),
      text = element_text(size = 20)
    )
    
    file_str = ''
    pathName = paste0(folder,gene,file_str,'.png')
    pathName = paste0(folder_feature,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  
  # Remove doubelts
  data_i_filtered_run = data_i_run
  data_i_filtered_run = data_i_filtered_run[,data_i_filtered_run@meta.data['predicted_doublet'] == F]
  
  data_i_filtered_run = NormalizeData(data_i_filtered_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/',folder_name,'/',sample_name,'/', 'pk_',pk_val,'/')
  dir.create(folder, recursive = T)
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Post_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Post_Doublet_Umap','Doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(   DimPlot(data_i_filtered_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = var))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E','CD4', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  folder_feature = paste0(folder,'Featureplots/Post_Doublet_Umap/' )
  dir.create(folder_feature,recursive = T)
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    plot = FeaturePlotFix(data_i_filtered_run, feature = gene, folder =folder_feature,
                          str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=24),
      legend.title=element_text(size=24),
      text = element_text(size = 20)
    )
    
    file_str = ''
    pathName = paste0(folder,gene,file_str,'.png')
    pathName = paste0(folder_feature,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
}