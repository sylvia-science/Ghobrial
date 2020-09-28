library(BiocSingular)
library(Seurat)
library(DoubletFinder)

library(scran)
library(scRNAseq)
library(singleCellTK)
library(dplyr)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)

source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
#metaData = metaData[metaData$Run== 1,]
metaData = metaData[metaData$`Sample Type` == 'PBMC',]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

folder = 'Soup_MT_scran_Doublet'

#pbmc.data <- Read10X(data.dir = '/home/sujwary/Downloads/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/')
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
i = 4
for (i in 1:nrow(metaData) ){
  
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder_output = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  path = paste0(folder_output,'/',sample_name,'.Robj')
  data_i_run = loadRData(path)
  
  #folder = paste0('/disk2/Projects/EloRD/Output/EmptyCells/',sample_name,'/')
  #br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv'))
  
  dbl.dens_cells = doubletCells(data_i_run@assays[["RNA"]]@counts)
  doublet_quantile = quantile(dbl.dens_cells,0.925) 
  scran_doublet = dbl.dens_cells >= doublet_quantile
  
  doublet_log =  log10(dbl.dens_cells+1)
  doublet_quantile_log =  quantile(doublet_log,0.925) 
  hist(doublet_log,breaks=20)
  
  scran_doublet = doublet_log > doublet_quantile_log
  folder_output = paste0('/disk2/Projects/EloRD/Output/',folder,'/',sample_name,'/')
  dir.create(folder_output, recursive = T)
  write.csv(scran_doublet, file=paste0(folder_output, 'Doublet','.csv'))
  #dbl.dens = doubletCluster(data_i_run@assays[["RNA"]]@counts, clusters= Idents(data_i_run))
  #data_i_run$doubletScore = log10(dbl.dens_cells+1)
  data_i_run$scran_doublet = scran_doublet
  print(100*sum(scran_doublet)/length(scran_doublet))
  #next
  
  dir.create(folder_output, recursive = T)
  
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  print('Plot')
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','scran_doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE,group.by ='scran_doublet' ))
  dev.off()
  
  
  # pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','doubletScore','.png')
  # png(file=pathName,width=1000, height=1000)
  # print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("doubletScore")))
  # dev.off()
  
  
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder_output = paste0('/disk2/Projects/EloRD/Output/',folder,'/',sample_name,'/','Featureplots/Pre_Doublet_Umap/')
    dir.create(folder_output,recursive = T)
    plot = FeaturePlotFix(data_i_run, feature = gene, folder =folder_output,
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
    pathName = paste0(folder,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  
  
  next
  # Remove doubelts
  data_i_filtered_run = data_i_run
  data_i_filtered_run = data_i_filtered_run[,data_i_filtered_run@meta.data[var] == 'Singlet']
  
  data_i_filtered_run = NormalizeData(data_i_filtered_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/','scran_Doublet/',sample_name,'/')
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
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/scran_Doublet/',sample_name,'/','Featureplots/Post_Doublet_Umap/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_run, feature = gene, folder =folder,
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
    pathName = paste0(folder,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
}
