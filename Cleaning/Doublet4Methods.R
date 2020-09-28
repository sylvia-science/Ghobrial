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
#metaData = metaData[metaData$Run== 1,]
metaData = metaData[metaData$`Sample Type` == 'PBMC',]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

#pbmc.data <- Read10X(data.dir = '/home/sujwary/Downloads/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/')
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
i = 1

# Soup + MT + Doublet
for (i in 3:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  #folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  #br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv')) # empty cells
  
  threshold=  sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  
  path = paste0('/disk2/Projects/EloRD/Output/Soup_MT_Scrublet/Umap Plots/',sample_name,'/', threshold,'/',sample_name,'_scrublet',threshold,'.csv')
  scrb = read.csv(path,header = T)
  scrb$predicted_doublet = ifelse(scrb$predicted_doublet=='True', T, F)
  scrb$Scrublet_Boolean = ifelse(scrb$Scrublet_Boolean=='True', T, F)
  
  
  pk_val = sampleParam$pk[sampleParam['Sample'] == sample_name]
  folder = paste0('/disk2/Projects/EloRD/Output/','Soup_MT_DoubletFinder','/',sample_name,'/','pk_',pk_val,'/')
  file=paste0(folder, 'Doublet',pk_val,'.csv')
  doublet_finder = read.csv(file)
  
  folder = paste0('/disk2/Projects/EloRD/Output/','Soup_MT_scran_Doublet/',sample_name,'/')
  scran_doublet = read.csv(paste0(folder, 'Doublet','.csv'))
  
  folder = paste0('/disk2/Projects/EloRD/Output/','Soup_MT_SCDS/',sample_name,'/')
  scds_doublet = read.csv(paste0(folder, 'Doublet','.csv'))
  
  
  doublet_summary = data.frame(matrix(ncol = 6, nrow = nrow((scds_doublet))))
  colnames(doublet_summary) = c('Cell',"scrublet","doublet_finder",
                                'scran_doublet','scds_doublet','Summary')
  doublet_summary$Cell = scrb$Cell
  doublet_summary$scrublet = scrb$Scrublet_Boolean
  doublet_summary$doublet_finder = doublet_finder[,2]
  doublet_summary$scran_doublet = scran_doublet[,2]
  doublet_summary$scds_doublet = scds_doublet[,2]
  doublet_summary$Summary = rowSums(doublet_summary[,2:5])
  #colSums(doublet_summary)

  #data_i_run$scrublet = scrb$Scrublet_Boolean
  #data_i_run$doublet_finder = doublet_finder[,2]
  #data_i_run$scran_doublet = scran_doublet[,2]
  #data_i_run$scds_doublet = scds_doublet[,2]

  #data_i_run$TrueDoublet = doublet_summary$Summary >= 2
  
 

  
  folder_name = 'Doublet4Methods'
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/',folder_name,'/',sample_name,'/')
  dir.create(folder, recursive = T)
  
  write.csv(doublet_summary, file = paste0(folder,'doublet_summary','.csv'),sep = ",")
  
  
  next
  pt.size  = 0.8
  print('Plot')
  
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  plot = DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE)
  plot = plot+ggtitle('Pre Doublet Umap') + 
    theme(plot.title = element_text(hjust = 0.5, size=22))
  print( plot )
  dev.off()
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap_','scrublet','.png')
  png(file=pathName,width=1000, height=1000)
  plot = DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'scrublet')
  plot = plot+ggtitle(paste0('scrublet ' , 100*sum(doublet_summary$scrublet)/length(doublet_summary$scrublet)) ) + 
    theme(plot.title = element_text(hjust = 0.5, size=22))
  print( plot  )
  dev.off()
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap_','doublet_finder','.png')
  png(file=pathName,width=1000, height=1000)
  plot = (   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'doublet_finder'))
  plot = plot+ggtitle(paste0('doubletfinder ' , 100*sum(doublet_summary$doublet_finder)/length(doublet_summary$doublet_finder)) ) + 
    theme(plot.title = element_text(hjust = 0.5, size=22))
  print( plot  )
  dev.off()
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap_','scran_doublet','.png')
  png(file=pathName,width=1000, height=1000)
  plot = (   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'scran_doublet'))
  plot = plot+ggtitle(paste0('scran doublet ' , 100*sum(data_i_run$scran_doublet)/length(data_i_run$scran_doublet)) ) + 
    theme(plot.title = element_text(hjust = 0.5, size=22))
  print( plot  )
  dev.off()
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap_','scds_doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'scds_doublet'))
  plot = plot+ggtitle(paste0('scds doublet ' , 100*sum(doublet_summary$scds_doublet)/length(doublet_summary$scds_doublet)) ) + 
    theme(plot.title = element_text(hjust = 0.5, size=22))
  print( plot  )
  dev.off()
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap_','TrueDoublet','.png')
  png(file=pathName,width=1000, height=1000)
  plot = (   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'TrueDoublet'))
  plot = plot+ggtitle(paste0('True Doublet ' , 100*sum(doublet_summary$scds_doublet)/length(doublet_summary$scds_doublet)) ) + 
    theme(plot.title = element_text(hjust = 0.5, size=22))
  print( plot  )
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  next
  
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
  data_i_filtered_run = data_i_filtered_run[,data_i_filtered_run@meta.data['TrueDoublet'] == F]
  
  data_i_filtered_run = NormalizeData(data_i_filtered_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/',folder_name,'/',sample_name,'/')
  dir.create(folder, recursive = T)
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Post_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Post_Doublet_Umap','TrueDoublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(   DimPlot(data_i_filtered_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = 'TrueDoublet'))
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
