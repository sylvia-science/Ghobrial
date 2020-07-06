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

for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
  
  colSum_list = colSums(data_i_raw ) # Needs to be from Matrix library
  keep = colSum_list >= 100
  data_i_filtered = data_i_raw[,keep]
  

  data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")
  
  
  data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  

  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/MT_Threshold/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  pt.size = 0.8
  pathName = paste0(folder,sample_name,'_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/MT_Threshold/',sample_name,'/','Featureplots/PreMT/')
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
  
  data_i_filtered_MT_run = data_i_filtered[, data_i_filtered$percent.mt < percent_mt]
  
  data_i_filtered_MT_run = data_i_filtered_MT_run
  data_i_filtered_MT_run = NormalizeData(data_i_filtered_MT_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_MT_run = FindVariableFeatures(data_i_filtered_MT_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_MT_run = ScaleData(data_i_filtered_MT_run)
  data_i_filtered_MT_run = RunPCA(data_i_filtered_MT_run,npcs = 30)
  data_i_filtered_MT_run = FindNeighbors(data_i_filtered_MT_run, dims = 1:30)
  data_i_filtered_MT_run = FindClusters(data_i_filtered_MT_run)
  data_i_filtered_MT_run = RunUMAP(data_i_filtered_MT_run, dims = 1:30)
  
  pt.size = 0.8
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/MT_Threshold/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_MT_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_MT_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_MT_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/MT_Threshold/',sample_name,'/','Featureplots/PostMT/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_MT_run, feature = gene, folder =folder,
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
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder)
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  save(data_i_filtered_run_soup,file= path)
  
  
  
}
####################################3
## Apply to soup corrected

for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)

  data_i_filtered_MT_run = data_soup[, data_soup$percent.mt < percent_mt]
  
  data_i_filtered_MT_run = data_i_filtered_MT_run
  data_i_filtered_MT_run = NormalizeData(data_i_filtered_MT_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_MT_run = FindVariableFeatures(data_i_filtered_MT_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_MT_run = ScaleData(data_i_filtered_MT_run)
  data_i_filtered_MT_run = RunPCA(data_i_filtered_MT_run,npcs = 30)
  data_i_filtered_MT_run = FindNeighbors(data_i_filtered_MT_run, dims = 1:30)
  data_i_filtered_MT_run = FindClusters(data_i_filtered_MT_run)
  data_i_filtered_MT_run = RunUMAP(data_i_filtered_MT_run, dims = 1:30)
  
  pt.size  = 0.8
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_Threshold/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_MT_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_MT_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_MT_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_Threshold/',sample_name,'/','Featureplots/PostMT/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_MT_run, feature = gene, folder =folder,
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
