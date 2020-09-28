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
i = 4

#######################################
# Soup + MT
for (i in 3:nrow(metaData) ){
#for (i in 1 ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  

  folder = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'.Robj')
  data_i_run = loadRData(path)

  
  folder_name = 'Soup_MT_DoubletFinder'
  
  # sweep.res.list = paramSweep_v3(data_i_run, PCs = 1:10, sct = FALSE)
  # sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  # pK = find.pK(sweep.stats)
  # 
  # folder = paste0('/disk2/Projects/EloRD/Output/',folder_name,'/',sample_name,'/')
  # dir.create(folder,recursive = T)
  # 
  # pathName = paste0(folder,sample_name,'_Sweep_PK','','.png')
  # png(file=pathName,width=800, height=500)
  # plot = ggplot(data=pK, aes(x=pK, y=BCmetric, group = 1)) +
  #    geom_line()+
  #    geom_point()
  #  print(plot)
  #  dev.off()
  #next

  #pk_val = pK
  pk_val = sampleParam$pk[sampleParam['Sample'] == sample_name]
  print(pk_val)
  nExp_poi <- round(0.075*length(colnames(data_i_run)))
  pn = 0.25
  #pK = 0.09 #Optimal pK values can be determined using mean-variance-normalized bimodality coefficient.
  data_i_run = doubletFinder_v3(data_i_run, PCs = 1:30, pN = pn, pK = pk_val, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  var = paste0('DF.classifications_',pn,'_',pk_val,'_',nExp_poi)
  
  
  folder = paste0('/disk2/Projects/EloRD/Output/',folder_name,'/',sample_name,'/','pk_',pk_val,'/')
  dir.create(folder, recursive = T)
  doublet_data =data_i_run@meta.data[var]
  doublet_data = ifelse(doublet_data=='Doublet', T, F)
  write.csv(doublet_data, file=paste0(folder, 'Doublet',pk_val,'.csv'))
  
  #next
  
  
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','Doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = var))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  folder_feature = paste0(folder,'Featureplots/Pre_Doublet_Umap/' )
  dir.create(folder_feature,recursive = T)
  
  gene_list = c('CD3D', 'CD3G', 'CD3E','CD4', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
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

##########################################
# Soup + Doublet
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv'))
  
  data_soup$emptyProb = br_e_sorted$LogProb
  data_soup$is_cell = br_e_sorted$is_cell

  data_i_run = data_soup

  #next
  nExp_poi <- round(0.075*length(colnames(data_i_run)))
  pn = 0.25
  pK = 0.09
  data_i_run = doubletFinder_v3(data_i_run, PCs = 1:10, pN = pn, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  var = paste0('DF.classifications_',pn,'_',pK,'_',nExp_poi)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/','Soup_DoubletFinder/',sample_name,'/')
  dir.create(folder, recursive = T)
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','Doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(   DimPlot(data_i_run, pt.size = pt.size,reduction = "umap",label = FALSE, group.by = var))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_DoubletFinder/',sample_name,'/','Featureplots/Pre_Doublet_Umap/')
    dir.create(folder,recursive = T)
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
    pathName = paste0(folder,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  

  
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
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/','Soup_DoubletFinder/',sample_name,'/')
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
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_DoubletFinder/',sample_name,'/','Featureplots/Post_Doublet_Umap/')
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

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
gt.calls <- pbmc@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(colnames(pbmc)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

DimPlot(pbmc, group.by = 'DF.classifications_0.25_0.09_202')
