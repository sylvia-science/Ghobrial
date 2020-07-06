library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)


i = 19


### Scrublet
for (i in 1:nrow(metaData) ){
  
  sample_name = metaData$Sample[i]
  print(sample_name)
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/Umap Plots/',sample_name,'/',Scrublet_threshold,'/',sample_name,'_scrublet',Scrublet_threshold,'.csv')
  scrb = read.csv(path,header = T)
  
  
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i_filtered = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_filtered = CreateSeuratObject(counts = data_i_filtered, project = "BM", min.cells = 3, min.features = 1)
  
  colSum_list = colSums(data_i_filtered ) # Needs to be from Matrix library
  keep = colSum_list >= 100
  data_i_filtered = data_i_filtered[,keep]
  
  
  
  data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  data_i_filtered_run$Scrublet_Boolean = as.logical(as.character(toupper(scrb$Scrublet_Boolean)))
  data_i_filtered_run$doublet_scores = as.numeric(scrb$doublet_scores)
  data_i_filtered_run$predicted_doublet = as.logical(as.character(toupper(scrb$predicted_doublet)))
  
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/Umap Plots/',sample_name,'/',Scrublet_threshold,'/')
  dir.create(folder,recursive = T)
  
  pathName = paste0(folder,sample_name,'_Umap',Scrublet_threshold,'.png')
  png(file=pathName,width=1000, height=1000)
  print(DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  #browser()
  pathName = paste0(folder,sample_name,'_Scrublet_Boolean_',Scrublet_threshold,'.png')
  png(file=pathName,width=1000, height=1000)
  print(DimPlot(data_i_filtered_run,pt.size = ifelse(data_i_filtered_run$Scrublet_Boolean == T, 2, 0.5), cols = c('navajowhite2','tomato2'), reduction = "umap",label = FALSE, group.by  = 'Scrublet_Boolean'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_doublet_scores_',Scrublet_threshold,'.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, features = 'doublet_scores'))
  dev.off()
  
  if (!all((scrb$predicted_doublet) == 0)){
    pathName = paste0(folder,sample_name,'_predicted_doublet_',Scrublet_threshold,'.png')
    png(file=pathName,width=1000, height=1000)
    print(DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE, group.by  = 'predicted_doublet'))
    dev.off()
  }
  data_i_filtered_run[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered_run, pattern = "^MT-")
  
  pathName = paste0(folder,sample_name,'_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  print(pathName)
  # 
  next
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')

  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/Umap Plots/',sample_name,'/',Scrublet_threshold,'/','Featureplots/')
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

  # Remove doublets
  data_i_filtered_run = data_i_filtered_run[,!(data_i_filtered_run$Scrublet_Boolean)]

  path = paste0(folder,sample_name,'_scrublet',Scrublet_threshold,'.Robj')
  save(data_i_filtered_run,file= path)


  
}

i = 1
# Soup, empty drops, MT
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  
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
  
  data_matrix_filtered = data_i_filtered_run@assays[["RNA"]]@counts
  data_matrix_raw = data_i_raw@assays[["RNA"]]@counts
  cluster_IDs = factor(as.character(Idents(data_i_filtered_run)))
  
  sc = SoupChannel(data_matrix_raw, data_matrix_filtered)
  sc = setClusters(sc, cluster_IDs)
  sc = setDR(sc, data_i_filtered_run@reductions[["umap"]]@cell.embeddings)
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes, HB = HBGenes))
  
  
  sc = calculateContaminationFraction(sc, list(IG = igGenes, HB = HBGenes), useToEst = useToEst)
  #head(sc$metaData)
  #quantile(sc$metaData$rho)
  
  out = adjustCounts(sc)
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_Empty_MT/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'Umap','_SoupAdjust','.png')
  png(file=pathName,width=1000, height=1000)
  print(plotChangeMap(sc, out, geneSet = HBGenes))
  dev.off()
  
  
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_Empty_MT/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'PreSoup_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'PreSoup_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'PreSoup_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  print(pathName)
  
  # Put soup data back into filtered run
  data_i_filtered_run_soup = data_i_filtered_run
  
  data_i_filtered_run_soup@assays[["RNA"]]@counts = out
  
  
  
  data_i_filtered_run_soup = NormalizeData(data_i_filtered_run_soup, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run_soup = FindVariableFeatures(data_i_filtered_run_soup, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run_soup = ScaleData(data_i_filtered_run_soup)
  data_i_filtered_run_soup = RunPCA(data_i_filtered_run_soup,npcs = 30)
  data_i_filtered_run_soup = FindNeighbors(data_i_filtered_run_soup, dims = 1:30)
  data_i_filtered_run_soup = FindClusters(data_i_filtered_run_soup)
  data_i_filtered_run_soup = RunUMAP(data_i_filtered_run_soup, dims = 1:30)
  data_i_filtered_run_soup[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered_run, pattern = "^MT-")
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder)
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  save(data_i_filtered_run_soup,file= path)
  
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run_soup,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run_soup,pt.size = 0.5, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run_soup,pt.size = 0.5, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  print(pathName)
  # 
  
  
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/','Featureplots/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_run_soup, feature = gene, folder =folder,
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


# Soup, MT, Empty
for (i in 17:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  br_e_sorted = read.csv(file = paste0(folder,'br_e_sorted_NFeatures100','.csv'))
  
  data_i_filtered = data_soup
  data_i_filtered$emptyProb = br_e_sorted$LogProb
  data_i_filtered$is_cell = br_e_sorted$is_cell
  
  data_i_filtered = data_i_filtered[,data_i_filtered$is_cell]
  data_i_filtered = data_i_filtered[, data_i_filtered$percent.mt < percent_mt]
  
  
  data_i_filtered_run = data_i_filtered
  data_i_filtered_run = NormalizeData(data_i_filtered_run, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  pt.size  = 0.8
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_Empty_MT/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_Post_Empty_MT_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_Empty_MT_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_Empty_MT_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_Empty_MT/',sample_name,'/','Featureplots/PostMT/')
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
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_Empty_MT/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted_Empty_MT.Robj')
  save(data_i_filtered_run,file= path)
  
  
}
