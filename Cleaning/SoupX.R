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
#metaData = metaData[metaData$Run== 1,]
metaData = metaData[metaData$`Sample Type` == 'PBMC',]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

i = 39

# Soup + MT + Normal threshold

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')

#sample_list = c('GL3404BM')
i = 1
run = T
for (i in 3:nrow(sampleParam)){
  sample_name =  sampleParam$Sample[i]
  #sample_name = metaData$Sample[i]
  #sample_name = sample_list[i]
  #sample_name = 'GL1160BM'
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  #RNA_features_min = sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
  #RNA_features_max = sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name]
  
  #filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
  
  colSum_list = colSums(data_i_raw ) # Needs to be from Matrix library
  keep = colSum_list >= 100
  data_i_filtered = data_i_raw[,keep]
  
  data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")
  
  percent_MT_list = data_i_filtered$percent.mt
  
  
  #folder = paste0('/home/sujwary/Desktop/scRNA/Output/tmp/Hist','/')
  #path = paste0(folder,sample_name,'_hist_MT.png')
  #png(file=path,width=1000, height=1000)
  #hist(percent_MT_list)
  #print( plot)
  #dev.off()
  
  ### Scran norm?
  data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  #data_i_filtered_run = ScranNorm(data_i_filtered)
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
  
  igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
              "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
  igGenes = igGenes[igGenes %in% sc[["toc"]]@Dimnames[[1]]] # usetoEst will break if not all genes are present
  HBGenes = c('HBB','HBA2')
  HBGenes = HBGenes[HBGenes %in% sc[["toc"]]@Dimnames[[1]]] 
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes, HB = HBGenes))
  
  
  sc = calculateContaminationFraction(sc, list(IG = igGenes, HB = HBGenes), useToEst = useToEst)
  
  out = adjustCounts(sc)
  print('Plot')
  
  
  # Put soup data back into filtered run
  data_i_filtered_soup = data_i_filtered
  
  data_i_filtered_soup@assays[["RNA"]]@counts = out
  
  data_i_filtered_soup = data_i_filtered_soup[, data_i_filtered_soup$percent.mt < 15]
  
  data_i_filtered_soup_run = NormalizeData(data_i_filtered_soup, normalization.method = "LogNormalize", scale.factor = 10000)
  #data_i_filtered_run = ScranNorm(data_i_filtered)
  data_i_filtered_soup_run = FindVariableFeatures(data_i_filtered_soup_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_soup_run = ScaleData(data_i_filtered_soup_run)
  data_i_filtered_soup_run = RunPCA(data_i_filtered_soup_run,npcs = 30)
  data_i_filtered_soup_run = FindNeighbors(data_i_filtered_soup_run, dims = 1:30)
  data_i_filtered_soup_run = FindClusters(data_i_filtered_soup_run)
  data_i_filtered_soup_run = RunUMAP(data_i_filtered_soup_run, dims = 1:30)
  
  
  
  folder = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  dir.create(folder,recursive = T)
  path = paste0(folder,'/',sample_name,'.Robj')
  save(data_i_filtered_soup_run,file= path)
  
  write(colnames(data_i_filtered_soup@assays[["RNA"]]@counts ), file = paste0(folder,sample_name,'_colnames.txt'))
  write(rownames(data_i_filtered_soup@assays[["RNA"]]@counts ), file = paste0(folder,sample_name,'_rownames.txt'))
  writeMM(data_i_filtered_soup@assays[["RNA"]]@counts , file = paste0(folder,sample_name,'_matrix.txt'))
  
  next
  data_i_filtered_soup = data_i_filtered_soup[, data_i_filtered_soup$percent.mt < percent_mt]
  data_i_filtered_soup = data_i_filtered_soup[, data_i_filtered_soup$nFeature_RNA > RNA_features_min]
  data_i_filtered_soup = data_i_filtered_soup[, data_i_filtered_soup$nFeature_RNA < RNA_features_max]
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_nFeature/',sample_name,'/')
  dir.create(folder,recursive = T)
  path = paste0(folder,'/',sample_name,'.Robj')
  save(data_i_filtered_run_soup,file= path)
}


# Load Soup + MT

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')

#sample_list = c('GL3404BM')
i = 1
for (i in 1:nrow(sampleParam)){
  sample_name =  sampleParam$Sample[i]
  #sample_name = sample_list[i]
  #sample_name = 'GL1420BM'
  print(sample_name)
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'.Robj')
  data_soup = loadRData(path)
  
  
  
  #next 
  var_feature = data_soup@assays[["RNA"]]@var.features
  
  data_matrix = data_soup@assays[["RNA"]]@counts
  write.csv(var_feature, file = paste0(folder,sample_name,'_var_feature.csv'))
  
  write(colnames(data_matrix), file = paste0(folder,sample_name,'_colnames.txt'))
  write(rownames(data_matrix), file = paste0(folder,sample_name,'_rownames.txt'))
  writeMM(data_matrix, file = paste0(folder,sample_name,'_matrix.txt'))
  
  
}



######################################################
# Just soup pre and post plots
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
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Feature100/',sample_name,'/')
  dir.create(folder,recursive = T)
  path = paste0(folder,'/',sample_name,'_Feature100.Robj')
  save(data_i_filtered_run,file= path)

  
  #next
  
  data_matrix_filtered = data_i_filtered_run@assays[["RNA"]]@counts
  data_matrix_raw = data_i_raw@assays[["RNA"]]@counts
  cluster_IDs = factor(as.character(Idents(data_i_filtered_run)))
  
  
  sc = SoupChannel(data_matrix_raw, data_matrix_filtered)
  sc = setClusters(sc, cluster_IDs)
  sc = setDR(sc, data_i_filtered_run@reductions[["umap"]]@cell.embeddings)
  
  igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
              "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
  igGenes = igGenes[igGenes %in% sc[["toc"]]@Dimnames[[1]]] # usetoEst will break if not all genes are present
  HBGenes = c('HBB','HBA2')
  HBGenes = HBGenes[HBGenes %in% sc[["toc"]]@Dimnames[[1]]] 
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes))
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_PreSoup_igGenes','','.png')
  png(file=pathName,width=1000, height=1000)
  print( plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst))
  dev.off()
  
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = HBGenes))
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_PreSoup_HBGenes','','.png')
  png(file=pathName,width=1000, height=1000)
  print( plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst))
  dev.off()
  
  

  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder, recursive = T)
  
  print('Plot')
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_PreSoup_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_PreSoup_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_PreSoup_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  next
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes, HB = HBGenes))
  
  
  sc = calculateContaminationFraction(sc, list(IG = igGenes, HB = HBGenes), useToEst = useToEst)
  
  out = adjustCounts(sc)
  print('Plot')
  
  pathName = paste0(folder,sample_name,'_PreSoup_Umap','_SoupAdjust','.png')
  png(file=pathName,width=1000, height=1000)
  print(plotChangeMap(sc, out, geneSet = HBGenes))
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
  pathName = paste0(folder,sample_name,'_Post_Soup_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run_soup,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_Soup_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run_soup,pt.size = 0.5, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_Soup_Umap','_nFeature_RNA','.png')
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
############################3
# Compare pre/post soup
i = 4
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Feature100/',sample_name,'/')
  dir.create(folder,recursive = T)
  path = paste0(folder,'/',sample_name,'_Feature100.Robj')
  data_i = loadRData(path)
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Pre_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  folder_featureplot = paste0(folder,'Featureplots/Pre/')
  dir.create(folder_featureplot,recursive = T)
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    plot = FeaturePlotFix(data_i, feature = gene, folder =folder,
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
    pathName = paste0(folder_featureplot,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Post_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_soup,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_soup,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_soup,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  folder_featureplot = paste0(folder,'Featureplots/Post/')
  dir.create(folder_featureplot,recursive = T)
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    plot = FeaturePlotFix(data_soup, feature = gene, folder =folder,
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
    pathName = paste0(folder_featureplot,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  
}



#####################################
i = 4
# Soup loading previous data filtering for MT
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
  #data_i_filtered = data_soup
  
  data_i_run = data_soup

  
  boxplot_df = data.frame(matrix(ncol = 3,nrow = length( data_i_run$is_cell)))
  colnames(boxplot_df) = c("Is_Cell",'Percent_MT', 'count')
  
  boxplot_df$Is_Cell = factor(data_i_run$is_cell)
  boxplot_df$Percent_MT = data_i_run$percent.mt
  boxplot_df$count = data_i_run$nCount_RNA
   
  boxplot_df = boxplot_df[!is.na(boxplot_df$Is_Cell),]
  
  boxplot_df_MT_15 = boxplot_df[ boxplot_df$Percent_MT < 15,]
  
  summary(boxplot_df_MT_15$count[ boxplot_df_MT_15$Is_Cell == F])
  
  boxplot_summary = data.frame(matrix(ncol = 4,nrow = 1))
  colnames(boxplot_summary) = c("Is_Cell_T",'Is_Cell_F', 'Is_Cell_MT15_T','Is_Cell_MT15_F')
  
  boxplot_summary$is_cell_F = sum(boxplot_df$Is_Cell == F)
  boxplot_summary$is_cell_T =  sum(boxplot_df$Is_Cell == T)
  
  boxplot_summary$Is_Cell_MT15_F  = sum(boxplot_df_MT_15$Is_Cell == F)
  boxplot_summary$Is_Cell_MT15_T  = sum(boxplot_df_MT_15$Is_Cell == T)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  dir.create(folder)
  
  write.csv(boxplot_summary, file=paste0(folder,'Boxplot/boxplot_summary'))

  
  # dir.create(paste0(folder,'Boxplot/'))
  # pathName =  paste0(folder,'Boxplot/', 'Boxplot_MTVsIs_cell_ylim15','','.png')
  # png(file=pathName,width=600, height=600)
  # plot = ggplot(boxplot_df_MT_15, aes(x = Is_Cell, y = Percent_MT)) +
  #   geom_boxplot() +
  #   theme_classic() +ylim(c(0,15))
  # plot = plot + theme(
  #   axis.title.x = element_text(color="black", size=24 ),
  #   axis.title.y = element_text(color="black", size=24),
  #   axis.text= element_text(color="black", size=24),
  #   legend.text=element_text(size=18),
  #   legend.title=element_text(size=18))
  # 
  # print(plot)
  # dev.off()
  # 
  # pathName =  paste0(folder,'Boxplot/', 'Boxplot_MTVsIs_cell_ylim100','','.png')
  # png(file=pathName,width=600, height=600)
  # plot = ggplot(boxplot_df, aes(x = Is_Cell, y = Percent_MT)) +
  #   geom_boxplot() +
  #   theme_classic() +ylim(c(0,100))
  # plot = plot + theme(
  #   axis.title.x = element_text(color="black", size=24 ),
  #   axis.title.y = element_text(color="black", size=24),
  #   axis.text= element_text(color="black", size=24),
  #   legend.text=element_text(size=18),
  #   legend.title=element_text(size=18))
  # 
  # print(plot)
  # dev.off()
  
  #next
  
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Pre_MT_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Pre_MT_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Pre_MT_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  folder_featureplot = paste0(folder,'Featureplots/PreMT/')
  dir.create(folder_featureplot,recursive = T)
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
    pathName = paste0(folder_featureplot,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  data_i_run = data_i_run[,data_i_run$percent.mt < 15]
  
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Post_MT_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  folder_featureplot = paste0(folder,'Featureplots/PostMT/')
  dir.create(folder_featureplot,recursive = T)
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
    pathName = paste0(folder_featureplot,gene,'','.png')
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
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  
  
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

