
###################
# Functions
###################
run_pipeline = function(filename,folder,sample_name,filter){
  
  if (filter){
    folder <- paste0(folder,sample_name,'/Filtered/')
  }
  else{
    folder <- paste0(folder,sample_name,'/Unfiltered/')
  }
  print(folder)
  print(sample_name)
  makeFolders(folder,sample_name)
  
  
  data = load_data(filename)
  nFeature_RNA_list <- list(200,2500)
  percent_mt <- 15
  data = quality_control(data,filter,nFeature_RNA_list,percent_mt,folder,sample_name)
  
  #Get Variable Genes
  norm_var <- 10000
  nfeatures_var <- 2000
  return_list= gene_var(data,norm_var,nfeatures_var)
  data = return_list$data
  top10 = return_list$top10
  
  # Scale
  data <- ScaleData(data, features = rownames(data)) 
  
  # PCA
  data <- RunPCA(data, features = VariableFeatures(object = data))
  visualize_PCA(data,folder,sample_name)
  
  data = visualize_dim(data)
  JackStrawPlot(data, dims = 1:15)
  ElbowPlot(data)
  
  # Cluster with Umap
  resolution_val<- 1.4
  data = getCluster(data,resolution_val)
  plot = DimPlot(data, reduction = "umap")
  
  pathName <- paste0(folder,'Cluster/Cluster.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  # Find Cluster Biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #BM_filtered_markers  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(data, features = top10$gene)
  
  
  # Get gene Descriptions
  #gene_desc_top10 =  get_gene_desc(top10)
  #gene_desc_top10
  return(data)
  
}

################
## Load Data
################

load_data <- function(filename) {
  
  data <- Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  
  data <- CreateSeuratObject(counts = data, project = "BM", min.cells = 3, min.features = 200)
  return (data)
}

##################
## Quality Control
##################

quality_control <- function(data,filter,nFeature_RNA_list,percent_mt,folder,sample_name){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  # Calculating percent Mitochondrial 
  print(nFeature_RNA_list)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  
  if (filter == TRUE){
    print('hi')
    print(nFeature_RNA_list)
    data <- subset(data, subset =  nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15) # For some reason this no longer works with variables
    
  }
  
  # Visualize QC metrics as a violin plot
  plot = VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)
  
  
  pathName <- paste0(folder,'QC Metrics/violin.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  # Visualize Feature-Feature Relationships
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  plot = CombinePlots(plots = list(plot1, plot2))
  
  pathName <- paste0(folder,'QC Metrics/scatter.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  return (data)
}


###################################
## Identify Highly Variable Genes
###################################

gene_var = function(data, norm_var, nfeatures_var){
  #Normalize
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = norm_var)
  
  # Identification of highly variable features
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_var) 
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(data), 10)
  top10
  return_list <- list("data" = data, "top10" = top10)
  
  return (return_list)
}




#######################################
## PCA Results
# Examine and visualize PCA results a few different ways
#######################################
visualize_PCA = function(data,folder,sample_name){
  print(data[["pca"]], dims = 1:5, nfeatures = 5)
  
  
  
  plot = VizDimLoadings(data, dims = 1:2, reduction = "pca")
  pathName <- paste0(folder,'PCA/DimLoading.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  plot = DimPlot(data, reduction = "pca")
  pathName <- paste0(folder,'PCA/DimPlot.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  plot = DimHeatmap(data, dims = 1:6, cells = 500, balanced = TRUE)
  pathName <- paste0(folder,'PCA/DimHeatMap.png')
  png(file=pathName)
  print(plot)
  dev.off()
}




## Get Visualize Dimension Data
visualize_dim = function(data){
  # Determine the 'dimensionality' of the dataset
  
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # computation time
  invisible(data <- JackStraw(data, num.replicate = 100))
  invisible(data <- ScoreJackStraw(data, dims = 1:20))
  # The JackStrawPlot function provides a visualization tool for comparing \
  # the distribution of p-values for each PC with a uniform distribution 
  # (dashed line). 'Significant' PCs will show a strong enrichment of 
  # features with low p-values (solid curve above the dashed line). In this # case it appears that there is a sharp drop-off in significance after the #first 10-12 PCs.
  
  
  
  return(data)
}


## Cluster
# Cluster 
getCluster = function(data,resolution_val){
  data <- FindNeighbors(data, dims = 1:10)
  data <- FindClusters(data, resolution = resolution_val)
  head(Idents(data), 5)
  data <- RunUMAP(data, dims = 1:10)
}



## Make subfolders
makeFolders = function(folder,sample_name){
  print('here')
  print(folder)
  print(sample_name)
  
  pathName <- paste0(folder,'QC Metrics')
  print(pathName)
  dir.create( pathName, recursive = TRUE)
  
  pathName <- paste0(folder,'PCA')
  dir.create( pathName, recursive = TRUE)
  
  pathName <- paste0(folder,'Cluster')
  dir.create( pathName, recursive = TRUE)
  
  pathName <- paste0(folder,'Cell Type')
  dir.create( pathName, recursive = TRUE)
}


## Get gene descriptions
# Get gene descriptions
get_gene_desc = function(top10){
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  attributes = listAttributes(ensembl)
  
  genedesc_pheno <- getBM(attributes=c('external_gene_name','phenotype_description'), filters = 'external_gene_name', values = top10['gene'], mart =ensembl)
  genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = top10['gene'], mart =ensembl)
  
  # Concatonate phenotype_description into one long description for each gene
  # VERY inefficient
  
  
  genedesc$phenotype_description <- ''
  for(i in seq_len(nrow(genedesc_pheno))){
    #print('i')
    #print(i)
    idx = which(genedesc$external_gene_name == genedesc_pheno$external_gene_name[i])
    val = genedesc_pheno$phenotype_description[i]
    #print(idx)
    #print(val)
    if (val!='' & genedesc['phenotype_description'][idx,] == ''){
      genedesc['phenotype_description'][idx,] = paste(genedesc['phenotype_description'][idx,],val)
    }
    else if (val!='')
    {
      genedesc['phenotype_description'][idx,] = paste(genedesc['phenotype_description'][idx,],',',val)
    }
  }
  
  tmp1 = genedesc['external_gene_name']
  tmp2 = top10['gene']
  
  result = setdiff(tmp2$gene,tmp1$external_gene_name) # Not all of the gene names might be found, so the description list can be shorter than the original
  
  if(length(result) > 0){
    genedesc[nrow(genedesc) + 1,] = list(result,vector("list", length(result)))
  }
  
  #genedesc
  #top10$gene
  index <- match(genedesc$external_gene_name,top10$gene)
  genedesc$index = index
  
  
  #genedesc['order'] = 
  #vector1 = genedesc$external_gene_name
  
  
  result = genedesc[
    with(genedesc, order(index)),
    ]
  #top10['transcript_biotype']= result['transcript_biotype']
  top10['description']= result['description']
  top10['phenotype_description']= result['phenotype_description']
  
  return(top10)
}




get_cellType = function(data,folder,sample_name,filter){
  
  sample <- data
  cell_list <- list(
    'bcell_activated',
    'bcell_plasma',
    'bcell_memory',
    'bcell_marginal',
    'bcell_follicular',
    'bcell_regulatory',
    'tcell_general',
    'tcell_activated',
    'tcell_effector',
    'tcell_regulatory',
    'tcell_exhausted',
    'tcell_helper_1',
    'tcell_helper_2',
    'tcell_naive',
    'tcell_cytotoxic',
    'tcell_memory',
    'tcell_helper_17',
    'monocyte_inflammatory',
    'monocyte_resident',
    'monocyte_CD14',
    'monocyte_FCGR3A',
    'macrophages',
    'nk_cell',
    'DC',
    'pDC',
    'megakaryocyte',
    'erythrocyte',
    'neutrophil',
    'eosionophil',
    'basophil',
    'mast',
    'hsc'
    
  )
  feature_list <- list(
    list('CD19','IL2RA','CD30'),
    list('CD27','CD38','SDC1','SLAMF7','IL6'),
    list('MS4A1','CD27','CD40','CD80','PDCD1LG2', 'CXCR3','CXCR4','CXCR5','CXCR6'),
    list('CD1A','CR2','CD37','NOTCH2'),
    list('CR2','CD22','FCER2'),
    list('CD1A','CD5','CR2','CD24','TLR4','IL10','TGFB1'),
    list('CD3D','CD3E','CD3G','PTPRC'),
    list('CD69','IL2RA'),
    list('CD3D','B3GAT1','PDCD1','FAS','CCR7'),
    list('CD4','IL2RA','FOXP3','SMAD3','STAT5A','STAT5B','IL10'),
    list('CD3D','PDCD1','FAS','CCR7','SELL','LAG3','HAVCR2','TIGIT','ENTPD1'),
    list('CD4','CXCR3','TBX21','STAT1','STAT6','IFNG'),
    list('CD4','CCR4','PTGDR2','GATA3','STAT5','STAT6','IL4','IL5','IL13'),
    list('IL7R','S100A4','CD3D','SELL','CD27'),
    list('CD8A', 'CD8B'),
    list('IL7R','CD3D','CCR7','SELL'),
    list('CD4','CCR6','RORC','RORA','STAT3','IL17A','IL17F'),
    list('CCR2'),
    list('CXCR1'),
    list('CD14', 'FCN1', 'LYZ'),
    list('CD14','FCGR3A','MS4A7'),
    list('CD68','CCR5','TFRC','ITGAM','FCGR1A','CSF1R','MRC1','CD163'),
    list('NKG7','GNLY','NCR1','KLRD1','NCAM1'),
    list('FCER1A', 'ITGAX', 'CD83', 'THBD','CD209','CD1C', 'LYZ'),
    list('IL3RA','CLEC4C','NRP1'),
    list('PPBP', 'ITGA2B','GP9', 'GP1BA', 'ITGB3'),
    list('GYPA', 'BLVRB', 'HBB', 'HBA1'),
    list('FUT4','ITGAM','FCGR3A','ITGB2','FCGR2A','CD44','CD55'),
    list('IL5RA','CCR3','EMR1','ADGRE1','SIGLEC8'),
    list('IL3RA', 'ENPP3'),
    list('KIT','CD33','TPSAB1', 'CD9'),
    list('CD34','THY1', 'CDK6')
    
  )
  
  if (filter){
    folder <- paste0(folder,sample_name,'/Filtered/')
  }
  else{
    folder <- paste0(folder,sample_name,'/Unfiltered/')
  }
  
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  for(i in seq_len(length(feature_list))){
    cell_type = cell_list[i]
    x <- unlist(feature_list[i])
    gene_list = (x[x %in% markers$gene])  
    # FeaturePlot gives error if no genes are found, so we have to check if they are included in the highly variable genes
    if (length(gene_list > 0)){
      print(paste0( cell_type,': Found'))
      print(x)
      plot = print( FeaturePlot(sample, features = c(x)))
      
      pathName <- paste0(folder,'Cell Type/',cell_type,'.png')
      png(file=pathName)
      print(plot)
      dev.off()
    }
    else {
      print(paste0( cell_type,': Not Found'))
      print(x)
    }
    print('')
    
  }
  
}

