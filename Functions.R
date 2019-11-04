
###################
# Functions
###################
run_pipeline = function(filename,folder_input,sample_name,sampleParam,filter,regress_TF){
  print(sample_name)
  
  cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  # if (cluster_IDs == 'tmp'){
  #   file_str = ''
  # }else
  # {
  #   file_str = '_label'
  # }
  file_str = ''
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
  metaData <- read_excel(filename_metaData)
  sampleParam <- read_excel(filename_sampleParam)
  data = load_data(filename)
  
  cell_features_file <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
  cell_features <- read_excel(cell_features_file)
  
  # QC
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  print('hello')
  data = quality_control(data,filter,nFeature_RNA_list,percent_mt,folder_input,sample_name)
  
  #Get Variable Genes
  nfeatures_val <- sampleParam$nfeatures_val[sampleParam['Sample'] == sample_name]
  return_list= gene_var(data,nfeatures_val)
  data = return_list$data
  top10 = return_list$top10
  
  #Score for Cell Cycle gene expression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  
  
  # Scale
  if (regress_TF){ # Check if we need to regress out
    # CHANGE vars.to.regress WITH VALUES FROM PARAM
    data <- ScaleData(data, features = rownames(data), vars.to.regress = c("nCount_RNA", "percent.mt"))
  }else{
    data <- ScaleData(data, features = rownames(data)) 
  }
  # PCA
  PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  data <- RunPCA(data, features = VariableFeatures(object = data), npcs = PCA_dim)
  visualize_PCA(data,folder_input,sample_name,PCA_dim)
  
  data = visualize_dim(data,PCA_dim)
  #JackStrawPlot(data, dims = 1:PCA_dim)
  
  print(folder_input)
  pathName <- paste0(folder_input,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  
  
  # Cluster with Umap
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  print(paste0('Resolution' ,': ', resolution_val))
  data <- getCluster (data,resolution_val, PCA_dim)
  
  # Name cells
  #data = label_cells(data,folder_input,sample_name,sampleParam,resolution_val,filter,regress_TF, file_str)
  
  
  
  # Find Cluster Biomarkers
  print('Find Cluster Biomarkers')
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  
  top20$Cell = NA
  # Add known markers to top20
  #browser()
  
  for (i in 1:nrow(top20)){
    marker_rows = grep(top20$gene[i], cell_features$Markers, value=TRUE)
    marker_idx = which( cell_features$Markers %in% marker_rows)
    #browser()
    if (length(marker_idx) > 0){
      #browser()
      cell_list = cell_features$Cell[marker_idx]
      cell_list = paste(cell_list, sep="", collapse=", ") 
      top20$Cell[i] = cell_list
    }else if (length(marker_idx) == 0 ){
      top20$Cell[i] = ''
    }
  }
  
  
  
  write.csv(top20, file = paste0(folder_input,'Top20Features',file_str,'.csv'),row.names=FALSE)
  #browser()
  # Plot Umap
  pathName <- paste0(folder_input,paste0('Cluster/ClusterUmap',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  # Visualize clustering
  # Cluster Metrics
  pathName <- paste0(folder_input,'Cluster/ClusterMetrics.png')
  png(file=pathName,width=600, height=350)
  print(FeaturePlot(data, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  # Cluster Heatmap
  plot = DoHeatmap(data, features = top10$gene)
  pathName <- paste0(folder_input,paste0('Cluster/HeatMap',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(plot)
  dev.off()
  
  
  # Get gene Descriptions
  #gene_desc_top10 =  get_gene_desc(top10)
  #gene_desc_top10
  #########################################################################################
  
  
  return(data)
  
}

################
## Load Data
################

load_data <- function(filename) {
  print(filename)
  data <- Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  
  data <- CreateSeuratObject(counts = data, project = "BM", min.cells = 3, min.features = 200)
  return (data)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
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
    
    print(nFeature_RNA_list)
    
    print(nFeature_RNA_list[2])
    nFeature_RNA_min = as.numeric(nFeature_RNA_list[1])
    nFeature_RNA_max = as.numeric(nFeature_RNA_list[2])
    #browser()
    nFeature_RNA_tmp = 100
    
    # Can't use subset function because doesn't work with variables
    # This is a workaround
    expr <- FetchData(object = data, vars = 'nFeature_RNA')
    data = data[, which(x = expr > nFeature_RNA_min & expr < nFeature_RNA_max)]
    
    expr <- FetchData(object = data, vars = 'percent.mt')
    data = data[, which(x = expr < percent_mt)]
    
  }
  
  # Visualize QC metrics as a violin plot
  plot = VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4)
  pathName <- paste0(folder,'QC Metrics/violin.png')
  print(pathName)
  png(file=pathName,width=600, height=600)
  print(plot)
  dev.off()
  
  # Visualize Feature-Feature Relationships
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") + xlim(0, 30000) + ylim(0, 100)
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + xlim(0, 30000)  + ylim(0, 5000)
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

gene_var = function(data, nfeatures_val){
  #Normalize
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  # Identification of highly variable features
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val*2) 
  
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
visualize_PCA = function(data,folder,sample_name,PCA_dim){
  print('Visualize PCA')
  print(data[["pca"]], dims = 1:5, nfeatures = 5)
  
  plot = VizDimLoadings(data, dims = 1:2, reduction = "pca")
  pathName <- paste0(folder,'PCA/DimLoading.png')
  png(file=pathName,width=600, height=350)
  print(plot)
  dev.off()
  
  ## ADD DOTPLOT
  pathName <- paste0(folder,'PCA/DimHeatMap1_6.png')
  png(file=pathName,width=2000, height=1000, res=300)
  print(DimHeatmap(data, dims = 1:6, cells = 500, balanced = TRUE))
  dev.off()
  
  if (PCA_dim >=12){
    pathName <- paste0(folder,'PCA/DimHeatMap7_12.png')
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 7:12, cells = 500, balanced = TRUE))
    dev.off()
  }else{
    pathName <- paste0(folder,'PCA/DimHeatMap7_12.png')
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 7:PCA_dim, cells = 500, balanced = TRUE))
    dev.off()
  }
  ## ADD DOTPLOT
  for (x in 1:(PCA_dim -1)){
    y <- x+1
    pathName <- paste0(folder,'PCA/PCA',x,'_',y,'.png')
    png(file=pathName,width=600, height=350)
    print(DimPlot(data, dims = c(x,y), reduction = "pca",pt.size = 2))
    dev.off()
  }
  
}




## Get Visualize Dimension Data
visualize_dim = function(data,PCA_dim){
  # Determine the 'dimensionality' of the dataset
  
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # computation time
  #invisible(data <- JackStraw(data, num.replicate = 100))
  #invisible(data <- ScoreJackStraw(data, dims = 1:PCA_dim))

  
  return(data)
}


## Cluster
# Cluster 
getCluster = function(data,resolution_val, PCA_dim){
  print('Get Cluster')
  data <- FindNeighbors(data, dims = 1:PCA_dim)
  data <- FindClusters(data, resolution = resolution_val) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
  head(Idents(data), 5)
  #data <- RunTSNE(data, dims = 1:PCA_dim)
  data <- RunUMAP(data, dims = 1:PCA_dim)
}



## Make subfolders
makeFolders = function(folder,sample_name,filter,regress_TF,makeFolder_TF){
  print('Make Folders')
  if (filter){
    folder <- paste0(folder,sample_name,'/Filtered/')
  }else if (filter == FALSE){
    folder <- paste0(folder,sample_name,'/Unfiltered/')
  }

  if (regress_TF){
    subfolder = 'Regress/' # Save files in regression folder
  }else{
    subfolder = 'No_Regress/' # Save files in No regression folder
  }
  folder_input = paste0(folder,subfolder)
  print(folder)
  print(folder_input)
  if (makeFolder_TF){
    pathName <- paste0(folder_input,'QC Metrics')
    dir.create( pathName, recursive = TRUE)
    
    pathName <- paste0(folder_input,'PCA')
    dir.create( pathName, recursive = TRUE)
    
    pathName <- paste0(folder_input,'Cluster')
    dir.create( pathName, recursive = TRUE)
    
    pathName <- paste0(folder_input,'Cell Type')
    dir.create( pathName, recursive = TRUE)
  }
  return (folder_input)
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




get_cellType = function(data,data_orig,folder,sample_name){
  
  print(folder)

  sample <- data
  
  cell_features_file <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
  cell_features <- read_excel(cell_features_file)
  
  # cell_list <- list(
  #   'bcell_activated',
  #   'bcell_plasma',
  #   'bcell_memory',
  #   'bcell_marginal',
  #   'bcell_follicular',
  #   'bcell_regulatory',
  #   'tcell_general',
  #   'tcell_activated',
  #   'tcell_effector',
  #   'tcell_regulatory',
  #   'tcell_exhausted',
  #   'tcell_helper_1',
  #   'tcell_helper_2',
  #   'tcell_naive',
  #   'tcell_cytotoxic',
  #   'tcell_memory',
  #   'tcell_helper_17',
  #   'monocyte_inflammatory',
  #   'monocyte_resident',
  #   'monocyte_CD14',
  #   'monocyte_FCGR3A',
  #   'macrophages',
  #   'nk_cell',
  #   'DC', # Add MZB1
  #   'pDC',
  #   'megakaryocyte',
  #   'erythrocyte',
  #   'neutrophil',
  #   'eosionophil',
  #   'basophil',
  #   'mast',
  #   'hsc'
  #   
  # )
  
  feature_list_old <- list(
    list('CD19','IL2RA','CD30'),
    list('CD27','CD38','SDC1','SLAMF7','IL6','CD138','TNFRSF17'),
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
    list('FCER1A', 'ITGAX', 'CD83', 'THBD','CD209','CD1C', 'LYZ', 'MZB1'),
    list('IL3RA','CLEC4C','NRP1'),
    list('PPBP', 'ITGA2B','GP9', 'GP1BA', 'ITGB3'),
    list('GYPA', 'BLVRB', 'HBB', 'HBA1'),
    list('FUT4','ITGAM','FCGR3A','ITGB2','FCGR2A','CD44','CD55'),
    list('IL5RA','CCR3','EMR1','ADGRE1','SIGLEC8'),
    list('IL3RA', 'ENPP3'),
    list('KIT','CD33','TPSAB1', 'CD9'),
    list('CD34','THY1', 'CDK6')

  )

  # Match cells with features
  # cell_features <- data.frame(cell=character(),
  #                  features=character(),
  #                  stringsAsFactors=FALSE)
  # 
  # #cell_features[,'cell'] = cell_list
  # 
  # 
  # for(i in seq_len(length(cell_list))){
  #   cell_features[i,'cell'] = cell_list[i]
  #   cell_features[i,'features'] = paste(unlist(feature_list[i]), collapse=", ")
  # 
  # }
  cell_list = cell_features$Cell
  feature_list = cell_features$Markers
  
  all_markers  = data_orig@assays[['RNA']]
  all_markers = all_markers@data@Dimnames[[1]]
  #all_markers2 <- FindAllMarkers(data)
  for(i in seq_len(length(feature_list))){
    cell_type = cell_list[i]
    x_old <- unlist(feature_list_old[i])
    

    x = unlist(strsplit(feature_list[i], ",")) 
    x = gsub("\\s", "", x)  
    
    gene_list = (x[x %in% all_markers])  
    # FeaturePlot gives error if no genes are found, so we have to check if they are included in the highly variable genes
    if (length(gene_list > 0)){
      print(paste0( cell_type,': Found'))
      print(x)
      
      # Add dotplot
      subtitle_str = paste0(cell_type,': ', toString(x))
      plot = FeaturePlot(sample, features = c(x))
      plot = plot + labs(subtitle=subtitle_str) + theme(plot.subtitle = element_text(hjust=0.5, size=16))
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
#################################
## Label Clusters
#################################

label_cells = function(data,sample_name,sampleParam){
  
  
  if (regress_TF){
    cluster_IDs <- sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  }else{
    cluster_IDs <- sampleParam$Cluster_IDs_pre_regress[sampleParam['Sample'] == sample_name]
  }
  
  
  cluster_IDs = unlist(strsplit(cluster_IDs, ", ")) # Remember to always put space after ,
  print(cluster_IDs)
  print(levels(data))
  
  if (length(cluster_IDs) == length(levels(data))){
    names(cluster_IDs) <- levels(data)
    data <- RenameIdents(data, cluster_IDs)
    label = '_label'
  }else{
    print('Cluster ID length does not match number of clusters')
  }
  
  
  return ((data))
}

###################
## Diff Expression
###################

diff_exp = function(data_input, cell_type){
  #Add active.ident to the metadata and convert to dataframe
  data_input@meta.data$active.ident <- data_input@active.ident
  data_permute <- data.frame(data_input@meta.data)
  
  #Calculate difference in cell proportion, following permutation of sample labels 10,000 times
  permuted_differences <- replicate(10000, {
    all <- sample(data_input$orig.ident)
    newPre <- which(all == "data_pre")
    newPost <- which(all == "data_post")
    diff <- (sum(data_permute$active.ident[newPost] %in% cell_type)/length(newPost)) 
    - (sum(data_permute$active.ident[newPre] %in% cell_type)/length(newPre))
    return(diff)
  })
  
  #What's the proportion of differences that were larger than the observed difference?
  obsdiff <- (sum(data_permute$active.ident[which(data_permute$orig.ident == "data_post")] 
                  %in% cell_type)/sum(data_permute$orig.ident == "data_post")) 
  - (sum(data_permute$active.ident[which(data_permute$orig.ident == "data_pre")] 
         %in% cell_type)/sum(data_permute$orig.ident == "data_pre"))
  
  
  output = (sum(abs(permuted_differences) > abs(obsdiff)) + 1) / (length(permuted_differences) + 1)
  return (output)
  
}

diff_exp_helper = function(data, folder_pre,sample_name,cell_type_list){
  
  # Load and merge original data
  
  df_expr <- data.frame(matrix(ncol = 2, nrow = length(cell_type_list)))
  colnames(df_expr) <- c("cell_type", "DE")
  df_expr$cell_type = cell_type_list
  for (i in 1:length(cell_type_list))
  {
    diff_exp_result = diff_exp(data,cell_type_list[i] )
    print(paste0(cell_type_list[i], ' diff_exp_result: ', diff_exp_result))
    df_expr$DE[i] = diff_exp_result
  
    
  }
  return (df_expr)
  
}


diffExpGene = function(data, folder_pre,sample_name,cell_type_list){
  
  # Load and merge original data
  
  cell_marker_list = list()

  for (i in 1:length(cell_type_list))
  {
    
    celltype.condition = data@meta.data[["celltype.condition"]]
    pre_sum = sum(celltype.condition == paste0(cell_type_list[i],'_pre'))
    post_sum = sum(celltype.condition == paste0(cell_type_list[i],'_post'))
    
    if (pre_sum > 3 & post_sum > 3){
      
      cell_markers <- FindMarkers(data, ident.1 = paste0(cell_type_list[i],'_post'), ident.2 = paste0(cell_type_list[i],'_pre'), verbose = FALSE)
      #browser()
      cell_markers['Cell'] =cell_type_list[i]
      setDT(cell_markers, keep.rownames = 'Gene')
      
      
      cell_marker_list = rbind(cell_marker_list,cell_markers)
    }
    
    
  }
  return (cell_marker_list)
  
}
