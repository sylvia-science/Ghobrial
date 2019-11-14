
###################
# Functions
###################
run_pipeline = function(filename,folder_input,sample_name,sampleParam,filter,regress_TF){
  print(sample_name)
  
  cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]

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
  
  
  pathName <- paste0(folder_input,'PCA/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data,ndims = PCA_dim))
  dev.off()
  
  
  
  # Cluster with Umap
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  print(paste0('Resolution' ,': ', resolution_val))
  print(folder_input)
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

quality_control <- function(data,filter,nFeature_RNA_list,percent_mt,sample_name){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  # Calculating percent Mitochondrial 
  print(nFeature_RNA_list)
  browser()
  data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-")
  
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
# TO DO: Run PCA with 30 dim first, save images, then run with param values
#######################################
visualize_PCA = function(data,folder,sample_name,PCA_dim){
  print('Visualize PCA')
  print(data[["pca"]], dims = 1:5, nfeatures = 5)
  
  pathName <- paste0(folder,'PCA/DimLoading.png')
  png(file=pathName,width=600, height=350)
  print(VizDimLoadings(data, dims = 1:2, reduction = "pca"))
  dev.off()
  
  
  ## JackStrawPlot
  # TO DO: Check if JackStrawPlot already exists with more dimensions
  if (FALSE){
    data <- JackStraw(data, num.replicate = 100)
    data <- ScoreJackStraw(data, dims = 1:PCA_dim)
  
    pathName <- paste0(folder,'PCA/jackstraw_',PCA_dim,'.png')
    png(file=pathName,width=600, height=350)
    print(JackStrawPlot(data, dims = 1:PCA_dim))
    dev.off()
  }
  ## ADD DOTPLOT
  pathName <- paste0(folder,'PCA/DimHeatMap1_6.png')
  png(file=pathName,width=2000, height=1000, res=300)
  print(DimHeatmap(data, dims = 1:(min(6,PCA_dim)), cells = 500, balanced = TRUE))
  dev.off()
  
  if (PCA_dim <=6){
    
    pathName <- paste0(folder,paste0('PCA/DimHeatMap1_',PCA_dim,'.png'))
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 1:PCA_dim, cells = 500, balanced = TRUE))
    dev.off()
  }else if (PCA_dim <=10){
    
    pathName <- paste0(folder,paste0('PCA/DimHeatMap1_',5,'.png'))
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 1:5, cells = 500, balanced = TRUE))
    dev.off()
    
    pathName <- paste0(folder,paste0('PCA/DimHeatMap6_',PCA_dim,'.png'))
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 6:PCA_dim, cells = 500, balanced = TRUE))
    dev.off()
  }else if (PCA_dim <=20){
    
    pathName <- paste0(folder,paste0('PCA/DimHeatMap1_',6,'.png'))
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 1:6, cells = 500, balanced = TRUE))
    dev.off()
    
    pathName <- paste0(folder,paste0('PCA/DimHeatMap6_',12,'.png'))
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = 7:12, cells = 500, balanced = TRUE))
    dev.off()
  }
  

  ## ADD DOTPLOT
  ## TO DO: Make clusters labeled by cell cycle
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  print(PCA_dim)
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
## TO DO: Make subfolder of cluster with res and pcadim
## paste0('_Res',resolution_val, '_PCA',PCA_dim)
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
  #print(folder)
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
    
    pathName <- paste0(folder_input,'DE')
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
  for(i in seq_len(length(feature_list))){
    cell_type = cell_list[i]

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

label_cells = function(data,cluster_IDs){
  print(cluster_IDs)
  print('')
  print('Label Cells')
  # if (regress_TF){
  #   cluster_IDs <- sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name]
  # }else{
  #   cluster_IDs <- sampleParam$Cluster_IDs_pre_regress[sampleParam['Sample'] == sample_name]
  # }
  
  
  cluster_IDs = unlist(strsplit(cluster_IDs, ", ")) # Remember to always put space after , (FIX THIS)
  print(cluster_IDs)
  print(levels(data))
  if (length(cluster_IDs) == length(levels(data))){
    names(cluster_IDs) = levels(data)
    data = RenameIdents(data, cluster_IDs)
    label = '_label'
  }else{
    print('Cluster ID length does not match number of clusters')
  }
  
  
  return ((data))
}

cellMarkers = function(data,cell_features){
  data$Cell = NA
  for (i in 1:nrow(data)){
    marker_rows = grep(data$gene[i], cell_features$Markers, value=TRUE)
    marker_idx = which( cell_features$Markers %in% marker_rows)
    #browser()
    if (length(marker_idx) > 0){
      #browser()
      cell_list = cell_features$Cell[marker_idx]
      cell_list = paste(cell_list, sep="", collapse=", ") 
      data$Cell[i] = cell_list
    }else if (length(marker_idx) == 0 ){
      data$Cell[i] = ''
    }
  }
  return (data)
}
