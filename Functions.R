
###################
# Functions
###################


################
## Load Data
################

load_data <- function(filename) {
  print(filename)
  data = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  
  data = CreateSeuratObject(counts = data, project = "BM", min.cells = 3, min.features = 200)
  return (data)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


###################################
## Identify Highly Variable Genes
###################################

gene_var = function(data, nfeatures_val){
  #Normalize
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  # Identification of highly variable features
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures_val) 
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(data), 10)
  top10
  return_list <- list("data" = data, "top10" = top10)
  
  return (return_list)
}




##############################
## Cluster
##############################
getCluster = function(data,resolution_val, PCA_dim){
  print('Get Cluster')
  data = FindNeighbors(data, dims = 1:PCA_dim)
  data = FindClusters(data, resolution = resolution_val) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
  head(Idents(data), 5)
  #data <- RunTSNE(data, dims = 1:PCA_dim)
  data = RunUMAP(data, dims = 1:PCA_dim)
  return (data)
}


##############################
## Make subfolders
##############################
makeFolders = function(folder,sample_name,filter,regress_TF,makeFolder_TF,nFeature_RNA_list = NA, percent_mt = NA){
  print('Make Folders')
  
  if (filter){
    folder = paste0(folder, sample_name,'/', 'nFeature',nFeature_RNA_list[1],'_',nFeature_RNA_list[2],'_MT',percent_mt,'/filtered/')
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
    
    pathName <- paste0(folder_input,'Stats')
    dir.create( pathName, recursive = TRUE)
  }
  return (folder_input)
}


#################################
## Label Clusters
#################################

label_cells = function(data,cluster_IDs){
  print(cluster_IDs)
  print('')
  print('Label Cells')
  
  cluster_IDs = unlist(strsplit(cluster_IDs, ",")) 
  cluster_IDs = trimws(cluster_IDs, which = c("both"), whitespace = "[ \t\r\n]")
  #cluster_IDs = trimws(cluster_IDs, which = c("both"), whitespace = " \t\n\r\v\f")
  #cluster_IDs = gsub("\\W", "", cluster_IDs)
  
  
  print(cluster_IDs)
  print(levels(data))
  if (length(cluster_IDs) == length(levels(data))){
    names(cluster_IDs) = levels(data)
    data = RenameIdents(data, cluster_IDs)
  }else{
    print('Cluster ID length does not match number of clusters')
  }
  
  
  return ((data))
}

#################################
## Add known markers colomn
#################################

cellMarkers = function(data,cell_features){
  #browser()
  data$Cell = NA
  for (i in 1:nrow(data)){
    marker_idx = c()
    
    #print(i)
    #print(data$gene[i])
    #browser()
    for (j in 1:nrow(cell_features)){

      gene_list = cell_features$Markers[j]
      gene_list = gsub("[[:space:]]", "", gene_list)
      
      gene_list = as.list(strsplit(gene_list, ",")[[1]])
      gene_list = trimws(gene_list, which = c("both"), whitespace = " \t\n\r\v\f")
      gene_list = gsub("\\W", "", gene_list)
      if (data$gene[i] == "CALM2"){
        #print(cell_features$Cell[j])
        #browser()
      }
      if (any(gene_list == data$gene[i])){
        
        marker_idx = c(marker_idx,j)
      }
    }
    # THIS WILL ALSO FIND MATCHING SUBSTRINGS
    # marker_rows = grep(data$gene[i], cell_features$Markers, value=TRUE)
    # marker_idx = which( cell_features$Markers %in% marker_rows)
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

clusterStats = function(data){
  #browser()
  cluster_num = table(Idents(data))
  cluster_num = as.data.frame(cluster_num)
  cluster_num$Percent = 100*cluster_num$Freq/sum(cluster_num$Freq)
  names(cluster_num) = c("Cluster", "Num",'Percent')
  cluster_num = cluster_num[order( levels(cluster_num$Cluster)),]
  return(cluster_num)
}

############################################################
# Convert list of individual genes into one list of csv
############################################################
ConvertCellMarkers = function(cellMarker){
  output = data.frame(matrix(ncol = 4, nrow = 1))
  colnames(output) = c("Paper",'Plot_marker', "Cell", "Markers")
  
  
  output$Plot_marker = TRUE
  output$Markers = paste( unlist(cellMarker$'Gene'), collapse=', ')
  output$Cell = cellMarker$'..Cell'[1]
  
  return(output)
}

SummariseSampleByCell = function(data){
  #browser()
  sample_list = unique(data$sample)
  celltype_list = unique(Idents(data))
  output = data.frame(matrix(ncol = length(celltype_list), nrow = length(sample_list)))
  colnames(output) = celltype_list
  rownames(output) = sample_list
  
  
  for (sample in sample_list){
    for (celltype in celltype_list){
      if (length(data$sample == sample & Idents(data) == celltype)  > 0){
      #browser()
      output[sample,celltype] = sum(data$sample == sample & Idents(data) == celltype)
      }else{
        output[sample,celltype]  = 0  
      }
    }
  }
  
  
  
  
  return(output)
}

getCellMarkers = function(folder){
  #browser()
  cell_features_file = paste0(folder,'Data/Cell_IDS.xlsx')
  cell_features = read_excel(cell_features_file)
  
  # cell_features_file = paste0(folder, '/Cell Markers/bcell_naive.csv')
  # bcell_naive = read.csv(cell_features_file)
  # bcell_naive = ConvertCellMarkers(bcell_naive)
  # 
  # cell_features_file = paste0(folder, '/Cell Markers/tcell_CD4_memory_TREG.csv') 
  # tcell_CD4_memory_TREG = read.csv(cell_features_file)
  # tcell_CD4_memory_TREG = ConvertCellMarkers(tcell_CD4_memory_TREG)
  # 
  # cell_features_file = paste0(folder, '/Cell Markers/tcell_CD4_naive.csv')  
  # tcell_CD4_naive = read.csv(cell_features_file)
  # tcell_CD4_naive = ConvertCellMarkers(tcell_CD4_naive)
  
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_naive_activated.csv'
  # tcell_CD4_naive_activated = read.csv(cell_features_file)
  # tcell_CD4_naive_activated = ConvertCellMarkers(tcell_CD4_naive_activated)
  # 
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_naive_TREG.csv'
  # tcell_CD4_naive_TREG = read.csv(cell_features_file)
  # tcell_CD4_naive_TREG = ConvertCellMarkers(tcell_CD4_naive_TREG)
  # 
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TFH.csv'
  # tcell_CD4_TFH = read.csv(cell_features_file)
  # tcell_CD4_TFH = ConvertCellMarkers(tcell_CD4_TFH)
  # 
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH1.csv'
  # tcell_CD4_TH1 = read.csv(cell_features_file)
  # tcell_CD4_TH1 = ConvertCellMarkers(tcell_CD4_TH1)
  # 
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH2.csv'
  # tcell_CD4_TH2 = read.csv(cell_features_file)
  # tcell_CD4_TH2 = ConvertCellMarkers(tcell_CD4_TH2)
  # 
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH17.csv'
  # tcell_CD4_TH17 = read.csv(cell_features_file)
  # tcell_CD4_TH17 = ConvertCellMarkers(tcell_CD4_TH17)
  # 
  # cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH117.csv'
  # tcell_CD4_TH117 = read.csv(cell_features_file)
  # tcell_CD4_TH117 = ConvertCellMarkers(tcell_CD4_TH117)
  
  
  
  # cell_features = rbind(cell_features,
  #                       bcell_naive,
  #                       tcell_CD4_memory_TREG,tcell_CD4_naive, 
  #                       tcell_CD4_naive_activated,
  #                       tcell_CD4_naive_TREG, tcell_CD4_TFH,
  #                       tcell_CD4_TH1, tcell_CD4_TH2,
  #                       tcell_CD4_TH17, tcell_CD4_TH117)
  
  
  cell_features = cell_features[cell_features$Plot_marker == 1,]
  #browser()
  return(cell_features)
}

SaveAsMatrix = function(data,folder_base_output){
  data_df = data.frame(data@assays$RNA@data) 
  
  data_df = t(data_df)
  data_df = data.frame(data_df)
  
  split_var =  data@meta.data$split_var
  data_df$category = split_var
  data_df$ident =   levels(data@active.ident)
  setDT(data_df)
  
  fwrite(data_df, paste0(folder_base_output ,"df.csv"))
  
  write.csv(data_df, paste0(folder_base_output ,"df.csv"))
}

ScranNorm = function(data){
  #browser()
  counts =  data@assays[['RNA']]@counts
  
  data_i_sc =  SingleCellExperiment(list(counts=counts))
  set.seed(100)
  clusters = quickCluster(data_i_sc)
  data_i_sc = computeSumFactors(data_i_sc, clusters=clusters)
  data_i_sc  = logNormCounts(data_i_sc, size_factors = sizeFactors(data_i_sc),log = T,center_size_factors = T )
  
  #colnames(data_i_sc) <- sub('[.]', '_', make.names(colnames(data_i_sc), unique=TRUE))
  
  data_out = data

  colnames(data_i_sc@assays@data@listData[["logcounts"]]) = colnames(counts)
  # DOUBLE CHECK
  data_out@assays[["RNA"]]@data = as.matrix(data_i_sc@assays@data@listData[["logcounts"]])
  data_out@assays[["RNA"]]@scale.data = as.matrix(data_i_sc@assays@data@listData[["logcounts"]]) # need to find out how to put an empty placeholder here
  size_factor_list = sizeFactors(data_i_sc)
  names(size_factor_list) = colnames(data_out)
  data_out$sizeFactors = size_factor_list
  
  return(data_out)
}

RemoveNan = function(data, slot){
  browser()
  # Remove cells with Nan values
  if (slot == 'data'){
    tmp = data[[data@active.assay]]@data
    tmp = as.matrix(tmp)
    cell_name = names(which(is.na(colSums(tmp))))
    cell_val = colnames(data[[data@active.assay]]@data) != cell_name
    if (length(cell_name) !=0){
      mask = !(colnames(data[[data@active.assay]]@data) %in% cell_name)
      data[[data@active.assay]]@data = data[[data@active.assay]]@data[,mask]
    }
    
    return(data)
  }else if (slot == 'scale.data'){
    tmp = data[[data@active.assay]]@scale.data
    tmp = as.matrix(tmp)
    cell_name = names(which(is.na(colSums(tmp))))
   
    
    if (length(cell_name) !=0){
      mask = !(colnames(data[[data@active.assay]]@scale.data) %in% cell_name)
      data[[data@active.assay]]@scale.data = data[[data@active.assay]]@scale.data[,mask]
    }
    
    #data_integrate_test = SubsetData(data, cells = cell_val)
    
    return(data)
  }
  
}

renameCells = function(data, idents,newident){
  #browser()
  newIdents = as.character(Idents(data))
  newIdents[newIdents %in% idents] = newident
  print(unique(newIdents))
  Idents(data) = newIdents
  return(data)
  
}
