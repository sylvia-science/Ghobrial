
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
}


##############################
## Make subfolders
##############################
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
  cluster_num = table(data)
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
  output$Cell = cellMarker$'ï..Cell'[1]
  
  return(output)
}

getCellMarkers = function(){
  #browser()
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell_IDS.xlsx'
  cell_features = read_excel(cell_features_file)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/bcell_naive.csv'
  bcell_naive = read.csv(cell_features_file)
  bcell_naive = ConvertCellMarkers(bcell_naive)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_memory_TREG.csv'
  tcell_CD4_memory_TREG = read.csv(cell_features_file)
  tcell_CD4_memory_TREG = ConvertCellMarkers(tcell_CD4_memory_TREG)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_naive.csv'
  tcell_CD4_naive = read.csv(cell_features_file)
  tcell_CD4_naive = ConvertCellMarkers(tcell_CD4_naive)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_naive_activated.csv'
  tcell_CD4_naive_activated = read.csv(cell_features_file)
  tcell_CD4_naive_activated = ConvertCellMarkers(tcell_CD4_naive_activated)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_naive_TREG.csv'
  tcell_CD4_naive_TREG = read.csv(cell_features_file)
  tcell_CD4_naive_TREG = ConvertCellMarkers(tcell_CD4_naive_TREG)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TFH.csv'
  tcell_CD4_TFH = read.csv(cell_features_file)
  tcell_CD4_TFH = ConvertCellMarkers(tcell_CD4_TFH)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH1.csv'
  tcell_CD4_TH1 = read.csv(cell_features_file)
  tcell_CD4_TH1 = ConvertCellMarkers(tcell_CD4_TH1)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH2.csv'
  tcell_CD4_TH2 = read.csv(cell_features_file)
  tcell_CD4_TH2 = ConvertCellMarkers(tcell_CD4_TH2)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH17.csv'
  tcell_CD4_TH17 = read.csv(cell_features_file)
  tcell_CD4_TH17 = ConvertCellMarkers(tcell_CD4_TH17)
  
  cell_features_file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Cell Markers/tcell_CD4_TH117.csv'
  tcell_CD4_TH117 = read.csv(cell_features_file)
  tcell_CD4_TH117 = ConvertCellMarkers(tcell_CD4_TH117)
  
  
  
  cell_features = rbind(cell_features,
                        bcell_naive,
                        tcell_CD4_memory_TREG,tcell_CD4_naive, 
                        tcell_CD4_naive_activated,
                        tcell_CD4_naive_TREG, tcell_CD4_TFH,
                        tcell_CD4_TH1, tcell_CD4_TH2,
                        tcell_CD4_TH17, tcell_CD4_TH117)
  
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


