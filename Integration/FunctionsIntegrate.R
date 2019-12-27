run_pipeline_integrate = function(folder_input,folder_output,sample_name_list,sampleParam){
  print(sample_name_list)
  
  
  #cluster_IDs = sampleParam$Cluster_IDs_post_regress[sampleParam['Sample'] == sample_name_list[1]]
  # if (cluster_IDs == 'tmp'){
  #   file_str = ''
  # }else
  # {
  #   file_str = '_label'
  # }
  file_str = ''
  #folder_str = substr(sample_name_list[1], nchar(sample_name_list[1])-2+1, nchar(sample_name_list[1])) 
  # Load data
  filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
  metaData <- read_excel(filename_metaData)
  
  ########################
  ## Known Cell Markers
  ########################
  cell_features = getCellMarkers()
  
  # Pre
  print(paste0(folder_input[1],'data.Robj'))
  sample_name_pre = sample_name_list[1]
  data_pre <- loadRData(paste0(folder_input[1],'data.Robj'))
  data_pre@meta.data$data_pre_renamed <- data_pre@active.ident
  
  # Post
  print(paste0(folder_input[2],'data.Robj'))
  sample_name_post = sample_name_list[2]
  data_post = loadRData(paste0(folder_input[2],'data.Robj'))
  data_post@meta.data$data_post_renamed <- data_post@active.ident

  
  #Calculate integration anchors

  data_pre@meta.data$orig.ident <- "data_pre"
  data_post@meta.data$orig.ident <- "data_post"
  
  data_pre@meta.data$data_pre_renamed <- data_pre@active.ident
  data_post@meta.data$data_post_renamed <- data_post@active.ident
  
  
  anchors <- FindIntegrationAnchors(object.list = list(data_pre,data_post),k.filter =100)
  #Integrate samples
  data_integrated <- IntegrateData(anchorset = anchors)
  
  #Score for cell cycle genes
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  data_integrated <- CellCycleScoring(data_integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  data_integrated <- ScaleData(data_integrated, features = rownames(data_integrated))
  
  PCA_dim<- sampleParam$PCA_dim[sampleParam['Sample'] == sample_name_pre]
  print(paste0('PCA: ', PCA_dim))
  data_integrated <- RunPCA(data_integrated, features = VariableFeatures(object =data_integrated), npcs = PCA_dim)

  ######################
  ## Start Plotting
  ######################
  
  #Visualize PCA results
  visualize_PCA(data_integrated,folder_output[1],PCA_dim)
  
  print(folder_pre)
  pathName <- paste0(folder_output[1],'PCA/elbow',file_str,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data_integrated,ndims = 15))
  dev.off()
  
  # Jackstraw
  data_integrated <- JackStraw(data_integrated, num.replicate = 100)
  data_integrated <- ScoreJackStraw(data_integrated, dims = 1:PCA_dim)
  
  
  # Find Neighbors and Clusters
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name_list[1]]
  data_integrated <- getCluster (data_integrated,resolution_val, PCA_dim)
  
  # Find Cluster Biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(data_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  all_markers =  markers %>% group_by(cluster)
  

  
  all_markers$Cell = NA
  # Add known markers to all
  #browser()
  
  for (i in 1:nrow(all_markers)){
    marker_rows = grep(all_markers$gene[i], cell_features$Markers, value=TRUE)
    marker_idx = which( cell_features$Markers %in% marker_rows)
    #browser()
    if (length(marker_idx) > 0){
      #browser()
      cell_list = cell_features$Cell[marker_idx]
      cell_list = paste(cell_list, sep="", collapse=", ") 
      all_markers$Cell[i] = cell_list
    }
  }
  Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name_list[1]]

  write.csv(all_markers, file = paste0(folder_output[1],'AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)


  # Visualize clustering
  filepath_cluster = paste0( folder_input[1], 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data_integrated, reduction = "umap",label = TRUE,pt.size = 1))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data_integrated, features = top10$gene))
  dev.off()
  
  pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  print(FeaturePlot(data_integrated, features = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")))
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_split', file_str,'.png'))  
  png(file=pathName,width=600, height=350)
  print(DimPlot(data_integrated, label=T, repel=F, reduction = "umap", split.by = "orig.ident"))
  dev.off()
  
  return(data_integrated)
  
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
  - (sum(data_permute$active.ident[which(data_permute$orig.ident == "data_pre")] %in% cell_type)
     /sum(data_permute$orig.ident == "data_pre"))
  
  
  pval = (sum(abs(permuted_differences) > abs(obsdiff)) + 1) / (length(permuted_differences) + 1)
  return (c(obsdiff,pval))
  
}

diff_exp_pval_helper = function(data, folder_pre,sample_name,cell_type_list){
  
  print('Get Percent Differences')
  # Load and merge original data
  
  df_expr <- data.frame(matrix(ncol = 3, nrow = length(cell_type_list)))
  colnames(df_expr) <- c("cell_type","Percent_Diff", "P_Val")
  df_expr$cell_type = cell_type_list
  for (i in 1:length(cell_type_list))
  {
    diff_exp_result = diff_exp(data,cell_type_list[i] )
    print(paste0(cell_type_list[i], ' Percent Diff: ', diff_exp_result[1]))
    print(paste0(cell_type_list[i], ' Pval: ', diff_exp_result[2]))
    df_expr$'Percent_Diff'[i] = 100*diff_exp_result[1]
    df_expr$'P_Val'[i] = diff_exp_result[2]
    
    
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
    
    # Make sure there's enough cells
    if (pre_sum > 3 & post_sum > 3){
      
      # positive logFC values indicate that the feature is more highly expressed in post.

      cell_markers <- FindMarkers(data, ident.1 = paste0(cell_type_list[i],'_post'), ident.2 = paste0(cell_type_list[i],'_pre'), verbose = FALSE)
      #browser()
      cell_markers['Cell'] =cell_type_list[i]
      setDT(cell_markers, keep.rownames = 'Gene')
      
      
      cell_marker_list = rbind(cell_marker_list,cell_markers)
    }
    
    
  }
  return (cell_marker_list)
  
}

diffExpPermute = function(data,cell_features){
  perm = permutations(n = length(levels(data)), r = 2, v = 0:length(levels(data)))

  Features_mvn_output = data.frame(matrix(ncol = 9, nrow = 0))
  names(Features_mvn_output) = c("p_val","avg_logFC", "pct.1", "pct.2", "p_val_adj","gene", "Cell","ident_1","ident_2")
  Features_mvn_df = vector(mode = "list", length = length(levels(data)))
  
  level_list = levels(data)
  
  for (level in 1:length(levels(data)) ){
    Features_mvn_df[[level]] = Features_mvn_output
  }
  
  for (i in 1:nrow(perm)){
    
    ident1 = level_list[perm[i,1] + 1]
    ident2 = level_list[perm[i,2] + 1]
    
    print(perm[i,])
    Features_mvn = FindMarkers(data, ident.1 = ident1, ident.2 = ident2
                               ,min.pct = 0.2, logfc.threshold = 0.2, only.pos = TRUE)
    Features_mvn$gene = rownames(Features_mvn)
    rownames(Features_mvn) = NULL
    Features_mvn = cellMarkers(Features_mvn,cell_features)
    
    data_level = levels(data)
    Features_mvn$ident_1 = data_level[perm[i,1] + 1]
    Features_mvn$ident_2 = data_level[perm[i,2] + 1]
    
    ## Filter high p values
    Features_mvn = Features_mvn[Features_mvn$p_val_adj < 0.1,]
    
    Features_mvn_df[[perm[i,1] + 1]] = rbind(Features_mvn_df[[perm[i,1] + 1]], Features_mvn)
  }
  
  return (Features_mvn_df)
}

getStatSummary = function(patient_list, sample_Integrate_pairs,folder_base_input){
  for (i in 1:nrow(data)){
    patient = patient_list 
    print('')
    print(patient)
    pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
    
    sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
    folder_pre = makeFolders(folder_base_input,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
    
    file = paste0(folder_pre,'DE/DE_Patient',patient,'.csv')
    
    # TO DO
    output = 0
    
    return (output)
  }
}

