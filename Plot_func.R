
###################
# Plot Functions
###################

################################
## Quality Control
################################

quality_control <- function(data,folder,filter,nFeature_RNA_list,percent_mt,sample_name){
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  # Calculating percent Mitochondrial 
  #browser()
  print(nFeature_RNA_list)
  data[["percent.mt"]] = PercentageFeatureSet(data, pattern = "^MT-")
  
  if (filter == TRUE){
    
    print(nFeature_RNA_list)
    
    print(nFeature_RNA_list[2])
    nFeature_RNA_min = as.numeric(nFeature_RNA_list[1])
    nFeature_RNA_max = as.numeric(nFeature_RNA_list[2])
    #\browser()
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


#######################################
## PCA Results
# Examine and visualize PCA results a few different ways
# TO DO: Run PCA with 30 dim first, save images, then run with param values
#######################################
visualize_PCA = function(data,folder,PCA_dim, reduction = 'pca'){
  print('Visualize PCA')
  #print(data[["pca"]], dims = 1:5, nfeatures = 5)
  dir.create( paste0(folder,'PCA'), recursive = TRUE)
  pathName <- paste0(folder,'PCA/DimLoading.png')
  png(file=pathName,width=600, height=350)
  print(VizDimLoadings(data, dims = 1:2, reduction =reduction))
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
  
  for (i in 1:PCA_dim ){
    
    pathName <- paste0(folder,paste0('PCA/',reduction,'_','DimHeatMap_',i,'.png'))
    png(file=pathName,width=2000, height=1000, res=300)
    print(DimHeatmap(data, dims = i,reduction = reduction, cells = 500, balanced = TRUE))
    dev.off()
  }
  
  
  ## ADD DOTPLOT
  ## TO DO: Make clusters labeled by cell cycle
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  #data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  print(PCA_dim)
  for (x in 1:(PCA_dim -1)){
    
    y <- x+1
    pathName <- paste0(folder,'PCA/',reduction,x,'_',y,'.png')
    png(file=pathName,width=600, height=350)
    print(DimPlot(data, dims = c(x,y), reduction = reduction,pt.size = 2))
    dev.off()
  }
  
  
}

#########################
## Make FeaturePlot
#########################
PlotKnownMarkers = function(data, 
                            folder, 
                            cell_features,
                            plotType ='FeaturePlot',
                            str = '', split_group = NA, plotTogether = T){
  
  #browser()
  
  all_markers  = rownames(data@assays[["RNA"]])
  #all_markers = rownames(data)
  
  #browser()
  dir.create( folder, recursive = TRUE)
  
  feature_list = as.character(cell_features$Markers)
  
  if ( plotType == 'FeaturePlotFix' & plotTogether == T){
    gene_list = c()
    for(i in 1:length(feature_list) ){
      x = unlist(strsplit(feature_list[i], ",")) 
      #x = gsub("\\s", "", x)  whitespace <- " \t\n\r\v\f"
      x = trimws(x, which = c("both"), whitespace = " \t\n\r\v\f")
      x = gsub("\\W", "", x)
      gene_list = c(gene_list,x)
    }
    
    gene_list = (gene_list[gene_list %in% all_markers])
    gene_list = unique(gene_list)
    
    for (j in 1:length(gene_list)){
      gene = gene_list[j]
      print(gene)
      #browser()
      plot = FeaturePlotFix(data, feature = gene,folder = '',str = '', markerSize = 1,
                            split = F, gene_TF = TRUE,title = '',saveTF = FALSE) 
      pathName = paste0(folder,'gene','_',gene,str,'.png')
      png(filename = pathName,width=1000, height=1000, res=100)
      print(plot)
      dev.off()
      remove(plot)
      
    }
    return 
  }
  
  for(i in 1:length(feature_list) ){
    #browser()
    cell_type = cell_features$Cell[i]
    
    
    x = unlist(strsplit(feature_list[i], ",")) 
    #x = gsub("\\s", "", x)  whitespace <- " \t\n\r\v\f"
    x = trimws(x, which = c("both"), whitespace = " \t\n\r\v\f")
    x = gsub("\\W", "", x)


    gene_list = (x[x %in% all_markers])  
    # FeaturePlot gives error if no genes are found, so we have to check if they are included in the highly variable gene list
    gene_list = unique(gene_list)
      
    if (plotType == 'HeatMap' & length(gene_list) > 0){
      #browser()
      print(cell_type)
      
      plot  = DoHeatmap(object = data, features = gene_list,assay = 'RNA', slot = "data",
                      group.by = "ident", label = T) +
                      ggtitle('' ) + 
                      theme(plot.title = element_text(hjust = 0.5)) +
                      theme(plot.title = element_text(size=24))
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        text = element_text(size = 20))
      
      
      
      pathName = paste0(folder,cell_type,str,'_heatmap','.png')
      if (length(gene_list) < 10){
        height = 1000
      }else{
        height = 40*length(gene_list)
      }
      png(file=pathName,width=3500, height=height)
      print(plot)
      dev.off()
      
    }else if (plotType == 'Violin' & length(gene_list) > 0){
      #browser()
      print(cell_type)
      plot  = StackedVlnPlot(obj = data, features = gene_list) +
        ggtitle('' ) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(size=24))
      plot = plot + theme(
        axis.text= element_text(color="black", size=24))
      
      pathName = paste0(folder,cell_type,str,'_violin','.png')
      
      png(file=pathName,width=3500, height=1000,res=100)
      print(plot)
      dev.off()
      
    }else if (length(gene_list) > 0 & plotType == 'FeaturePlotFix'){
      #browser()
      print(paste0( cell_type,': Found'))
      plot_list = vector('list',length(gene_list))
      for (j in 1:length(gene_list)){
        gene = gene_list[j]
        print(gene)
        #browser()
        plot = FeaturePlotFix(data, feature = gene,folder = '',str = '', markerSize = 1,
                              split = F, gene_TF = TRUE,title = '',saveTF = FALSE) 
        pathName = paste0(folder,cell_type,'_',gene,str,'.png')
        png(filename = pathName,width=1000, height=1000, res=100)
        print(plot)
        dev.off()
        remove(plot)
        
        
        }
      }else {
      print(paste0( cell_type,': Not Found'))
    }
    
    ##
    
   
    
    print('')
    
  } 
  ####################################  

  
  
  
  
  
  
}

FeaturePlotFix = function(data, feature,folder,str, split,markerSize = 2, gene_TF,title = '',saveTF = TRUE){
  #browser()
  dir.create(folder,recursive = TRUE)
  
  data_umap = FetchData(data, vars = c("ident", "orig.ident","UMAP_1", "UMAP_2"))
  data_umap = FetchData(data, vars = c("ident", "orig.ident","UMAP_1", "UMAP_2"))
  if (data@active.assay == 'RNA'){
    data_df = data.frame(data@assays$RNA@data) 
  }else if (data@active.assay == 'integrated'){
    data_df = data.frame(data@assays$RNA@data) 
    
  }else if (data@active.assay == 'SCT'){
    data_df = data.frame(data@assays$SCT@data) 
  }
  if (gene_TF){
    gene <- data_df[rownames(data_df) == feature,]
    gene <- data.frame(t(gene))
    if (ncol(gene) > 0){
      data_umap$g = gene[,1]
    }
  }else{
    
    data_matrix = data[[]]
    data_matrix = select(data_matrix, feature)
    data_umap$g = data_matrix[,1]
  }
  #browser()
  if (split){

    
    
  }else{
   
    color_low = 'grey'
    fontSize = 24
    
    g_post = ggplot(data_umap, aes(UMAP_1, UMAP_2)) +  
      geom_point(aes(colour = data_umap$g), size = markerSize) +
      scale_color_gradient(low = color_low, high = "blue", name = "") +
      #geom_text_repel(data=umap_label_pre, aes(label=ident, x, y))  +
      ggtitle(title,  subtitle = feature) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.subtitle = element_text(hjust = 0.5))+
      theme(plot.title = element_text(size=fontSize)) + 
      theme(plot.subtitle = element_text(size=fontSize)) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #browser()
    if (saveTF){
      pathName <- paste0(folder,str,'_',feature,'.png')
      png(file=pathName, width=2000, height=1000, res=100)
      print(g_post)
      dev.off()
    }
    
    return (g_post)
  }
  #browser()
  
}


FeaturePlotFix_Scran = function(data, feature,folder,str, split,markerSize = 2, gene_TF,title = '',saveTF = TRUE){
  #browser()
  dir.create(folder,recursive = TRUE)
  
  data_umap = FetchData(data, vars = c("ident", "orig.ident","umap_1", "umap_2"))
  if (data@active.assay == 'RNA'){
    data_df = data.frame(data@assays$RNA@data) 
  }else if (data@active.assay == 'integrated'){
    data_df = data.frame(data@assays$RNA@data) 
    
  }else if (data@active.assay == 'SCT'){
    data_df = data.frame(data@assays$SCT@data) 
  }
  if (gene_TF){
    gene <- data_df[rownames(data_df) == feature,]
    gene <- data.frame(t(gene))
    if (ncol(gene) > 0){
      data_umap$g = gene[,1]
    }
  }else{
    
    data_matrix = data[[]]
    data_matrix = select(data_matrix, feature)
    data_umap$g = data_matrix[,1]
  }
  #browser()
  if (split){
    
    
    
  }else{
    
    color_low = 'grey'
    fontSize = 24
    
    g_post = ggplot(data_umap, aes(umap_1, umap_2)) +  
      geom_point(aes(colour = data_umap$g), size = markerSize) +
      scale_color_gradient(low = color_low, high = "blue", name = "") +
      #geom_text_repel(data=umap_label_pre, aes(label=ident, x, y))  +
      ggtitle(title,  subtitle = feature) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.subtitle = element_text(hjust = 0.5))+
      theme(plot.title = element_text(size=fontSize)) + 
      theme(plot.subtitle = element_text(size=fontSize)) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    #browser()
    if (saveTF){
      pathName <- paste0(folder,str,'_',feature,'.png')
      png(file=pathName, width=2000, height=1000, res=100)
      print(g_post)
      dev.off()
    }
    
    return (g_post)
  }
  #browser()
  
}


plotHLA = function(data,folder_base_output,str){
  #################
  ## dCD14+ genes
  #################
  #browser()
  folder_output_analysis =  paste0( folder_base_output, 'Analysis/FeaturePlot_HLA/',str,'/' )
  dir.create( folder_output_analysis, recursive = TRUE)
  
  Features_dCD14_MonoVsn_label_Patient10 = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/GL1080BM/Filtered/Regress/DE/nVsm/Features_dCD14+ MonoVsn_label_Patient10.csv'    
  Features_dCD14_MonoVsn_label_Patient10  = read.csv(Features_dCD14_MonoVsn_label_Patient10)
  Features_dCD14_MonoVsn_label_Patient10 = Features_dCD14_MonoVsn_label_Patient10[Features_dCD14_MonoVsn_label_Patient10$ident_2 == 'CD14+ Mono', ]
  Features_dCD14_MonoVsn_label_Patient10 = Features_dCD14_MonoVsn_label_Patient10[Features_dCD14_MonoVsn_label_Patient10$p_val_adj < 0.05, ]
  Features_dCD14_MonoVsn_label_Patient10 = levels(Features_dCD14_MonoVsn_label_Patient10$gene)
  
  DEgeneHLA = c(Features_dCD14_MonoVsn_label_Patient10)
  DEgeneHLA = DEgeneHLA[grep('HLA', DEgeneHLA)]
  
  for (i in 1:(length(DEgeneHLA))){
    feature = as.character(DEgeneHLA[i])
    print(feature)
    FeaturePlotFix(data, feature,folder_output_analysis,str, split = FALSE, gene_TF = TRUE)
  }
  
  data_score = AddModuleScore(object = data, features = list(DEgeneHLA),nbin = 25, ctrl = 7, name = 'HLA_Feature')
  score_HLA_orig = data_score[[]]
  score_HLA = subset(score_HLA_orig, select = c(orig.ident,HLA_Feature1))
  FeaturePlotFix(data_score, 'HLA_Feature1',folder_output_analysis,'DexaT', split = FALSE, gene_TF = FALSE)
  
}
###################################################################

DoHeatMapHelper = function(data_run,folder_base_output,folder_heatMap = NA, Features = NA,
                           ident1,ident2,celltype_list,category_list, split_var = 'split_var'
                           ,sortby,cellPair,gene_num,str = '',Feature_str = '', Ident_order = NA){
  
  #browser()
  print(ident1)
  print(ident2)
  print(unique(cluster_IDs))
  #browser()
  folder_name = paste0('DoHeatmap/',ident1, '_', ident2)
  print(folder_name)
  
  
  if (is.nan(folder_heatMap)){
    
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    
  }
  if (sortby == 'avg_logFC'){
    
    Features = Features[Features$p_val_adj < 0.05,]
    Features = Features[order(Features$avg_logFC),]
    gene_list = rownames(Features)
    gene_list = c(gene_list[1:gene_num[1]], tail(gene_list, n=gene_num[2]))
    gene_list = unique(gene_list)
  }else if(sortby == 'p_val_adj') {
    Features = Features[Features$p_val_adj < 0.05,]
    Features = Features[order(Features$p_val_adj),]
    gene_list = rownames(Features)
    gene_list = gene_list[1:60]
  }
  
  folder_heatMap = gsub("CD14+ Mono","CD14M",folder_heatMap,fixed = T)
  folder_heatMap = gsub("CD16+ Mono","CD16M",folder_heatMap,fixed = T)
  folder_heatMap  = gsub('CD8+ T Cell','CD8T',folder_heatMap, fixed = T)
  folder_heatMap = gsub("baseline","B",folder_heatMap,fixed = T)
  folder_heatMap = gsub("Good Response","GR",folder_heatMap,fixed = T)
  folder_heatMap = gsub("Poor Response","PR",folder_heatMap,fixed = T)
  
  
  
  pathName <- paste0(folder_heatMap,paste0('Markers','.csv')) 
  dir.create( folder_heatMap, recursive = TRUE)
  print(pathName)
  write.csv(Features, file = pathName,row.names = TRUE)
  print('gene_list')
  print(gene_list)
  if (all(is.na(gene_list))){
    print('No genes found')
    return()
  }
  
  data_run_label = label_cells(data_run, cluster_IDs)
  data_run_label_subset_cell = subset(data_run_label, idents = celltype_list)
  Idents(data_run_label_subset_cell) = paste0(Idents(data_run_label_subset_cell),' ', data_run_label_subset_cell@meta.data[,split_var])
  #Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run_label@meta.data[,split_var])
  
  #data_run_label = label_cells(data_run, cluster_IDs)
  celltype_list = c(ident1, ident2)
  ident1 = gsub("CD14+ Mono","CD14M",ident1,fixed = T)
  ident1  = gsub('CD8+ T Cell','CD8T',ident1, fixed = T)
  ident1  = gsub('CD16+ Mono','CD16M',ident1, fixed = T)
  ident2 = gsub('CD14+ Mono','CD14M',ident2, fixed = T)
  ident2  = gsub('CD8+ T Cell','CD8T',ident2, fixed = T)
  ident2  = gsub('CD16+ Mono','CD16M',ident2, fixed = T)
  #browser()
  
  
  if (!cellPair){
    if (sum( Idents(data_run_label_subset_cell) ==celltype_list[1] ) > 3 
        && sum(Idents(data_run_label_subset_cell) == celltype_list[2] ) > 3){
      
      data_run_label_subset = subset(data_run_label_subset_cell, idents = celltype_list)
      #browser()
      #Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run_label@meta.data[,split_var])
      filename = paste0(ident1, '_', ident2,' Split_',split_var,str,' Subset','.png')
      filename = gsub("Baseline","B",filename,fixed = T)
      filename = gsub("Good Response","GR",filename,fixed = T)
      filename = gsub("Poor Response","PR",filename,fixed = T)
      pathName = as.character(paste0(folder_heatMap,filename) )
      
      
      
      print(pathName)
      #browser()
      plot  = DoHeatmap(object = data_run_label_subset, features = gene_list,
                        group.by = "ident",label = T) +
        ggtitle('' ) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(size=24))
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=18 ),
        axis.title.y = element_text(color="black", size=18),
        axis.text= element_text(color="black", size=18),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        text = element_text(size = 20)
        
      )
      
      png(file=pathName,width=1000, height=1500)
      print(plot)
      dev.off()
      
    }else{
      print('Cell type not found')
      #return()
    }
    
    
    
  }else{
    
    #Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run_label@meta.data[,split_var])
    #data_run_label = subset(data_run_label, idents = celltype_list)
    
  }
  #browser()
  filename = paste0(ident1, '_', ident2,' Split_',split_var,str,'.png')
  filename = gsub("baseline","B",filename,fixed = T)
  filename = gsub("Good Response","GR",filename,fixed = T)
  filename = gsub("Poor Response","PR",filename,fixed = T)
  pathName = as.character(paste0(folder_heatMap,filename) )
  #pathName = paste0(folder_heatMap,paste0(ident1, '_', ident2,' Split_',split_var,str,' 2','.png')) 
  print(pathName)
  #browser()
  
  if (!is.na(Ident_order)){
    Idents(data_run_label_subset_cell) = factor(Idents(data_run_label_subset_cell) , levels = Ident_order)
    #data_run_label_subset_cell = subset(data_run_label_subset_cell, idents = Ident_order)
  }else{
    
  }
  
  plot  = DoHeatmap(object = data_run_label_subset_cell, features = gene_list,
                    group.by = "ident",label = T) +
    ggtitle('' ) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=24))
  
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=18 ),
    axis.title.y = element_text(color="black", size=18),
    axis.text= element_text(color="black", size=18),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20)
    
  )
  png(file=pathName,width=1000, height=1500)
  print(plot)
  dev.off()
  
  #browser()
  
  
  
}

###################

CompareCellNum = function(data,folder,split_var,metaData){
  browser()
  
  
  dir.create( folder, recursive = TRUE)
  
  category_list = unique(data@meta.data[,split_var])
  
  cell_type_list = sort(as.character(unique(Idents(data))))
  stats_category = data.frame(matrix(ncol = length(cell_type_list) + 1, nrow = 0))
  colnames(stats_category) = c(cell_type_list, 'category')
  
  cnt = 1
  for (category in category_list){ 
    
    data_subset_category = SubsetData(object = data, cells = data@meta.data[,split_var] == category )
    sample_list = sort(unique(data_subset_category@meta.data[,'sample_name'] ))
    #cell_type_list = sort(as.character(unique(Idents(data_subset_category))))
    
    stats_df = data.frame(matrix(ncol = length(cell_type_list) + 1, nrow = length(sample_list)))
    colnames(stats_df) = c(cell_type_list, 'category')
    rownames(stats_df) = sample_list
    #browser()
    
    for (sample in sample_list){
      
      #browser()
      data_subset_sample = SubsetData(object = data_subset_category, cells = data_subset_category$sample_name == sample )
      
      stats = clusterStats(data_subset_sample)
      stats_df[sample,1:length(cell_type_list)] = stats$Percent
      stats_df$category = category
      #stats_df$cell_num = length(data_subset_sample@assays[["integrated"]]@data@p)
      #browser()
    }
    #browser()
    stats_category = rbind(stats_category,stats_df)
    cnt = cnt + 1
  }
  
  browser()
  cell_list = c("T Cell", "CD8+ T Cell", "HSC","NK","Mature B Cell", 
                "Immature B Cell", "CD16+ Mono", "pDC",'Erythrocyte','CD14+CD16+ Mono','Pre B Cell','DC')
  cell_percent = c(100, 50, 25, 30,20,5,10,5,20,10, 10,15)
  cell_type_list = unique(Idents(data))
  for (cell in cell_type_list){
    #browser()
    print(cell)
    #browser()
    if (length(category_list) == 3){
      if (max(stats_category[,cell]) > 50 ||max(stats_category[,cell])  > 50 || max(stats_category[,cell])  > 50){
        max_val = 100
      }else if (max(stats_category[[1]][,cell]) > 30 ||max(stats_category[[2]][,cell])  > 30 || max(stats_category[[3]][,cell])  > 30){
        max_val = 50
      }else {
        max_val = 30
      }
      if (cell %in% cell_list){
        names(cell_percent) <- cell_list
        
        max_val = cell_percent[[cell]]
      }
      print(max_val)
      
      #browser()
      
      stats_summary = rbind(stats_category[[1]][, c(cell,'category')],stats_category[[2]][, c(cell,'category')],stats_category[[3]][, c(cell,'category')])
      stats_summary$category = factor(stats_summary$category , levels = category_list)
      stats_summary$Sample = rownames(stats_summary)
      stats_summary =merge(stats_summary, metaData, by = 'Sample')
      #browser()
      stats_summary$Response = gsub("Minimal Response then Progression", "TMP", stats_summary$Response )
      stats_summary$Response = gsub("Minimal Response", "TMP", stats_summary$Response )
      stats_summary$Response = gsub("TMP", "MR",stats_summary$Response )
      
      stats_summary$splitvar_response = paste0(stats_summary$category, ' ', stats_summary$Response)
      stats_summary$splitvar_response = gsub("NBM NBM", "NBM", stats_summary$splitvar_response)
      splitvar_response_factor = c('NBM', 'Pre MR','Post MR','Pre VGPR','Post VGPR')
      splitvar_response_factor = unique(stats_summary$splitvar_response)
      stats_summary$splitvar_response = factor(stats_summary$splitvar_response ,splitvar_response_factor)
      
      stats_summary$response_dexa = paste0(stats_summary$Response, ' ', stats_summary$"Dexa or not", ' Dexa')
      stats_summary$response_dexa = gsub("NBM NBM Dexa", "NBM", stats_summary$response_dexa)
      
      ###############
      #browser()
      # Just pre and post
      category1 = 'C9D1'
      category2 = 'EOT'
      stats_summary_line = stats_summary[(stats_summary$category == category1 | stats_summary$category == category2),]
      stats_summary_line$category = as.character(stats_summary_line$category)
      category_levels = c('Pre', 'Post')
      category_levels = unique(stats_summary_line$category)
      stats_summary_line$category = factor(stats_summary_line$category, levels  = category_levels)
      stats_summary_line$Patient = stats_summary_line$`Patient Number`
      
      # With lines
      pathName =  paste0(folder,'Boxplot Lines/','boxplot_', cell,'_lines','.png')
      png(file=pathName,width=600, height=600)
      plot = ggplot(stats_summary_line, aes(x = category, y = !!ensym(cell))) +
        geom_boxplot() +
        coord_cartesian(ylim = c(0, max_val))+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        geom_point(color="black", size=2) +
        geom_line(aes(group=Patient, color=stats_summary_line$Response,alpha = Response ), size = 1) +
        #scale_colour_manual(values=c(category1="blue",category2="red"))+
        scale_alpha_manual(values = c( 0.9, 0.9))+
        guides(alpha = FALSE)+
        labs(color = "Response")+
        theme_classic()
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
        
      )
      
      print(plot)
      dev.off()
      # Without Lines
      
      pathName =  paste0(folder,'Boxplot/','boxplot_prepost', cell,'','.png')
      png(file=pathName,width=600, height=600)
      plot = ggplot(stats_summary_line, aes(x = category, y = !!ensym(cell))) +
        geom_boxplot() +
        coord_cartesian(ylim = c(0, max_val))+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        geom_point(color="black", size=2) +
        labs(color = "Response")+
        theme_classic()
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
        
      )
      
      print(plot)
      dev.off()
      
      ##########################
      #browser()
      # Split by response
      stats_summary_line$splitvar_response = paste0(stats_summary_line$category, ' ', stats_summary_line$Response)
      splitvar_response_factor = c('C9D1 Poor Response', 'EOT Poor Response','C9D1 Good Response', 'EOT Good Response')
      stats_summary_line$splitvar_response = factor(stats_summary_line$splitvar_response ,splitvar_response_factor)
      
      pathName =  paste0(folder,'Boxplot Lines Response/','boxplot_splitvar_response_', cell,'_lines','.png')
      png(file=pathName,width=1200, height=600)
      plot = ggplot(stats_summary_line, aes(x = splitvar_response, y = !!ensym(cell))) +
        geom_boxplot() +
        coord_cartesian(ylim = c(0, max_val))+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        geom_point(color="black", size=2) +
        geom_line(aes(group=Patient, color=stats_summary_line$Response,alpha = Response), size = 1) +
        #scale_colour_manual(values=c(MR="blue",VGPR="red"))+
        scale_alpha_manual(values = c(0.9,0.9))+
        guides(alpha = FALSE)+
        labs(color = "Response")+
        theme_classic()
      
      plot = plot + theme(
        #plot.title = element_text(hjust = 0.5,color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=20),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
      )
      print(plot)
      dev.off()
      
      category1 = 'baseline'
      category2 = 'C9D1'
      stats_summary$Patient =stats_summary$`Patient Number`
      stats_summary_line = stats_summary[(stats_summary$category == category1 | stats_summary$category == category2),]
      stats_summary_line$splitvar_response = paste0(stats_summary_line$category, ' ', stats_summary_line$Response)
      splitvar_response_factor = c('baseline Poor Response', 'C9D1 Poor Response','baseline Good Response','C9D1 Good Response')
      stats_summary_line$splitvar_response = factor(stats_summary_line$splitvar_response ,splitvar_response_factor)
      
      pathName =  paste0(folder,'Boxplot Lines Response Baseline/','boxplot_splitvar_response_', cell,'_lines','.png')
      png(file=pathName,width=1200, height=600)
      plot = ggplot(stats_summary_line, aes(x = splitvar_response, y = !!ensym(cell))) +
        geom_boxplot() +
        coord_cartesian(ylim = c(0, max_val))+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        geom_point(color="black", size=2) +
        geom_line(aes(group=Patient, color=stats_summary_line$Response,alpha = Response), size = 1) +
        #scale_colour_manual(values=c(MR="blue",VGPR="red"))+
        scale_alpha_manual(values = c(0.9,0.9))+
        guides(alpha = FALSE)+
        labs(color = "Response")+
        theme_classic()
      
      plot = plot + theme(
        #plot.title = element_text(hjust = 0.5,color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=20),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
      )
      print(plot)
      dev.off()
      
      #######################
      #browser()
      # Split by response with dexa color
      pathName =  paste0(folder,'Boxplot Lines Response/','boxplot_splitvar_response_dexa_', cell,'_lines','.png')
      png(file=pathName,width=900, height=600)
      plot = ggplot(stats_summary_line, aes(x = splitvar_response, y = !!ensym(cell))) +
        geom_boxplot() +
        coord_cartesian(ylim = c(0, max_val))+
        ggtitle(cell)+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        #geom_point(colour="black", size=2) +
        #geom_line(aes(group=Patient, colour=stats_summary_line$"Dexa or not"),alpha = 0.4, size = 1) +
        guides(alpha = FALSE)+
        labs(color = "Dexa")+
        theme_classic()
      
      plot = plot + theme(
        #plot.title = element_text(hjust = 0.5,color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
      )
      print(plot)
      dev.off()
      
      
      
      ############
      pathName =  paste0(folder,'Boxplot/','boxplot_', cell,'.png')
      print(pathName)
      x_name = 'category'
      png(file=pathName,width=600, height=600)
      plot = ggplot(stats_summary, aes(x = !!ensym(x_name), y = !!ensym(cell))) + 
        geom_boxplot()+
        coord_cartesian(ylim = c(0, max_val))+
        ggtitle('')+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        theme_classic()
      
      plot = plot + geom_jitter(shape=16, position=position_jitter(0.2))
      
      plot = plot + theme(
        #plot.title = element_text(hjust = 0.5,color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
      )
      
      #plot = plot + theme(plot.title = element_text(hjust = 0.5))
      print(plot)
      dev.off()
      
      ######################
      #browser()
      
      stats_summary_input = stats_summary[stats_summary$splitvar_response!= 'NBM',]
      pathName =  paste0(folder,'Boxplot Lines Response/','boxplot_', cell,'SplitByResponse','.png')
      x_name = 'splitvar_response'
      png(file=pathName,width=900, height=600)
      plot = ggplot(stats_summary_input, aes(x = !!ensym(x_name), y = !!ensym(cell))) + 
        geom_boxplot()+
        coord_cartesian(ylim = c(0, max_val))+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        theme_classic()
      # Box plot with dot plot
      plot = plot + geom_jitter(shape=16, position=position_jitter(0.2), color="black", size=2)
      
      plot = plot + theme(
        #plot.title = element_text(hjust = 0.5,color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
      )
      
      plot = plot + theme(plot.title = element_text(hjust = 0.5))
      print(plot)
      dev.off()
      # ########################################
      
      
      
      # Pre = stats_summary[stats_summary$category == 'Pre',][,cell]
      # Post = stats_summary[stats_summary$category == 'Post',][,cell]
      # NBM = stats_summary[stats_summary$category == 'NBM',][,cell]
      # 
      # wilcox_PrePost = wilcox.test(x = Pre,y =  Post, paired = FALSE)
      # wilcox_NBMPre = wilcox.test(x = NBM,y =  Pre, paired = FALSE)
      # wilcox_NBMPost = wilcox.test(x = NBM,y =  Post, paired = FALSE)
      # 
      # wilcox_output = data.frame(matrix(ncol = 2, nrow = 3))
      # colnames(wilcox_output) = c('Compare', 'Pval')
      # wilcox_output[1,] = c('Pre Post',  wilcox_PrePost$p.value)
      # wilcox_output[2,] = c('NBM Pre',  wilcox_NBMPre$p.value)
      # wilcox_output[3,] = c('NBM Post',  wilcox_NBMPost$p.value)
      
      
      #write.csv(wilcox_output, file = paste0(folder,'wilcox_PrePost_',cell,'.csv'),row.names = FALSE)
    }else if (length(category_list) == 2){
      #browser()
      if (max(stats_category[[1]][,cell]) > 50 ||max(stats_category[[2]][,cell])  > 50){
        max_val = 100
      }else{
        max_val = 50
      }
      
      stats_summary = rbind(stats_category[[1]][, c(cell,'category')],stats_category[[2]][, c(cell,'category')])
      stats_summary$category = factor(stats_summary$category , levels = c('Pre','Post'))
      #stats_summary =merge(stats_summary, sampleParam, by = , sort = TRUE)
      #browser()
      
      pathName =  paste0(folder,'boxplot_', cell,'.png')
      x_name = 'category'
      png(file=pathName,width=600, height=600)
      plot = ggplot(stats_summary, aes(x = !!ensym(x_name), y = !!ensym(cell))) + 
        geom_boxplot()+
        coord_cartesian(ylim = c(0, max_val))+
        ggtitle(cell)+
        xlab("") + ylab("%")
      
      plot = plot + theme(
        plot.title = element_text(color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24, face="bold"),
        axis.title.y = element_text(color="black", size=24, face="bold"),
        axis.text=element_text(size=24)
      )
      
      plot = plot + theme(plot.title = element_text(hjust = 0.5))
      
      print(plot)
      dev.off()
      
      #browser()
      pathName =  paste0(folder,'boxplot_', cell,'_diff','.png')
      x_name = 'category'
      png(file=pathName,width=600, height=600)
      plot = ggplot(stats_summary, aes(x = !!ensym(x_name), y = !!ensym(cell))) + 
        geom_boxplot()+
        coord_cartesian(ylim = c(0, max_val))+
        ggtitle(cell)+
        xlab("") + ylab("%")
      
      plot = plot + theme(
        plot.title = element_text(color="black", size=24, face="bold.italic"),
        axis.title.x = element_text(color="black", size=24, face="bold"),
        axis.title.y = element_text(color="black", size=24, face="bold"),
        axis.text=element_text(size=24)
      )
      
      plot = plot + theme(plot.title = element_text(hjust = 0.5))
      
      print(plot)
      dev.off()
    }else if (length(category_list) == 4){
      if (max(stats_category[,cell]) > 50){
        max_val = 100
      }else if (max(stats_category[,cell]) > 30){
        max_val = 50
      }else {
        max_val = 30
      }
      
      if (cell %in% cell_list){
        
        names(cell_percent) <- cell_list
        
        max_val = cell_percent[[cell]]
      }
      print(max_val)
      
      
      
      stats_summary = stats_category[, c(cell,'category')]
      
      stats_summary$Sample = rownames(stats_summary)
      stats_summary =merge(stats_summary, metaData, by = 'Sample')
      #browser()
      stats_summary$Response = gsub("Minimal Response then Progression", "TMP", stats_summary$Response )
      stats_summary$Response = gsub("Minimal Response", "TMP", stats_summary$Response )
      stats_summary$Response = gsub("TMP", "MR",stats_summary$Response )
      
      stats_summary$splitvar_response = paste0(stats_summary$category, ' ', stats_summary$Response)
      stats_summary$splitvar_response = gsub("NBM NBM", "NBM", stats_summary$splitvar_response)
      splitvar_response_factor = c('NBM', 'baseline Poor Response','C9D1 Poor Response','EOT Poor Response',
                                   'baseline Good Response','C9D1 Good Response','EOT Good Response')
      splitvar_response_factor = unique(stats_summary$splitvar_response)
      stats_summary$splitvar_response = factor(stats_summary$splitvar_response ,splitvar_response_factor)
      
      stats_summary$response_dexa = paste0(stats_summary$Response, ' ', stats_summary$"Dexa or not", ' Dexa')
      stats_summary$response_dexa = gsub("NBM NBM Dexa", "NBM", stats_summary$response_dexa)
      
      ###############
      #browser()
      # Just pre and post
      category1 = 'baseline'
      category2 = 'C9D1'
      category3 = 'EOT'
      stats_summary_line = stats_summary[(stats_summary$category == category1 | 
                                            stats_summary$category == category2 | 
                                            stats_summary$category == category3 ),]
      stats_summary_line$category = as.character(stats_summary_line$category)
      category_levels = c('baseline', 'C9D1','EOT')
      category_levels = unique(stats_summary_line$category)
      stats_summary_line$category = factor(stats_summary_line$category, levels  = category_levels)
      stats_summary_line$Patient = stats_summary_line$`Patient Number`
      
      # With lines
      pathName =  paste0(folder,'Boxplot Lines/','boxplot_', cell,'_lines','.png')
      png(file=pathName,width=800, height=600)
      plot = ggplot(stats_summary_line, aes(x = category, y = !!ensym(cell))) +
        geom_boxplot() +
        coord_cartesian(ylim = c(0, max_val))+
        xlab("") + ylab(paste0(cell, ' Proportion'))+
        geom_point(color="black", size=2) +
        geom_line(aes(group=Patient, color=stats_summary_line$Response,alpha = Response ), size = 1) +
        #scale_colour_manual(values=c(category1="blue",category2="red"))+
        scale_alpha_manual(values = c( 0.9, 0.9,0.9))+
        guides(alpha = FALSE)+
        labs(color = "Response")+
        theme_classic()
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18)
        
      )
      
      print(plot)
      dev.off()
    }
    
  }
  
  #browser()
  patient_list = unique((data$Patient))
  patient_list = c()
  # Split by Patient
  for (patient in patient_list){
    print(patient)
    max_val = 60
    stats_summary = rbind(stats_category[[1]],stats_category[[2]],stats_category[[3]])
    stats_summary$category = factor(stats_summary$category , levels = c('NBM','Pre','Post'))
    stats_summary$Sample = rownames(stats_summary)
    stats_summary =merge(stats_summary, metaData, by = 'Sample')
    stats_summary_patient = stats_summary[stats_summary$`Patient Number` == patient,]
    
    stats_summary_trans = data.frame(r1=names(stats_summary_patient), t(stats_summary_patient))
    stats_summary_trans_cells = stats_summary_trans[stats_summary_trans$r1 %in% cell_type_list,]
    
    cat = stats_summary_trans[stats_summary_trans$r1 == 'category',]
    stats_summary_trans_cells$category = cat[1,2]
    stats_summary_trans_cat1 = stats_summary_trans_cells[,c(2,4)]
    stats_summary_trans_cat1$patient = 1
    stats_summary_trans_cat1$cell = rownames(stats_summary_trans_cat1)
    stats_summary_trans_cells$category = cat[1,3]
    stats_summary_trans_cat2 = stats_summary_trans_cells[,c(3,4)]
    stats_summary_trans_cat2$patient = 2
    stats_summary_trans_cat2$cell = rownames(stats_summary_trans_cat2)
    colnames(stats_summary_trans_cat1) = c('percent','category','patient','cell')
    colnames(stats_summary_trans_cat2) = c('percent','category','patient','cell')
    
    stats_summary_final = rbind(stats_summary_trans_cat1, stats_summary_trans_cat2)
    stats_summary_final$percent = as.numeric(levels(stats_summary_final$percent))[stats_summary_final$percent]
    
    folder_name = paste0(folder,'Boxplot Lines Patient/')
    dir.create( folder_name, recursive = TRUE)
    pathName =  paste0(folder_name,'boxplot_lines_patient', patient,'_',stats_summary_patient$Response[1],'','.png')
    png(file=pathName,width=600, height=800)
    plot = ggplot(stats_summary_final, aes(x = category, y = percent)) +
      coord_cartesian(ylim = c(0, max_val))+
      ggtitle(paste0('Patient ', patient, ' ', stats_summary_patient$Response[1]))+
      xlab("") + ylab("%")+
      geom_point(colour="black", size=2) +
      geom_text(aes(label=cell),hjust=0, vjust=0) +
      geom_line(aes(group=cell, colour=stats_summary_final$"cell"),alpha = 0.4, size = 1.5) +
      theme_classic()
    plot = plot + theme(plot.title = element_text(hjust = 0.5))
    print(plot)
    dev.off()
  }
  
  browser()
  
  
}


TopNHeatmap = function(data,markers,filepath_cluster,PCA_dim,resolution_val,num_markers,file_str){
  browser()
  topN = markers %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
  topN$gene = as.character(topN$gene)
  plot = DoHeatmap(data, features = topN$gene)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text=element_text(size=18)
  )
  
  #browser()
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=1000, height=1500)
  print(plot)
  dev.off()
  
  return()
}


modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = 24, angle = 0), 
          axis.text.y = element_text(size = 24), 
          plot.margin = plot.margin, ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(size=24), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


StackedVlnPlotHelper = function(data,gene_list,folder_heatMap,filename){
  #browser()
  gene_list = gene_list[!is.na(gene_list)]
  if (length(gene_list) > 0){
    plot  = StackedVlnPlot(obj = data, features = gene_list ) +
      ggtitle('' ) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=24))
    plot = plot + theme(
      axis.text= element_text(color="black", size=32))
    
    
    folder_heatMap = gsub('CD14+ Mono','CD14M',folder_heatMap, fixed = T)
    folder_heatMap  = gsub('CD8+ T Cell','CD8T',folder_heatMap, fixed = T)
    folder_heatMap  = gsub('CD16+ Mono','CD16M',folder_heatMap, fixed = T)
    
    folder_heatMap = gsub("baseline","B",folder_heatMap,fixed = T)
    folder_heatMap = gsub("Good Response","GR",folder_heatMap,fixed = T)
    folder_heatMap = gsub("Poor Response","PR",folder_heatMap,fixed = T)
    
    dir.create(folder_heatMap,recursive  = T)
    
    pathName = as.character(paste0(folder_heatMap,filename) )
    pathName = gsub('CD14+ Mono','CD14M',pathName, fixed = T)
    pathName  = gsub('CD8+ T Cell','CD8T',pathName, fixed = T)
    pathName  = gsub('CD16+ Mono','CD16M',pathName, fixed = T)
    
    pathName = gsub("Baseline","B",pathName,fixed = T)
    pathName = gsub("Good Response","GR",pathName,fixed = T)
    pathName = gsub("Poor Response","PR",pathName,fixed = T)
    
    print(pathName)
    png(file=pathName,width=2000, height=3000,res=100)
    print(plot)
    dev.off()
  }
}