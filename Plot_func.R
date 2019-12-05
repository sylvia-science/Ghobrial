
###################
# Plot Functions
###################

################################
## Quality Control
################################

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

#########################
## Make FeaturePlot
#########################
MakeFeaturePlot = function(data,data_orig, folder, cell_features,split){
  
  print(folder)
  
  sample <- data
  
  if (is.na(cell_features)){
    
    cell_features = getCellMarkers()
    cell_features = cell_features[cell_features$Plot_marker == TRUE,]
    
  } 
  
  cell_list = cell_features$Cell
  feature_list = cell_features$Markers
  
  all_markers  = data_orig@assays[['RNA']]
  all_markers = all_markers@data@Dimnames[[1]]
  for(i in seq_len(length(feature_list))){
    cell_type = cell_list[i]
    
    
    x = unlist(strsplit(feature_list[i], ",")) 
    x = gsub("\\s", "", x)  
    
    # Take first 9 markers if longer than 9 genes
    if (length(x) > 12){
      x = x[1:12]
    }
    gene_list = (x[x %in% all_markers])  
    # FeaturePlot gives error if no genes are found, so we have to check if they are included in the highly variable gene list
    if (length(gene_list > 0)){
      #print(paste0( cell_type,': Found'))
      #print(x)
      
      # TO DO Add dotplot
      subtitle_str = paste0(cell_type,': ', toString(x))
      if (split){
        plot = FeaturePlot(sample, features = (x),split.by = "orig.ident")
      }else{
        plot = FeaturePlot(sample, features = (x))
      }
      plot = plot + labs(subtitle=subtitle_str) + theme(plot.subtitle = element_text(hjust=0.5, size=16))
      pathName <- paste0(folder,cell_type,'.png')
      png(file=pathName)
      print(plot)
      dev.off()
    }
    else {
      #print(paste0( cell_type,': Not Found'))
      #print(x)
    }
    print('')
    
  }
  
}

FeaturePlotFix = function(data, feature,folder,str, split, gene_TF){
  # Feature: 1 gene name
  if (feature == "HLA-DQA2"){
    #browser()
  }
  data_umap = FetchData(data, vars = c("ident", "orig.ident","UMAP_1", "UMAP_2"))
  data_df = data.frame(data@assays$RNA@data) 
  
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
  
  
    if (split){
      class_umap_data_pre = subset(data_umap, data_umap$orig.ident == "data_pre",)
      class_umap_data_post = subset(data_umap, data_umap$orig.ident == "data_post",)
    
      #data_umap_cell = subset(data_umap, data_umap$ident == str,)
      #max_val = ceiling(max(data_umap_cell$g))
      max_val = 4
      min_val = 0
      # Adding cluster label to center of cluster on UMAP
    
      umap_label <- FetchData(data, vars = c("ident", "orig.ident","UMAP_1", "UMAP_2"))  %>%
      group_by(orig.ident, ident) %>% 
      summarise(x=mean(UMAP_1), y = mean(UMAP_2))
      
      umap_label_pre <- subset(umap_label, umap_label$orig.ident == "data_pre",)
    
      fontSize = 24
      markerSize = 1
      g_pre = ggplot(class_umap_data_pre, aes(UMAP_1, UMAP_2)) +  
        geom_point(aes(colour = class_umap_data_pre$g), size = markerSize) + 
        scale_color_gradient(low = "grey", high = "red", name = "", limits = c(min_val,max_val) ) +
        ggtitle("Pre",  subtitle = feature) +
        theme(plot.title = element_text(hjust = 0.5))+
        theme(plot.subtitle = element_text(hjust = 0.5))+
        theme(plot.title = element_text(size=fontSize)) + 
        theme(plot.subtitle = element_text(size=fontSize)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      
      # Adding cluster label to center of cluster on UMAP
      
      umap_label <- FetchData(data, vars = c("ident", "orig.ident","UMAP_1", "UMAP_2"))  %>%
        group_by(orig.ident, ident) %>% 
        summarise(x=mean(UMAP_1), y = mean(UMAP_2))
      
      umap_label_post <- subset(umap_label, umap_label$orig.ident == "data_post",)
      
      
      g_post = ggplot(class_umap_data_post, aes(UMAP_1, UMAP_2)) +  
        geom_point(aes(colour = class_umap_data_post$g), size=markerSize) + 
        scale_color_gradient(low = "grey", high = "red", name = "" , limits = c(min_val,max_val) ) +
        ggtitle("Post",  subtitle = feature) +
        theme(plot.title = element_text(hjust = 0.5))+
        theme(plot.subtitle = element_text(hjust = 0.5))+
        theme(plot.title = element_text(size=fontSize)) + 
        theme(plot.subtitle = element_text(size=fontSize)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      g_output = CombinePlots(plots = list(g_pre, g_post))
      
      pathName <- paste0(folder,str,'_',feature,'.png')
      png(file=pathName, width=2000, height=1000, res=100)
      print(g_output)
      dev.off()
      
    }
  #browser()
  
}

