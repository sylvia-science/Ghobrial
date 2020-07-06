library(scRNAseq)
library(Seurat)
library(dplyr)
library(ggplot2)
library(OneR )
library(cowplot)
library(Hmisc)
library(dplyr)
library(readxl)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')

filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam = read_excel(filename_sampleParam)

filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluseter_id_param = read_excel(filename)


sample_name = metaData$Sample[i]

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')
#sample_list = c('GL1497BM', 'GL1160BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')
sample_list = c('GL1290BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM')

#folder_input ='Soup_MT_nFeature'
folder_input ='Soup_MT_C100'
run = F
i = 6
# Soup + MT
for (i in 1:length(sample_list) ){
#for (i in 4 ){
  sample_name = sample_list[i]
  print(sample_name)

  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  RNA_features_min = sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
  RNA_features_max = sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name]
  
  if ( folder_input == 'Soup_MT_nFeature'){
    cluster_IDs = cluseter_id_param$Cluster_Ids_nFeature_Sctransform[cluseter_id_param$Sample == sample_name]
  }else{
    cluster_IDs = cluseter_id_param$Cluster_Ids_Sctransform[cluseter_id_param$Sample == sample_name]
    
  }
  
  
  if (run){
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/',folder_input,'/',sample_name,'/')
    path = paste0(folder,'/',sample_name,'_data.Robj')
    data_i_run = loadRData(path)
    
    # filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
    # data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
    # data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
    # 
    # colSum_list = colSums(data_i_raw ) # Needs to be from Matrix library
    # keep = colSum_list >= 100
    # data_i_run = data_i_raw[,keep]
    # 
    # data_i_run[["percent.mt"]] = PercentageFeatureSet(data_i_run, pattern = "^MT-")
    data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
    #data_i_run = data_i_run[, data_i_run$nFeature_RNA > RNA_features_min]
    #data_i_run = data_i_run[, data_i_run$nFeature_RNA < RNA_features_max]
    # 
    
    if ( folder_input == 'Soup_MT_nFeature'){
      data_i_run = data_i_run[, data_i_run$nFeature_RNA > RNA_features_min]
      data_i_run = data_i_run[, data_i_run$nFeature_RNA < RNA_features_max]
    }
    
    
    data_i_filtered_run = SCTransform(data_i_run, verbose = T, do.center = F)
     
    data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
    data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
    data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
    data_i_filtered_run = FindClusters(data_i_filtered_run,resolution = 1.2)
  
    
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/',folder_input,'/Sctransform/',sample_name,'/')
    dir.create(folder,recursive = T)
    path = paste0(folder,'/',sample_name,'_centerF.Robj')
    save(data_i_filtered_run,file= path)
  }else{
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/',folder_input,'/Sctransform/',sample_name,'/')
    path = paste0(folder,'/',sample_name,'_centerF.Robj')
    data_i_filtered_run = loadRData(path)
  }
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/',folder_input,'/Sctransform/',sample_name,'/')
  
  
  data_i_filtered_run = label_cells(data_i_filtered_run, cluster_IDs)
  
  if(ncol(data_i_filtered_run)< 200){
    pt.size = 3
  }else{
    pt.size = 1
  }
  print('Plot')
  
  
  
  #plot the total UMI count per cell on the x axis (aka colSums) and normalized expression on the y axis. 
  # This will result in 6 plots, so use cowplot to put all six of these in three rows and two columns (2 plots per row).

  rowMeans_list = rowMeans(data_i_filtered_run@assays[["SCT"]]@data)
  bin_list = as.numeric(cut_number(rowMeans_list,6))
  
  
  colSums_list = colSums(data_i_filtered_run@assays[["RNA"]]@counts)
  
  plot_list = vector('list', 6)
  plot_list_10000 = vector('list', 6)
  plot_list_noLog =  vector('list', 6)
  for (bin in 1:6){
    #bin = bin + 1
    print(bin)
    y_list = as.matrix(data_i_filtered_run@assays[["SCT"]]@data[bin_list == bin,])
    y_list_nolog = exp((y_list - 1))
    y_list_sd = colSds(y_list)
    df = as.data.frame(matrix(ncol = 2, nrow = length(colSums_list)))
    colnames(df) = c('colsum','exp')
    df$colsum = (colSums_list) # Unormalized Mean per cell
    df$exp = colMeans(y_list) # Norm Mean per cell
    df$exp_nolog = colMeans(y_list_nolog)
    df$ID = Idents(data_i_filtered_run[bin_list == bin,])
    df$SD = y_list_sd
    
    plot = ggplot(df, aes(colsum, exp)) +
      geom_point(aes(colour = df$ID), size = 2)+
      stat_smooth(aes(x=colsum,y=exp), method="loess", se=F, color="tomato2") +
      theme(text = element_text(size=20))#+ xlim(c(0,10000))
    #print(plot)
    plot_list[[bin]]= plot
    
    plot = ggplot(df, aes(colsum, exp_nolog)) +
      geom_point(aes(colour = df$ID), size = 2)+
      #stat_smooth(aes(x=colsum,y=exp), method="loess", se=F, color="tomato2") +
      theme(text = element_text(size=20))#+ xlim(c(0,10000))
    #print(plot)
    plot_list_noLog[[bin]]= plot

    
  }
  
  str = ''
  pathName = paste0(folder,sample_name,'_Sctransform_NormVsGeneSum',str,'.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list[[1]], plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Sctransform_NormVsGeneSum_nolog',str,'.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list_noLog[[1]], plot_list_noLog[[2]],plot_list_noLog[[3]],
    plot_list_noLog[[4]],plot_list_noLog[[5]],plot_list_noLog[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  
  #next
  
  # plot housekeeping gene variance across all cells on the y axis and total UMI count per cell on the x axis
  
  gene_list = read.csv('/home/sujwary/Desktop/scRNA/Data/HousekeepingGenes.csv')
  gene_list = as.character(gene_list[gene_list$Gene %in% rownames(data_i_filtered_run),1])

  data_subset = as.matrix(data_i_filtered_run@assays[["SCT"]]@data[gene_list,])
  data_subset = exp((data_subset - 1))
  var_list = colVars((data_subset))
  mean_list = colMeans(data_subset)
  UMI_per_cell = colSums(as.matrix(data_i_filtered_run@assays[["RNA"]]@counts))
  
  df = data.frame(matrix(ncol = 3, nrow = length(var_list)))
  colnames(df) = c("LibrarySize",'variance','mean')
  df$'LibrarySize' = UMI_per_cell
  df$var_list = var_list
  df$mean = mean_list
  
  pathName = paste0(folder,sample_name,'_Sctransform_varVsCount_nolog','.png')
  png(file=pathName,width=1000, height=1000)
  plot = ggplot(df, aes(LibrarySize, var_list)) +
    geom_point(size = 4)+
    theme(text = element_text(size=20))
  plot = plot +xlab("Library Size") + ylab("Variance")
  print(plot)
  dev.off()
  
  
  
  #next
  
  data_i_filtered_noLabel = FindClusters(data_i_filtered_run,resolution = 1.4)
  
  pathName = paste0(folder,sample_name,'_Sctransform_Umap','_noLabel','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_noLabel,pt.size = pt.size, reduction = "umap",label = T))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Sctransform_Umap','_Label','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap",label = T))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Sctransform_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Sctransform_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  next
  
  PCA_dim = 30
  resolution_val = 0.8
  num_markers = 10
  markers = FindAllMarkers(data_i_filtered_run, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  top10 = markers %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
  all_markers =  markers %>% group_by(cluster)
  
  pathName = paste0(folder,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,'.png'))
  png(file=pathName,width=500, height=500)
  print(DoHeatmap(data_i_filtered_run, features = top10$gene))
  dev.off()
  
  
  write.csv(all_markers, file = paste0(folder,'Features','.csv'),row.names=FALSE)
  
  #next
  
  if(ncol(data_i_filtered_run)< 200){
    pt.size = 6
  }else{
    pt.size = 3
  }
  
  folder_feature = paste0(folder,'Featureplots/' )
  dir.create(folder_feature,recursive = T)
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17',' PRTN3', 'NAMPT')
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    plot = FeaturePlotFix(data_i_filtered_run, feature = gene, folder =folder,
                          str = '',split = F, markerSize = pt.size,gene_TF = TRUE,title = '',saveTF = FALSE)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=24),
      legend.title=element_text(size=24),
      text = element_text(size = 20)
    )
    
    file_str = ''
    pathName = paste0(folder_feature,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
}


