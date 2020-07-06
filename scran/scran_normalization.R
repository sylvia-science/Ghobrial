library(scRNAseq)
library(scran)
library(Seurat)
library(scater)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readxl)
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
library(Matrix)

filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]



sample_name = metaData$Sample[i]

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')

sample_list = c('GL1290BM','GL1497BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM','GL1502BM', 'GL3417BM', 'GL2653BM')

sample_list = c('GL1290BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM')

sample_list = c('GL1502BM')
#folder_input ='Soup_MT_nFeature'
folder_input ='Soup_MT_C100'

sample_list = c('GL1080BM','GL1478BM','GL1160BM','GL1320BM', 'GL2257BM','GL2941BM','GL3563BM','NBM10CD138N','NBM11CD138N')
sample_list = c('GL1797BM')


sample_list = c('GL1080BM', 
'GL1110BM',
'GL1478BM',
'GL1320BM',
'GL2257BM',
'GL2941BM',
'GL1502BM' ,
'GL3417BM' ,
'GL3404BM' ,
'GL3563BM' ,
'NBM10CD138N', 
'NBM11CD138N' ,
'NBM12CD138N')


i = 12
# Soup + MT
run = F
for (i in 1:nrow(sampleParam)){
#for (i in 20:35){
#for (i in 4 ){
  sample_name =   sampleParam$Sample[i]
  #sample_name ='GL1374BM'
  #sample_name =  sample_list[i]#sampleParam$Sample[i]
  print(sample_name)
  resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  RNA_features_min = sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
  RNA_features_max = sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name]

  if ( folder_input == 'Soup_MT_nFeature'){
    cluster_IDs = cluster_id_param$Cluster_Ids_nFeature_Scran[cluster_id_param$Sample == sample_name]
  }else{
    cluster_IDs = cluster_id_param$Cluster_Ids_Scran[cluster_id_param$Sample == sample_name]
    
  }
  
  if (run){
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/',folder_input,'/',sample_name,'/')
    path = paste0(folder,'/',sample_name,'.Robj')
    data_i_run = loadRData(path)
   
    data_i_run = data_i_run[, data_i_run$percent.mt < percent_mt]
    #data_i_run = data_i_run[, data_i_run$nFeature_RNA > RNA_features_min]
    
    if ( folder_input == 'Soup_MT_nFeature'){
      data_i_run = data_i_run[, data_i_run$nFeature_RNA > RNA_features_min]
      data_i_run = data_i_run[, data_i_run$nFeature_RNA < RNA_features_max]
    }
    
    
    
    
    #all.equal( as.matrix(data_i@assays[["RNA"]]@data),as.matrix(data_i_sc@assays@data))
    
    #data_i@assays[["RNA"]]@data[c('CD3D','HBB','FCGR3A','CD14','IGKC','NKG7'),1:10]
    #data_i_sc@assays@data$logcounts[c('CD3D','HBB','FCGR3A','CD14','IGKC','NKG7'),1:10]
    
  
    # lot = F
    # xlim of 2000 to 5000
    data_i_filtered_run = ScranNorm(data_i_run)
    data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
    data_i_filtered_run = ScaleData(data_i_filtered_run)
    data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
    data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
    data_i_filtered_run = FindClusters(data_i_filtered_run,resolution = resolution_val)
    data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
    
    
    
    data_i_filtered_run$oldLabel = Idents(data_i_filtered_run)
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/',folder_input,'/Scran/',sample_name,'/')
    dir.create(folder,recursive = T)
    path = paste0(folder,'/',sample_name,'.Robj')
    
    data_i_filtered_run$sample = sample_name
    if (!file.exists(path)){
      save(data_i_filtered_run,file= path)
    }
    
  }else{
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/',folder_input,'/Scran/',sample_name,'/')
    path = paste0(folder,'/',sample_name,'.Robj')
    data_i_filtered_run = loadRData(path)
  }
  data_i_filtered_run$sample = sample_name
  data_i_filtered_run = addMetaData(data_i_filtered_run, metaData)
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/',folder_input,'/Scran/',sample_name,'/')
  
  data_i_filtered_run$oldLabel = Idents(data_i_filtered_run)
  data_i_filtered_run = label_cells(data_i_filtered_run, cluster_IDs)
  CellIdents = Idents(data_i_filtered_run)
  print(CellIdents[1:10])
  #browser()
  
  write.csv(CellIdents, file = paste0(folder,'cellIdents.csv'))
  next
  data_matrix = data_i_filtered_run@assays[["RNA"]]@data
  var_feature = data_i_filtered_run@assays[["RNA"]]@var.features
  write.csv(var_feature, file = paste0(folder,sample_name,'_var_feature.csv'))
  
  write.csv(colnames(data_matrix), file = paste0(folder,sample_name,'_colnames.csv'))
  write.csv(rownames(data_matrix), file = paste0(folder,sample_name,'_rownames.csv'))
  write.csv(data_matrix, file = paste0(folder,sample_name,'_NormData.csv'))
  write.csv(data_i_filtered_run@assays[["RNA"]]@var.features, file = paste0(folder,sample_name,'_varFeatures.csv'))
  
  next
  if(ncol(data_i_filtered_run)< 200){
    pt.size = 3
  }else{
    pt.size = 1
  }
  print('Plot')
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )

    
  data_i_filtered_run$split_var = ''
  resolution_val = 1.4
  cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
  #groupBy_list = c('sample','Diagnosis','10X kit','Treatment','Batch')
  groupBy_list = c('oldLabel')
  splitBy_list = NA
  plotAll(data_i_filtered_run, folder = folder,
          sample_name,sampleParam = NA,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = F, keepOldLabels = T, 
          groupBy = groupBy_list, splitBy = splitBy_list,
          PCA_dim = 30,resolution_val = resolution_val)
  
  next
  PlotKnownMarkers(data_i_filtered_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                   plotType ='FeaturePlotFix' , str = '', plotTogether = F)
  
  
  next
  ###################################
  rowMeans_list = rowMeans(data_i_filtered_run@assays[["RNA"]]@data)
  bin_list = as.numeric(cut_number(rowMeans_list,6))
  
  
  colSums_list = colSums(data_i_filtered_run@assays[["RNA"]]@counts)
  
  plot_list = vector('list', 6)
  plot_list_flat = vector('list', 6)
  plot_list_sizeFactor = vector('list', 6)
  plot_list_sizeFactor_log = vector('list', 6)
  plot_list_librarySize = vector('list', 6)
  plot_list_librarySize_zscore = vector('list', 6)
  plot_list_sizeFactor_zscore = vector('list', 6)
  
  #ylim_list = c(0.02,0.1,0.25,0.5,0.5,5)
  for (bin in 1:6){
    #bin = bin + 1
    print(bin)
    
    y_list_count = data_i_filtered_run@assays[["RNA"]]@counts[bin_list == bin,]
    y_list_count_stdev = colSds(as.matrix(y_list_count))
    y_list_count_sf =as.matrix( t(t(y_list_count)/data_i_filtered_run$sizeFactors))
    y_list_count_sf_zscore = t((t(y_list_count_sf) - colMedians(y_list_count_sf))/colSds(y_list_count_sf))
    
    y_list_count_sf_flat = as.matrix(y_list_count_sf)
    dim(y_list_count_sf_flat) = NULL
    
    y_list_count_librarySizeNorm = as.matrix(t(t(y_list_count)/(colSums_list)))
    y_list_count_librarySizeNorm_zscore = t((t(y_list_count_librarySizeNorm) - colMedians(y_list_count_librarySizeNorm))/colSds(y_list_count_librarySizeNorm))
    y_list = data_i_filtered_run@assays[["RNA"]]@data[bin_list == bin,]
    
    df = as.data.frame(matrix(ncol = 6, nrow = length(colSums_list)))
    colnames(df) = c('librarySize','sizeFactorNorm','sizeFactorNorm_zscore','librarySizeNorm','librarySizeNorm_zscore','ID')
    df$librarySize = colSums_list
    df$exp = colMeans(y_list) # Norm Mean per cell
    df$SD = y_list_count_stdev
    df$sizeFactorNorm = colMeans(as.matrix(y_list_count_sf))
    df$sizeFactorNorm_zscore = colMeans((y_list_count_sf_zscore))
    #df$log1p = log(df$exp_sizeFactors + 1)
    df$librarySizeNorm = colMeans(y_list_count_librarySizeNorm)
    df$librarySizeNorm_zscore = colMeans(y_list_count_librarySizeNorm_zscore)
    
    df$ID = Idents(data_i_filtered_run[bin_list == bin,])
    
    plot = ggplot(df, aes(librarySize, exp)) +
      geom_point(aes(colour = df$ID), size = 2)+
      stat_smooth(aes(x=librarySize,y=exp), method="loess", se=F, color="tomato2") +
      theme(text = element_text(size=20))
    plot_list[[bin]]= plot
    
    
    plot = ggplot(df, aes(librarySize, sizeFactorNorm)) +
      geom_point(aes(colour = df$ID), size = 2)+
      theme(text = element_text(size=20))
    plot_list_sizeFactor[[bin]]= plot
    
    plot = ggplot(df, aes(librarySize, sizeFactorNorm_zscore)) +
      geom_point(aes(colour = df$ID), size = 2)+
      theme(text = element_text(size=20))
    plot_list_sizeFactor_zscore[[bin]]= plot
    
    plot = ggplot(df, aes(librarySize, librarySizeNorm)) +
      geom_point(aes(colour = df$ID), size = 2)+
      theme(text = element_text(size=20))
    plot_list_librarySize[[bin]]= plot
    
    plot = ggplot(df, aes(librarySize, librarySizeNorm_zscore)) +
      geom_point(aes(colour = df$ID), size = 2)+
      theme(text = element_text(size=20))
    plot_list_librarySize_zscore[[bin]]= plot

    
    
    
  }
  #librarySize_flat = rep(colSums_list, times = nrow(y_list_count_sf)) 
  #ID_flat =  Idents(data_i_filtered_run[bin_list == bin,])
  #ID_flat = rep(ID_flat, times = nrow(y_list_count_sf) )
  #df_flat <- data.frame(cbind(y_list_count_sf_flat, librarySize_flat,ID_flat), f)
  #colnames(df_flat) <- c("sizeFactorNorm_flat","librarySize_flat",'ID')
  
  
  #y_list_count_t = log(y_list_count_t + 1)
  #y_list_count_librarySizeNorm = (prop.table(as.matrix(y_list_count),2))
  #y_list_count_librarySizeNorm = colMeans(t(t(y_list_count)/colSums(y_list_count)))

  pathName = paste0(folder,sample_name,'_Scran_NormVsGeneSum','.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list[[1]], plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  
  
  
  pathName = paste0(folder,sample_name,'_Scran_NormVsGeneSum_SizeFactor','.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list_sizeFactor[[1]], plot_list_sizeFactor[[2]],plot_list_sizeFactor[[3]],
    plot_list_sizeFactor[[4]],plot_list_sizeFactor[[5]],plot_list_sizeFactor[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Scran_NormVsGeneSum_SizeFactor_zscore','.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list_sizeFactor_zscore[[1]], plot_list_sizeFactor_zscore[[2]],plot_list_sizeFactor_zscore[[3]],
    plot_list_sizeFactor_zscore[[4]],plot_list_sizeFactor_zscore[[5]],plot_list_sizeFactor_zscore[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Scran_NormVsGeneSum_librarySizeNorm','.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list_librarySize[[1]], plot_list_librarySize[[2]],plot_list_librarySize[[3]],
    plot_list_librarySize[[4]],plot_list_librarySize[[5]],plot_list_librarySize[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Scran_NormVsGeneSum_librarySizeNorm_zscore','.png')
  png(file=pathName,width=1500, height=1000,res = 100)
  plot = plot_grid(
    plot_list_librarySize_zscore[[1]], plot_list_librarySize_zscore[[2]],plot_list_librarySize_zscore[[3]],
    plot_list_librarySize_zscore[[4]],plot_list_librarySize_zscore[[5]],plot_list_librarySize_zscore[[6]],
    labels = "AUTO", ncol = 2)
  print(plot)
  dev.off()
  

  
  ######
  
  gene_list = read.csv('/home/sujwary/Desktop/scRNA/Data/HousekeepingGenes.csv')
  gene_list = as.character(gene_list[gene_list$Gene %in% rownames(data_i_filtered_run),1])
  
  data_subset = (data_i_filtered_run@assays[["RNA"]]@counts[gene_list,])
  data_subset_sf =as.matrix( t(t(as.matrix(data_subset))/data_i_filtered_run$sizeFactors))
  
  data_subset_library = as.matrix(t(t(as.matrix(data_subset))/(colSums_list)))
  
  var_list_sf = colVars((data_subset_sf))
  mean_list_sf = colMeans(data_subset_sf)
  var_list_library = colVars((data_subset_library))
  mean_list_library = colMeans(data_subset_library)
  UMI_per_cell = colSums(as.matrix(data_i_filtered_run@assays[["RNA"]]@counts))
  
  df = data.frame(matrix(ncol = 3, nrow = length(var_list_library)))
  colnames(df) = c("library Size",'var_list_sf','Variance of house-keeping genes')
  df$'LibrarySize' = UMI_per_cell
  df$var_list_sf = var_list_sf
  df$var_list_library = var_list_library

  
  pathName = paste0(folder,sample_name,'_scran_varVsCount_sf','.png')
  png(file=pathName,width=1000, height=1000)
  plot = ggplot(df, aes(LibrarySize, var_list_sf)) +
    geom_point(size = 4)+
    theme(text = element_text(size=20))
  plot = plot +xlab("Library Size") + ylab("Variance")
  print(plot)
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_varVsCount_librarysize','.png')
  png(file=pathName,width=1000, height=1000)
  plot = ggplot(df, aes(LibrarySize, var_list_library)) +
    geom_point(size = 4)+
    theme(text = element_text(size=20))
  plot = plot +xlab("Library Size") + ylab("Variance")
  print(plot)
  dev.off()

  
  ####################
  #next
  
  
}


