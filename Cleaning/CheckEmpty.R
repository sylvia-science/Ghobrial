

library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
library(DoubletFinder )
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)
i = 4

# Doublet
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/C100_Soup_MT/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_data.Robj')
  data_i_run = loadRData(path)
  
  data_i_run$nCount_RNA_log = log(data_i_run$nCount_RNA)
  data_i_run$nFeature_RNA_log = log(data_i_run$nFeature_RNA)
  data_i_run$SmallFeature = data_i_run$nFeature_RNA <= 200
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/','C100_Soup_MT/',sample_name,'/')
  dir.create(folder, recursive = T)
  
  
  
  #num_markers = 10
  #markers = FindAllMarkers(data_i_run, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #markers  %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
  # Plotting the top 10 markers for each cluster.
  #top10 = markers %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
  #all_markers =  markers %>% group_by(cluster)
  
  #all_markers = all_markers[order(-all_markers$avg_logFC),]
  
  
  # Add known markers
  #all_markers = cellMarkers(all_markers,cell_features)
  #write.csv(all_markers, file = paste0(filepath_cluster,'AllFeatures',file_str,'_Patient',Patient_num,'.csv'),row.names=FALSE)
  #write.csv(all_markers, file = paste(folder, 'markers','.csv'),row.names=FALSE)
  
  #cluster = '10'
  #markers = FindMarkers(data_i_run, ident.1 = cluster)
  #markers = markers[markers$p_val_adj < 0.05,]
  #markers = markers[order(-markers$avg_logFC),]
  
  #path = 
  #write.csv(markers_all, file = path,row.names = FALSE)
  
  print('Plot')
  pathName = paste0(folder,sample_name,'_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()

  pathName = paste0(folder,sample_name,'_Umap_SmallFeature','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE, group.by = 'SmallFeature'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap_is_cell','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE, group.by = 'is_cell'))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_nCount_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nCount_RNA'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_nCount_RNA_log','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nCount_RNA_log'))
  dev.off()
  
  
  pathName = paste0(folder,sample_name,'_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_nFeature_RNA_log','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA_log'))
  dev.off()
  
  gene_list = c( 'PTPRC','CD53', 'MME', 'VNN2','VNN3')
  gene_list = c('CD164', 'CXCR4')
  #next
  gene_list = 'NAMPT'
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'PTPRC', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO','PRTN3', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17', 'NAMPT')
  
  folder_feature = paste0(folder,'Featureplots/Umap/')
  
  #gene_list = rownames(markers)[1:50]
  #folder_feature = paste0(folder,'Featureplots/Umap/Cluster ', cluster, 'Top Markers/')
  dir.create(folder_feature,recursive = T)
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    
    plot = FeaturePlotFix(data_i_run, feature = gene, folder =folder,
                          str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE)
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
