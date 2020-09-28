library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)

library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
#metaData = metaData[metaData$Run== 1,]
metaData = metaData[metaData$`Sample Type` == 'PBMC',]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

# Soup + MT
for (i in 1:nrow(metaData) ){
  sample_name = metaData$Sample[i]
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder_output = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  path = paste0(folder_output,'/',sample_name,'.Robj')
  data_i_run = loadRData(path)
  
  data_i_run_sc =  as.SingleCellExperiment(data_i_run)
  
  #sce = cxds(data_i_run_sc,retRes = TRUE)
  #sce = bcds(sce,retRes = TRUE,verb=TRUE)
  sce = cxds_bcds_hybrid(data_i_run_sc)
  

  
  doublet_log =  log10(sce$hybrid_score+1)
  doublet_quantile_log =  quantile(doublet_log,0.925) 
  hist(doublet_log,breaks=20)
  
  folder_output = paste0('/disk2/Projects/EloRD/Output/','Soup_MT_SCDS/',sample_name,'/')
  dir.create(folder_output, recursive = T)
  scds_doublet = doublet_log > doublet_quantile_log
  write.csv(scds_doublet, file=paste0(folder_output, 'Doublet','.csv'))
  
  #next
  data_i_run$cxds_score = sce$cxds_score
  data_i_run$bcds_score = sce$bcds_score
  data_i_run$hybrid_score = sce$hybrid_score
  data_i_run$scds_doublet = scds_doublet
  #par(mfcol=c(1,3))
  #boxplot(sce$cxds_score   ~ sce$doublet_true_labels, main="cxds")
  #boxplot(sce$bcds_score   ~ sce$doublet_true_labels, main="bcds")
  boxplot(sce$hybrid_score, main="hybrid")
  
  
  pt.size  = 0.8
  print('Plot')
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','_scds_doublet','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_run,pt.size = pt.size, reduction = "umap",label = FALSE, group.by= 'scds_doublet' ))
  dev.off()
  
  
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','hybrid_score','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("hybrid_score")))
  dev.off()
  
  
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, features = c("percent.mt")))
  dev.off()
  
  
  pathName = paste0(folder_output,sample_name,'_Pre_Doublet_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_run,pt.size =pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder_output = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/','Featureplots/Pre_Doublet_Umap/')
    dir.create(folder_output,recursive = T)
    plot = FeaturePlotFix(data_i_run, feature = gene, folder =folder_output,
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
    pathName = paste0(folder_output,gene,file_str,'.png')
    pathName = paste0(folder_output,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
}
