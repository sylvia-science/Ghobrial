# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/AnalyzeP10P12.R')
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
library(matrixStats)
require(gridExtra)

#library(biomaRt)
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

folder_base_input <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'
folder_base_output <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_sampleParam_integrate <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters.xlsx'
sampleParam_integrate <- read_excel(filename_sampleParam_integrate)

filename_sample_Integrate_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Integrate_pairs.xlsx'
sample_Integrate_pairs <- read_excel(filename_sample_Integrate_pairs)

#filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
#metaData <- read_excel(filename_metaData)
sample_type = 'BM'

cell_type_list = c('T Cell','T Cell Cytotoxic', 'Monocyte CD14','Monocyte FCGR3A','NK','B Cell', 'DC')

#######################################################

## get DE HLA Genes from Patient 10

patient = 10

pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$`Patient Number` == patient, ]
sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
data = loadRData(paste0(folder_pre,'data.Robj'))


Features_CD14_MonoVsn_label_Patient10 = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/GL1080BM/Filtered/Regress/DE/nVsm/Features_CD14+ MonoVsn_label_Patient10.csv'    
Features_CD14_MonoVsn_label_Patient10  = read.csv(Features_CD14_MonoVsn_label_Patient10)
Features_CD14_MonoVsn_label_Patient10 = Features_CD14_MonoVsn_label_Patient10[Features_CD14_MonoVsn_label_Patient10$ident_2 == 'dCD14+ Mono', ]
Features_CD14_MonoVsn_label_Patient10 = Features_CD14_MonoVsn_label_Patient10[Features_CD14_MonoVsn_label_Patient10$p_val_adj < 0.05, ]
Features_CD14_MonoVsn_label_Patient10 = levels(Features_CD14_MonoVsn_label_Patient10$gene)

Features_dCD14_MonoVsn_label_Patient10 = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/GL1080BM/Filtered/Regress/DE/nVsm/Features_dCD14+ MonoVsn_label_Patient10.csv'    
Features_dCD14_MonoVsn_label_Patient10  = read.csv(Features_dCD14_MonoVsn_label_Patient10)
Features_dCD14_MonoVsn_label_Patient10 = Features_dCD14_MonoVsn_label_Patient10[Features_dCD14_MonoVsn_label_Patient10$ident_2 == 'CD14+ Mono', ]
Features_dCD14_MonoVsn_label_Patient10 = Features_dCD14_MonoVsn_label_Patient10[Features_dCD14_MonoVsn_label_Patient10$p_val_adj < 0.05, ]
Features_dCD14_MonoVsn_label_Patient10 = levels(Features_dCD14_MonoVsn_label_Patient10$gene)

DEgeneHLA = c(Features_dCD14_MonoVsn_label_Patient10)
DEgeneHLA = DEgeneHLA[grep('HLA', DEgeneHLA)]
#DEgeneHLA = "HLA-DMB"

####################
## Plot
####################

patient = 10

pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$`Patient Number` == patient, ]
sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
data = loadRData(paste0(folder_pre,'data.Robj'))

PCA_dim = sampleParam_integrate$PCA_dim[sampleParam_integrate['Sample'] == sample_name_pre]
resolution_val<- sampleParam_integrate$resolution_val[sampleParam_integrate['Sample'] == sample_name_pre]
data = getCluster(data,resolution_val, PCA_dim)

cluster_IDs = sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Sample'] == sample_name_pre]
data = label_cells(data,cluster_IDs)
cluster_IDs = unlist(strsplit(cluster_IDs, ",")) 
cluster_IDs = trimws(cluster_IDs, which = c("both"), whitespace = "[ \t\r\n]")

plot = FALSE
if (plot){
  
  
  
  
  filepath_cluster = paste0( folder_pre, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_NoLabel','.png'))
  png(file=pathName,width=600, height=350, res = 100)
  print(DimPlot(data, reduction = "umap",label = FALSE, pt.size = 1)+ NoLegend())
  dev.off()
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_split','_NoLabel','.png'))  
  png(file=pathName,width=600, height=350)
  print(DimPlot(data, label=FALSE, repel=F, reduction = "umap", split.by = "orig.ident")+ NoLegend())
  dev.off()
  
  
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers  %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
  top = markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
  
  if ('dCD14+ Mono' %in% cluster_IDs){
    gene_monocyte = top[top$cluster == 'CD14+ Mono' | top$cluster ==	'dCD14+ Mono' | top$cluster == 'DC',]
    data_monocyte = subset(data, idents = c("CD14+ Mono", "dCD14+ Mono", "DC"))
  }else{
    gene_monocyte = top[top$cluster == 'CD14+ Mono'| top$cluster == 'DC',]
    data_monocyte = subset(data, idents = c("CD14+ Mono", "DC"))
  }
  
  
  pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,'_Monocyte','.png'))
  png(file=pathName,width=1000, height=1200)
  print(DoHeatmap(data_monocyte, features = gene_monocyte$gene) )
dev.off()
}
#######################################
## Score
#######################################
if ('dCD14+ Mono' %in% cluster_IDs){
  data_monocyte = subset(data, idents = c("CD14+ Mono", "dCD14+ Mono"))
}else{
  data_monocyte = subset(data, idents = c("CD14+ Mono"))
}
DEgeneHLA_list = list(DEgeneHLA)
score_HLA_orig = AddModuleScore(object = data_monocyte, features = DEgeneHLA_list, ctrl = 8, name = 'HLA_Feature')
score_HLA_orig = score_HLA_orig[[]]

score_HLA = subset(score_HLA_orig, select = c(orig.ident,HLA_Feature1))
#score_HLA = score_HLA_orig
#score_HLA = score_HLA[ , grepl( "HLA_Feature" , names( score_HLA ) ) ]
#colnames(score_HLA)=DEgeneHLA

score_HLA_output = score_HLA
#score_HLA_output$orig.ident = score_HLA_orig$orig.ident
#score_HLA_output$old.ident = score_HLA_orig$old.ident
#score_HLA_output$seurat_clusters = score_HLA_orig$seurat_clusters


#score_HLA_summary = colMeans(score_HLA)
#score_HLA_summary = data.frame(as.list(score_HLA_summary))
#score_HLA_SD = colSds(data.matrix(score_HLA))
#score_HLA_summary[2,] = score_HLA_SD
#row.names(score_HLA_summary) = c('Mean','SD')


score_HLA_summary = data.frame(matrix(ncol = 2, nrow = 2))
colnames(score_HLA_summary) = c('Score Pre', 'Score Post')
row.names(score_HLA_summary) = c('Mean','SD')
score_HLA_summary$'Score Pre'[1] = mean(score_HLA$HLA_Feature1[score_HLA$orig.ident == 'data_pre' ])
score_HLA_summary$'Score Post'[1] = mean(score_HLA$HLA_Feature1[score_HLA$orig.ident == 'data_post' ])
score_HLA_summary$'Score Pre'[2] = sd(score_HLA$HLA_Feature1[score_HLA$orig.ident == 'data_pre' ])
score_HLA_summary$'Score Post'[2] = sd(score_HLA$HLA_Feature1[score_HLA$orig.ident == 'data_post' ])

filepath_analysis = paste0( folder_pre, 'Analysis/' )
#dir.create( filepath_analysis, recursive = TRUE)


#write.csv(score_HLA_orig, file = paste0(filepath_analysis,'score_HLA_patient',patient,'.csv'),row.names=TRUE)
#write.csv(score_HLA_summary, file = paste0(filepath_analysis,'score_HLA_summary_patient',patient,'.csv'),row.names=TRUE)

# view cell cycle scores and phase assignments
#head(object)

#browser()

patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all

for(patient in patient_list){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  folder_output =  paste0( folder_pre, 'Analysis/FeaturePlot_HLA/' )
  dir.create( folder_output, recursive = TRUE)
  
  PCA_dim = sampleParam_integrate$PCA_dim[sampleParam_integrate['Sample'] == sample_name_pre]
  resolution_val<- sampleParam_integrate$resolution_val[sampleParam_integrate['Sample'] == sample_name_pre]
  cluster_IDs = sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Patient Number'] == patient]
  
  data = loadRData(paste0(folder_pre,'data.Robj'))
  
  data = getCluster(data,resolution_val, PCA_dim)
  
  data = label_cells(data,cluster_IDs)
  
  
  #browser()
  for (i in 1:(length(DEgeneHLA))){
    feature = as.character(DEgeneHLA[i])
    cell_str = paste0('patient_',patient)
    print(feature)
    #FeaturePlotFix(data, feature,folder_output,cell_str, split = TRUE, gene_TF = TRUE)
  }
  
  data = AddModuleScore(object = data, features = list(DEgeneHLA),nbin = 25, ctrl = 7, name = 'HLA_Feature')
  score_HLA_orig = data[[]]
  
  score_HLA = subset(score_HLA_orig, select = c(orig.ident,HLA_Feature1))
  cell_str = paste0('patient_',patient)
  FeaturePlotFix(score_HLA, 'HLA_Feature1',folder_output,cell_str, split = TRUE, gene_TF = FALSE)
  
  

}


##############
## Individual
##############


# folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'
# 
# filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
# sampleParam <- read_excel(filename_sampleParam)
# 
# filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
# metaData <- read_excel(filename_metaData)
# 
# 
# 
# for(i in 1:nrow(metaData)){
#   sample_name <- metaData$Sample[i]
#   filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
# 
#   print('Starting Plotting')
#   print(paste0('Sample: ', sample_name))
# 
#       filter = sampleParam$filter_TF[sampleParam['Sample'] == sample_name]
#       regress_TF = sampleParam$regress_TF[sampleParam['Sample'] == sample_name]
#       
#       folder = makeFolders(folder_base,sample_name,filter,regress_TF,makeFolder_TF = TRUE)
#       folder_output =  paste0( folder, 'Analysis/FeaturePlot_HLA/' )
#       dir.create( folder_output, recursive = TRUE)
#       
#       print(folder_output)
#       
#       data = loadRData(paste0(folder,'data.Robj'))
# 
#       patient = NA
#       
#       
#       data = AddModuleScore(object = data, features = list(DEgeneHLA),nbin = 25, ctrl = 7, name = 'HLA_Feature')
#       score_HLA_orig = data[[]]
#       
#       score_HLA = subset(score_HLA_orig, select = c(orig.ident,HLA_Feature1))
#       
#       cell_str = paste0('patient_',patient)
#       #FeaturePlotFix(data, 'HLA_Feature1',folder_output,cell_str, split = TRUE, gene_TF = FALSE)
# 
# }
# 


