# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/MainIntegrateAll_RunPipline.R')
rm(list = ls())
gc()

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')

source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PipelineIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PlotFunctionIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/IntegrateAll_ClusterUmap.R')


library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)

require(gridExtra)
require(data.table) 
integrate_merge = 'Integrate'
ConvertCategorical = 'ConvCatF'

# Cell type if subsetting
celltype = 'Mono'
#celltype = 'T Cell'
#celltype = 'NK'
celltype = ''
# Variables to regress
regress_var = sort(c("Patient","dexa","kit","Response"))
regress_var = sort(c('dexa','kit',"nCount_RNA", "percent.mt"))
regress_var = sort(c('kit','SampleType'))
#regress_var = sort(c("PatientLast","dexa",'kit'))
#regress_var = sort(c("dexa","kit"))
regress_var = sort(c(""))
regress_var = sort(c("kit"))

#regress_var = sort(c("dexa","kit",'Patient'))
#regress_var = sort(c(""))

#regress_var = c('Patient')
#regress_var = ''
str = paste0('Reg',paste(regress_var, collapse = '_'))

# Set to Clean if cell types have been removed
clean = '/3TimePoints_'
clean = '/'
saveClean = FALSE

# Set to rpca is rpca was used to integrate 
rpca = '_rpca'
rpca = ''
if (integrate_merge == 'Merge'){
  rpca = ''
}


print(paste0('integrate_merge:', integrate_merge))

if (integrate_merge == 'Integrate' || integrate_merge == 'Merge'){
  filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
  filename_metaData <- paste0('/home/sujwary/Desktop/scRNA/Data/EloRD Meta','','.xlsx')
  sampleParam <- read_excel(filename_sampleParam)
  metaData = read_excel(filename_metaData)
}else{
  filename_sampleParam <- '/home/sujwary/Desktop/scRNA/Data/Data/sample_parameters.xlsx'
  sampleParam <- read_excel(filename_sampleParam)
}

print(filename_sampleParam)

sample_type = 'NBMHCL_MT15'
sample_type = 'Tcell_MergeInt_Sub'
sample_type = 'PrePostEOTNBM'
sample_type = 'PrePostEOTNBM_MT15'

if (celltype == ''){
  sample_name = paste0(integrate_merge,'_',sample_type,'_',str)
}else{
  sample_name = paste0(integrate_merge,'_',sample_type,'_',str,'_',celltype)
}
folder_base = '/home/sujwary/Desktop/scRNA/'
folder_base_output = paste0('/home/sujwary/Desktop/scRNA/Output/',
                            integrate_merge ,' All/',sample_type,'',clean,'',str,'/')


PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name] 

cell_features = getCellMarkers(folder_base)

print(PCA_dim)

############################################
RunPipeline = T
if (RunPipeline){
  print('Run')
  #celltype = 'T Cell'
  if (celltype != ''){
    folder_base_output = paste0(folder_base_output, 'SubsetIntegrate/',celltype,'/')
    path = paste0(folder_base_output, 'data_',integrate_merge,'_',sample_type,'.Robj')
    regress_var = c('')
  }else{
    path = paste0(folder_base,'Output/',integrate_merge ,' All/',sample_type,'/data_',integrate_merge,'_',sample_type,'.Robj')
    
  }
  
  str_data_input = '_features2000'
  
  #path = paste0('/home/sujwary/Desktop/scRNA/Output/',integrate_merge ,' All/',sample_type
  #              ,'/data',str_data_input,'.Robj')
  data_integrate = loadRData(path)
  #browser()
  #folder_base_output = paste0(folder_base_output,str,'/')
  #folder_base_output = paste0(folder_base_output, ConvertCategorical,'/')
  pathName = paste0(folder_base_output,'PCA')
  dir.create( pathName, recursive = TRUE)
  
  pathName = paste0(folder_base_output,'Cluster')
  dir.create( pathName, recursive = TRUE)
  
  data_integrate$Patient = data_integrate$'Patient Number'
  data_integrate$kit = data_integrate$'10X kit'
  

  
  #browser()
  ##################################
  #regress_var = sort(c(""))
  cell_features = getCellMarkers(folder_base)
  data_run = PipelineIntegrateAll(data_integrate,
                                  sample_name,
                                  folder_base_output,
                                  sampleParam,
                                  integrate_merge,regress_var =regress_var, 
                                  markersTF = FALSE,
                                  cell_features,
                                  ConvertCategorical = F)

  browser()
  

  # Renaming variables to be simpler
  
  data_run@meta.data$Treatment = gsub("baseline", "Baseline", data_run@meta.data$Treatment)
  data_run$split_var = data_run$Treatment
  
  path = paste0(folder_base_output,
                '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
                '.Robj')
  save(data_run,file= path)

  
  
  
  
  browser()
  groupBy_list = c('dexa','sample_name','Treatment','kit','Patient')
  #groupBy_list = c('dexa','sample_name', 'Response','Treatment','kit','Patient','Cell_Label','SampleType')
  
  plotAll(data_run, folder = folder_base_output,
          sample_name,sampleParam,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF = F, markersTF =T, keepOldLabels = F, groupBy = groupBy_list, 
          PCA_dim = NA,resolution_val = NA)
  
  
  filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  cell_features = getCellMarkers(folder_base)
  PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/HeatMap/'), cell_features = cell_features,
                   plotType ='HeatMap' , str = '')
  PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/Violin/'), cell_features = cell_features,
                   plotType ='Violin' , str = '')
  PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                   plotType ='FeaturePlotFix' , str = '')
  
  #browser()
  ###################################
}else{
  # Load data and plot
  
  if (celltype != ''){
    
    folder_base_output = paste0(folder_base_output, 'SubsetIntegrate/',celltype,'/')
    path = paste0(folder_base_output,
                  '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
                  '.Robj')
  }else{
    path = paste0(folder_base_output,
                  '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
                  '.Robj')

  }
  
  #PCA_dim = 10
  #resolution_val =
  filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )


  data_run = loadRData(path) # The main data we're working with
  data_run = FindClusters(data_run, resolution = resolution_val)
  

  #data_run@meta.data$split_var = paste(data_run@meta.data$orig.ident,data_run@meta.data$dexa)
  #data_run@meta.data$split_var = gsub("data_post Yes", "Post D", data_run@meta.data$split_var)
  #data_run@meta.data$split_var = gsub("data_pre Yes", "Pre D", data_run@meta.data$split_var)
  #data_run@meta.data$split_var = gsub("data_post No", "Post ND", data_run@meta.data$split_var)
  #data_run@meta.data$split_var = gsub("data_pre No", "Pre ND", data_run@meta.data$split_var)
  #data_run@meta.data$split_var = gsub("data_NBM NBM", "NBM", data_run@meta.data$split_var)
  
  # the split_var metadata will contain relevent info for our project
  data_run@meta.data$split_var =  data_run@meta.data$orig.ident
  # Renaming variables to be simpler
  data_run@meta.data$split_var = gsub("data_baseline", "Baseline", data_run@meta.data$split_var)
  data_run@meta.data$split_var = gsub("data_EOT", "EOT", data_run@meta.data$split_var)
  data_run@meta.data$split_var = gsub("data_C9D1", "C9D1", data_run@meta.data$split_var)
  data_run@meta.data$split_var = gsub("data_NBM", "NBM", data_run@meta.data$split_var)
  
  data_run@meta.data$Response = gsub("Minimal Response then Progression", "MRP", data_run@meta.data$Response)
  data_run@meta.data$Response = gsub("Minimal Response", "MR", data_run@meta.data$Response)

  data_run@meta.data$Response[is.na(data_run@meta.data$Response)] = 'NBM'
  data_run@meta.data$kit[(data_run@meta.data$Response) == 'NA'] = 'Microwell-seq'
  data_run$Best_Overall_Response = ''
  data_run$Current_Status = ''
  data_run$Response = ''
  for (sample in unique(data_run$sample_name)){
    data_run$Response[data_run$sample_name == sample] = metaData[metaData$Sample == sample,]$General_Response
    data_run$Best_Overall_Response[data_run$sample_name == sample] = metaData[metaData$Sample == sample,]$Best_Overall_Response
    data_run$Current_Status[data_run$sample_name == sample] = metaData[metaData$Sample == sample,]$Current_Status
  }
  
  sample_list = unique(data_run$sample_name)
  SampleNum = data.frame(matrix(ncol = 2, nrow = length(sample_list)))
  colnames(SampleNum) = c('sample','Cells')

  for (i in 1:length(sample_list)){
    sample = sample_list[i]
    print (sample)
    data_subset = data_run[,data_run$sample_name == sample]
    SampleNum$sample[i] = sample
    SampleNum$Cells[i]  = ncol(data_subset)
  }
  
  #SampleNum = SampleNum[order(SampleNum$Cells),]
  SampleNum = SampleNum[order(SampleNum$Cells),]
  
  write.csv(SampleNum, file=paste0(filepath_cluster,'/Stats/','SampleNum.csv'), row.names = F)
  
  umap = data_run@reductions[["umap"]]@cell.embeddings
  write.csv(umap, file=paste0(filepath_cluster,'Umap.csv'))
  
  cell_list = c('0', '2', '3', '4', '5', '6', '7', '9', '10', '11', '13', '22', '28', '37', '38')
  data_run_tcell = data_run[,Idents(data_run) %in% cell_list]
  
  data_run_NK = data_run[,Idents(data_run) %in% c('8','17')]
  
  cell_list = c('12', '14', '15', '16', '23', '25', '39', '41')
  data_run_mono = data_run[,Idents(data_run) %in% cell_list]
  
  data_matrix = data_run@assays[["RNA"]]@counts
  #data_matrix = data_matrix[rownames(data_matrix) %in% data_run@assays[["integrated"]]@var.features,]
  #data_matrix = NormalizeData(data_matrix, normalization.method = "LogNormalize", scale.factor = 10000)
  write.table(data_matrix, file='/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_MT15_All_AllFeatures.tsv', quote=FALSE, sep='\t')
  
  data_matrix = data_run_tcell@assays[["RNA"]]@counts
  data_matrix = data_matrix[rownames(data_matrix) %in% data_run@assays[["integrated"]]@var.features,]
  data_matrix = NormalizeData(data_matrix, normalization.method = "LogNormalize", scale.factor = 10000)
  write.table(data_matrix, file='/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_MT15_TCell.tsv', quote=FALSE, sep='\t')
  
  
  data_matrix = data_run_NK@assays[["RNA"]]@counts
  data_matrix = data_matrix[rownames(data_matrix) %in% data_run@assays[["integrated"]]@var.features,]
  data_matrix = NormalizeData(data_matrix, normalization.method = "LogNormalize", scale.factor = 10000)
  write.table(data_matrix, file='/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_MT15_NK.tsv', quote=FALSE, sep='\t')
  
  data_matrix = data_run_mono@assays[["RNA"]]@counts
  data_matrix = data_matrix[rownames(data_matrix) %in% data_run@assays[["integrated"]]@var.features,]
  data_matrix = NormalizeData(data_matrix, normalization.method = "LogNormalize", scale.factor = 10000)
  write.table(data_matrix, file='/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_MT15_Mono.tsv', quote=FALSE, sep='\t')
  
  
  # Label data
  data_run_label = label_cells(data_run,cluster_IDs)
  cell_list =  c('NK','mNK','Tgd','CD8+ T Cell','Treg','TEM','TEMRA','TCM','TSCM','U1', 'U2','U3')
  data_run_label_clean = data_run_label[,Idents(data_run_label) %in% cell_list]
  Idents(data_run_label_clean) = factor(Idents(data_run_label_clean) , levels = cell_list)
  
  
  cell_list =  c('NK','CD8+ T Cell','T Cell')
  data_run_label_clean = data_run_label[,Idents(data_run_label) %in% cell_list]
  Idents(data_run_label_clean) = factor(Idents(data_run_label_clean) , levels = cell_list)
  

  idents_unique = sort(unique(Idents(data_run)))
  for (i in idents_unique){
    
    data_subset = subset(data_run, idents = i)
    print(i)
    print( summary(data.frame(data_subset$Cell_Label)))
  }
  
  data_run_label_orig = data_run_label[,data_run_label$kit != 'Microwell-seq']
  for (i in (unique(Idents(data_run_label)))){
    
    print('')
    data_subset = subset(data_run_label, idents = i)
    print(i)
    data_summary = summary(data.frame(data_subset$Cell_Label))
    print(data_summary )
    data_subset_orig = data_subset[,data_subset$kit != 'Microwell-seq']
    
    sample_names = unique(data_subset_orig$sample_name)
    print('Number of cells')
    print(length(data_subset_orig$sample_name))
    
    print('Percentage')
    print(ncol(data_subset_orig)/ncol(data_run_label_orig)*100)
    
    #print('SM')
    #sum(str_count(sample_names, "GL"))
    #print('NBM')
    #sum(str_count(sample_names, "NBM"))
    print('Num Samples')
    print(length(sample_names))
  }
  
  browser()
  # Begin plotting data
  cell_features = getCellMarkers(folder_base)
  groupBy_list = c('dexa','sample_name', 'Response','split_var','kit','Patient')
  cell_list =  c("CD8+ T Cell","mNK","NK","TCM","TEM","TEMRA","Tgd","Treg","TSCM" )
  data_run_label_clean = data_run_label[,Idents(data_run_label) %in% cell_list]
  plotAll(data_run, folder = folder_base_output,
          sample_name,sampleParam,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = T, keepOldLabels = F, groupBy = groupBy_list, 
          PCA_dim = NA,resolution_val = NA)
  plotAll(data_run, folder = folder_base_output,
          sample_name = sample_name,sampleParam = sampleParam,
          cell_features = cell_features,
          label_TF = T,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF = F, markersTF = T, keepOldLabels = F, groupBy = groupBy_list)
  plotAll(data_run_label_clean, folder = folder_base_output,
          sample_name = sample_name,sampleParam = sampleParam,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF = F, markersTF = T, keepOldLabels = T, groupBy = groupBy_list, str = '_Clean')
  
  print('Done Plotting')
  browser()
  
############################
  
  celltype = data.frame(matrix(ncol = 4, 0))
  colnames(stats_summary_line) = c("Sample",'Cluster','Num','Percent')
  
  sample_list = unique(data_run_label$sample_name)
  celltype_iterate = as.character(unique(Idents(data_run_label)))
  patient_list = unique(data_run_label$Patient)
  patient = patient_list[1]
  celltype = celltype_iterate[1]
  
  stats_summary_line = data.frame(matrix(ncol = 4, 0))
  colnames(stats_summary_line) = c("Sample",'Cluster','Num','Percent')
  
  
  for (sample in sample_list){
    
    data_run_input = data_run_label
    data_run_input = data_run_input[,data_run_input$sample_name == sample]
    cluster_num = clusterStats(data_run_input)
    
    cluster_num$Sample = sample
    stats_summary_line = rbind(stats_summary_line,cluster_num)
    
  }
  stats_summary_line =merge(stats_summary_line, metaData, by = 'Sample')
  
  stats_summary_line_baseline = stats_summary_line[stats_summary_line$Treatment == 'baseline'| stats_summary_line$Treatment == 'NBM',]
  stats_summary_celltype = data.frame(matrix(ncol = 4, nrow = 46 ))
  colnames(stats_summary_celltype) = c('Patient','NBM','Cluster','Proportion')
  stats_summary_celltype$Patient = stats_summary_line_baseline$`Patient Number`
  stats_summary_celltype$NBM = stats_summary_line_baseline$`Treatment`
  stats_summary_celltype$Cluster = stats_summary_line_baseline$Cluster
  stats_summary_celltype$Proportion = stats_summary_line_baseline$Percent
  

  write.csdevtools::install_github("constantAmateur/SoupX")
v(stats_summary_celltype,'/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM/Regkit/SubsetIntegrate/NK/Analysis/stats_summary_celltype.csv')
  
  stats_summary_celltype = data.frame(matrix(ncol = 2, nrow = length(patient_list) ))
  colnames(stats_summary_celltype) = c('Patient','Time')
  rownames(stats_summary_celltype) = patient_list
  
  for (patient in patient_list){
    rowdata = stats_summary_line[stats_summary_line$`Patient Number` == patient,]
    for(i in 1:nrow(rowdata)){
      print(i)
      patient = (rowdata$`Patient Number`)[i]
      celltype = rowdata$Cluster[i]
      time = rowdata$Treatment[i]
      stats_summary_celltype[patient,'Patient'] = patient
      colname = paste0(celltype,' ', time)
      stats_summary_celltype[[colname]] = NA
      
    }
    
    
  }
  
  for (patient in patient_list){
    rowdata = stats_summary_line[stats_summary_line$`Patient Number` == patient,]
    for(i in 1:nrow(rowdata)){
      print(i)
      patient = (rowdata$`Patient Number`)[i]
      celltype = rowdata$Cluster[i]
      print(celltype)
      time = rowdata$Treatment[i]
      stats_summary_celltype[patient,'Patient'] = patient
      colname = paste0(celltype,' ', time)
      stats_summary_celltype[patient,colname] = rowdata$Percent[i]
      
    }
  }
  
  stats_summary_celltype = stats_summary_celltype[ , order(colnames(stats_summary_celltype))]
  write.csv(stats_summary_celltype, file = paste0(filepath_cluster,'Stats/stats_summary_celltype.csv'),row.names = T)
  
  stats_summary_celltype
  
  
  ####################
  # Plot known markers
  ####################
  #data_run = getCluster (data_run,resolution_val, PCA_dim)
  cell_features = getCellMarkers(folder_base)
  PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/HeatMap/'), cell_features = cell_features,
                   plotType ='HeatMap' , str = '')
  PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/Violin/'), cell_features = cell_features,
                   plotType ='Violin' , str = '')
  PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                   plotType ='FeaturePlotFix' , str = '')
  
  idents_unique = unique(Idents(data_run))
  
  data_run_subset_small = subset(data_run, idents = round(length(idents_unique)/2):(length(idents_unique) -1))
  data_run_subset_large = subset(data_run, idents = 0:round(length(idents_unique)/2))
  
  inc = round(length(idents_unique)/4)
  data_run_subset_small = subset(data_run, idents = (inc*3 + 1):(length(idents_unique)-1))
  data_run_subset_medium = subset(data_run, idents = (inc*2 + 1):(inc*3))
  data_run_subset_large = subset(data_run, idents = (inc-1):(inc*2))
  data_run_subset_extralarge = subset(data_run, idents = 0:(inc-2))
  PlotKnownMarkers(data_run_subset_small, folder = paste0(filepath_cluster,'Cell Type/HeatMap/Small Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F, str = '')
  PlotKnownMarkers(data_run_subset_medium, folder = paste0(filepath_cluster,'Cell Type/HeatMap/Medium Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F, str = '')
  PlotKnownMarkers(data_run_subset_large, folder = paste0(filepath_cluster,'Cell Type/HeatMap/Large Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F, str = '')
  PlotKnownMarkers(data_run_subset_extralarge, folder = paste0(filepath_cluster,'Cell Type/HeatMap/Extra Large Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F,str = '')
  
  data_run_subset = subset(data_run, idents = c('0','1','2','3','4','5','6','7','9','17','20'))
  data_run_subset = data_run_label_clean
  pathName = paste0(folder_base_output,
                     paste0('HeatMapTCellMarkers_Clean','.png'))
  png(file=pathName,width=2000, height=2500, res = 100)
  plot = DoHeatmap(object = data_run_subset, features = gene_list,assay = 'RNA', slot = "data",
                   group.by = "ident", label = T)
  plot = plot + theme(
    axis.text= element_text(color="black", size=38))
  print(plot)
  dev.off()
  
  
  pathName <- paste0(folder_base_output,
                     paste0('ViolinTCellMarkers8_14','.png'))
  png(file=pathName,width=1000, height=1000, res = 100)
  plot  = StackedVlnPlot(obj = data_run_subset, features = gene_list ) +
    ggtitle('' ) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=24))
  plot = plot + theme(
    axis.text= element_text(color="black", size=24))
  print(plot)
  dev.off()

  gene_list = c("PTPRC", "HNRNPLL", "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","SELL","CCR7","CD27",
                "CD28", "CD69", "IL2RA", "IL2RB", "IL2RG", "CD38", "FAS", "IL7R", "KLRG1", "ZEB1","ZEB2",
                "ZEB2", "PRF1", "GNLY", "NKG7","FCGR3A", "ITGAL", "CX3CR1", "B3GAT1", "BCL2", "MCL1", "LEF1",
                "TRDC", "TRGC1", "TRGC2", "TRAV10", "KLRB1","LAMP1","TRAC",'TCF7',"FOXP3",'IL10', 'TGFB1', 
                'CTLA4', 'TNFRSF18', 'LTB','NOSIP','NTDP1','GZMA','GZMB','GZMK','GZMH','GZMM',
                'CCL3','IFNG','KLRD1','ITGAM','HAVCR2','LAG3','PDCD1','TIGIT','TBX21','ADGRG1',
                'NCAM1','SRGN',"HLA-DRA" , "HLA-DRB5", "HLA-DRB1","ENTB1","TRDC") 

  gene_list = c( "CD3D","CD3E","CD3G","CD4","CD8A","CD8B",'NCAM1','CX3CR1', 'SELL', 'CXCR4','GZMA','GZMB','GZMK','GZMH','GZMM',
                'NKG7', 'PRF1', 'GNLY', 'FCGR3A', 'ITGAL', 'KLRG1', 'KLRB1',
                'TRDC', 'ZEB1', 'ZEB2', 'IL7R','CD27','TRDC','TRGC1','TRGC2','TCF7')



  ##################
  #data_run_NKGenes = data_run[,Idents(data_run) != '6']
  #VariableFeatures(data_run_NKGenes) = gene_list
  #file_str = '_NKGenes_removeTgd_label'
  #data_run_NKGenes = ScaleData(data_run_NKGenes)
  #data_run_NKGenes = RunPCA(data_run_NKGenes, npcs = 10)
  #data_run_NKGenes = FindNeighbors(data_run_NKGenes, dims = 1:10)
  #data_run_NKGenes = FindClusters(data_run_NKGenes, resolution = 0.5) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.)
  #data_run_NKGenes = RunUMAP(data_run_NKGenes, dims = 1:10)
  
  #data_run_NKGenes = data_run_NKGenes[,Idents(data_run_NKGenes)!= '4']
  #data_run_NKGenes = ScaleData(data_run_NKGenes)
  #data_run_NKGenes = RunPCA(data_run_NKGenes, npcs = 10)
  #data_run_NKGenes = FindNeighbors(data_run_NKGenes, dims = 1:10)
  #data_run_NKGenes = FindClusters(data_run_NKGenes, resolution = 2) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.)
  #data_run_NKGenes = RunUMAP(data_run_NKGenes, dims = 1:10)
  
  #data_run_NKGenes = label_cells(data_run_NKGenes,'CX3CR1+GZMK+, CXCR4+GZMK-, CXCR4+GZMK+, SELL+GZMK+,  CX3CR1+ GZMK-, CXCR4+GZMK-')
  #pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',10,'_res',2,file_str,'.png'))
  #png(file=pathName,width=1000, height=1000, res = 100)
  #print(DimPlot(data_run_NKGenes,pt.size = 0.8, reduction = "umap",label = TRUE))
  #dev.off()
  

 
  ###################
  file_str = ''
  gene_list = unique(gene_list)
  folder = paste0(filepath_cluster,'Cell Type/HeatMap Gene/') 
  dir.create(folder,recursive = T)
  pathName = paste0(folder,'NKGenes','_heatmap',file_str,'.png')
  pathName = paste0(folder,'LymphoGenes','_heatmap','',file_str,'.png')
  
  plot  = DoHeatmap(object = data_run, features = gene_list,assay = 'RNA', slot = "data",
                    group.by = "ident", label = T) +
    ggtitle('' ) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=24))
  
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=12),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20))
  
  
  
  
  png(file=pathName,width=1000, height=1000,res=100)
  print(plot)
  dev.off()
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/T Cell/') 
    folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/NK/') 
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_run, feature = gene, folder =folder,
                          str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE) 
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=24),
      legend.title=element_text(size=24),
      text = element_text(size = 20)
    )
    
    
    pathName = paste0(folder,gene,file_str,'.png')
    pathName = paste0(folder,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  ##################################
  pathName <- paste0(folder_base_output,
                     paste0('FeaturePlotTCellMarkersSmall','.png'))
  png(file=pathName,width=2000, height=2500, res = 100)
  plot = DoHeatmap(object = data_run_subset_small, features = gene_list,assay = 'RNA', slot = "data",
                   group.by = "ident", label = T)
  plot = plot + theme(
    axis.text= element_text(color="black", size=42))
  print(plot)
  dev.off()
  
  
  PlotKnownMarkers(data_run_subset_large, paste0(filepath_cluster,'Cell Type/FeaturePlot/Large Clusters/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = F,featurePlotFix = F,str = '')
  PlotKnownMarkers(data_run_subset_small, paste0(filepath_cluster,'Cell Type/FeaturePlot/Small Clusters/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = F,featurePlotFix = F,str = '')
  
  PlotKnownMarkers(data_run,data_run, paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = F,featurePlotFix = TRUE,str = '_FIX')
  
  
  data_run_label = label_cells(data_run,cluster_IDs)
  
  
  PlotKnownMarkers(data_run_label,data_run_label, folder = paste0(filepath_cluster,'Cell Type/HeatMapLabel/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F,featurePlotFix = TRUE, str = '')
  
  PlotKnownMarkers(data_run_label,data_run_label, folder = paste0(filepath_cluster,'Cell Type/ViolinLabel/'), cell_features = cell_features,
                   plotType ='Violin' ,split = F,featurePlotFix = TRUE, str = '')
  
  idents_unique = levels(unique(Idents(data_run_label)))
  
  data_run_subset_small = subset(data_run_label, idents = idents_unique[round(length(idents_unique)/2):(length(idents_unique) -1)])
  data_run_subset_large = subset(data_run_label, idents = idents_unique[0:round(length(idents_unique)/2)])
  
  inc = round(length(idents_unique)/3)
  data_run_subset_small = subset(data_run_label, idents = idents_unique[(inc*2):(length(idents_unique)-1)])
  data_run_subset_medium = subset(data_run_label, idents = idents_unique[inc:(inc*2)])
  data_run_subset_large = subset(data_run_label, idents = idents_unique[0:(inc-1)])
  
  PlotKnownMarkers(data_run_subset_small, folder = paste0(filepath_cluster,'Cell Type/HeatMapLabel/Small Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F,featurePlotFix = TRUE, str = '')
  PlotKnownMarkers(data_run_subset_medium, folder = paste0(filepath_cluster,'Cell Type/HeatMapLabel/Medium Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F,featurePlotFix = TRUE, str = '')
  PlotKnownMarkers(data_run_subset_large, folder = paste0(filepath_cluster,'Cell Type/HeatMapLabel/Large Clusters/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = F,featurePlotFix = TRUE, str = '')
  
  ############
  
  data_run_subset = data_run[,Idents(data_run) %in% c('1','4', '18', '34', '4', '3', '0')]
  num_markers = 100
  file_str = ''
  file = paste0(filepath_cluster,'Features',file_str,'.csv')
  markers = read.csv(file)
  markers = markers[markers$cluster == '1',]
  
  str = '_Cluster1'
  markers_small = markers[markers$cluster %in% round(length(idents_unique)/2):(length(idents_unique) -1) ,]
  markers_large = markers[markers$cluster %in% 1:round(length(idents_unique)/2) ,]
  TopNHeatmap(data_run_subset,markers = markers,filepath_cluster,PCA_dim,resolution_val,num_markers,file_str = paste0('_',num_markers,str))
  TopNHeatmap(data_run_subset_small,markers = markers_small,filepath_cluster,PCA_dim,resolution_val,num_markers,file_str = paste0('smallCluster_',num_markers,str))
  TopNHeatmap(data_run_subset_large,markers = markers_large,filepath_cluster,PCA_dim,resolution_val,num_markers,file_str = paste0('largeCluster_',num_markers,str))
  
  file_str = '_label'
  file = paste0(filepath_cluster,'Features',file_str,'.csv')
  markers = read.csv(file)
  TopNHeatmap(data_run_label,markers = markers,filepath_cluster,PCA_dim,resolution_val,num_markers,file_str = paste0('label_',num_markers))
  
  
  ## Plot Markers from DE between each cluster and all other clusters combined
  ## The files are saved from plotAll when markersTF == T
  
  
  #
  # Saving subset of data if desired
  if (saveClean){
    # cluster_IDs_list = unlist(strsplit(cluster_IDs, ",")) 
    # cluster_IDs_list = trimws(cluster_IDs_list, which = c("both"), whitespace = "[ \t\r\n]")
    # #cluster_IDs_new = cluster_IDs_list[cluster_IDs_list!='Erythrocyte' & cluster_IDs_list!='Ignore','cluster_IDs_new']
    # cluster_IDs_new = cluster_IDs_list[cluster_IDs_list!='CD14+ Mono']
    # data_run_clean = subset(data_run_label, idents = cluster_IDs_new)
    # data_run_clean = FindVariableFeatures(data_run_clean, selection.method = "vst", nfeatures = 2000)
    # # 
    # path = paste0(folder_base_output,
    #               '/data_run_clean_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,'.Robj')
    # save(data_run_clean,file= path)
    
    
    cluster_IDs_list = unlist(strsplit(cluster_IDs, ",")) 
    cluster_IDs_list = trimws(cluster_IDs_list, which = c("both"), whitespace = "[ \t\r\n]")
    #cluster_IDs_new = cluster_IDs_list[cluster_IDs_list!='Erythrocyte' & cluster_IDs_list!='Ignore','cluster_IDs_new']
    remove_cluster_list = c('CD14+ Mono','Mast','CXCL9 High','Basal','DC')
    cluster_IDs_new = cluster_IDs_list[!(cluster_IDs_list %in% remove_cluster_list)]
    
    
    patient_complete_list = c('5','6','12','16','20','30','31','40','51')
    data_run_clean = data_run_label[,data_run$Patient %in% patient_complete_list]
    data_run_clean = subset(data_run_clean, idents = cluster_IDs_new)
    
    data_run_clean = subset(data_run_label, idents = c('T Cell','CD8+ T Cell','NK'))
    path = paste0(folder_base_output,
                  '/data_run','_clean','_PCAdim',PCA_dim,'_','features2000','.Robj')
    save(data_run_clean,file= path)
    
    
  }
  
  
  ################
  ## DE Between conditions
  ################
  sortby = 'avg_logFC'
  gene_num = c(40,40)
  category_list = c('Pre ND','Pre D', 'Post ND', 'Post D')
  cluster_ids_new = gsub("CD14+ Mono1", "CD14+ Mono", cluster_ids_new)

  celltype_iterate = c('NK', 'T Cell','CD8+ T Cell','CD14+ Mono')
  #celltype_iterate = c('CD14+ Mono')
  
  # 
  split_var = 'split_var'
  # Split by pre/post + response
  # Compare markers using pre/post response
  celltype_iterate = c('NK', 'T Cell','CD8+ T Cell','CD14+ Mono','CD16+ Mono')
  #celltype_iterate = c('NK', 'T Cell','CD14+ Mono')
  #celltype_iterate = c('CD8+ T Cell')
  #celltype_iterate = c('NK', 'T Cell','CD8+ T Cell','CD14+ Mono')
  #data_run_remove10 = data_run[,data_run$Patient !='10']
  #data_run_remove21 = data_run[,data_run$Patient !='21']
  # Go through each cell type and compare conditions in that cell type
  for (celltype in celltype_iterate){
    print(celltype)
    celltype_list = c(celltype)
    
    
    #########################
    data_run_label = label_cells(data_run, cluster_ids_new)
    #data_run_label = RenameIdents(object = data_run_label, 'dCD14+ Mono' = 'CD14+ Mono')
    #data_run_label = data_run
    
    data_run_label@meta.data$split_var = gsub("Pre ND", "Pre", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Pre D", "Pre", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Post ND", "Post", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Post D", "Post", data_run_label@meta.data$split_var)
    data_run_label@meta.data$Response = gsub("MRP", "TMP", data_run_label@meta.data$Response)
    data_run_label@meta.data$Response = gsub("MR", "TMP", data_run_label@meta.data$Response)
    data_run_label@meta.data$Response = gsub("TMP", "MR", data_run_label@meta.data$Response)
    
    data_run_label$response_split_var = paste0(data_run_label$split_var ,' ',data_run_label$Response )
    
    data_run_input = data_run_label
    Idents(data_run_input) = paste0(Idents(data_run_input),' ', 
                                    data_run_input@meta.data[,'response_split_var'])
    
    
    ident1 = paste0(celltype, ' ', 'baseline',' ','Poor Response')
    ident2 = paste0(celltype, ' ', 'baseline',' ','Good Response')
    
    folder_name = paste0('DoHeatmap/Response/',ident1, '_', ident2)
    
    folder_name = paste0('',ident1, '_', ident2)
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    #folder_heatMap = paste0(folder_base_output,'', folder_name,'/')
    Features = FindMarkers(data_run_input, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    cluster_ids_new = cluster_IDs
    data_run_input = data_run_label
    DoHeatMapHelper(data_run_input,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = 'response_split_var',
                    cluster_ids_new,cellPair = FALSE,gene_num,str = '')
    
    
    
    
    ################################        
    category_list = unique(data_run_label$response_split_var)
    for (category in category_list){
      data_run_input = data_run_label
      data_run_input = SubsetData(object = data_run_input, cells = data_run_input$response_split_var == category )
      
      DoHeatMapHelper(data_run_input,folder_base_output,folder_heatMap, Features,
                      ident1,ident2,celltype_list,category_list,sortby,split_var = 'sample_name',
                      cluster_ids_new,cellPair = FALSE,gene_num,str = paste0('_',category))
      
      Idents(data_run_input) = paste0(Idents(data_run_input),' ', 
                                      data_run_input@meta.data[,'response_split_var'])
      #stats = clusterStats(data_run_input)
      #write.csv(stats, file = paste0(folder_heatMap,'Stats_',category,'.csv'),row.names = FALSE)
      
    }
    
    data_run_input = data_run_label
    Idents(data_run_input) = paste0(Idents(data_run_input),' ', 
                                    data_run_input@meta.data[,'split_var'])
    
    ident1 = paste0(celltype, ' ', 'C9D1')
    ident2 = paste0(celltype, ' ', 'EOT')
    
    folder_name = paste0('DoHeatmap/Response/',ident1, '_', ident2)
    folder_name = paste0('',ident1, '_', ident2)
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    #folder_heatMap = paste0(folder_base_output,'', folder_name,'/')
    Features = FindMarkers(data_run_input, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    
    data_run_input = data_run_label
    data_run_input$response_split_var[data_run_input$response_split_var == 'NBM NBM'] = 'NBM'
    DoHeatMapHelper(data_run_input,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = 'split_var',
                    cluster_ids_new,cellPair = FALSE,gene_num,str = '', Ident_order = NA)
    
    Ident_order = c('NBM','Pre MR','Pre VGPR','Post MR', 'Post VGPR')
    Ident_order = paste0(celltype,' ', Ident_order)
    # DoHeatMapHelper(data_run_input,folder_base_output,folder_heatMap, Features,
    #                 ident1,ident2,celltype_list,category_list,sortby,split_var = 'response_split_var',
    #                 cluster_ids_new,cellPair = FALSE,gene_num,str = '', Ident_order = Ident_order)
    
    category_list = unique(data_run_label$response_split_var)
    for (category in category_list){
      data_run_input = data_run_label
      data_run_input = SubsetData(object = data_run_input, cells = data_run_input$response_split_var == category )
      
      DoHeatMapHelper(data_run_input,folder_base_output,folder_heatMap, Features,
                      ident1,ident2,celltype_list,category_list,sortby,split_var = 'sample_name',
                      cluster_ids_new,cellPair = FALSE,gene_num,str = paste0('_',category))
      Idents(data_run_input) = paste0(Idents(data_run_input),' ', 
                                      data_run_input@meta.data[,'response_split_var'])
      stats = clusterStats(data_run_input)
      #write.csv(stats, file = paste0(folder_heatMap,'Stats_',category,'.csv'),row.names = FALSE)
    }
    
  }
  
  ############################################################
  ## Check which clusters have only one sample
  
  ignore_cluster_list = c()
  for (cluster in unique(data_run@active.ident)){
    print('Cluster')
    print(cluster)
    data_subset = subset(data_run, idents = cluster)
    patient_list = unique(data_subset$Patient)
    sample_list = unique(data_subset$sample_name)
    if (length(sample_list) == 1){
      print('Patient List')
      print(sample_list)
      print('Cluster with only 1 patient')
      ignore_cluster_list = c(ignore_cluster_list,cluster)
    }
    print('')
    
  }
  
  
  #####
  # Score for CD14 markers in Oksana paper
  m2_markers = c('HLA-DPB1','HLA-DPA1','TAGLN2','HLA-DRB5','F13A1','LIPA','HLA-DQA1')
  m7_markers = c('ZFP36L1','PSME2','IFITM2','PLAC8','APOBEC3A','TNFSF10','LY6E','ISG15','IFITM3','IFI44L')
  
  data_run <- AddModuleScore(
    object = data_run,
    features = m2_markers,
    name = 'm2_markers',
    assay = 'RNA'
    
  )
  
  data_run <- AddModuleScore(
    object = data_run,
    features = m7_markers,
    name = 'm7_markers',
    assay = 'RNA'
  )
  
  print(FeaturePlot(data_run,pt.size = 0.5, features = c("m2", "m7")))
  
  ## 
  # Permute for p values
  ##
  perm_df = data.frame(matrix(ncol = 2, nrow = length( unique(cluster_IDs_list))))
  colnames(perm_df) = c("PercentDiff",'p_val')
  
  
  cluster_IDs_list = unlist(strsplit(cluster_IDs, ",")) 
  cluster_IDs_list = trimws(cluster_IDs_list, which = c("both"), whitespace = "[ \t\r\n]")
  cluster_IDs_list = unique(Idents(data_run_label))
  rownames(perm_df) = cluster_IDs_list
  for (j in 1:length(unique(cluster_IDs_list))){
    cell = unique(cluster_IDs_list)[j]
    print(cell)
    diff_result = diff_exp (data_run_label, cell)
    print(diff_result)
    perm_df[j,] = diff_result
  }
  
  write.csv(perm_df, file = paste0(filepath_cluster,'Stats/perm_pval.csv'),row.names = T)
 
  
  
  #browser()
  ####
  # cell_features = getCellMarkers()
  # ident1 = '4'
  # ident2 = '5'
  # Features_mvn = FindMarkers(data_run, ident.1 = ident1, ident.2 = ident2
  #                            ,min.pct = 0.1, logfc.threshold = 0.01, only.pos = TRUE)
  # Features_mvn$gene = rownames(Features_mvn)
  # rownames(Features_mvn) = NULL
  # Features_mvn = cellMarkers(Features_mvn,cell_features)
  # filepath_mvn = paste0( filepath_cluster, 'DE/nVsm/')
  # filepath = paste0(filepath_mvn
  #                   ,'Features_',ident1,'Vs',ident2
  #                   ,'.csv')
  # write.csv(Features_mvn, file = filepath,row.names = FALSE)
  # 
  
  #data_run_label@meta.data$split_var = gsub("Pre ND", "Pre", data_run_label@meta.data$split_var)
  #data_run_label@meta.data$split_var = gsub("Pre D", "Pre", data_run_label@meta.data$split_var)
  #data_run_label@meta.data$split_var = gsub("Post ND", "Post", data_run_label@meta.data$split_var)
 


  # 

  #PlotKnownMarkers(data_run_label,data_run, paste0(filepath_cluster,'Cell Type/'), cell_features = NA,
  #                plotType ='HeatMap' ,split = FALSE,featurePlotFix = TRUE, str = 'label')
  
  #browser()
  
  #cell_features = c('CD34','CDK6','CD24', 'CD9','MZB1', 'GYPA', 'HBB','MZB1', 'TIGIT', 'LAG3','CD8A','CD8B')
  
  cell_features = c('FCER1A','ITGAX','CD83','THBD','CD209','CD1C','LYZ','MZB1')
  cell_features = c('IL3RA','CLEC4C','NRP1', 'MZB1')
  cell_features = c('CD19','MS4A1')
  
  
  cell_features = c('CD27','CD38','SDC1','SLAMF7','IL6','CD138','TNFRSF17') # pDC
  cell_features = c('CD21','CD5','CD10','PAX5','CD24','CD34','CD93', 'CD23'
                    ,'CD19','CD21','MS4A1','CD27','IL10','CD1D','TGFB1','TGFB2','EBI3','PRDM1','IGHM','CD1D','CD4') # B Cell 

  cell_features = c('CD10')
  
  cell_features = c('FCGR3A')
  #cell_features = gene_list[13:22]
  
  celltype = 'NK'
  data_run_subset = data_run[,Idents(data_run) == celltype]
  
  data_run_subset = data_run
  data_run_subset@meta.data$Response = gsub("MRP", "TMP", data_run_subset@meta.data$Response)
  data_run_subset@meta.data$Response = gsub("MR", "TMP", data_run_subset@meta.data$Response)
  data_run_subset@meta.data$Response = gsub("TMP", "MR", data_run_subset@meta.data$Response)
  data_run_subset = data_run_subset[,data_run_subset@meta.data$Response != 'NBM']
  Idents(data_run_subset) = paste0(Idents(data_run_subset), ' ',data_run_subset@meta.data$Response)
  
  Idents(data_run_subset) = factor(Idents(data_run_subset) , levels = paste0(celltype,c(' MR',' VGPR')))
  
  data_run_subset_VGPR = data_run_subset[,data_run_subset$Response == 'VGPR']
  data_run_subset_MR = data_run_subset[,data_run_subset$Response == 'MR']
  
  PlotKnownMarkers(data_run_subset_VGPR,data_run_subset_VGPR, paste0(filepath_cluster,'Cell Type/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = TRUE, str = paste0('_VGPR'))
  PlotKnownMarkers(data_run_subset_MR,data_run_subset_MR, paste0(filepath_cluster,'Cell Type/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = TRUE, str = paste0('_MR'))
  
  PlotKnownMarkers(data_run_subset,data_run_subset, paste0(filepath_cluster,'Cell Type/'), cell_features = cell_features,
                   plotType ='HeatMap' ,split = FALSE,featurePlotFix = TRUE, str = paste0(celltype,' ', cell_features))

  data_pre = data_run[,data_run@meta.data$orig.ident == "data_pre"]
  data_post = data_run[,data_run@meta.data$orig.ident == "data_post"]
  PlotKnownMarkers(data_pre,data_pre, paste0(filepath_cluster,'Cell Type/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = TRUE, str = '_Pre')
  PlotKnownMarkers(data_post,data_post, paste0(filepath_cluster,'Cell Type/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = TRUE, str = '_Post')
  
  PlotKnownMarkers(data_run,data_run, paste0(filepath_cluster,'Cell Type/'), cell_features = cell_features,
                   plotType ='FeaturePlot' ,split = FALSE,featurePlotFix = TRUE, str = '_Post')
  

  
  #data_run = getCluster (data_run,resolution_val, PCA_dim)
  data_run_label = label_cells(data_run,cluster_IDs)
  #browser()
  ###############################################################################
  
  for (celltype in celltype_iterate){
    celltype_list = c(celltype)
    data_run_label = label_cells(data_run, cluster_ids_new)

    
    Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run_label@meta.data[,'split_var'])
  
    ident1 = paste0(celltype, ' ', 'Post ND')
    ident2 = paste0(celltype, ' ', 'Post D')
    folder_name = paste0('DoHeatmap/',ident1, '_', ident2)
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    data_run_label = label_cells(data_run, cluster_ids_new)

    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = split_var,
                    cluster_ids_new,cellPair = FALSE,gene_num,str = '')
    
    ident1 = paste0(celltype, ' ', 'Pre ND')
    ident2 = paste0(celltype, ' ', 'Post ND')
    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = split_var,
                    cluster_ids_new,cellPair = TRUE,gene_num,str = '')
    
    ident1 = paste0(celltype, ' ', 'Pre D')
    ident2 = paste0(celltype, ' ', 'Post D')
    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = split_var,
                    cluster_ids_new,cellPair = TRUE,gene_num,str = '')
    
    ident1 = paste0(celltype, ' ', 'Post ND')
    ident2 = paste0(celltype, ' ', 'Post D')
    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = split_var,
                    cluster_ids_new,cellPair = TRUE,gene_num,str = '')
    
  }
  
  # Split by pre/post + response
  # Use D/ND markers
  split_var = 'response_split_var'
  for (celltype in celltype_iterate){
    celltype_list = c(celltype)
    data_run_label = label_cells(data_run, cluster_ids_new)
    
    
    Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run_label@meta.data[,'split_var'])
    
    ident1 = paste0(celltype, ' ', 'Pre D')
    ident2 = paste0(celltype, ' ', 'Post D')
    folder_name = paste0('DoHeatmap/Response/',ident1, '_', ident2)
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    data_run_label = label_cells(data_run, cluster_ids_new)
    data_run_label@meta.data$split_var = gsub("Pre ND", "Pre", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Pre D", "Pre", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Post ND", "Post", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Post D", "Post", data_run_label@meta.data$split_var)
    data_run_label$response_split_var = paste0(data_run_label$split_var ,' ',data_run_label$Response )
    
    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = split_var,
                    cluster_ids_new,cellPair = FALSE,gene_num,str = '')
    
    
  }
  # Split by pre/post + response
  # Don't Use D/ND markers
  #data_run_label = RenameIdents(object = data_run_label, 'dCD14+ Mono' = 'CD14+ Mono')
  for (celltype in celltype_iterate){
    print(celltype)
    celltype_list = c(celltype)
    #data_run_label = label_cells(data_run, cluster_ids_new)
    
    #
    
    data_run_label@meta.data$split_var = gsub("Pre ND", "Pre", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Pre D", "Pre", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Post ND", "Post", data_run_label@meta.data$split_var)
    data_run_label@meta.data$split_var = gsub("Post D", "Post", data_run_label@meta.data$split_var)
    
    Idents(data_run_label) = paste0(Idents(data_run_label),' ', 
                                    data_run_label@meta.data[,'split_var'])

    ident1 = paste0(celltype, ' ', 'Pre')
    ident2 = paste0(celltype, ' ', 'Post')
    folder_name = paste0('DoHeatmap/Response/',ident1, '_', ident2)
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                           ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
    data_run_label = label_cells(data_run, cluster_ids_new)
    data_run_label$split_var = gsub("Pre ND", "Pre", data_run_label$split_var)
    data_run_label$split_var = gsub("Pre D", "Pre", data_run_label$split_var)
    data_run_label$split_var = gsub("Post ND", "Post", data_run_label$split_var)
    data_run_label$split_var = gsub("Post D", "Post", data_run_label$split_var)
    
    data_run_label$response_split_var = paste0(data_run_label$split_var ,' ',data_run_label$Response )
    
    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = 'response_split_var',
                    cluster_ids_new,cellPair = FALSE,gene_num,str = '')

    DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                    ident1,ident2,celltype_list,category_list,sortby,split_var = 'split_var',
                    cluster_ids_new,cellPair = FALSE,gene_num,str = '')
    
  }
  #############################

  ################################
  # Make heatmaps with specific genes
  celltype = 'CD14+ Mono'
  celltype = 'NK'
  celltype = 'CD16+ Mono'
  celltype = 'CD8+ T Cell'
  celltype = 'T Cell'
  celltype_list = c(celltype)
  
  sortby = 'avg_logFC'
  gene_num = c(40,40)
  

  folder_name = paste0('DoHeatmap/Response/',celltype,' Genes/')
  folder_name = paste0('HeatMap/',celltype, '')
  folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
  
  data_run_label = label_cells(data_run, cluster_IDs)
  data_run_label = RenameIdents(object = data_run_label, 'dCD14+ Mono' = 'CD14+ Mono')
  
  data_run_label@meta.data$split_var = gsub("Pre ND", "Pre", data_run_label@meta.data$split_var)
  data_run_label@meta.data$split_var = gsub("Pre D", "Pre", data_run_label@meta.data$split_var)
  data_run_label@meta.data$split_var = gsub("Post ND", "Post", data_run_label@meta.data$split_var)
  data_run_label@meta.data$split_var = gsub("Post D", "Post", data_run_label@meta.data$split_var)
  data_run_label@meta.data$Response = gsub("MRP", "TMP", data_run_label@meta.data$Response)
  data_run_label@meta.data$Response = gsub("MR", "TMP", data_run_label@meta.data$Response)
  data_run_label@meta.data$Response = gsub("TMP", "MR", data_run_label@meta.data$Response)
  
  data_run_label$response_split_var = paste0(data_run_label$split_var ,' ',data_run_label$Response )
  
  ident1 = paste0(celltype, ' ', 'Pre MR')
  ident2 = paste0(celltype, ' ', 'Pre VGPR')
  
  #ident1 = paste0(celltype, ' ', 'Pre')
  #ident2 = paste0(celltype, ' ', 'Post')
  
  Ident_order = c('Pre MR','Post MR','Pre VGPR', 'Post VGPR')
  Ident_order = paste0(celltype,' ', Ident_order)
  category_list = unique(data_run_label$response_split_var)
  
  gene_list =  c('FCER1G','CD302','ILF2','MYCBP2') # Mono pre VGPR Pre MR
  gene_list =  c('IRF8','AHNAK','DDAH2','PIK3R1','SKAP1','RIOK3','DUSP2') # NK pre MR Pre VGPR
  gene_list = c('NKG7','CMTM6') # NK Pre Post
  gene_list = c('PA2G4') # NK MR Pre Post VGPR Pre Post
  gene_list = c('FOSB', 'KLF2', 'PPP4R3A', 'NKG7', 'CD58', 'SERB1') # CD8+ T Cells Pre MR Pre VGPR
  #gene_list = c('HIPK2', 'ZFR') # CD16+
  gene_list = c('FOSB', 'IQGAP1', 'HLA-DPA1', 'ITGB1', 'TRAM1')
  
  Features=  data.frame(matrix(ncol = 0, nrow = length(gene_list)))
  rownames(Features) =gene_list
  Features$p_val_adj = integer( length(gene_list))
  Features$avg_logFC = integer( length(gene_list))
  split_var = 'response_split_var'
  #split_var = 'split_var'
  DoHeatMapHelper(data_run_label,folder_base_output,folder_heatMap, Features,
                  ident1,ident2,celltype_list,category_list,sortby,split_var = split_var,
                  cluster_IDs,cellPair = FALSE,gene_num,str = '', Ident_order = Ident_order)
  

  split_var = 'response_split_var'
  data_run_label_subset_cell = subset(data_run_label, idents = celltype_list)
  data_run_label_subset_cell =data_run_label_subset_cell[,data_run_label_subset_cell$response_split_var != 'NBM NBM']
  Idents(data_run_label_subset_cell) = data_run_label_subset_cell@meta.data[,split_var]

  ident_list = c('Pre MR','Post MR','Pre VGPR','Post VGPR')
  ident_list = c('Pre MR','Pre VGPR')
  data_run_label_subset = subset(data_run_label_subset_cell, idents = ident_list)
  Idents(data_run_label_subset) = factor(Idents(data_run_label_subset) , levels = ident_list)

  str = ''
  
  if (all(ident_list == c('Pre MR','Post MR','Pre VGPR','Post VGPR'))){
    color_list = c('Blue','Blue','Red','Red')
  }else{
    color_list = c('Blue','Red')
  }
  for (gene in gene_list){
    folder_name = paste0('ViolinPlots/',celltype)
    folder = paste0(folder_base_output,'Analysis/', folder_name,'/')
    dir.create(folder,recursive = TRUE)
    
    pathName = as.character(paste0(folder,paste0('Violin ',celltype,' Split_',split_var,str,gene,'.png')) )
    png(file=pathName,width=900, height=600)
    plot= VlnPlot(data_run_label_subset, gene, pt.size = 0.1,cols =color_list)
    plot = plot + theme(
      title =element_text(size=24, color="black"),
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18),
      legend.position = 'none'
    )
    
    print(plot)
    dev.off()
    
    folder_name = paste0('Ridge Plots/',celltype)
    folder = paste0(folder_base_output,'Analysis/', folder_name,'/')
    dir.create(folder,recursive = TRUE)
    pathName = as.character(paste0(folder,paste0('Ridge ',celltype,' Split_',split_var,str,gene,'.png')) )
    png(file=pathName,width=900, height=600)
    plot= RidgePlot(data_run_label_subset, features = gene,cols = color_list)
    plot = plot + theme(
      title =element_text(size=24, color="black"),
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18),
      legend.position = 'none'
    )
    
    print(plot)
    dev.off()
    
    
    
  }
  ######################################33
  # CompareCellNum
  data_run_label = data_run
  data_run_label = label_cells(data_run, cluster_ids_new)
  data_run_label = RenameIdents(object = data_run_label, 'dCD14+ Mono' = 'CD14+ Mono')
  
  data_run_label@meta.data$split_var = gsub("Pre ND", "Pre", data_run_label@meta.data$split_var)
  data_run_label@meta.data$split_var = gsub("Pre D", "Pre", data_run_label@meta.data$split_var)
  data_run_label@meta.data$split_var = gsub("Post ND", "Post", data_run_label@meta.data$split_var)
  data_run_label@meta.data$split_var = gsub("Post D", "Post", data_run_label@meta.data$split_var)
  data_run_label@meta.data$Response = gsub("MRP", "TMP", data_run_label@meta.data$Response)
  data_run_label@meta.data$Response = gsub("MR", "TMP", data_run_label@meta.data$Response)
  data_run_label@meta.data$Response = gsub("TMP", "MR", data_run_label@meta.data$Response)
  
  data_run_label$response_split_var = paste0(data_run_label$split_var ,' ',data_run_label$Response )
  folder_stats = paste0(folder_base_output, 'Analysis/Stats/')
  CompareCellNum(data_run_label,folder_stats,split_var = 'split_var',metaData = metaData)
    
  
  stats = data.frame(table(data_run_label$sample_name))
  names(stats)[names(stats) == "Var1"] <- "Sample"
  stats =merge(stats, metaData, by = 'Sample')
  
  pathName =  paste0(folder_stats,'boxplot_', 'Cell Num','.png')
  x_name = 'Treatment'
  y_name = 'Freq'
  png(file=pathName,width=600, height=600)
  plot = ggplot(stats, aes(x = !!ensym(x_name), y = !!ensym(y_name), fill = stats$`10X kit`)) + 
    geom_boxplot()+
    coord_cartesian(ylim = c(0, 3000))+
    ggtitle('Cell Num')+
    xlab("") + ylab("%") 
  # Box plot with dot plot
  plot = plot + geom_jitter(shape=16, position=position_jitter(0.2))
  
  plot = plot + theme(
    plot.title = element_text(color="black", size=24, face="bold.italic"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=24),
  )
  plot = plot + theme(plot.title = element_text(hjust = 0.5))
  print(plot)
  dev.off()
  
  
  plot = VlnPlot(data_run_label, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 4,group.by = 'split_var', point.size.use = 0.0001)
  pathName <- paste0(folder_stats,'violin.png')
  print(pathName)
  png(file=pathName,width=1000, height=600)
  print(plot)
  dev.off()
  
  
  browser()

  celltype = 'NK1'
  celltype_list = c('NK1','NK2')
  ident1 = paste0(celltype, ' ', 'Post D')
  ident2 = paste0(celltype, ' ', 'Post D')
  DoHeatMapHelper(data_run,folder_base_output,folder_heatMap = NA, Features = NA, 
                  ident1,ident2,celltype_list,category_list,sortby,cluster_IDs,cellPair = FALSE)
  
 
  
  ################################################# 
 
  data_run_label = label_cells(data_run,cluster_ids_new)
  
  Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run@meta.data$split_var)
  
  ident1 = paste0('T Cell', ' ', 'Post ND')
  ident2 = paste0('T Cell', ' ', 'Post D')
  Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                         ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)

  Features = Features[Features$p_val_adj < 0.05,]
  folder_name = paste0('Feature Plots/',ident1, '_', ident2 )
  
  Features = Features[order(Features$avg_logFC),]
  gene_list = rownames(Features)
  gene_list = c(gene_list[1:20], tail(gene_list, n=20))
  FeaturePlot_GeneList(data_run_label,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
  
  
  ##############
  
  ident1 = paste0('CD14+ Mono', ' ', 'Pre ND')
  ident2 = paste0('CD14+ Mono', ' ', 'Post ND')
  Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                         ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
  
  Features = Features[order(Features$avg_logFC),]
  Features = Features[Features$p_val_adj < 0.05,]
  gene_list = rownames(Features)
  gene_list = gene_list[1:50]
  folder_name = 'Dexa CD14+ ND Gene Feature Plots 12-17-2019'
  FeaturePlot_GeneList(data_run,gene_list,folder_base_output,folder_name,sample_type, FeaturePlotFix = TRUE)
  
#################################################
  
  
  path = paste0(folder_base_output,'DE/','T','/')
  dir.create( path, recursive = TRUE)
  path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
  print(path)
  write.csv(Features, file = path,row.names=TRUE)
  
  
  browser()
  #PlotKnownMarkers(data_run,data_run, paste0(folder_base_output,'Cell Type/'), cell_features = NA,split = FALSE)
  # Label:
  
  #IntegrateAll_ClusterUmap(data_run,sample_type,folder_base_output,PCA_dim,resolution_val,label = FALSE)
  #IntegrateAll_ClusterUmap(data_run_label,sample_type,folder_base_output,PCA_dim,resolution_val,label = TRUE)
  
  
  #########################
  data_run_label = label_cells(data_run,cluster_ids_new)
  filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_splitAll', '','.png'))  
  png(file=pathName,width=2600, height=500,res = 100)
  print(DimPlot(data_run_label, label=T, repel=F, reduction = "umap", split.by = "split_var"))
  dev.off()
  ###############################
  
  #browser()
  

  
  #data_run = label_cells(data_run,cluster_IDs)
  
  
  # 
  # data = SubsetData(object = data_run, cells = data_run$split_var == 'Post D' )
  # plot = DimPlot(data, label=T, repel=F, reduction = "umap",pt.size = 1) + 
  #   ggtitle('Post D' ) + 
  #   theme(plot.title = element_text(hjust = 0.5))
  # pathName <- paste0(filepath_cluster,paste0('ClusterUmap',resolution_val,'_splitAll_PostD', '','.png'))  
  # png(file=pathName,width=1000, height=500,res = 100)
  # print(plot)
  # dev.off()
  # 

  
  #browser()
  
  ####################

  print('Get Markers')
  data_run_label = label_cells(data_run,gsub("dCD14+ Mono", "CD14+ Mono", cluster_IDs))
  Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run@meta.data$split_var)
  
  cell_list = c('T Cell','Mono1','Mono2','Mono3', 'NK')
  category_list = c('Pre ND','Pre D', 'Post ND', 'Post D')
  
  for (cat_i in category_list){
    for (cat_j in category_list){
      for (cell in cell_list){
        if (cat_i != cat_j){
          ident1 = paste0(cell, ' ', cat_i)
          ident2 = paste0(cell, ' ', cat_j)
          Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                                              ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
          
          path = paste0(folder_base_output,'DE/',cell,'/')
          dir.create( path, recursive = TRUE)
          path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
          print(path)
          write.csv(Features, file = path,row.names=TRUE)
        }
      }
    }
    
  }

  

  ######################################33
  data = SubsetData(object = data_run_label, cells = (data_run$dexa == "Yes" && data_run$orig.ident == "data_post"))
  data_run_label = label_cells(data_run,cluster_IDs)
  gene_list = c('GZMA','GZMB','GZMH','GZMK','FCER1G' ,'CXCR4','KLRF1', 
                'KLRB1', 'KLRD1', 'KLRC1', 'KLRG1', 'IL2RB', 'IL2RG', 'TSC22D3', 'NR4A2', 
                'EVL', 'IFITM2', 'TNFAIP3','TGFB1','NFKBIA', 'GNLY', 'NKG7','FCGR3A', 'CCND3'
                , 'LTB','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4',
                'GIMAP7','FTH1','THBS1','CCR2','FCGR1A','HLA-DRB5','HLA-DQB1')
  
  gene_list = c('KLRK1','ITGAX','CX3CR1','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4',
                'GIMAP7','FTH1','THBS1','CCR2','FCGR1A','HLA-DRB5','HLA-DQB1')
  folder_name = 'Dexa Gene Feature Plots'

  FeaturePlot_GeneList(data_run_label,folder_base_output,folder_name)
    
    
  
  
  
  # 
  # gene_list = c('GZMA','GZMB','GZMH','GZMK','FCER1G' ,'CXCR4','KLRF1', 
  #               'KLRB1', 'KLRD1', 'KLRC1', 'KLRG1', 'IL2RB', 'IL2RG', 'TSC22D3', 'NR4A2', 
  #               'EVL', 'IFITM2', 'TNFAIP3','TGFB1','NFKBIA', 'GNLY', 'NKG7','FCGR3A', 'CCND3',
  #               'LTB','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4','GIMAP7',
  #               'KLRK1','ITGAX','CX3CR1','RGS1', 'CXCR4','TSC22D3','JUN', 'JUNB','JUND', 'FOS', 'GIMAP4',
  #               'GIMAP7','FTH1','THBS1','CCR2','FCGR1A','HLA-DRB5','HLA-DQB1','IL2RB')
  # 
###############################################################################



  ###############################################################################
  # HeatMap Compare Cell Types
  #gene_list = c('CXCR4','NFKBIA', 'TNFAIP3', 'NR4A2', 'TGFB1', 'RGS1','TSC22D3') 
  category_list = c('Pre ND','Pre D', 'Post ND', 'Post D')

  data_run_label = label_cells(data_run,cluster_IDs)
  
    
  
  
  Features = FindMarkers(data_run_label, ident.1 = ident1, ident.2 = ident2
                         ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
  Features = Features[Features$p_val_adj < 0.05,]
  Features = Features[order(Features$avg_logFC),]
  gene_list = rownames(Features)
  gene_list = c(gene_list[1:20], tail(gene_list, n=20))
  
  for (category in category_list){
    folder_name = paste0('DoHeatmap/',ident1, '_', ident2)
    folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
    dir.create( folder_featureplot, recursive = TRUE)
    data = SubsetData(object = data_run_label, cells = data_run_label$split_var == category )
    print(category)
    print(data)
    data = subset(data, idents =celltype_list)
    
    plot  = DoHeatmap(object = data, features = gene_list,
      group.by = "ident") +
      ggtitle(category ) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=24))
    
    pathName <- paste0(folder_featureplot,paste0(category,'.png'))  
    png(file=pathName,width=600, height=600)
    print(plot)
    dev.off()
    
    pathName <- paste0(folder_featureplot,paste0(ident1, '_', ident2,' Markers','.csv'))  
    
    write.csv(Features, file = pathName,row.names = TRUE)
    
    
  }  
  
  folder_name = paste0('DoHeatmap/',ident1, '_', ident2)
  folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
  data_run_label = label_cells(data_run, cluster_IDs)
  Idents(data_run_label) = paste0(Idents(data_run_label),' ', data_run@meta.data$split_var)
  celltype_list = c(ident1,ident2)
  data = subset(data_run_label, idents = celltype_list)
  
  
  plot  = DoHeatmap(object = data, features = gene_list,
                    group.by = "ident") +
    ggtitle('' ) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=24))
  
  pathName <- paste0(folder_featureplot,paste0(ident1, '_', ident2,'.png'))  
  png(file=pathName,width=600, height=1000)
  print(plot)
  dev.off()
  #################################################   
  # HeatMap
  gene_list = c('HLA-DQB1', 'HLA-DRB5', 'FCGR1A', 'CCR2', 'CX3CR1') 
  for (category in category_list){
    folder_name = paste0('Dexa Gene DoHeatmap ',category)
    folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
    dir.create( folder_featureplot, recursive = TRUE)
    data = SubsetData(object = data_run, cells = data_run$split_var == category )
    data = subset(data, idents = c('CD14+ Mono'))
    
    plot  = DoHeatmap(object = data, features = gene_list,
                      group.by = "ident", disp.min = 0, disp.max = 2.5) +
      ggtitle(category ) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=24))
    pathName <- paste0(folder_featureplot,paste0(category,' All Mono','.png'))  
    png(file=pathName,width=600, height=600)
    print(plot)
    dev.off()
  }    
  
  
  ###############################
  ## Save Data in matrix
  ###############################
  
  #SaveAsMatrix(data_run,folder_base_output)
  
  ###############################
  ## Plot samples seperately
  ###############################

  sample_list = unique(data_run$sample_name)
  data_run_label = label_cells(data_run,cluster_IDs)
  
  
  for (sample in sample_list){
    folder_output = paste0(folder_base_output,'Samples Seperate/',sample,'/')
    print(folder_output)
    pathName <- paste0(folder_output,'Cluster')
    dir.create( pathName, recursive = TRUE)
  
    pathName <- paste0(folder_output,'DE')
    dir.create( pathName, recursive = TRUE)
    
    pathName <- paste0(folder_output,'Stats')
    dir.create( pathName, recursive = TRUE)
    
    pathName <- paste0(folder_output,'PCA')
    dir.create( pathName, recursive = TRUE)
    
    data = SubsetData(object = data_run_label, cells = data_run_label$sample_name == sample )
    
    
    plotAll(data,folder_output,sample_name,sampleParam,label_TF = TRUE,integrate_TF = FALSE,  DE_perm_TF = FALSE,clusterTF = FALSE)
  
    PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample]
    resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample]
    
    #data_run_label = getCluster (data,resolution_val, PCA_dim)
    #data_run_label = data
    
    #cluster_IDs_sample <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample]
    #browser()
    #data_run_label = label_cells(data_run_label,cluster_IDs)
    sample_type = ''
    folder_name  = 'Analysis/Feature Plots'
    #FeaturePlot_GeneList(gene_list,folder_output,folder_name,sample_type, FeaturePlotFix = TRUE)
    data = subset(data, idents = celltype)
   
    
  }
  
  data_run_label = label_cells(data_run,cluster_IDs)
  
  celltype = 'T Cell'
  ident1 = paste0(celltype, ' ', 'Post ND')
  ident2 = paste0(celltype, ' ', 'Post D')
  folder_name = paste0('DoHeatmap/',ident1, '_', ident2)
  folder_featureplot = paste0(folder_base_output,'Analysis/', folder_name,'/')
  filename = paste0(folder_featureplot, ident1, '_', ident2,' Markers.csv')
  Features = read.csv(filename)
  Features = Features[Features$p_val_adj < 0.05,]
  Features = Features[order(Features$avg_logFC),]
  gene_list = as.character(Features$X)
  gene_list = c(gene_list[1:10], tail(gene_list, n=50))
  
  
  category_list = unique(data_run_label$split_var)
  folder_output = paste0(folder_base_output,'Samples Seperate/')
  for (category in category_list){
    
    data = SubsetData(object = data_run_label, cells = data_run$split_var == category )
    data = subset(data, idents = celltype)
    plot  = DoHeatmap(object = data, features = gene_list,
      group.by = "sample_name") +
      ggtitle(celltype ) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=24))
    
    folder_name = paste0('DoHeatmap/',ident1, '_', ident2)
    folder_heatMap = paste0(folder_output,'Analysis/', folder_name,'/')
    dir.create( folder_heatMap, recursive = TRUE)
    pathName <- paste0(folder_heatMap,paste0(ident1, '_', ident2,'_',category,'.png'))  
    png(file=pathName,width=1500, height=1000)
    print(plot)
    dev.off()
  }
  

  
}

