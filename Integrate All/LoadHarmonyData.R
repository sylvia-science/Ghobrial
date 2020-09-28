
LoadHarmonyData = function(metaData, sampleParam_combine){ 
  
  library(Seurat)
  
  
  source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
  source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
  source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
  source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
  source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')
  source('/home/sujwary/Desktop/scRNA/Code/Visualization/PlotCellPhoneDB.R')
  

  metaData = metaData[metaData$Run== 1,]
  
  filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
  sampleParam_combine <- read_excel(filename_sampleParam)
  
  downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
  downsample  = NA #downsample$x
  
  
  sample_type = 'Harmony_AllSamples_Sample_Kit'
  PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
  resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
  cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 
  
  
  patient_list = c(12, 16, 20)
  
  i = 1
  
  filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
  Samples_runs = read_excel(filename_testIntRun)
  
  folder_name = 'AllSamples'
  sample_list = metaData$Sample
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                  '/Batch_Sample_Kit/','/')
  dir.create(folder,recursive = T)
  

  path = paste0(folder,'data_run','.Robj')
  data_merge_run = loadRData(path)
  tmp = data_merge_run@meta.data[paste0('RNA_snn_res.', resolution_val)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_merge_run) = tmp
  data_merge_run_label = label_cells(data_merge_run,cluster_IDs)
  
  celltype = 'Mono_DC'
  cell_list = c('CD14+ Mono','CD16+ Mono','DC')
  resolution_val_subset = 1.6
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  
  Ident_main = colnames(data_merge_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_merge_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_merge_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_merge_run_label) = newIdents2
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
  
  #(0,48)
  #data_merge_run = FindClusters(data_merge_run,resolution = resolution_val)
  return_list <- list("data" = data_merge_run_label, "filepath_cluster" = filepath_cluster)
  
  return(return_list)
}