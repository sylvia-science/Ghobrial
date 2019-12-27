# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/main_integrate.R')
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)

require(gridExtra)

#library(biomaRt)
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

run = TRUE

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


print('Start')
patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all
#patient_list = c(16, 51, 6, 40) # v3

patient_list_dexaF = c(5,12, 16)
patient_list_dexaT = c(10,20,34,28,21,31,51,6,40)

redo = c( 20, 34, 28, 21, 31, 40)
for(patient in 28){ # Patient numbers 
  
  filter = TRUE
  regress_TF = TRUE
  
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  sample_name_post = pair_list[[paste0('Sample Post ', sample_type)]]
  
  sample_name_list = c(sample_name_pre,sample_name_post)
  
  print('Running')
  print(folder_base_input)
  folder_pre = makeFolders(folder_base_input,sample_name_pre,filter,regress_TF,TRUE)

  folder_post = makeFolders(folder_base_input,sample_name_post,filter,regress_TF,FALSE)
  folder_input = c(folder_pre,folder_post)
  
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter,regress_TF,TRUE)
  folder_post = makeFolders(folder_base_output,sample_name_post,filter,regress_TF,FALSE)
  folder_output = c(folder_pre,folder_post)
  
  
  filter_pre = sampleParam$filter_TF[sampleParam['Sample'] == sample_name_pre]
  filter_post = sampleParam$filter_TF[sampleParam['Sample'] == sample_name_post]
  regress_pre = sampleParam$regress_TF[sampleParam['Sample'] == sample_name_pre]
  regress_post = sampleParam$regress_TF[sampleParam['Sample'] == sample_name_post]
  
  folder_orig_pre = makeFolders(folder_base_input,sample_name_pre,filter = filter_pre,regress_TF = regress_pre,FALSE)
  data_orig_pre = loadRData(paste0(folder_orig_pre,'data.Robj' ))
  folder_orig_post = makeFolders(folder_base_input,sample_name_post,filter = filter_post,regress_TF = regress_post,FALSE)
  data_orig_post = loadRData(paste0(folder_orig_post,'data.Robj'))
  
  cluster_IDs_pre <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name_pre]
  cluster_IDs_post <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name_post]
  #browser()
  data_orig_pre = label_cells(data_orig_pre,cluster_IDs_pre)
  data_orig_post = label_cells(data_orig_post,cluster_IDs_post)
  
  data_integrated_orig = merge(data_orig_pre, data_orig_post)
  
  
  if (run == TRUE){
    
    print('running')
    
    data_integrated = run_pipeline_integrate(folder_input,folder_output,sample_name_list,sampleParam_integrate)
    #browser()  
    
    PCA_dim =  sampleParam_integrate$PCA_dim[sampleParam_integrate['Sample'] == sample_name_pre]
    save(data_integrated,file=paste0(folder_output,'data_','PCA_dim',PCA_dim,patient,'.Robj'))
    browser()
    param_output =  sampleParam_integrate[ sampleParam_integrate$Sample == sample_name_pre, ]
    write.csv(param_output, file = paste0(folder_output[1],'parameters.csv'),row.names=FALSE)
    
    MakeFeaturePlot(data_integrated,data_integrated_orig, cell_features = NA, split = FALSE)
    #data_integrated = label_cells(data_integrated,cluster_IDs)
    
    
    
    print('done!')
    
  }else{
    
    
    print('Starting Plotting')
    print(paste0('Sample: ', sample_name_pre))
    print(paste0('folder: ', folder_output[1]))
    
    data = loadRData(paste0(folder_output[1],'data.Robj'))
    
    #print(FeaturePlot(data, features))
    
    plotAll(data,folder_output[1],sample_name_pre,sampleParam_integrate, 
            label_TF = FALSE, integrate_TF = TRUE, DE_perm_TF = FALSE)
    #plotAll(data,folder_output[1],sample_name_pre,sampleParam_integrate, 
           #label_TF =TRUE, integrate_TF = TRUE, DE_perm_TF = FALSE)

    PCA_dim = sampleParam_integrate$PCA_dim[sampleParam_integrate['Sample'] == sample_name_pre]
    resolution_val<- sampleParam_integrate$resolution_val[sampleParam_integrate['Sample'] == sample_name_pre]
    data = getCluster(data,resolution_val, PCA_dim)
    cluster_IDs = sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Sample'] == sample_name_pre]
    data = label_cells(data,cluster_IDs)

    df_expr_pval = diff_exp_pval_helper(data, folder_pre,sample_name_pre,cell_type_list)
    write.csv(df_expr_pval, file = paste0(folder_output[1],'DE/DE_pval_Patient',patient,'.csv'),row.names=FALSE)
    
    ## Get DE Genes
    data@meta.data$condition = gsub("^.*_","",data@meta.data$orig.ident)
    data$celltype.condition = paste(Idents(data), data$condition, sep = "_")
    data$celltype = Idents(data)
    Idents(data) = "celltype.condition"
    df_expr_gene = diffExpGene(data, folder_pre,sample_name_pre,cell_type_list)
    df_expr_gene = df_expr_gene[df_expr_gene$p_val_adj < 0.05,] # Should save all but only use those <0.05
    write.csv(df_expr_gene, file = paste0(folder_output[1],'DE/DE_Patient',patient,'.csv'),row.names=FALSE)

  }
  
}
