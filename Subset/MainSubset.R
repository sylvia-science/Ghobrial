
# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Subset/MainSubset.R')
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(stringr)
library(data.table)


#library(biomaRt)
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/Functions_integrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Subset/runPipelineSubset.R')


folder_base_output <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_sampleParam_integrate <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters.xlsx'
sampleParam_integrate <- read_excel(filename_sampleParam_integrate)


filename_sampleParam_subset <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters_Tcell.xlsx'
sampleParam_integrate_subset <- read_excel(filename_sampleParam_subset)

filename_sample_Integrate_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Integrate_pairs.xlsx'
sample_Integrate_pairs <- read_excel(filename_sample_Integrate_pairs)


sample_type = 'BM'
#sample_type = 'PB'
cell_type = 'T Cell'

print('Start')
patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all
#patient_list = c(16, 51, 6, 40) # v3
#patient_list = c(34, 28, 21, 31, 16, 51, 6, 40) # all

patient_list_dexaF = c(5,12, 16)
patient_list_dexaT = c(10,20,34,28,21,31,51,6,40)

run = TRUE
if (run){
  for(patient in patient_list){ # Patient numbers 
    
  
    pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
    
    sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  
  
    folder_output = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,TRUE)
    data <- loadRData(paste0(folder_output,'data.Robj'))
    
    folder_output = paste0(folder_output,cell_type,'/')
  
    dir.create( paste0(folder_output,'PCA'), recursive = TRUE)
    dir.create( paste0(folder_output,'Cluster'), recursive = TRUE)
    dir.create( paste0(folder_output,'Cell Type'), recursive = TRUE)
    dir.create( paste0(folder_output,'DE'), recursive = TRUE)
    
    print('Starting Plotting')
    print(paste0('Sample: ', sample_name_pre))
      
    print(paste0('folder: ', folder_output))
      
    cluster_IDs = sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Sample'] == sample_name_pre]
    data = label_cells(data,cluster_IDs)
      
    ##
    ## Cluster T cells seperately and add back to original data
    ##
      
    data_subset_orig = subset(data, idents = cell_type)
    data_subset = runPipelineSubset(data_subset_orig,folder_output,sample_name_pre,sampleParam_integrate_subset,filter = TRUE,regress_TF = TRUE)
    save(data_subset,file=paste0(folder_output,'data.Robj'))
    
    #get_cellType(data_subset,data_subset_orig,folder_output,sample_name_pre) 
      
    ##
    ## Diff Exp
    ##
    levels_list = levels(data_subset)
    
    
    #df_expr_gene = diffExpGene(data_subset, folder_output,sample_name_pre,NA)
    #write.csv(df_expr_gene, file = paste0(folder_output,'diffExprGene_Patient',patient,'.csv'),row.names=FALSE)
  }
  
}else{
  
  print('Starting Plotting')
  print(paste0('Sample: ', sample_name_pre))
  
  print(paste0('folder: ', folder_output))
  
  data <- loadRData(paste0(folder_output,'data.Robj'))
  
  ## MAKE SURE plotAll EQUAL TO plotAll_integrate
  plotAll(data,folder_output,sample_name_pre,sampleParam_integrate_subset,label_TF = FALSE)
  plotAll(data,folder_output,sample_name_pre,sampleParam_integrate_subset, label_TF = TRUE)
  
  cluster_IDs <- sampleParam_integrate[[paste0('Cluster_IDs_',sample_type)]][sampleParam_integrate['Sample'] == sample_name_pre]
  data = label_cells(data,cluster_IDs)
  
  
  df_expr_pval = diff_exp_pval_helper(data, folder_pre,sample_name_pre,cell_type_list)
  write.csv(df_expr_pval, file = paste0(folder_output,'diffExpr_pval_Patient',patient,'.csv'),row.names=FALSE)
  
  data@meta.data$condition = gsub("^.*_","",data@meta.data$orig.ident)
  data$celltype.condition = paste(Idents(data), data$condition, sep = "_")
  data$celltype = Idents(data)
  Idents(data) = "celltype.condition"
  df_expr_gene = diffExpGene(data, folder_pre,sample_name_pre,cell_type_list)
  write.csv(df_expr_gene, file = paste0(folder_output,'diffExprGene_Patient',patient,'.csv'),row.names=FALSE)
  
}