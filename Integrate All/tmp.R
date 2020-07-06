# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/MainIntegrateAll_RunPipline.R')
rm(list = ls())
gc()

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/PlotAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/PipelineIntegrateAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/PlotFunctionIntegrateAll.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integrate All/IntegrateAll_ClusterUmap.R')


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

regress_var = sort(c("Patient","dexa","kit","Response"))
str = paste0('Reg',paste(regress_var, collapse = '_'))
clean = 'Clean_'
#clean = ''
rpca = '_rpca_TEST'
if (integrate_merge == 'Merge'){
  rpca = ''
}


sample_type_list = c('PreNBM','PostNBM','NBM')

sample_type_list = c('PrePostNBM')
#sample_type_list = c('PrePostNBM_filterF')

print(paste0('integrate_merge:', integrate_merge))

if (integrate_merge == 'Integrate' || integrate_merge == 'Merge'){
  filename_sampleParam <- paste0('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_','Combine','_parameters.xlsx')
  sampleParam <- read_excel(filename_sampleParam)
}else{
  filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
  sampleParam <- read_excel(filename_sampleParam)
}

print(filename_sampleParam)
i = 1


sample_type = sample_type_list[i]

folder_base_output = paste0('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/',
                            integrate_merge ,' All/',sample_type,'/')

sample_name = paste0(integrate_merge, '_',sample_type)
PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]

path = paste0(folder_base_output,'data','_',integrate_merge,rpca,'_',sample_type,'.Robj')
print(path)
dataPreNorm = loadRData(path)

dataPreNorm$Patient = dataPreNorm$'Patient Number'
dataPreNorm$kit = dataPreNorm$'10X kit'

data = NormalizeData(dataPreNorm, normalization.method = "LogNormalize", scale.factor = 10000)

#########################
regress_var = sort(c("")) 
data_RegNone = ScaleData(data)
data_RegNone_NoNorm = ScaleData(dataPreNorm)

###############################################3
regress_var = sort(c("Patient"))

if (regress_var!= ''){    
  
  ## Convsert variables to numeric for regression
  for (var in regress_var){
    data@meta.data[,var] = as.character(factor(data@meta.data[,var]))
    data_var  = unique(data@meta.data[,var])
    cnt = 1
    data@meta.data[,paste0(var,'_numeric')] =  NA
    for (level in data_var){
      data@meta.data[,paste0(var,'_numeric')][data@meta.data[,var] == level] = as.numeric(cnt)
      cnt = cnt + 1
    }
    
  }
  data_RegPatient = ScaleData(data, vars.to.regress = paste0(regress_var))
  data_RegPatient_Numeric = ScaleData(data, vars.to.regress = paste0(regress_var,'_numeric'))
}else{
  #data = ScaleData(data)
}
############################################
regress_var = sort(c("dexa"))

if (regress_var!= ''){    
  
  ## Convert variables to numeric for regression
  for (var in regress_var){
    data@meta.data[,var] = as.character(factor(data@meta.data[,var]))
    data_var  = unique(data@meta.data[,var])
    cnt = 1
    data@meta.data[,paste0(var,'_numeric')] =  NA
    for (level in data_var){
      data@meta.data[,paste0(var,'_numeric')][data@meta.data[,var] == level] = as.numeric(cnt)
      cnt = cnt + 1
    }
    
  }
  #data_Regdexa = ScaleData(data, vars.to.regress = paste0(regress_var))
  data_Regdexa2 = ScaleData(data, vars.to.regress = paste0(regress_var))
  #data_Regdexa_Numeric = ScaleData(data, vars.to.regress = paste0(regress_var,'_numeric'))
}else{
  #data = ScaleData(data)
}

###############################

regress_var = sort(c("Patient","dexa"))

if (regress_var!= ''){    
  
  ## Convert variables to numeric for regression
  for (var in regress_var){
    data@meta.data[,var] = as.character(factor(data@meta.data[,var]))
    data_var  = unique(data@meta.data[,var])
    cnt = 1
    data@meta.data[,paste0(var,'_numeric')] =  NA
    for (level in data_var){
      data@meta.data[,paste0(var,'_numeric')][data@meta.data[,var] == level] = as.numeric(cnt)
      cnt = cnt + 1
    }
    
  }
  data_RegdexaPatient = ScaleData(data, vars.to.regress = paste0(regress_var))
  data_RegdexaPatient_Numeric = ScaleData(data, vars.to.regress = paste0(regress_var,'_numeric'))
  
}else{
  #data = ScaleData(data)
}

#################################################
regress_var = sort(c("dexa"))

if (regress_var!= ''){    
  
  ## Convert variables to numeric for regression
  for (var in regress_var){
    data@meta.data[,var] = as.character(factor(data@meta.data[,var]))
    data_var  = unique(data@meta.data[,var])
    cnt = 1
    data@meta.data[,paste0(var,'_numeric')] =  NA
    for (level in data_var){
      data@meta.data[,paste0(var,'_numeric')][data@meta.data[,var] == level] = as.numeric(cnt)
      cnt = cnt + 1
    }
    
  }
  data_Regdexa_noNorm = ScaleData(dataPreNorm, vars.to.regress = paste0(regress_var))
  data_Regdexa_Numeric_noNorm = ScaleData(dataPreNorm, vars.to.regress = paste0(regress_var,'_numeric'))
}else{
  #data = ScaleData(data)
}

data_RegNone
data_RegPatient
data_RegPatient_Numeric
data_Regdexa 
data_Regdexa_Numeric
data_RegdexaPatient
data_RegdexaPatient_Numeric

case1 = data_RegNone@assays[["integrated"]]@scale.data == data_RegPatient@assays[["integrated"]]@scale.data
case2 = data_RegNone@assays[["integrated"]]@scale.data == data_RegPatient_Numeric@assays[["integrated"]]@scale.data
case3 = data_RegNone@assays[["integrated"]]@scale.data == data_Regdexa@assays[["integrated"]]@scale.data
case4 = data_RegNone@assays[["integrated"]]@scale.data == data_Regdexa_Numeric@assays[["integrated"]]@scale.data

case5 = data_Regdexa@assays[["integrated"]]@scale.data == data_RegdexaPatient@assays[["integrated"]]@scale.data

case6 = data_Regdexa@assays[["integrated"]]@scale.data == data_Regdexa_Numeric@assays[["integrated"]]@scale.data

case7 = data_RegNone_NoNorm@assays[["integrated"]]@scale.data == data_Regdexa_Numeric_noNorm@assays[["integrated"]]@scale.data

tmp1 = data_Regdexa@assays[["integrated"]]@scale.data
tmp2 = data_Regdexa2@assays[["integrated"]]@scale.data


tmp0 = data_RegNone@assays[["integrated"]]@scale.data
tmp1 = data_Regdexa_cluster@assays[["integrated"]]@scale.data
tmp2 = data_RegdexaPatient_cluster@assays[["integrated"]]@scale.data

PCA_dim = 20
data_Regdexa_cluster = RunPCA(data_Regdexa, npcs = PCA_dim)
data_Regdexa_cluster = FindNeighbors(data_Regdexa_cluster, dims = 1:PCA_dim)
data_Regdexa_cluster = FindClusters(data_Regdexa_cluster, resolution = resolution_val) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.)
data_Regdexa_cluster = RunUMAP(data_Regdexa_cluster, dims = 1:PCA_dim)

PCA_dim = 20
data_RegdexaPatient_cluster = RunPCA(data_RegdexaPatient, npcs = PCA_dim)
data_RegdexaPatient_cluster = FindNeighbors(data_RegdexaPatient_cluster, dims = 1:PCA_dim)
data_RegdexaPatient_cluster = FindClusters(data_RegdexaPatient_cluster, resolution = resolution_val) # Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.)
data_RegdexaPatient_cluster = RunUMAP(data_RegdexaPatient_cluster, dims = 1:PCA_dim)


print(DimPlot(data_Regdexa_cluster,pt.size = 0.5, reduction = "umap",label = TRUE))
print(DimPlot(data_RegdexaPatient_cluster,pt.size = 0.5, reduction = "umap",label = TRUE))



path1 ='C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate All/PrePostNBM/RegPatient/ConvCatF/data_run_Integrate_rpca_PCAdim20_PrePostNBM.Robj'
path2 ='C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate All/PrePostNBM/Regdexa_Patient/ConvCatF/data_run_Integrate_rpca_PCAdim20_PrePostNBM.Robj'

data1 = loadRData(path1)
data2 = loadRData(path2)


tmp1 = data1@assays[["integrated"]]@scale.data
tmp2 = data2@assays[["integrated"]]@scale.data
