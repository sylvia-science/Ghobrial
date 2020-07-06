
library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
library(sc)
library(scater)
library(dplyr)
library(scran)
library(reshape2)

library(tidyverse)
library(rstatix)
library(ggpubr)


source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')

source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PipelineIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PlotFunctionIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/IntegrateAll_ClusterUmap.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')

source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/TestIntegration/SaveScranSeuratMergeInt_Helper.R')
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)
sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

sample_type = 'Integrate_AllSamples'
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 


filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

folder = 'AllSamples'

sample_list = metaData$Sample

folder_output = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Seurat/Merge/',folder,'/','DetectionRate','/')
dir.create()

data_list = vector(mode = "list",length = length(sample_list))
data_list_noSoup = vector(mode = "list",length = length(sample_list))
#for (i in 1:nrow(sampleParam)){
for (i in 1:length(sample_list)){
  
  
  #sample_name = sampleParam$Sample[i]
  sample_name = sample_list[i]
  
  
  # No Soup
  print(sample_name)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
  
  colSum_list = colSums(data_i_raw ) # Needs to be from Matrix library
  keep = colSum_list >= 100
  data_i_filtered = data_i_raw[,keep]
  data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")
  data_i_filtered = data_i_filtered[, data_i_filtered$percent.mt < percent_mt]
  data_i_filtered$sample = sample_name
  data_list_noSoup[[i]] = data_i_filtered
  next
  # soup
  print(sample_name)
  folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name , '/')
  data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
  data_i$sample = sample_name
  path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
  cellIdents = read.csv(path,sep = ',',row.names = 1)
  cellIdents$x = paste0(cellIdents$x, ' S', i)
  data_i$CellType = cellIdents$x
  data_list[[i]] = data_i
}
data_merge = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = T)
for (i in 3:length(data_list)){
  print(i)
  data_merge = merge(x =  data_merge,y = data_list[[i]],merge.data = T)
}

data_merge_noSoup = merge(x =  data_list_noSoup[[1]],y = data_list_noSoup[[2]], merge.data = T)
for (i in 3:length(data_list)){
  print(i)
  data_merge_noSoup = merge(x =  data_merge_noSoup,y = data_list_noSoup[[i]],merge.data = T)
}


metaData = read_excel(filename_metaData)
data_merge = addMetaData(data_merge, metaData)
data_merge = load_emptyDrops(data_merge)
data_merge = load_Doublets(data_merge)
data_merge = load_CellLabel(data_merge)
data_merge$GeneralCellType = str_match(data_merge$CellType, "(^.+)\\s")[, 2]
data_merge$kit = data_merge$`10X kit`
data_merge$split_var = ''

metaData = read_excel(filename_metaData)
data_merge_noSoup = addMetaData(data_merge_noSoup, metaData)
data_merge_noSoup = load_emptyDrops(data_merge_noSoup)
data_merge_noSoup = load_Doublets(data_merge_noSoup)
data_merge_noSoup = load_CellLabel(data_merge_noSoup)
data_merge_noSoup$GeneralCellType = str_match(data_merge_noSoup$CellType, "(^.+)\\s")[, 2]
data_merge_noSoup$kit = data_merge_noSoup$`10X kit`
data_merge_noSoup$split_var = ''

path = paste0(folder_output,'data_raw','.Robj')
save(data_merge,file= path)
path = paste0(folder_output,'data_raw_noSoup','.Robj')
save(data_merge_noSoup,file= path)

path = paste0(folder_output,'data_raw','.Robj')
data_merge = loadRData(path)

path = paste0(folder_output,'data_raw_noSoup','.Robj')
data_merge_noSoup = loadRData(path)


## DR
data = as.data.frame(data_merge@assays[["RNA"]]@counts)
DR = apply(data,2, function(x) sum(x > 0)/nrow(data))
DR = as.data.frame(DR)
DR$Treatment = data_merge$Treatment
DR$Treatment = factor(DR$Treatment, levels = c('baseline','C9D1','EOT','NBM'))
DR$kit = data_merge$kit
DR$kit = factor(DR$kit, levels = c('v2','v3','v3_1'))


## DR no Soup
data_noSoup = as.data.frame(data_merge_noSoup@assays[["RNA"]]@counts)
DR_noSoup = apply(data,2, function(x) sum(x > 0)/nrow(data_noSoup))
DR_noSoup = as.data.frame(DR_noSoup)
DR_noSoup$DR = DR_noSoup$DR_noSoup
DR_noSoup$Treatment = data_merge_noSoup$Treatment
DR_noSoup$Treatment = factor(DR_noSoup$Treatment, levels = c('baseline','C9D1','EOT','NBM'))
DR_noSoup$kit = data_merge_noSoup$kit
DR_noSoup$kit = factor(DR_noSoup$kit, levels = c('v2','v3','v3_1'))

change = (DR$DR - DR_noSoup$DR)/( DR_noSoup$DR)

plotDR(DR)
plotDR(DR_noSoup)

plotDR = function(DR){
  browser()
  folder_output = '/home/sujwary/Desktop/scRNA/Output/MiscPlots/DR/'
  
  pwc <- DR %>% pairwise_t_test(DR ~ Treatment, p.adjust.method = "bonferroni")
  pwc <- pwc %>% add_xy_position(x = "Treatment")
  
  pathName <- paste0(folder_output,'DR_Treatment' ,'.png')
  png(file=pathName,width=500, height=500,res = 100)
  plot = ggboxplot(DR, x = "Treatment", y = "DR") +
    stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0, step.increase = 0.1) +
    labs(
      caption = get_pwc_label(pwc)
    )
  print(plot)
  dev.off()
  
  pwc_kit <- DR %>%
    pairwise_t_test(DR ~ kit, p.adjust.method = "bonferroni")
  pwc_kit <- pwc_kit %>% add_xy_position(x = "kit")
  pathName <- paste0(folder_output,'DR_kit' ,'.png')
  png(file=pathName,width=500, height=500,res = 100)
  plot = ggboxplot(DR, x = "kit", y = "DR") +
    stat_pvalue_manual(pwc_kit, label = "p.adj", tip.length = 0, step.increase = 0.1) +
    labs(
      caption = get_pwc_label(pwc_kit)
    )
  print(plot)
  dev.off()
  
  patient_list = unique(data_merge$`Patient Number`)
  for (patient in patient_list){
    print(patient)
    DR_patient = DR[data_merge$'Patient Number' == patient,]
    
    pathName <- paste0(folder_output,'DR_patient',patient ,'.png')
    png(file=pathName,width=500, height=500,res = 100)
    
    pwc <- DR_patient %>%
      pairwise_t_test(DR ~ Treatment, p.adjust.method = "bonferroni")
    
    
    DR_patient$Treatment = factor(DR_patient$Treatment, levels = c('baseline','C9D1','EOT','NBM'))
    pwc <- pwc %>% add_xy_position(x = "Treatment")
    plot = ggboxplot(DR_patient, x = "Treatment", y = "DR") +
      stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0, step.increase = 0.1) +
      labs(
        caption = get_pwc_label(pwc)
      )
    print(plot)
    dev.off()
  }
}
