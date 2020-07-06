library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)

base = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM/Regdexa_kit_Patient/SubsetIntegrate/T Cell/'
path = paste0(base,'anchors_IntegratePrePostEOTNBM_T Cell_features2000.Robj')
anchors = loadRData(path)

data_integrate_final = IntegrateData(anchorset = anchors, k.weight = 100, dims = 1:20)
numfeatures = 2000
path = paste0(base,'data','_features',numfeatures,'.Robj')
save(data_integrate_final,file=path)
