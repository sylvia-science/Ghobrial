library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
library(hdf5r)
require(gridExtra)
#library(biomaRt)
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

clean = '/'
rpca = ''
integrate_merge = 'Integrate'


sample_type = 'PrePostEOTNBM'
PCA_dim = 30
resolution_val = 2

folder_base_output = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM/RegkitClean/'
filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',30,'/res',2,'/' )


path = paste0(folder_base_output,
              '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
              '.Robj')
data_old=loadRData(path)
##################
sample_type = 'PrePostEOTNBM_MT15'
PCA_dim = 30
resolution_val = 2

folder_base_output = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM_MT15/Regkit/'
filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',30,'/res',2,'/' )


path = paste0(folder_base_output,
              '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
              '.Robj')

data_run = loadRData(path)




cellnames_all = data_run@assays[["integrated"]]@data@Dimnames[[2]] 
cellnames_all = gsub("\\-.*","",cellnames_all) 
cellnames_all = paste0(cellnames_all, '_', data_run@meta.data[["sample_name"]])

cellnames_old = data_old@assays[["integrated"]]@data@Dimnames[[2]] 
cellnames_old = gsub("\\-.*","",cellnames_old) 
cellnames_old = paste0(cellnames_old, '_', data_old@meta.data[["sample_name"]])

Idents_old = as.character(Idents(data_old))

data_run$oldlabels = ''

idx_list = c()
for (i in 1:length(cellnames_all )){
  
  if(cellnames_all[i] %in% cellnames_old){
    print(i)
    #browser()
    idx = which(cellnames_all[i] == cellnames_old)
    data_run$oldlabels[i] = Idents_old[idx]
    idx_list = c(idx_list,idx)
  }
}

data_run$oldlabels = 'New Cells'
data_run$oldlabels[which(cellnames_all %in% cellnames_old)] = as.character(Idents_old[which(cellnames_old %in% cellnames_all)])
data_run_clean = data_run[,data_run$oldlabels != 'New Cells']
data_run_clean = data_run

group = 'oldlabels'
pathName = paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_GroupBy',group,'.png'))


plot = DimPlot(data_run_clean,pt.size = 1, reduction = "umap",label = T,group.by  = group,label.size = 6)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
)
png(file=pathName,width=1500, height=1000)
print(plot)

dev.off()

#####
