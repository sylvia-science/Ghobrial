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
sample_type = 'PrePostEOTNBM_MT15'
PCA_dim = 30
resolution_val = 2

folder_base_output = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM_MT15/Regkit/'
filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',30,'/res',2,'/' )


path = paste0(folder_base_output,
              '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
              '.Robj')

data_run_res2 = loadRData(path) 
data_nick = loadRData('/home/sujwary/Desktop/scRNA/Nick data/data.Tcells.091919.0.Robj')

scores_nick = read.table(file = '/home/sujwary/Desktop/scRNA/Nick data/T-cells.log_data.var.noX.txt_H.txt', sep = '\t', header = TRUE)
cellnames_nick = data_nick@data@Dimnames[[2]]
cellnames_nick = gsub("\\-.*","",cellnames_nick) 
cellnames_nick = paste0(cellnames_nick, '_', data_nick@meta.data[["sample"]])
length(cellnames_nick)
length(unique(cellnames_nick))

cellnames_all = data_run_res2@assays[["integrated"]]@data@Dimnames[[2]] 
cellnames_all = gsub("\\-.*","",cellnames_all) 
cellnames_all = paste0(cellnames_all, '_', data_run_res2@meta.data[["sample_name"]])


Idents_nick = as.character(data_nick@meta.data[["cell.type"]])
length(cellnames_nick)
length(Idents_nick)
idx_cellnames_in_cellnames_nick = which(cellnames_all %in%  cellnames_nick)
idx_cellnames_nick_in_cellnames = which(cellnames_nick %in%  cellnames_all)

cellnames = cellnames_all[idx_cellnames_in_cellnames_nick]
cellnames_nick_in_all = cellnames_nick[idx_cellnames_nick_in_cellnames]

cellname_matrix = data.frame()

cellname_matrix$cellname_original = cellnames
cellname_matrix$cellname_original2 = cellnames

print(length(cellnames))
print(length(cellnames_nick_in_all))
all(sort(cellnames) == sort(cellnames_nick_in_all))

cellnames_all[idx_cellnames_in_cellnames_nick[1]]
cellnames_nick[which(cellnames_nick_in_all %in% cellnames_all)[1]]

i_list = c()
idx_list = c()
for (i in 1:length(cellnames_all )){
  if(cellnames_all[i] %in% cellnames_nick){
     print(i)
     #browser()
    
     idx = which(cellnames_nick == cellnames_all[i] )
     data_run_res2$nicklabels[i] = Idents_nick[idx]
     i_list = c(i_list,i)
     idx_list = c(idx_list,idx)
   }
}


tmp1 = which(cellnames_all %in% cellnames_nick)
tmp2 = which(cellnames_nick %in% cellnames_all )
i_list %in% tmp2
idx_list %in% tmp1


data_run_res2$nicklabels = ''

data_run_res2$nicklabels[which(cellnames_all %in% cellnames_nick)] = Idents_nick[which(cellnames_nick %in% cellnames_all )]
data_run_res2$nick_sylvia_labels = paste0(Idents(data_run_res2), '_', data_run_res2$nicklabels)


data_run_res2_clean = data_run_res2[,data_run_res2$nicklabels != '']
group = 'nicklabels'
pathName = paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=1000, height=1000)
plot = DimPlot(data_run_res2_clean,pt.size = 2, reduction = "umap",label = T,group.by  = group,label.size = 6)

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
)
print(plot)

dev.off()

####################################

gene_list = c("PTPRC", "HNRNPLL", "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","SELL","CCR7","CD27",
              "CD28", "CD69", "IL2RA", "IL2RB", "IL2RG", "CD38", "FAS", "IL7R", "KLRG1", "ZEB1",
              "ZEB2", "PRF1", "GNLY", "NKG7","FCGR3A", "ITGAL", "CX3CR1", "B3GAT1", "BCL2", "MCL1", "LEF1",
              "TRDC", "TRGC1", "TRGC2", "TRAV10", "KLRB1","LAMP1","TRAC",'TCF7',"FOXP3",'IL10', 'TGFB1', 
              'CTLA4', 'TNFRSF18', 'LTB','NOSIP','NTDP1','GZMA','GZMB','GZMK','GZMH','GZMM',
              'CCL3','IFNG','KLRD1','ITGAM','HAVCR2','LAG3','PDCD1','TIGIT','TBX21','ADGRG1',
              'NCAM1','SRGN',"HLA-DRA" , "HLA-DRB5", "HLA-DRB1",'ENTB1') 

pathName = paste0(filepath_cluster,'/HeatMap/','Heatmap_nick', '_TCellGenes','.png')
plot  = DoHeatmap(object = data_run_res2_clean, features = gene_list,assay = 'RNA', slot = "data",
                  group.by = "nicklabels", label = T) +
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

pathName = paste0(filepath_cluster,'/HeatMap/','Heatmap_nick_sylvia', '_TCellGenes','.png')

plot  = DoHeatmap(object = data_run_res2_clean, features = gene_list,assay = 'RNA', slot = "data",
                  group.by = "nick_sylvia_labels", label = T) +
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


png(file=pathName,width=3000, height=1000,res=100)
print(plot)
dev.off()



data_run_res2_clean_small = data_run_res2_clean[,data_run_res2_clean$nicklabels %in% c('CD4 Cytotoxic',"CD8 Cytotoxic","Helper 1","Helper 2","IFN-producing" ,'gamma-delta','NKT' )]


pathName = paste0(filepath_cluster,'/HeatMap/','Heatmap_nick_small', '_TCellGenes','.png')
plot  = DoHeatmap(object = data_run_res2_clean_small, features = gene_list,assay = 'RNA', slot = "data",
                  group.by = "nicklabels", label = T) +
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





data_run_res2_clean_supersmall = data_run_res2_clean[,data_run_res2_clean$nicklabels %in% c("IFN-producing" ,'Helper 2','NKT' )]


pathName = paste0(filepath_cluster,'/HeatMap/','Heatmap_nick_supersmall', '_TCellGenes','.png')
plot  = DoHeatmap(object = data_run_res2_clean_supersmall, features = gene_list,assay = 'RNA', slot = "data",
                  group.by = "nicklabels", label = T) +
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

