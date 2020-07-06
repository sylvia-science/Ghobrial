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

print(paste0('integrate_merge:', integrate_merge))

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
filename_metaData <- paste0('/home/sujwary/Desktop/scRNA/Data/EloRD Meta','','.xlsx')
sampleParam <- read_excel(filename_sampleParam)
metaData = read_excel(filename_metaData)


#############################################

# Res = 2.0
celltype = ''

regress_var = sort(c("kit"))
str = paste0('Reg',paste(regress_var, collapse = '_'))

# Set to Clean if cell types have been removed
clean = '/'
saveClean = FALSE

# Set to rpca is rpca was used to integrate 
rpca = ''
print(filename_sampleParam)

sample_type = 'PrePostEOTNBM'
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

filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )

data_run_res2 = loadRData(path) # The main data we're working with
data_run_res2 = FindClusters(data_run_res2, resolution = resolution_val)

data_run_res2 = label_cells(data_run_res2, cluster_IDs)


data_run_res2@meta.data$split_var =  data_run_res2@meta.data$orig.ident
data_run_res2@meta.data$split_var = gsub("data_baseline", "Baseline", data_run_res2@meta.data$split_var)
data_run_res2@meta.data$split_var = gsub("data_EOT", "EOT", data_run_res2@meta.data$split_var)
data_run_res2@meta.data$split_var = gsub("data_C9D1", "C9D1", data_run_res2@meta.data$split_var)
data_run_res2@meta.data$split_var = gsub("data_NBM", "NBM", data_run_res2@meta.data$split_var)
data_run_res2@meta.data$Response[is.na(data_run_res2@meta.data$Response)] = 'NBM'
data_run_res2$Best_Overall_Response = ''
data_run_res2$Current_Status = ''
for (sample in unique(data_run_res2$sample_name)){
  data_run_res2$Response[data_run_res2$sample_name == sample] = metaData[metaData$Sample == sample,]$General_Response
  data_run_res2$Best_Overall_Response[data_run_res2$sample_name == sample] = metaData[metaData$Sample == sample,]$Best_Overall_Response
  data_run_res2$Current_Status[data_run_res2$sample_name == sample] = metaData[metaData$Sample == sample,]$Current_Status
}


# Label data

data_run_res3_2 = FindClusters(data_run_res2, resolution = 3.2)
#####################################################
celltype = 'TCell_NK'

regress_var = sort(c("kit"))
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



print(filename_sampleParam)

sample_type = 'PrePostEOTNBM'
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

filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )

data_run_subset_TCell = loadRData(path) # The main data we're working with
data_run_subset_TCell = FindClusters(data_run_subset_TCell, resolution = resolution_val)
# the split_var metadata will contain relevent info for our project
data_run_subset_TCell@meta.data$split_var =  data_run_subset_TCell@meta.data$orig.ident
data_run_subset_TCell@meta.data$Response[is.na(data_run_subset_TCell@meta.data$Response)] = 'NBM'
data_run_subset_TCell$Best_Overall_Response = ''
data_run_subset_TCell$Current_Status = ''
for (sample in unique(data_run_subset_TCell$sample_name)){
  data_run_subset_TCell$Response[data_run_subset_TCell$sample_name == sample] = metaData[metaData$Sample == sample,]$General_Response
  data_run_subset_TCell$Best_Overall_Response[data_run_subset_TCell$sample_name == sample] = metaData[metaData$Sample == sample,]$Best_Overall_Response
  data_run_subset_TCell$Current_Status[data_run_subset_TCell$sample_name == sample] = metaData[metaData$Sample == sample,]$Current_Status
}
# Label data
data_run_TCell = label_cells(data_run_subset_TCell,cluster_IDs)
cell_list =  c('NK','mNK','Tgd','CD8+ T Cell','Treg','TEM','TEMRA','TCM','TSCM','U1', 'U2','U3')
data_run_TCell = data_run_TCell[,Idents(data_run_TCell) %in% cell_list]
Idents(data_run_TCell) = factor(Idents(data_run_TCell) , levels = cell_list)

## Load NK subset
############################################################################
sample_name = 'Integrate_PrePostEOTNBM_Regkit_NK'
PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name] 

path = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM/Regkit/SubsetIntegrate/NK//data_run_Integrate_PCAdim30_PrePostEOTNBM.Robj'
data_run_NK = loadRData(path) # The main data we're working with
data_run_NK = FindClusters(data_run_NK, resolution = resolution_val)
data_run_NK = label_cells(data_run_NK,cluster_IDs)
##################################################################
# Get original cell names
###########################################################

#Use subsetted run to get mNK

#Use 2.3 run for TEMRA. TEMRAs are 10 GNLY FCGR3A


cell_list =  c('mNK')
data_run_TCell = data_run_TCell[,Idents(data_run_TCell) %in% cell_list]
Idents(data_run_TCell) = factor(Idents(data_run_TCell) , levels = cell_list)


data_run_res2$cell_names = data_run_res2@assays[["integrated"]]@data@Dimnames[[2]]
data_run_res2$cell_names = gsub("_.*", "", data_run_res2$cell_names)
data_run_res2$cell_sample = paste0(data_run_res2$cell_names,'_',data_run_res2$sample_name)

data_run_res3_2$cell_names = data_run_res3_2@assays[["integrated"]]@data@Dimnames[[2]]
data_run_res3_2$cell_names = gsub("_.*", "", data_run_res3_2$cell_names)
data_run_res3_2$cell_sample = paste0(data_run_res3_2$cell_names,'_',data_run_res3_2$sample_name)

data_run_TCell$cell_names = data_run_TCell@assays[["integrated"]]@data@Dimnames[[2]]
data_run_TCell$cell_names = gsub("_.*", "", data_run_TCell$cell_names)
data_run_TCell$cell_sample = paste0(data_run_TCell$cell_names,'_',data_run_TCell$sample_name)


data_run_NK$cell_names = data_run_NK@assays[["integrated"]]@data@Dimnames[[2]]
data_run_NK$cell_names = gsub("_.*", "", data_run_NK$cell_names)
data_run_NK$cell_sample = paste0(data_run_NK$cell_names,'_',data_run_NK$sample_name)


data_run_TEMRA = data_run_res3_2[,Idents(data_run_res3_2) %in% '9']
data_run_Neutro = data_run_res3_2[,Idents(data_run_res3_2) %in% '40']




idx_maininMNK = which(data_run_res2$cell_sample %in%  data_run_TCell$cell_sample )
idx_maininNeutro = which(data_run_res2$cell_sample %in%  data_run_Neutro$cell_sample )
idx_maininTEMRA = which(data_run_res2$cell_sample %in%  data_run_TEMRA$cell_sample )

#Use original run to get rest

data_run_res2$newlabels = as.character(Idents(data_run_res2))

data_run_res2$newlabels[idx_maininMNK] = 'NK'
data_run_res2$newlabels[idx_maininNeutro] = 'Neutrophil'
data_run_res2$newlabels[idx_maininTEMRA] = 'TEMRA'
Idents(data_run_res2) = data_run_res2$newlabels
data_run_res2 = RenameIdents(object = data_run_res2, 'mNK' = "NK CD56dim", 'NK' = "NK CD56br")


celltype = ''

regress_var = sort(c("kit"))
str = paste0('Reg',paste(regress_var, collapse = '_'))

clean = '/'
rpca = ''

sample_type = 'PrePostEOTNBM'
sample_name = paste0(integrate_merge,'_',sample_type,'_',str)
folder_base_output = '/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM/RegkitClean/'
filepath_cluster = paste0( folder_base_output, 'Cluster/', 'PCA',30,'/res',2,'/' )

path = paste0(folder_base_output,
              '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
              '.Robj')
#save(data_run_res2,file= path)

data_run_res2 = loadRData(path) 
##############################################################


gene_list = c("PTPRC", "HNRNPLL", "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","SELL","CCR7","CD27",
              "CD28", "CD69", "IL2RA", "IL2RB", "IL2RG", "CD38", "FAS", "IL7R", "KLRG1", "ZEB1",
              "ZEB2", "PRF1", "GNLY", "NKG7","FCGR3A", "ITGAL", "CX3CR1", "B3GAT1", "BCL2", "MCL1", "LEF1",
              "TRDC", "TRGC1", "TRGC2", "TRAV10", "KLRB1","LAMP1","TRAC",'TCF7',"FOXP3",'IL10', 'TGFB1', 
              'CTLA4', 'TNFRSF18', 'LTB','NOSIP','NTDP1','GZMA','GZMB','GZMK','GZMH','GZMM',
              'CCL3','IFNG','KLRD1','ITGAM','HAVCR2','LAG3','PDCD1','TIGIT','TBX21','ADGRG1',
              'NCAM1','SRGN',"HLA-DRA" , "HLA-DRB5", "HLA-DRB1") #ENTB1

gene_list = c('PTPRC', 'CD14', 'FCN1', 'LYZ', 'HLA-DRA', 'HLA-DRB1', 'S100A4', 'CD68', 'CD74', 'FTH1', 'FTL',
              'CD302', 'FCGR3A' , 'CCR5', 'ITGAM', 'ITGAX', 'CD68', 'TFRC') #macro
gene_list = c('MPO', 'PRTN3', 'ELANE', 'CD34') # GMPC
gene_list = c('NCAM1','CD27', 'SELL', 'TCF7') # NK CD56
gene_list = c('MPO', 'PRTN3', 'ELANE') # Neutro

gene_list = c('CD68','FUT4','ITGAM','FCGR3A','ITGB2','FCGR2A',
              'CD44','CD55', 'SRGN', 'MNDA', 'CLINT1', 'S100A12', 
              'CD117', 'CD34', 'SOX4', 'MPO', 'ELANE', 'PRTN3', 'MMP8')
gene_list = c('PRTN3','MMP8', 'CAMP', 'NGAL')
gene_list = c('CEACAM8','CGM6','NCA95','CD66B')
gene_list =c('CD115', 'CD135', 'CD64','LY6G6D')
gene_list = c('CD74', 'FTH1', 'FTL')
gene_list = c('PTPRC', 'CD14', 'FCN1', 'LYZ', 'HLA-DRA',' HLA-DRB1', 'S100A4', 'CD68', 'CD74', 'FTH1', 'FTL','CD302')
gene_list = c('CD34', 'MME', 'CD19', 'CD2', 'CD14', 'FUT4', 'NCAM1', 'GYPA', 'CD38', 'MS4A1', 'IGHM', 'DNTT') #CLP
gene_list = c('FUT4','CD38', 'MS4A1', 'IGHM', 'DNTT')

gene_list = c('SLAMF7', 'FCGR3A')

gene_list = c('CXCR4')
for (j in 1:length(gene_list)){
  gene = gene_list[j]
  print(gene)
  #browser()
  folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/') 
  plot = FeaturePlotFix(data_run_res2, feature = gene, folder =folder,
                        str = '',split = F,markerSize =2, gene_TF = TRUE,title = '',saveTF = FALSE) 
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20)
  )
  
  
  pathName = paste0(folder,gene,'','.png')
  png(filename = pathName,width=2000, height=2000)
  print(plot)
  dev.off()
  remove(plot)
  
}

############################

file_str = '_2markers'
cell_list =  c("NK CD56br","NK CD56dim","CD8+ T Cell",'TEMRA','Tgd')
data_run_label_clean = data_run_res2[,Idents(data_run_res2) %in% cell_list]
Idents(data_run_label_clean) = factor(Idents(data_run_label_clean) , levels = cell_list)

pathName <- paste0(filepath_cluster,
                   paste0('HeatMap/HeatMapTCellMarkers_Clean_cyto',file_str,'.png'))
png(file=pathName,width=2000, height=3000, res = 100)
plot = DoHeatmap(object = data_run_label_clean, features = gene_list,assay = 'RNA', slot = "data",
                 group.by = "ident", label = T,size = 12)
plot = plot + theme(
  axis.text= element_text(color="black", size=25))
print(plot)
dev.off()

cell_list =  c("Erythrocyte","CD14+ Mono","CD16+ Mono",'CD14+CD16+ Mono','B Cell',
               'HSC','pDC','Macrophage','DC','Pre B Cell','Pro B Cell','Neutrophil','Plasma Cell')
data_run_label_clean = data_run_res2[,Idents(data_run_res2) %in% cell_list]
Idents(data_run_label_clean) = factor(Idents(data_run_label_clean) , levels = cell_list)


pathName <- paste0(filepath_cluster,
                   paste0('HeatMap/HeatMapTCellMarkers_Clean_other',file_str,'.png'))
png(file=pathName,width=2000, height=3000, res = 100)
plot = DoHeatmap(object = data_run_label_clean, features = gene_list,assay = 'RNA', slot = "data",
                 group.by = "ident", label = T,size = 12)
plot = plot + theme(
  axis.text= element_text(color="black", size=25))
print(plot)
dev.off()


cell_list =  c("TSCM",'TCM',"TTM",'TEM',"TEMRA","T reg",'T Naive')
data_run_label_clean = data_run_res2[,Idents(data_run_res2) %in% cell_list]
Idents(data_run_label_clean) = factor(Idents(data_run_label_clean) , levels = cell_list)

pathName <- paste0(filepath_cluster,
                   paste0('HeatMap/HeatMapTCellMarkers_Clean_Notcyto',file_str,'.png'))
png(file=pathName,width=2000, height=3000, res = 100)
plot = DoHeatmap(object = data_run_label_clean, features = gene_list,assay = 'RNA', slot = "data",
                 group.by = "ident", label = T,size = 12)
plot = plot + theme(
  axis.text= element_text(color="black", size=25))
print(plot)
dev.off()

#########################################
## Plot everything with clean data
#########################################

cell_list =  c("Erythrocyte","32","33",'36','Plasma Cell','MT Cell','T Cell')
data_run_label_clean = data_run_res2[,!(Idents(data_run_res2) %in% cell_list)]

cell_features = getCellMarkers(folder_base)
groupBy_list = c('dexa','sample_name', 'Response','split_var','kit','Patient')

plotAll(data_run_label_clean, folder = folder_base_output,
        sample_name = sample_name,sampleParam = sampleParam,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF = F, markersTF = T, keepOldLabels = T, groupBy = groupBy_list, str = '_Clean')

print('Done Plotting')



#############################
## Compare Individual samples
#############################

celltype_iterate = as.character(unique(Idents(data_run_res2)))
celltype_iterate = celltype_iterate[!(celltype_iterate %in% c('32','33','36'))]
celltype_iterate = c("mNK","NK","CD8+ T Cell",'TEMRA','Tgd','NK','CD8+ T Cell',
                     "TSCM",'TCM',"TTM",'TEM',"TEMRA","T Reg",'T Naive','Erythrocyte',
                     'CD14+ Mono', 'CD16+ Mono','CD14+CD16+ Mono','HSC','DC','Pre B Cell')
celltype_iterate = c('NK CD56br','NK CD56dim','Pre B Cell','Pro B Cell')
category_list = unique(data_run_res2$Patient)
category_list = c('30','31','34','40','51')
celltype = celltype_iterate[1]
category = category_list[1]
sortby = 'avg_logFC'
gene_num = c(40,40)
var = 'split_var'
#category_list = category_list[2:length(category_list)]

for (category in category_list){
  
  for (celltype in celltype_iterate){
    print(celltype)
    celltype_list = c(celltype)
    print(category)
    
    data_run_input = data_run_res2
    data_run_input = SubsetData(object = data_run_input, cells = data_run_input$Patient == category )
    
    
    Idents(data_run_input) = paste0(Idents(data_run_input),' ', data_run_input@meta.data[,var])
    
    ident1 = paste0(celltype, ' ', 'Baseline')
    ident2 = paste0(celltype, ' ', 'C9D1')
    
    folder_name = paste0('Patient ',category,'/',ident1, '_', ident2)
    folder_heatMap = paste0(folder_base_output,'Analysis/', folder_name,'/')
    #folder_heatMap = paste0(folder_base_output,'', folder_name,'/')
    
    if (sum( Idents(data_run_input) ==ident1 ) > 3 && sum(Idents(data_run_input) == ident2 ) > 3){
      
      folder_heatMap = gsub('CD14+ Mono','CD14M',folder_heatMap, fixed = T)
      folder_heatMap  = gsub('CD8+ T Cell','CD8T',folder_heatMap, fixed = T)
      folder_heatMap  = gsub('CD16+ Mono','CD16M',folder_heatMap, fixed = T)
      
      folder_heatMap = gsub("Baseline","B",folder_heatMap,fixed = T)
      folder_heatMap = gsub("Good Response","GR",folder_heatMap,fixed = T)
      folder_heatMap = gsub("Poor Response","PR",folder_heatMap,fixed = T)
      
      
      
      Features = FindMarkers(data_run_input, ident.1 = ident1, ident.2 = ident2
                             ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)
      
      pathName <- paste0(folder_heatMap,paste0('Markers','.csv')) 
      dir.create( folder_heatMap, recursive = TRUE)
      print(pathName)
      
      
      DoHeatMapHelper(data_run_res2,folder_base_output,folder_heatMap, Features,
                      ident1,ident2,celltype_list = celltype_list,
                      category_list,sortby,split_var = 'split_var',
                      cellPair = FALSE,gene_num,str = paste0('_',category))
      data_run_input = data_run_res2[,(data_run_res2$Patient == category | data_run_res2$Response == 'NBM')]
      data_run_input = SubsetData(object = data_run_input, cells = Idents(data_run_input) == celltype )
      
      Idents(data_run_input) = paste0(Idents(data_run_input), ' ', data_run_input$split_var)
      Idents(data_run_input)  = factor(Idents(data_run_input) , levels = paste0(celltype,' ',c('NBM','Baseline','C9D1','EOT') ))
      
      
      Features = Features[Features$p_val_adj < 0.05,]
      Features = Features[order(Features$avg_logFC),]
      write.csv(Features, file = pathName,row.names = TRUE)
      
      gene_list_pos = rev(rownames(Features[Features$avg_logFC > 0,]))
      gene_list_neg = rownames(Features[Features$avg_logFC < 0,])
      split_var = 'split_var'
      for (i in 0:ceiling(length(gene_list_pos)/10)){
        gene_list_posi = gene_list_pos[(i*10 + 1):(i*10 + 10)]
        filename = paste0(ident1, '_', ident2,' Split_',split_var,'_',category,'','_Pos',(i*10 + 10),'.png')
        StackedVlnPlotHelper(data_run_input,gene_list_posi,folder_heatMap,filename)
      }
      
      for (i in 0:ceiling(length(gene_list_neg)/10)){
        gene_list_negi = gene_list_neg[(i*10 + 1):(i*10 + 10)]
        filename = paste0(ident1, '_', ident2,' Split_',split_var,'_',category,'','_Neg',(i*10 + 10),'.png')
        StackedVlnPlotHelper(data_run_input,gene_list_negi,folder_heatMap,filename)
      }
      gene_list_pos10 = gene_list_pos[1:10]
      gene_list_pos20 = gene_list_pos[11:20]
      gene_list_neg10 = gene_list_neg[1:10]
      gene_list_neg20 = gene_list_neg[11:20]
      
    }
    #stats = clusterStats(data_run_input)
    #write.csv(stats, file = paste0(folder_heatMap,'Stats_',category,'.csv'),row.names = FALSE)
    
  }
  
 
  
}


## Plot umap of each sample
patient = 20
for (patient in unique(data_run_res2$Patient)){
  
  data_run_input = data_run_res2[,data_run_res2$Patient == patient]
  
  pathName <- paste0(filepath_cluster,
                     paste0('Split Patient Time/ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,' Patient',patient,'.png'))
  height = 1000*length(unique(data_run_input$split_var))
  png(file=pathName,width=height, height=1000, res = 100)
  print(DimPlot(data_run_input,pt.size = 0.5, reduction = "umap",label = TRUE, 
                split.by = 'split_var',label.size = 5))
  dev.off()
}


############################
## Make Line plots
############################

cell_list =  c("Erythrocyte","32","33",'36','Plasma Cell','MT Cell','T Cell')
data_run_label_clean = data_run_res2[,!(Idents(data_run_res2) %in% cell_list)]


sample_list = unique(data_run_label_clean$sample_name)
celltype_iterate = as.character(unique(Idents(data_run_label_clean)))
celltype = celltype_iterate[1]

stats_summary_line = data.frame(matrix(ncol = 4, 0))
colnames(stats_summary_line) = c("Sample",'Cluster','Num','Percent')


for (sample in sample_list){
  print(sample)
  data_run_input = data_run_label_clean
  data_run_input = data_run_input[,data_run_input$sample_name == sample]
  cluster_num = clusterStats(data_run_input)
  
  cluster_num$Sample = sample
  stats_summary_line = rbind(stats_summary_line,cluster_num)
  
}

stats_summary_line =merge(stats_summary_line, metaData, by = 'Sample')


stats_summary_line$Sample= gsub("NBM1CD45P", "NBM1", stats_summary_line$Sample)
stats_summary_line$Sample = gsub("NBM6CD138N", "NBM6", stats_summary_line$Sample)
stats_summary_line$Sample = gsub("NBM7CD138N", "NBM7", stats_summary_line$Sample)
stats_summary_line$Sample = gsub("NBM9CD138N", "NBM9", stats_summary_line$Sample)
stats_summary_line$Sample = gsub("NBM10CD138N", "NBM10", stats_summary_line$Sample)
stats_summary_line$Sample = gsub("NBM11CD138N", "NBM11", stats_summary_line$Sample)
stats_summary_line$Sample = gsub("NBM12CD138N", "NBM12", stats_summary_line$Sample)


stats_summary_long = data.frame(matrix(ncol = 4, nrow = nrow(stats_summary_line) ))
colnames(stats_summary_long) = c('Patient','NBM','Cluster','Proportion')
stats_summary_long$Patient = stats_summary_line$`Patient Number`
stats_summary_long$NBM = stats_summary_line$`Treatment`
stats_summary_long$Cluster = stats_summary_line$Cluster
stats_summary_long$Proportion = stats_summary_line$Percent
#stats_summary_long = stats_summary_long[stats_summary_long$NBM == 'baseline'| stats_summary_long$NBM == 'NBM',]

write.csv(stats_summary_long, file = paste0(filepath_cluster,'Stats/stats_summary_long.csv'),row.names = T)


all_cell_list = c('TSCM', 'NK CD56dim','Tgd','TCM','T Naive' ,
                  'NK CD56br', 'TEMRA','T-Reg', 'CD8+ T Cell', 'TEM','CD14+ Mono', 'B Cell' ,       
                  'CD16+ Mono ' ,'HSC' ,'pDC','Macrophage','DC','Pre B Cell','CD14+CD16+ Mono','Neutrophil' ,'Pro B Cell' )
stats_summary_long_lympho = stats_summary_long[stats_summary_long$Cluster %in% c('TSCM',),]
pdf("/home/sujwary/Desktop/scRNA/Output/Integrate All/PrePostEOTNBM/RegkitClean/Analysis/StackedBarplot.pdf")
ggplot(stats_summary_long, aes(fill=Cluster, y=Proportion, x=Patient)) + 
  geom_bar(position="fill", stat="identity")
dev.off()


stats_summary_short = data.frame(matrix(ncol = 2, nrow = length(patient_list) ))

colnames(stats_summary_short) = c('Patient','Time')
rownames(stats_summary_short) = patient_list

for (patient in patient_list){
  rowdata = stats_summary_line[stats_summary_line$`Patient Number` == patient,]
  for(i in 1:nrow(rowdata)){
    print(i)
    patient = (rowdata$`Patient Number`)[i]
    celltype = rowdata$Cluster[i]
    time = rowdata$Treatment[i]
    stats_summary_short[patient,'Patient'] = patient
    colname = paste0(celltype,' ', time)
    stats_summary_short[[colname]] = NA
    
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
    stats_summary_short[patient,'Patient'] = patient
    colname = paste0(celltype,' ', time)
    stats_summary_short[patient,colname] = rowdata$Percent[i]
    
  }
}

stats_summary_short = stats_summary_short[ , order(colnames(stats_summary_short))]
write.csv(stats_summary_short, file = paste0(filepath_cluster,'Stats/stats_summary_short.csv'),row.names = T)

stats_summary_short


patient  = 20

celltype_list = c('TSCM','NK CD56dim','Tgd','TCM','T Naive','NK CD56br',
                  'TEMRA','T-Reg','CD8+ T Cell','TEM','B Cell','Pre B Cell','Pro B Cell')
celltype_list = c('CD14+ Mono','CD16+ Mono',
                  'CD14+CD16+ Mono', 'Macrophage','DC','Neutrophil','pDC')

celltype_list = c('TSCM','NK CD56dim','Tgd','TCM','T Naive','NK CD56br',
                  'TEMRA','T-Reg','CD8+ T Cell','TEM','B Cell','Pre B Cell','Pro B Cell',
                  'CD14+ Mono','CD16+ Mono','CD14+CD16+ Mono', 'Macrophage','DC','Neutrophil')
celltype_list = c('HSC','Pre B Cell','Pro B Cell')
ymax_list = c(40,40,40,40,40,
              40,40,40,40,40,
              40,40,40,40,40,
              40,40,40,40,40,40)

ymax_list = c(5, 20,20,20,20,
              20,20,20,20,20,
              20,20,20,20,20,
              20,20,20,20,20,20)
path = paste0(folder_base_output,'Analysis/','Boxplots/Patient/')
dir.create(path,recursive = T)

for (i in 1:length(patient_list)){
  
  patient = patient_list[i]
  ymax= ymax_list[i]
  stats_summary_input = stats_summary_line[stats_summary_line$`Patient Number` == patient,]
  stats_summary_input = stats_summary_input[stats_summary_input$Cluster %in% celltype_list,]
  pathName =paste0(folder_base_output,'Analysis/','Boxplots/Patient/Patient',patient,'_All','.png') 
  #pathName =paste0(folder_base_output,'Analysis/','Boxplots/Patient/Patient',patient,'_Other','.png') 
  
  #pathName =paste0(folder_base_output,'Analysis/','Boxplots/Patient/Patient',patient,'_Lymphoid','.png') 
  pathName =paste0(folder_base_output,'Analysis/','Boxplots/Patient/Patient',patient,'_Myeloid','.png') 
  png(file=pathName,width=800, height=1000)
  plot = ggplot(stats_summary_input, aes(x = Treatment, y = Percent)) +
    coord_cartesian(ylim = c(0, ymax))+
    xlab("") + ylab(paste0('Proportion'))+
    geom_point(color="black", size=3,fill = NA) +
    geom_line(aes(group=`Cluster`, color=stats_summary_input$Cluster), size = 2,alpha = 0.5) +
    labs(colour = "Cell Type",title =paste0('Patient ',patient)) +
    theme_classic()
  
  
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=18),
    legend.title=element_text(size=18),
    plot.title = element_text(hjust = 0.5, size=24)
    
  )
  
  print(plot)
  dev.off()
  
  
  
  
  
  
}


cell_list =  c("TSCM",'Tgd','TCM',"TTM",'TEM',"TEMRA","T-Reg",'T Naive',
               'CD8+ T Cell' )
data_run_res2_Tcell = data_run_res2[,as.character(Idents(data_run_res2)) %in% cell_list]

umap_coord = data_run_res2_Tcell@reductions[["umap"]]@cell.embeddings
Idents_list = as.character(Idents(data_run_res2_Tcell) )

write.table(umap_coord, file=paste0(folder_base_output, 'Data/umap_coord','.tsv'), sep='\t')
write.table(Idents_list, file=paste0(folder_base_output, 'Data/Idents_list','.tsv'),sep='\t')


# get matrix for nmf
data_matrix = data_run_res2_Tcell@assays[["RNA"]]@counts
data_matrix = data_matrix[rownames(data_matrix) %in% data_run_res2@assays[["integrated"]]@var.features,]
data_matrix = NormalizeData(data_matrix, normalization.method = "LogNormalize", scale.factor = 10000)
write.table(data_matrix, file='/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_TCell.tsv', quote=FALSE, sep='\t')

#write.csv(data_matrix, '/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_TCell.csv')


cell_list = c('NK CD56dim','NK CD56br' )
data_run_res2_NK = data_run_res2[,Idents(data_run_res2) %in% cell_list]

path = paste0(folder_base_output,
              '/data_run_',integrate_merge,rpca,'_PCAdim',PCA_dim,'_',sample_type,
              '.Robj')
save(data_run_res2,file= path)

