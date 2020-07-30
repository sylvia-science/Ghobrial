# merge the samples

# select highly variable genes

# run PCA and then 

# run BKNN instead of findneighbors. 
library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(harmony)
library(ggplot2)
library(SoupX)
library(sc)
library(scater)
library(dplyr)
library(scran)
library(reshape2)
library(stringr)
library(edgeR)
library(SingleCellExperiment)
library(MAST)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')
source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')

celltype = 'T Cell'
cell_list = c('T Cell','CD8+ T Cell')
resolution_val_subset = 3.5
# 
celltype = 'Mono_DC'
cell_list = c('CD14+ Mono','CD16+ Mono','DC')
resolution_val_subset = 1.6
# # 
# celltype = 'NK'
# cell_list = c('NK')
# resolution_val_subset = 3


filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

#downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  = NA #downsample$x


sample_type = 'Harmony_AllSamples_Sample_Kit'
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 
cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
 
patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

#folder = 'Intra-v3_1'
#folder = 'Inter-version'
folder_name = 'AllSamples'
#folder_name = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample

folder_main = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                '/Batch_Sample_Kit/','/')
dir.create(folder_main,recursive = T)




run = F


if (run){
  path = paste0(folder_main,'data_run','.Robj')
  data_orig = loadRData(path)
  data_orig =label_cells(data_orig,cluster_IDs)
  data_orig = data_orig[,Idents(data_orig) %in% cell_list]
  cell_names_all = colnames(data_orig)
  cell_names_all = sub("_.*", "", cell_names_all)
  data_orig$cell_sample = paste0(cell_names_all,' ',data_orig$sample )
  
  data_list = vector(mode = "list",length = length(sample_list))
  data_list_norm = vector(mode = "list",length = length(sample_list))
  #for (i in 1:nrow(sampleParam)){
  for (i in 1:length(sample_list)){
    #sample_name = sampleParam$Sample[i]
    sample_name = sample_list[i]
    print(sample_name)
    folder_input = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/', sample_name , '/')
    data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
    data_i$sample = sample_name
    path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
    cellIdents = read.csv(path,sep = ',',row.names = 1)
    cellIdents$x = paste0(cellIdents$x, ' S', i)
    data_i$CellType = cellIdents
    data_i$GeneralCellType = str_match(data_i$CellType, "(^.+)\\s")[, 2]
    
    data_i = data_i[,data_i$nFeature_RNA > 200]
    data_i = load_emptyDrops(data_i)
    data_i = data_i[,data_i$is_cell]
    
    cell_names =colnames(data_i)
    cell_names = sub("_.*", "", cell_names)
    data_i$cell_sample = paste0(cell_names ,' ', data_i$sample)
    
    if (any(data_i$cell_sample %in% data_orig$cell_sample)){
      data_i = data_i[, data_i$cell_sample %in% data_orig$cell_sample]
  
      if (!is.na(downsample)){
        downsample = sub("_.*", "", downsample)
        cellnames = colnames(data_i)
        cellnames = sub("_.*", "", cellnames)
        data_i = data_i[,cellnames %in% downsample]
        #browser()
      }
      
      print(ncol(data_i))
      if (ncol(data_i) > 100){
        #if (T){
        
        #data_list[[i]] = data_i
        data_list_norm[[i]] = ScranNorm(data_i)
        
      }
    }
  }
  
  data_list_norm =data_list_norm[lengths(data_list_norm) != 0]
  data_merge = merge(x =  data_list_norm[[1]] ,y = data_list_norm[2:length(data_list_norm)], merge.data = T)
  
  #data_merge = merge(x =  data_list_norm[[1]],y = data_list_norm[[2]], merge.data = T)
  #for (i in 3:length(data_list_norm)){
  #  print(i)
  #  data_merge = merge(x =  data_merge,y = data_list_norm[[i]],merge.data = T)
  #}
  
  data_merge_subset = data_orig
  data_merge_subset = addMetaData(data_merge_subset, metaData)
  data_merge_subset = load_emptyDrops(data_merge_subset)
  data_merge_subset = load_Doublets(data_merge_subset)
  data_merge_subset = load_CellLabel(data_merge_subset)
  data_merge_subset$GeneralCellType = str_match(data_merge_subset$CellType, "(^.+)\\s")[, 2]
  data_merge_subset$kit = data_merge_subset$`10X kit`
  data_merge_subset$split_var = ''
  
  data_run_subset = FindVariableFeatures(data_merge_subset, selection.method = "vst", nfeatures = 2000)
  data_run_subset = ScaleData(data_run_subset)
  data_run_subset = RunPCA(data_run_subset,npcs = 40)
  
  data_run_subset = RunHarmony(data_run_subset,group.by.vars =  c("sample", "10X kit"),
                              dims.use = 1:30)
  
  data_run_subset = RunUMAP(reduction = "harmony",data_run_subset, dims = 1:30)
  data_run_subset = FindNeighbors(data_run_subset, reduction = "harmony", dims = 1:30)
  data_run_subset = FindClusters(data_run_subset,resolution = resolution_val_subset)
  
  
  folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
  dir.create(folder_subcluster,recursive = T)
  path = paste0(folder_subcluster,'data_run','.Robj')
  save(data_run_subset,file= path)
  
}else{
  folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  data_run_subset = FindClusters(data_run_subset,resolution = resolution_val_subset)
  
  
}
gene_name_vector= rownames(data_run_subset)
gene_name_vector[grep("HLA-DR", gene_name_vector)]
gene_list = 'LY6C1, LY6G6D, MMP9, TIMP1, LGALS3, S100A6, S100A11, ISG15, MS4A4C'
gene_list = unlist(strsplit(gene_list, ",")) 
gene_list = trimws(gene_list, which = c("both"), whitespace = "[ \t\r\n]")
gene_list [ !(gene_list %in% gene_name_vector )]

data_run_subset = addMetaData(data_run_subset, metaData)
data_run_subset = load_Doublets(data_run_subset)
data_run_subset = load_CellLabel(data_run_subset)
print(unique(data_run_subset$GeneralCellType))

data_run_subset$FeatureLessThan400 = data_run_subset$nFeature_RNA < 400

data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')

groupBy_list = c('sample','Diagnosis','kit',
                 'Treatment','Batch','LowCount',
                 'Doublet','GeneralCellType','FeatureLessThan400')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
splitBy_list = NA

filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',30,'/res',resolution_val_subset,'/' )


plotAll(data_run_subset, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val_subset)

data_run_subset_label = data_run_subset_label[,Idents(data_run_subset_label) != 14]
plotAll(data_run_subset_label, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = T, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = resolution_val_subset, str = '_label')

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
PlotKnownMarkers(data_run_subset, folder = paste0(filepath_cluster,'Cell Type/T Cell Genes/07-23-2020/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',plotTogether = F)

stats_category = SummariseSampleByCell(data_run_subset_label)
CompareCellNum(data_run_subset_label,stats_category,filepath_cluster,split_var = 'Treatment',metaData)
  
# Heatmap

newIdents = as.character(Idents(data_run_subset))

newIdents[newIdents == '8'] = '8_11'
newIdents[newIdents == '11'] = '8_11'
Idents(data_run_subset) = newIdents


gene_name_vector= rownames(data_run_subset)

gene_list = 'PTPRC, CD3D, CD3G, CD3E, CD4, CD8A, CD8B, IL2RA, FOXP3, CTLA4, 
  TIGIT, LAG3, SELL, CCR7, CD27, CD28, IL7R, NOSIP, LTB, CDK6, 
  LEF1, TCF7, BCL2, MYC, FCGR3A, CX3CR1, CXCR4, ADGRG1, GZMA, 
  GZMK, GZMH, GZMB, CCL3, CCL4, CCL5, PRF1, NKG7, RORA, KLRB1, 
  CD74, S100A4, S100A11, LGALS1, VIM, ANXA2, CD69, FAS, HNRNPLL, 
  GNLY, ITGB1, JUN, JUNB, JUND, FOS, FOSB, DUSP1, TNFAIP3, NFKBIA, 
  TSC22D3, CTSW, LDHB, ACTB, ACTG1, LIMS1, PASK, SOCS3, GATA3, 
  CCL4L2, CMC1, HLA-DRA, HLA-DRB1, HLA-DPB1, HLA-DPA1'
gene_list = unlist(strsplit(gene_list, ",")) 
gene_list = trimws(gene_list, which = c("both"), whitespace = "[ \t\r\n]")
gene_list [ !(gene_list %in% gene_name_vector )]

pathName <- paste0(filepath_cluster,paste0('HeatMap', '_genelist','.png'))
pdf(file=pathName, height=12, width=20)
#png(file=pathName)

plot = DoHeatmap(data_run_subset, features = gene_list)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()

########
## DE
########

newIdents = as.character(Idents(data_run_subset))
newIdents[newIdents == '8'] = '8_11'
newIdents[newIdents == '11'] = '8_11'
Idents(data_run_subset) = newIdents

# Seurat
ident1 = '8_11'
ident2 = '0'
ident1 = 'CD16+ Mono'
ident2 = 'SELL+ CD14+ Mono'
Features = FindMarkers(data_run_subset, ident.1 = ident1, ident.2 = ident2
                       ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)

path = paste0(filepath_cluster,'DE/')
dir.create( path, recursive = TRUE)
path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
print(path)
write.csv(Features, file = path,row.names=TRUE)

#

ident1 = 'CD16+ Mono'
ident2 = 'SELL+ CD14+ Mono'

celltype = 'SELL+ CD14+ Mono'
ident1 = paste0('baseline ',celltype)
ident2 = paste0('C9D1 ',celltype)

DE_input = data_run_subset_label
DE_input = renameCells(DE_input,idents = c('cDC1','cDC2'),newident = 'DC')
#DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
#                       newident = 'CD14+ Mono')
DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'

#DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
#                           DE_input$Best_Overall_Response, ' ', Idents(DE_input))
DE_input$DE_ident = paste0(DE_input$Treatment, ' ', Idents(DE_input))
#DE_input$DE_ident = paste0(Idents(DE_input))
#DE_input = DE_input[,DE_input$Treatment =='baseline']
DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
ncol(DE_input)
unique(DE_input$DE_ident)

DE_input_sce = as.SingleCellExperiment(DE_input)

data = as.data.frame(DE_input@assays[["RNA"]]@counts)
DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
DE_input$DR = DR
kit = factor(DE_input$kit)
DE_ident = factor(DE_input$DE_ident)
patient = factor(DE_input$`Patient Number`)

formula = ~ 0 + DE_ident + DR + kit
#formula = ~ 0 + DE_ident + DR + patient
design <- model.matrix(formula)
colnames(design) <- gsub("DE_ident", "", colnames(design))
colnames(design)


int1 = match(ident1,colnames(design) )
int2 = match(ident2,colnames(design) )
con <- integer(ncol(design))
con[int1] <- 1 
con[int2] <- -1 

subfolder = paste0(ident1,' Vs ',ident2)

# colnames(design) = gsub( '+','',colnames(design),fixed = T)
# colnames(design) = gsub( ' ','',colnames(design),fixed = T)
# colnames(design)
# con <- makeContrasts('CD16Mono - SELLCD14Mono', levels=design)

## EdgeR
library('org.Hs.eg.db')

result_edgeR = runEdgeR(DE_input,design, contrast = con, 
                        folder_output = filepath_cluster,
                        subfolder = subfolder)

dge_edgeR = result_edgeR[[1]]
fit_edgeR = result_edgeR[[2]]
qlf = result_edgeR[[3]]
res_edgeR = result_edgeR[[4]]


geneid = mapIds(org.Hs.eg.db, rownames(qlf), 'ENTREZID', 'SYMBOL')

go <- goana(qlf,geneid = geneid, species="Hs")
topGO(go, sort="up")
keg <- kegga(qlf,geneid = geneid, species="Hs")
topKEGG(keg, sort="up")

res_edgeR = res_edgeR[order(-res$logFC),]
res_edgeR[1:20,]
sum(res_edgeR$FDR < 0.05)
sum(res_edgeR$FDR > 0.05)
nrow(res_edgeR)
res_edgeR_sig = rownames(res_edgeR)
res_edgeR_sig = res_edgeR_sig[res_edgeR$FDR< 0.05]

base = paste0(filepath_cluster, 'DE/','EdgeR/')
plotEdgeR(dge_edgeR, DE_ident = DE_input$DE_ident,
          qlf,base,subfolder = subfolder)

path = paste0(base,ident1,' Vs ',ident2,'/',base,'qlf','.Robj')
qlf = loadRData(path)

# Voom-limma

result_VL= runVoomLimma(DE_input,design,
                        contrast = con, 
                        folder_output =filepath_cluster, subfolder= subfolder)

dge_VL = result_VL[[1]]
fit_VL = result_VL[[2]]
res_VL = result_VL[[3]]

res_VL[1:10,]
sum(res_VL$adj.P.Val < 0.05)
sum(res_VL$adj.P.Val > 0.05)
nrow(res_VL)

res_VL_sig = rownames(res_VL)
res_VL_sig = res_VL_sig[res_VL$adj.P.Val < 0.05]

plotMD(fit_VL)

## MAST
formula = ~ DE_ident + DR + kit
result_MAST = runMAST(DE_input,formula = formula, ident_list = DE_input$DE_ident, 
                      ident = ident2, 
                      folder_output = filepath_cluster, subfolder= subfolder)

dge_MAST = result_MAST[[1]]
summaryCond = result_MAST[[2]]
res_MAST = result_MAST[[3]]

sum(res_MAST$FDR < 0.05)
sum(res_MAST$FDR > 0.05)
nrow(res_MAST)

res_MAST_sig = res_MAST$primerid[res_MAST$FDR < 0.05]

## DEseq2
library(DESeq2)
library(tidyverse)

result_DESeq2 = runDESeq2(DE_input,design,contrast = con, 
                          folder_output = filepath_cluster, subfolder = subfolder)

dds = result_DESeq2[[1]]
res_DESeq2 = result_DESeq2[[2]]

plotDispEsts(dds)

res <- results(dds, 
               contrast = con,
               alpha = 0.05)

res <- lfcShrink(dds,type ='normal',contrast =  con,res=res)
res_DESeq2 = res_DESeq2[order(res_DESeq2$padj),]
sum(res_DESeq2$padj < 0.05, na.rm = T)
sum(res_DESeq2$padj > 0.05, na.rm = T)
nrow(res_DESeq2)
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_DESeq2_sig = res_DESeq2$gene[res_DESeq2$padj < 0.05 & !is.na(res_DESeq2$padj)]
#######################
doubletMethod = 'Doublet3Methods'
mean_doublets = 100*sum((data_run_subset@meta.data[,doubletMethod]))/ncol(data_run_subset)
title_str = paste0('% Doublets: ',specify_decimal(mean_doublets,2))
pathName <- paste0(filepath_cluster,paste0(doubletMethod,'_GeneralCellType','.png'))
png(file=pathName,width=1000, height=500,res = 100)
plot = DimPlot(data_run_subset,
               pt.size = ifelse(data_run_subset@meta.data[,doubletMethod] == T, 1, 0.5), cols = c('navajowhite2','tomato2'), 
               reduction = "umap",label = F, group.by  = doubletMethod) + 
  ggtitle(title_str)
plot = plot +  DimPlot(data_run_subset,pt.size = 0,
                       reduction = "umap",label = F, group.by  = 'GeneralCellType')
print(plot)
dev.off()

## Barplot for CD3+/CD14+ cluster by sample

cells  = colnames(data_run_subset)
cells = cells[Idents(data_run_subset) %in% c(8,17)]
samples = data_run_subset$sample[Idents(data_run_subset) %in% c(8,17)]

TcellMono_cells = data.frame(cbind(cells,samples))
TcellMono_cells$cell_sample = data_run_subset$cell_sample[Idents(data_run_subset) %in% c(8,17)]
write.csv(TcellMono_cells, file = paste0(filepath_cluster,'/Stats/','TcellMono_cells'))                                


# Stats


df = SummariseSampleByCell(data_run_subset_label)
df = table(data_run_subset_label$sample,Idents(data_run_subset_label) )
path = paste0(filepath_cluster,'Stats/SampleByCellType.csv')
write.csv(df, file = path,row.names=TRUE,col.names=TRUE)


# Matrix for NMF
path = paste0('/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit_',
              celltype,'.tsv')
data_matrix = data_run_subset@assays[["RNA"]]@data
data_matrix = data_matrix[rownames(data_matrix) %in% data_run_subset@assays[["RNA"]]@var.features,]
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')

# UmapCoord and metadata
sample = as.character(data_run_subset_label$sample)
ident = as.character(Idents(data_run_subset_label))
UmapCoord = data_run_subset_label@reductions[["umap"]]@cell.embeddings

output = cbind(sample,ident,UmapCoord)
path = paste0(filepath_cluster,'Stats/Umap.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)

