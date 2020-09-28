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
library(UpSetR)
library(SingleCellSignalR)
library(fgsea)
data(examplePathways)
data(exampleRanks)



source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')
source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')

celltype = 'T Cell'
cell_list = c('T Cell','CD8+ T Cell')
resolution_val_subset = 4
# 
celltype = 'Mono_DC'
cell_list = c('CD14+ Mono','CD16+ Mono','DC')
resolution_val_subset = 1.6
# # 
celltype = 'NK'
cell_list = c('NK')
resolution_val_subset = 3

# celltype = 'B_Cell'
# cell_list = c('B Cell','Pro B Cell','Pre B Cell')
# resolution_val_subset = 3.4

# celltype = 'Stem_Cell'
# cell_list = c('HSC','CMPC', 'GMPC','MDPC','42')
# resolution_val_subset = 3


filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]
#metaData = metaData[metaData$`Sample Type` == 'PBMC',]
#metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]


sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

#downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  = NA #downsample$x


sample_type = 'Harmony_AllSamples_Sample_Kit'
#sample_type = 'Harmony_PBMC_Sample_Kit'

folder_name = 'AllSamples'
#folder_name = 'PBMC'

PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 
cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
resolution_val_subset = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]

harmony_dim_subset = sampleParam_combine$harmony_dim[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]

patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

#folder = 'Intra-v3_1'
#folder = 'Inter-version'

#folder_name = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample

folder_main = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                '/Batch_Sample_Kit/','/')
dir.create(folder_main,recursive = T)
folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',30,'/res',resolution_val_subset,'/' )
dir.create(filepath_cluster,recursive = T)



run = F

folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')


if (run){
  path = paste0(folder_main,'data_run','.Robj')
  data_orig = loadRData(path)
  data_orig =label_cells(data_orig,cluster_IDs)
  data_orig = data_orig[,Idents(data_orig) %in% cell_list]
  plot = DimPlot(data_orig,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  print(plot)
  
  cell_names_all = colnames(data_orig)
  cell_names_all = sub("_.*", "", cell_names_all)
  data_orig$cell_sample = paste0(cell_names_all,' ',data_orig$sample )
  
  data_merge_subset = data_orig
  data_merge_subset = addMetaData(data_merge_subset, metaData)
  data_merge_subset = load_emptyDrops(data_merge_subset)
  data_merge_subset = load_Doublets(data_merge_subset)
  data_merge_subset = load_CellLabel(data_merge_subset)
  data_merge_subset$origIdent = Idents(data_orig)
  data_merge_subset$kit = data_merge_subset$`10X kit`
  data_merge_subset$split_var = ''
  
  data_run_subset = FindVariableFeatures(data_merge_subset, selection.method = "vst", nfeatures = 2000)
  data_run_subset = ScaleData(data_run_subset)
  data_run_subset = RunPCA(data_run_subset,npcs = 40)
  
  reduction = 'pca'
  visualize_PCA(data_run_subset,folder_subcluster,40,reduction)
  dir.create( paste0(folder_subcluster, reduction), recursive = TRUE)
  pathName <- paste0(folder_subcluster,reduction,'/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data_run_subset,ndims = PCA_dim))
  dev.off()
  
  data_run_subset = RunHarmony(data_run_subset,group.by.vars =  c("sample", "10X kit"),
                              dims.use = 1:harmony_dim_subset)
  
  reduction = 'harmony'
  visualize_PCA(data_run_subset,folder_subcluster,harmony_dim_subset,reduction)
  dir.create( paste0(folder_subcluster, reduction), recursive = TRUE)
  pathName <- paste0(folder_subcluster,reduction,'/elbow_',PCA_dim,'.png')
  png(file=pathName,width=600, height=350)
  print(ElbowPlot(data_run_subset,ndims = harmony_dim_subset))
  dev.off()
  
  
  data_run_subset = RunUMAP(reduction = "harmony",data_run_subset, dims = 1:harmony_dim_subset)
  data_run_subset = FindNeighbors(data_run_subset, reduction = "harmony", dims = 1:harmony_dim_subset)
  data_run_subset = FindClusters(data_run_subset,resolution = resolution_val_subset)
  data_run_subset$origIdent = Idents(data_orig)
  
  
  dir.create(folder_subcluster,recursive = T)
  path = paste0(folder_subcluster,'data_run','.Robj')
  save(data_run_subset,file= path)
  
}else{
  folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  #data_run_subset = FindClusters(data_run_subset,resolution = resolution_val_subset)
  plot = DimPlot(data_run_subset,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  print(plot)
  
}
gene_name_vector= rownames(data_run_subset)
gene_name_vector[grep("HLA-DR", gene_name_vector)]
gene_list = 'LY6C1, LY6G6D, MMP9, TIMP1, LGALS3, S100A6, S100A11, ISG15, MS4A4C'
gene_list = 'NCR2'
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

cell_features_file = paste0('/home/sujwary/Desktop/scRNA','/Data/Cell_IDS.xlsx')
cell_features = read_excel(cell_features_file)
plot = DimPlot(data_run_subset_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)


plotAll(data_run_subset, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF =F, 
        clusterTF =F, markersTF = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val_subset,
        pt.size = 2)

#data_run_subset_label = data_run_subset_label[,Idents(data_run_subset_label) != 14]
plotAll(data_run_subset_label, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = T, 
        clusterTF =F, markersTF = F,  
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = resolution_val_subset, pt.size = 2,
        str = '_label')


PlotKnownMarkers(data_run_subset, folder = paste0(filepath_cluster,'Cell Type/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1)

stats_category = SummariseSampleByCell(data_run_subset_label)
CompareCellNum(data_run_subset_label,stats_category,filepath_cluster,split_var = 'Treatment',metaData)
  
### Plot DE

pathways.hallmark <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/h.all.v7.2.symbols.gmt')
pathways.c2 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c2.all.v7.2.symbols.gmt')
pathways.c5 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c5.all.v7.2.symbols.gmt')
pathways.c7 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c7.all.v7.2.symbols.gmt')
pathways.c8 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c8.all.v7.2.symbols.gmt')
pathways = c(pathways.hallmark,pathways.c2,pathways.c5, pathways.c7,pathways.c8)


filepath = paste0(filepath_cluster
                  ,'Features',
                  '.csv')
DE = read.csv(file = filepath)
cluster_list = unique(DE$cluster)
str = ''
DE = DE[DE$p_val_adj < 0.05,]
DE = DE[!grepl("MT-", DE$gene),]
DE = DE[!grepl(" ?RP\\w+ ?", DE$gene),]

for (cluster in cluster_list){
  DE_cluster = DE[DE$cluster == cluster,]
  ranks <- DE_cluster$avg_logFC
  names(ranks) <- DE_cluster$gene
  
  DE_cluster = DE_cluster[DE_cluster$avg_logFC > 0 ,]
  
  
  
  gene_list = DE_cluster
  #gene_list = as.data.frame(as.character(gene_list[!is.na(gene_list)]))
  #colnames(gene_list) = 'Markers'
  
  if (nrow(gene_list) > 3){
    path = paste0(filepath_cluster,'DE/GSEA/')
    dir.create(path,recursive = T)
    path = paste0(path,cluster,'VsAll',str,'GSEA.csv')
    
    fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)
    fgseaRes  = fgseaRes[order(fgseaRes$padj, fgseaRes$pval),]
    fgseaRes$leadingEdge =  sapply( fgseaRes$leadingEdge , paste0, collapse=",")
    write.csv(fgseaRes, file = path,row.names=F)
    #next
    
    gene_list$Cell = cluster
    gene_list$Plot_marker = 1
    path = paste0(filepath_cluster,'DE/Cluster',cluster,'VsAll',str,'/')
    PlotKnownMarkers(data_run_subset, 
                     folder = path, 
                     cell_features = gene_list,
                     plotType ='FeaturePlotFix' , 
                     prefix_logFC = T, str = '',markerSize = 1)
  }
  
  
}


### Plot DE permutations
data_input = data_run_subset_label
cluster_list = sort(unique(Idents(data_input)))
cluster1 = 23
str = '_label'
str = ''
for (cluster1 in cluster_list){
  #
  cluster1 = gsub("/", "", cluster1, fixed = T)
  filepath_mvn = paste0( filepath_cluster, 'DE/nVsm/')
  filepath = paste0(filepath_mvn
                    ,'Features_',cluster1,'Vsn',str,
                    '.csv')
  DE_cluster = read.csv(file = filepath)
  
  for (cluster2 in cluster_list){
    #data_cluster  = data_run_subset[, Idents(data_run_subset) == cluster]
    #print(cluster)
    #print(summary(data_cluster$origIdent))
    cluster2 = gsub("/", "", cluster2, fixed = T)
    DE_cluster_m = DE_cluster[DE_cluster$ident_2 == cluster2,]
    DE_cluster_m = DE_cluster_m[order(-abs(DE_cluster_m$avg_logFC )),]
    DE_cluster_m = DE_cluster_m[DE_cluster_m$p_val_adj < 0.05 & DE_cluster_m$avg_logFC > 0 ,]
    DE_cluster_m = DE_cluster_m[!grepl("MT-", DE_cluster_m$gene),]
    DE_cluster_m = DE_cluster_m[!grepl(" ?RP\\w+ ?", DE_cluster_m$gene),]

    gene_list = DE_cluster_m$gene[1:20]
    gene_list = as.data.frame(as.character(gene_list[!is.na(gene_list)]))
    colnames(gene_list) = 'Markers'
    if (nrow(gene_list) > 0){
      gene_list$Cell = paste0(cluster1,'Vs',cluster2)
      gene_list$Plot_marker = 1
      path = paste0(filepath_mvn,'FeaturePlot/Cluster',cluster1,str,'/')
      PlotKnownMarkers(data_input, 
                       folder = path, 
                       cell_features = gene_list,
                       plotType ='FeaturePlotFix' , str = '',markerSize = 2)
    }
  }
  
  
}

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

celltype = 'CD16+ Mono'
ident1 = paste0('baseline ','GR ',celltype)
ident2 = paste0('C9D1 ','PR ',celltype)

DE_input = data_run_subset_label
DE_input = renameCells(DE_input,idents = c('cDC1','cDC2'),newident = 'DC')
#DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
#                       newident = 'CD14+ Mono')
DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'

DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
                          DE_input$Best_Overall_Response, ' ', Idents(DE_input))
#DE_input$DE_ident = paste0(DE_input$Treatment, ' ', Idents(DE_input))
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


keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                     min.count = 1,min.total.count=10, 
                     large.n = 10,min.prop = 0.1)


countpergene_mean <- rowMeans(DE_input@assays[["RNA"]]@counts)

#keep = countpergene_mean > 0.02

sum(keep)

subfolder = paste0(ident1,' Vs ',ident2)

# colnames(design) = gsub( '+','',colnames(design),fixed = T)
# colnames(design) = gsub( ' ','',colnames(design),fixed = T)
# colnames(design)
# con <- makeContrasts('CD16Mono - SELLCD14Mono', levels=design)

## EdgeR
library('org.Hs.eg.db')

result_edgeR = runEdgeR(DE_input,design, contrast = con, keep = keep,
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
                        contrast = con, keep = keep,
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
                      ident = ident2, keep = keep,
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
library(SciViews)
library(tibble)
result_DESeq2 = runDESeq2(DE_input,design,contrast = con, keep = keep,
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

res_DESeq2_sig = res_DESeq2$gene[res_DESeq2$padj < 0.05 & !is.na(res_DESeq2$padj)]


# Upset plot
folder_output = paste0(filepath_cluster,
                       '/DE/Test DE/LowThreshold/',subfolder,'/')
dir.create(folder_output,recursive = T)
allSigGenes_isSig = read.csv(file = paste0(folder_output,'allSigGenes_manual.csv'),row.names = 1)  
allSigGenes_isSig$x = as.character(allSigGenes_isSig$x)
SigGenes = allSigGenes_isSig$x[as.logical(allSigGenes_isSig$IsSig)]

allgenes = rownames(res_edgeR)
EdgeR_geneList =  allgenes %in% res_edgeR_sig 
MAST_geneList =  allgenes %in% res_MAST_sig 
VL_geneList =  allgenes %in% res_VL_sig 
DESeq2_geneList =  allgenes %in% res_DESeq2_sig 
SigGenes_geneList = allgenes %in% SigGenes
sigGeneData = data.frame(cbind(EdgeR_geneList,
                               MAST_geneList,
                               VL_geneList, 
                               DESeq2_geneList, SigGenes_geneList
                               ))

res_edgeR_sig_LowThreshold = read.csv(paste0(filepath_cluster,'DE/EdgeR/LowThreshold/',
                                             subfolder,'/','DE_EdgeR ',ident1,' Vs ', ident2,'.csv'))
res_MAST_sig_LowThreshold = read.csv(paste0(filepath_cluster,'DE/MAST/LowThreshold/',
                                             subfolder,'/fcHurdle','.csv'))
res_VL_sig_LowThreshold = read.csv(paste0(filepath_cluster,'DE/VoomLimma/LowThreshold/',
                                             subfolder,'/','DE_Voom ',ident1,' Vs ', ident2,'.csv'))
res_DESeq2_sig_LowThreshold = read.csv(paste0(filepath_cluster,'DE/DESeq2/LowThreshold/',
                                             subfolder,'/','DE_DESeq2 ',ident1,' Vs ', ident2,'.csv'))
length(res_edgeR_sig)
sum(res_edgeR_sig %in% res_edgeR_sig_LowThreshold$X[res_edgeR_sig_LowThreshold$FDR<0.05])
length(res_MAST_sig)
sum(res_MAST_sig %in% res_MAST_sig_LowThreshold$primerid[res_MAST_sig_LowThreshold$FDR<0.05])
length(res_VL_sig)
sum(res_VL_sig %in% res_VL_sig_LowThreshold$X[res_VL_sig_LowThreshold$adj.P.Val<0.05])
length(res_DESeq2_sig)
sum(res_DESeq2_sig %in% res_DESeq2_sig_LowThreshold$X[res_DESeq2_sig_LowThreshold$padj<0.05])

colnames(sigGeneData) = c('EdgeR','MAST','VL','DESeq2', 'SigGenes')
colnames(sigGeneData) = c('EdgeR','MAST','VL','DESeq2')
sigGeneData = sigGeneData*1

pathName <- paste0(filepath_cluster,'DE/UpsetSigGenes','.png')
png(file=pathName,width=1000, height=500,res = 100)
plot = upset(sigGeneData, sets =  c('EdgeR','MAST','DESeq2','VL','SigGenes'),
             point.size = 3.5, line.size = 2, order.by = "freq",
             mainbar.y.label = "Gene", sets.x.label = "Sig Genes per Method")
print(plot)
dev.off()



allSigGenes = sort(unique(c(res_edgeR_sig,
                 res_MAST_sig, 
                 res_DESeq2_sig)))

write.csv(allSigGenes, file = paste0(folder_output,
                                     'allSigGenes','.csv'))                                


# Venn
library(VennDiagram)

library(RColorBrewer)
myCol <- brewer.pal(5, "Pastel2")

venn.diagram(
  x = list(res_edgeR_sig, res_MAST_sig, res_DESeq2_sig,res_VL_sig,SigGenes),
  category.names = c("EdgeR" , "MAST" , "DESeq2",'VL','True Positive'), 
  filename = paste0(folder_output,'Venn.png'),
  imagetype="png",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol

)
gene_input = res_DESeq2_sig
sum((gene_input %in% SigGenes ))/length(gene_input)*100
sum(!(gene_input %in% SigGenes ))/length(gene_input)*100


#color_list = c('Blue','Red')
group_var = 'DE_ident'
DE_input = ScaleData(DE_input, features = rownames(DE_input))
for (gene in allSigGenes){
  print(gene)
  logFC = specify_decimal(res_edgeR[gene,]$logFC,2)
  
  title = paste0('EdgeR logFC: ', logFC)
  folder_name = paste0('ViolinPlots/',group_var,'/')
  folder = paste0(folder_output, folder_name,'/')
  dir.create(folder,recursive = TRUE)
  
  pathName = as.character(paste0(folder,'Violin_',gene,'_Split_',group_var,'.png') )
  png(file=pathName,width=900, height=600)
  plot= VlnPlot(DE_input, slot = 'scale.data', gene, pt.size = 0.1,group.by	 = group_var)
  plot = plot + ggtitle(title) +
    theme_classic() + theme(
    title =element_text(size=24, color="black"),
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=18),
    legend.title=element_text(size=18),
    legend.position = 'none'
  )
  
  print(plot)
  dev.off()
  
  #boxplot 
  folder_name = paste0('BoxPlots/',group_var,'/')
  folder = paste0(folder_output, folder_name,'/')
  dir.create(folder,recursive = TRUE)
  data = as.data.frame( DE_input@assays[["RNA"]]@scale.data[gene,])
  colnames(data) = 'count'
  data$DE_ident = DE_input$DE_ident
  pathName = as.character(paste0(folder,'Boxplot_ ',gene,'_Split_',group_var,'.png') )
  png(file=pathName,width=900, height=600)
  plot = ggplot(data, aes(x=DE_ident, y=count)) + 
    geom_boxplot()
  plot = plot + ggtitle(title) +
    theme_classic()  + theme(
    title =element_text(size=24, color="black"),
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=18),
    legend.title=element_text(size=18),
    legend.position = 'none'
  )
  print(plot)
  dev.off()
}
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
# metadata
sample = as.character(data_run_subset_label$sample)
ident = as.character(Idents(data_run_subset_label))
UmapCoord = data_run_subset_label@reductions[["umap"]]@cell.embeddings

output = cbind(sample,ident,UmapCoord)
#output = data_run_subset_label@meta.data
dir.create(paste0(filepath_cluster,'Data/'))
path = paste0(filepath_cluster,'Data/Umap.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)


## AUCell

exprMatrix = data_run_subset@assays[["RNA"]]@counts

library(AUCell)
library(GSEABase)
gmtFile <- paste0('/home/sujwary/Desktop/scRNA/Data/msigdb.v7.1.symbols.gmt')
geneSets <- getGmt(gmtFile)

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE) 

cells_assignment

path = paste0(folder_subcluster,'cells_AUC','.Robj')
save(cells_AUC,file= path)

cellAssigned <- cells_assignment$Oligodendrocyte_Cahoy$assignment

geneSetName <- rownames(cells_AUC)[grep("Oligodendrocyte_Cahoy", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]

miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),100)]
library(NMF)
aheatmap(miniAssigMat, scale="none", color="black", legend=FALSE)


## scVelo

data_input = data_run_subset
data_input = seurat[,!(Idents(data_input) %in% c(0,8,11,17,18,21,22))]
data_input$sample_cell = paste0(data_input$Sample,'_', colnames(data_input))

write.csv(unname(data_input$sample_cell), 
          file = paste0(filepath_cluster,"Data/cellID_obs.csv"))
write.csv(Embeddings(data_input, reduction = "umap"), 
          file = paste0(filepath_cluster,"Data/cell_embeddings.csv"))
write.csv(Idents(data_input), file = paste0(filepath_cluster,"Data/clusters.csv"))

library(loomR)

for(j in 1:ncol(data_input@meta.data)){
  if(is.factor(data_input@meta.data[,j]) == T){
    data_input@meta.data[,j] = as.character(data_input@meta.data[,j]) # Force the variable type to be character
    data_input@meta.data[,j][is.na(data_input@meta.data[,j])] <- "NA"
  }
  if(is.character(data_input@meta.data[,j]) == T){
    data_input@meta.data[,j][is.na(data_input@meta.data[,j])] <- "N.A"
  }
}
meta_data = data_input@meta.data
data_input@meta.data = meta_data[,c("orig.ident","nCount_RNA" ,"percent.mt" , "seurat_clusters")]
pfile <- as.loom(data_input, filename = paste0(filepath_cluster,"Data/data.loom"), verbose = FALSE)
 