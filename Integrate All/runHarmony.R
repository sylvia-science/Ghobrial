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
library(UpSetR)
library(edgeR)
library(SingleCellExperiment)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')
source('/home/sujwary/Desktop/scRNA/Code/Visualization/PlotCellPhoneDB.R')

filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
#metaData = metaData[metaData$Run== 1,]
metaData = metaData[metaData$`Sample Type` == 'PBMC',]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  = NA #downsample$x


sample_type = 'Harmony_AllSamples_Sample_Kit'
sample_type = 'Harmony_PBMC_Sample_Kit'

PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 


patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

#folder = 'Intra-v3_1'
#folder = 'Inter-version'
folder_name = 'AllSamples'
folder_name = 'PBMC'

#folder_name = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample

sample_list = sample_list[c(1, 3:8)]

folder = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                '/Batch_Sample_Kit/','/')
dir.create(folder,recursive = T)


run = F

if (run){
  data_list = vector(mode = "list",length = length(sample_list))
  data_list_norm = vector(mode = "list",length = length(sample_list))
  #for (i in 1:nrow(sampleParam)){
  for (i in 1:length(sample_list)){
    #sample_name = sampleParam$Sample[i]
    sample_name = sample_list[i]
    print(sample_name)
    folder_input = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/', sample_name , '/')
    data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
    data_i$sample = sample_name

    path = paste0('/disk2/Projects/EloRD/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
    cellIdents = read.csv(path,sep = ',',row.names = 1)
    cellIdents$x = paste0(cellIdents$x, ' S', i)
    data_i$CellType = cellIdents
    data_i = data_i[,data_i$nFeature_RNA > 200]
    data_i = load_emptyDrops(data_i)
    data_i = data_i[,data_i$is_cell]
    # remove empty drops
    # Remove MT > 15 already done
    if (!is.na(downsample)){
      downsample = sub("_.*", "", downsample)
      cellnames = colnames(data_i)
      cellnames = sub("_.*", "", cellnames)
      data_i = data_i[,cellnames %in% downsample]
      #browser()
    }
    
    print(ncol(data_i))
    if (ncol(data_i) > 100){
      data_list[[i]] = data_i
      data_list_norm[[i]] = ScranNorm(data_i)
      
    }
    
  }
  
  data_list_norm =data_list_norm[lengths(data_list_norm) != 0]
  data_merge = merge(x =  data_list_norm[[1]] ,y = data_list_norm[2:length(data_list_norm)], merge.data = T)
  
  data_merge = addMetaData(data_merge, metaData)
  data_merge = load_emptyDrops(data_merge)
  data_merge = load_Doublets(data_merge)
  data_merge = load_CellLabel(data_merge)
  data_merge$GeneralCellType = str_match(data_merge$CellType, "(^.+)\\s")[, 2]
  data_merge$kit = data_merge$`10X kit`
  data_merge$split_var = ''
  cell_names_all = sub("_.*", "", colnames(data_merge))
  data_merge$cell_sample = paste0(cell_names_all,' ',data_merge$sample )
  
  
  data_merge_run = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
  data_merge_run = ScaleData(data_merge_run)
  data_merge_run = RunPCA(data_merge_run,npcs = 40)

  data_merge_run = RunHarmony(data_merge_run,group.by.vars =  c("sample", "10X kit"),
                              dims.use = 1:30)
  
  data_merge_run = RunUMAP(reduction = "harmony",data_merge_run, dims = 1:30)
  data_merge_run = FindNeighbors(data_merge_run, reduction = "harmony", dims = 1:30)
  data_merge_run = FindClusters(data_merge_run,resolution = resolution_val)
  

  
  path = paste0(folder,'data_run','.Robj')
  save(data_merge_run,file= path)

}else{
  resolution_val = 3
  
  path = paste0(folder,'data_run','.Robj')
  data_merge_run = loadRData(path)
  tmp = data_merge_run@meta.data[paste0('RNA_snn_res.', resolution_val)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_merge_run) = tmp
  data_merge_run_label = label_cells(data_merge_run,cluster_IDs)
  

  celltype = 'Mono_DC'
  cell_list = c('CD14+ Mono','CD16+ Mono','DC')
  resolution_val_subset = 1.6
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  
  Ident_main = colnames(data_merge_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)

  newIdents = as.character(Idents(data_merge_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_merge_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_merge_run_label) = newIdents2
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
  
  #(0,48)
  data_merge_run = FindClusters(data_merge_run,resolution = resolution_val)
  

}


#metaData = read_excel(filename_metaData)
data_merge_run = addMetaData(data_merge_run, metaData)
#data_merge_run = load_emptyDrops(data_merge_run)
data_merge_run = load_Doublets(data_merge_run)
data_merge_run = load_CellLabel(data_merge_run)
print(unique(data_merge_run$GeneralCellType))
#data_merge_run$GeneralCellType = str_match(data_merge_run$CellType, "(^.+)\\s")[, 2]
#data_merge_run$kit = data_merge_run$`10X kit`
#data_merge_run$split_var = ''
#data_merge_run_label = label_cells(data_merge_run,cluster_IDs)
data_merge_run$FeatureLessThan400 = data_merge_run$nFeature_RNA < 400

#data_merge_run_label = data_merge_run_label[,data_merge_run_label$is_cell]
#data_merge_run_label = data_merge_run_label[,Idents(data_merge_run_label) != 'Remove']


plot = DimPlot(data_merge_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)

print(plot)

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')

groupBy_list = c('sample','Diagnosis','kit',
                 'Treatment','Batch','LowCount',
                 'Doublet','GeneralCellType','FeatureLessThan400')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
splitBy_list = NA


plotAll(data_merge_run, folder = folder, 
        sample_name,sampleParam = NA,
        cell_features = cell_features, plot_PCA = T,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = T,
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val)


plotAll(data_merge_run_label, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = resolution_val,str = '_label')

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
PlotKnownMarkers(data_merge_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',plotTogether = F)
print(DimPlot(data_merge_run, reduction = "umap", label = TRUE, pt.size = .1))

gene_list = c('IRF8','KLF4')
gene_list = c('FCGR2A','CSF1R','FLT3','NR4A1')

FeaturePlot_GeneList(data_merge_run,gene_list,
                     folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'),
                     FeaturePlotFix = T,str = 'gene_')


group = 'CellType'
pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=5000, height=1500)
plot = DimPlot(data_merge_run,pt.size = 1, reduction = "umap",label = FALSE,group.by  = group)

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
)


#plot = plot + scale_color_manual(values=color_list)

print(plot)

dev.off()



### Plot DE
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
  DE_cluster = DE_cluster[DE_cluster$avg_logFC > 0 ,]
  
  
  gene_list = DE_cluster$gene
  gene_list = as.data.frame(as.character(gene_list[!is.na(gene_list)]))
  colnames(gene_list) = 'Markers'
  if (nrow(gene_list) > 0){
    gene_list$Cell = cluster
    gene_list$Plot_marker = 1
    path = paste0(filepath_cluster,'FeaturePlot/Cluster',cluster,'VsAll',str,'/')
    PlotKnownMarkers(data_merge_run, 
                     folder = path, 
                     cell_features = gene_list,
                     plotType ='FeaturePlotFix' , str = '',plotTogether = F)
  }
  
  
}

#######################
## DE
#######################
# Detection rate
data = as.data.frame(data_merge_run_label@assays[["RNA"]]@counts)
DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
data_merge_run_label$DR = DR
kit = factor(data_merge_run_label$kit)

#  Identifies differentially expressed genes between two groups of cells 
#  using a hurdle model tailored to scRNA-seq data. 
# Utilizes the MAST package to run the DE testing.

# MAST
markers = FindMarkers(data_merge_run_label,  
                      ident.1 = "CD8+ T Cell", ident.2 = "T Cell", 
                      latent.vars = 'kit', test.use = "MAST")

##
# edgeR
data_merge_run_label_sce = as.SingleCellExperiment(data_merge_run_label)
## Convert to DGEList, calculate logCPMs
dge <- scran::convertTo(data_merge_run_label_sce, type = "edgeR")
plotMDS(dge)

design <- model.matrix(~kit+DR)

dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)
topTags(lrt)

##
y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- limma::lmFit(y, design)



# voom limma
dge <- DGEList(data_merge_run_label@assays[["RNA"]]@counts, 
               group = data_merge_run_label$kit)
dge <- calcNormFactors(dge)

design <- model.matrix(~ kit + DR)

vm <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

hist(tt$P.Value, 50)
hist(tt$adj.P.Val, 50)
#limma::plotMDS(dge, col = as.numeric(as.factor(data_merge_run_label$kit)), pch = 19)
plotMD(fit)

# Cluster comparisions using design matrix made from cell idetents
design <- model.matrix(~ 0 + ident, data = colData(data_merge_run_label_sce))
colnames(design) <- gsub("ident", "", colnames(design))
colnames(design)

y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- limma::lmFit(y, design)

## Perform pairwise comparisons
nclust <- length(unique(data_merge_run_label_sce$ident))
all.results <- all.pairs <- list()
counter <- 1

for (i in seq_len(nclust)) {
  for (j in seq_len(i - 1L)) {
    con <- integer(ncol(design))
    con[i] <- 1
    con[j] <- -1
    fit2 <- limma::contrasts.fit(fit, con)
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
    
    res <- limma::topTable(fit2, number = Inf, sort.by = "none")
    all.results[[counter]] <- res
    all.pairs[[counter]] <- colnames(design)[c(i, j)]
    counter <- counter + 1L
    
    ## Also filling the reverse comparison.
    res$logFC <- -res$logFC
    all.results[[counter]] <- res
    all.pairs[[counter]] <- colnames(design)[c(j, i)]
    counter <- counter + 1L
  }
}

## Combine results across all pairwise tests
all.pairs <- do.call(rbind, all.pairs)
combined <- scran::combineMarkers(all.results, all.pairs, 
                                  pval.field = "P.Value",
                                  pval.type = "any")
head(combined[["cluster1"]])

markers <- scran:::findMarkers(data_merge_run_label_sce, design=design)

#scran:::findMarkers()

######################
cluster_list = levels(unique(Idents(data_merge_run)))

for (cluster in cluster_list){
  data_subset = data_merge_run[,Idents(data_merge_run) == cluster]
  print(cluster)
  cell_percent =  100*table(data_subset$GeneralCellType)/ncol(data_subset)
  print(cell_percent)
  print(max(cell_percent))
  print('')
}

for (cluster in cluster_list){
  data_subset = data_merge_run[,Idents(data_merge_run) == cluster]
  print(cluster)
  sample_percent =  100*table(data_subset$sample)/ncol(data_subset)
  print(sample_percent)
  print(max(sample_percent))
  print('')
}

filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
k_num = 30
dim_num = 30
folder_output = paste0(filepath_cluster,'Entropy/')
dir.create(folder_output)
folder_output = paste0(folder_output,'Harmony')


#data_merge_run_label_downsample = data_merge_run_label[,cellNames]

compute_entropy_Seurat(data_merge_run, 
                       corrected_assay = 'harmony',folder = folder_output,
                       k_num = k_num, dim_num = dim_num)


file_entropy = paste0(folder_output,'_k',k_num ,'_entropy.csv')
entropy = read.csv(file_entropy,sep = ',')

entropy$Method = 'Harmony'

entropy$batch_entropy <- NULL
plotEntropy(entropy,folder_output)


##########################
sample_summary =table(data_merge_run$sample)
write.csv(sample_summary,file = paste0(filepath_cluster,'Stats/','celltype_summary','.csv'))

####################
## Doublets
####################
doubletMethod = 'Doublet3Methods'
mean_doublets = 100*sum(data_merge_run$Doublet3Methods)/ncol(data_merge_run)
title_str = paste0('% Doublets: ',specify_decimal(mean_doublets,2))
pathName <- paste0(filepath_cluster,paste0(doubletMethod,'_GeneralCellType','.png'))
png(file=pathName,width=1000, height=500,res = 100)
plot = DimPlot(data_merge_run,
              pt.size = ifelse(data_merge_run$Doublet3Methods == T, 1, 0.5), cols = c('navajowhite2','tomato2'), 
              reduction = "umap",label = F, group.by  = doubletMethod) + 
              ggtitle(title_str)
plot = plot +  DimPlot(data_merge_run,pt.size = 0,
               reduction = "umap",label = T, group.by  = 'GeneralCellType')
print(plot)
dev.off()


doublet_data = data.frame(cbind(data_merge_run$scrublet,
                                data_merge_run$doublet_finder, 
                                data_merge_run$scran_doublet, 
                                data_merge_run$scds_doublet))
colnames(doublet_data) = c('scrublet','doublet_finder','scran_doublet','scds_doublet')

doublet_data = doublet_data*1

pathName <- paste0(filepath_cluster,paste0('Stats/','doublet_upset','.pdf'))
pdf(file=pathName)
plot = upset(doublet_data, sets =  c('scrublet','doublet_finder','scran_doublet','scds_doublet'),
      point.size = 3.5, line.size = 2, order.by = "freq",
      mainbar.y.label = "Doublet Intersections", sets.x.label = "Doublets Per method")
print(plot)
dev.off()

##########################
##
##########################
pathName <- paste0(filepath_cluster,paste0('Empty_GeneralCellType','.png'))
png(file=pathName,width=1000, height=500,res = 100)
plot = DimPlot(data_merge_run,
               pt.size = ifelse(data_merge_run$is_cell == T, 0.5, 1), cols = c('tomato2','navajowhite2'), 
               reduction = "umap",label = F, group.by  = 'is_cell')
plot = plot +  DimPlot(data_merge_run,pt.size = 0,
                       reduction = "umap",label = T, group.by  = 'GeneralCellType')
print(plot)
dev.off()

downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
# Plot scatterplots of gene


gene_list = c('HBB','IGKC','IGLL5')

gene1 = 'HBB'
for (gene1 in gene_list){
  plotGeneScatter(data_merge_run,gene1, 'PTPRC')
  plotGeneScatter(data_merge_run,gene1, 'CD3D')
  plotGeneScatter(data_merge_run,gene1, 'CD14')
  plotGeneScatter(data_merge_run,gene1, 'FCGR3A')
  plotGeneScatter(data_merge_run,gene1, 'NKG7')
}


## Get T cell, mono stats
path = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/Mono_DC/Cluster/PCA30/res1.6/Stats/TcellMono_cells'
TcellMono_cells = read.csv(path)

unique(TcellMono_cells$samples)

orig_cell_labels = as.character(Idents(data_merge_run_label))

orig_cell_labels[data_merge_run_label$cell_sample %in% TcellMono_cells$cell_sample] = 'T_Mono'
orig_cell_labels[orig_cell_labels == 'CD8+ T Cell'] = 'T Cell'
orig_cell_labels[orig_cell_labels == 'CD14+ Mono' | orig_cell_labels == 'CD16+ Mono'] = 'Mono'

df = data.frame(cbind(orig_cell_labels,data_merge_run_label$sample))

colnames(df) = c('Cell','Sample')

df = df[df$Cell %in% c('T Cell','Mono','T_Mono'),]
df$Cell = as.character(df$Cell)
df$Sample = as.character(df$Sample)
#df  = df[df$Sample %in% unique(TcellMono_cells$samples),]
df_stats = t(table(df))

write.csv(df_stats,file = paste0(filepath_cluster,'Stats/','T_Mono_Stats.csv'))



path = paste0('/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit',
              '.tsv')
data_matrix = data_merge_run@assays[["RNA"]]@data
data_matrix = data_matrix[rownames(data_matrix) %in% data_merge_run@assays[["RNA"]]@var.features,]
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')

# Save all data
#data_merge_run_label = data_merge_run_label[,Idents(data_merge_run_label) %in% c('NK','CD14+ Mono')]
data_merge_run_label = label_cells(data_merge_run,cluster_IDs)
#data_merge_run_label = data_merge_run_label[1:1000,1:1000]

data_mere_run_label = data_main_label

data_merge_run_label_input = data_merge_run_label[,data_merge_run_label$Treatment == 'baseline']
library("org.Hs.eg.db") # remember to install it if you don't have it already

data_matrix = as.data.frame(data_merge_run_label_input@assays[["RNA"]]@data)
gene_ensemble= mapIds(org.Hs.eg.db, keys = rownames(data_matrix), keytype = "SYMBOL", column="ENSEMBL", multiVals = "asNA")
gene_ensemble = gene_ensemble[!is.na(gene_ensemble)]
gene_ensemble = gene_ensemble[!duplicated(unname(gene_ensemble)) ] # Keeps only first dup
data_matrix = data_matrix[names(gene_ensemble),]
rownames(data_matrix) = unname(gene_ensemble)

path = paste0(filepath_cluster,'/data/')
dir.create(path,recursive = T)

str = '_baseline'

write.csv(data_matrix, 
            file=paste0(path,'Harmony_AllSamples_Sample_Kit',str,'.csv'), 
            quote=FALSE, row.names = T)

write.table(data_matrix, 
            file=paste0(path,'Harmony_AllSamples_Sample_Kit',str,'.tsv'), 
            quote=FALSE, sep='\t', row.names = T, col.names=NA)

labels = as.character(Idents(data_merge_run_label_input))
write.csv(labels, file = paste0(path,'Labels.csv'), )

labels = as.data.frame(colnames(data_matrix))
colnames(labels) = 'Cell'
labels$cell_type = (as.character(Idents(data_merge_run_label_input)))
write.table(labels, file = paste0(path,'Labels',str,'.tsv'), quote=FALSE, sep='\t', row.names= F)
write.csv(labels, file = paste0(path,'Labels',str,'.csv'), quote=FALSE, row.names= F)


tmp1 = colnames(data_matrix)
tmp2 = labels$Cell

# metadata
sample = as.character(data_merge_run_label$sample)
ident = as.character(Idents(data_merge_run_label))
UmapCoord = data_merge_run_label@reductions[["umap"]]@cell.embeddings

output = cbind(sample,ident,UmapCoord)
output = data_merge_run_label@meta.data
dir.create(paste0(filepath_cluster,'data/'))
path = paste0(filepath_cluster,'Data/Umap.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)

# CellPhoneDB
output = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3/CellPhoneDB/out baseline/'
setwd(output)
dot_plot(selected_rows = NULL,
         selected_columns = c('CD8+ T Cell|IFN+ Mono'),
         filename = 'plot.pdf',
         width = 8,
         height = 10,
         means_path = paste0(output,'/means.txt'),
         pvalues_path = paste0(output,'/pvalues.txt'),
         means_separator = '\t',
         pvalues_separator = '\t',
         output_extension = '.pdf'
)
