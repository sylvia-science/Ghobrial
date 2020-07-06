library(DropletUtils)


library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
library(annotables)
library(tidyverse)
library(rtracklayer)
library("biomaRt")
library("dplyr")
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

foldername = 'Oct 10 '
foldername = 'Nov 9 '
filename_list = list.files(path = paste0("/home/sujwary/Desktop/scRNA/Data/SequencingBatches/",foldername,'Data',"/" ))
mult.mol.info = paste0('/home/sujwary/Desktop/scRNA/Data/SequencingBatches/',foldername,'Data/',filename_list)


s.out <- swappedDrops(mult.mol.info, min.frac=0.9)
length(s.out$cleaned)


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)

i = 39
sample_name = metaData$Sample[i]
print(sample_name)
idx = 13
filename_list[idx]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
count_matrix = s.out$cleaned[[idx]]
genes <-  rownames(count_matrix)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)
genes_df = as.data.frame(genes)
colnames(genes_df) = c('gene_id')


filename = '/home/sujwary/Desktop/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'

gtf <- rtracklayer::import(filename)
gtf_df=as.data.frame(gtf)
gene_list_orig = gtf_df$gene_name

cellranger = read.csv(filename, sep = '\t', header = TRUE)

filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)

# check if naming issue
 
sum(!(rownames(data_i_raw) %in% gene_list_orig))
rownames(data_i_raw) [!(rownames(data_i_raw) %in% gene_IDs$hgnc_symbol )]

#rownames(gene_IDs) = gene_IDs$ensembl_gene_id
#gene_IDs = gene_IDs[!duplicated(gene_IDs$ensembl_gene_id), ]


genes_df$HUGO <- gtf_df$gene_name[match(genes_df$gene_id, gtf_df$gene_id )]

dimnames(count_matrix) = list(genes_df$HUGO,colnames(count_matrix))


data_i_raw_Barswap_fixed = CreateSeuratObject(counts = count_matrix, project = "BM", min.cells = 3, min.features = 1)
nrow(data_i_raw_Barswap_fixed)
ncol(data_i_raw_Barswap_fixed)


#filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
#data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
#data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
#nrow(data_i_raw)
#ncol(data_i_raw)


colSum_list = colSums(data_i_raw_Barswap_fixed ) # Needs to be from Matrix library
keep = colSum_list >= 100
data_i_filtered = data_i_raw_Barswap_fixed[,keep]

data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")

data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
data_i_filtered_run = ScaleData(data_i_filtered_run)
data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
data_i_filtered_run = FindClusters(data_i_filtered_run)
data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)

folder = paste0("/home/sujwary/Desktop/scRNA/Data/SequencingBatches/",foldername,'Plots',"/",'NBM2CD45P','/' )
dir.create(folder)
pathName = paste0(folder,'Umap','','.png')
png(filename = pathName,width=2000, height=2000)
print(  DimPlot(data_i_filtered_run,pt.size = 3, reduction = "umap",label = FALSE))
print(plot)
dev.off()

folder_feature = paste0(folder,'Featureplots/' )
dir.create(folder_feature,recursive = T)
gene_list = c('CD3D', 'CD3G', 'CD3E', 'HBA2', 'HBB', 'IGKC', 'CD14', 'MS4A1', 'CDK6', 'FCGR3A')
for (j in 1:length(gene_list)){
  gene = gene_list[j]
  print(gene)
  #browser()
  plot = FeaturePlotFix(data_i_filtered_run, feature = gene, folder ='',
                        str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20)
  )
  
  file_str = ''
  pathName = paste0(folder_feature,gene,'','.png')
  png(filename = pathName,width=2000, height=2000)
  print(plot)
  dev.off()
  remove(plot)
  
}

##############################
data_matrix_filtered = data_i_filtered_run@assays[["RNA"]]@counts
data_matrix_raw = data_i_raw_Barswap_fixed@assays[["RNA"]]@counts
cluster_IDs = factor(as.character(Idents(data_i_filtered_run)))


sc = SoupChannel(data_matrix_raw, data_matrix_filtered)
sc = setClusters(sc, cluster_IDs)
sc = setDR(sc, data_i_filtered_run@reductions[["umap"]]@cell.embeddings)

igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
igGenes = igGenes[igGenes %in% sc[["toc"]]@Dimnames[[1]]] # usetoEst will break if not all genes are present
HBGenes = c('HBB','HBA2')
HBGenes = HBGenes[HBGenes %in% sc[["toc"]]@Dimnames[[1]]] 

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes))
print( plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst))



useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = HBGenes))
print( plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst))






print(  DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))



print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'percent.mt'))


print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap", features = 'nFeature_RNA'))


useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes, HB = HBGenes))
sc = calculateContaminationFraction(sc, list(IG = igGenes, HB = HBGenes), useToEst = useToEst)

out = adjustCounts(sc)
print('Plot')

pathName = paste0(folder,sample_name,'_PreSoup_Umap','_SoupAdjust','.png')
png(file=pathName,width=1000, height=1000)
print(plotChangeMap(sc, out, geneSet = HBGenes))
dev.off()


