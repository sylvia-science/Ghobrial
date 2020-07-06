library(Seurat)
library(SoupX)
library(ggplot2)
library(readxl)
library(DropletUtils)

threshold =0.3

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sample_name = metaData$Sample[4]
path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_scrublet',threshold,'.csv')
scrb = read.csv(path,header = T)


filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
data_i_raw = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3,min.features = 1)


# Load scrublet data
path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_scrublet',threshold,'.Robj')
load(path)
get(ls()[ls() != "fileName"])

DimPlot(data_i_filtered_run)


row_val = rownames(data_i_filtered_run)
col_val = colnames(data_i_filtered_run)

##
br.out = barcodeRanks(data_i_raw@assays[["RNA"]]@counts)
e.out = emptyDrops(data_i_raw@assays[["RNA"]]@counts)


br = data.frame(br.out@rownames,br.out$rank, br.out$total)
e = data.frame(e.out@rownames, e.out$FDR,e.out$LogProb)
colnames(br) = c('rownames','rank','total')
colnames(e) = c('rownames','FDR','LogProb')


e_sorted = e[match(e$rownames, col_val),]

br_e <- merge(br,e,by="rownames")
br_e$is_cell <- br_e$FDR <= 0.01

ggplot(br_e) + geom_point(aes(rank,total,color=is_cell)) + scale_y_log10() + scale_x_log10()



plot(log10(e.out$Total), -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")

br_e_sorted = br_e[match(as.character(rownames(data_i_filtered_run@meta.data)), as.character(br_e$rownames)   ),]


data_i_filtered_run$emptyProb = br_e_sorted$LogProb
data_i_filtered_run$is_cell = br_e_sorted$is_cell
data_i_filtered_run[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered_run, pattern = "^MT-")

markers1 = FindMarkers(data_i_filtered_run, ident.1 = 1)
  
DimPlot(data_i_filtered_run,pt.size = 0.8, reduction = "umap", label = T)

DimPlot(data_i_filtered_run,pt.size = 0.8, reduction = "umap", group.by = 'is_cell')

FeaturePlot(data_i_filtered_run,pt.size = 0.8, reduction = "umap", features = 'percent.mt')

FeaturePlot(data_i_filtered_run,pt.size = 0.8, reduction = "umap", features = 'LogProb')
FeaturePlot(data_i_filtered_run,pt.size = 0.8, reduction = "umap", features = 'HBA2')


data_i_filtered_run = NormalizeData(data_i_filtered_run, normalization.method = "LogNormalize", scale.factor = 10000)
data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
data_i_filtered_run = ScaleData(data_i_filtered_run)
data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
data_i_filtered_run = FindClusters(data_i_filtered_run)
data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)

data_matrix_filtered = data_i_filtered_run@assays[["RNA"]]@counts
data_matrix_raw = data_i_raw@assays[["RNA"]]@counts
cluster_IDs = factor(as.character(Idents(data_i_filtered_run)))

sc = SoupChannel(data_matrix_raw, data_matrix_filtered)
sc = setClusters(sc, cluster_IDs)
sc = setDR(sc, data_i_filtered_run@reductions[["umap"]]@cell.embeddings)

gene_list = c('IGKC','IGLC1','HBB','HBA2')
gene_list = gene_list[gene_list %in% sc[["toc"]]@Dimnames[[1]]]
for (gene in gene_list){
  print(gene)
  #png(filename = )
  print(plotMarkerMap(sc, geneSet = gene))
}

igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
igGenes = igGenes[igGenes %in% sc[["toc"]]@Dimnames[[1]]] # usetoEst will break if not all genes are present
HBGenes = c('HBB','HBA2')
HBGenes = HBGenes[HBGenes %in% sc[["toc"]]@Dimnames[[1]]] 

pdf()

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes))
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = HBGenes))
plotMarkerMap(sc,geneSet=HBGenes,useToEst=useToEst)

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes, HB = HBGenes))


sc = calculateContaminationFraction(sc, list(IG = igGenes, HB = HBGenes), useToEst = useToEst)
#head(sc$metaData)
#quantile(sc$metaData$rho)

out = adjustCounts(sc)
plotChangeMap(sc, out, geneSet = HBGenes)

FeaturePlot(data_i_filtered_run,pt.size = 0.8, reduction = "umap", features = 'HBA2')
# Put soup data back into filtered run
data_i_filtered_run_soup = data_i_filtered_run

data_i_filtered_run_soup@assays[["RNA"]]@counts = out

path = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'_Soup_adjusted.Robj')
save(data_i_filtered_run_soup,file= path)


data_i_filtered_run_soup = NormalizeData(data_i_filtered_run_soup, normalization.method = "LogNormalize", scale.factor = 10000)
data_i_filtered_run_soup = FindVariableFeatures(data_i_filtered_run_soup, selection.method = "vst", nfeatures = 2000)
data_i_filtered_run_soup = ScaleData(data_i_filtered_run_soup)
data_i_filtered_run_soup = RunPCA(data_i_filtered_run_soup,npcs = 30)
data_i_filtered_run_soup = FindNeighbors(data_i_filtered_run_soup, dims = 1:30)
data_i_filtered_run_soup = FindClusters(data_i_filtered_run_soup)
data_i_filtered_run_soup = RunUMAP(data_i_filtered_run_soup, dims = 1:30)

DimPlot(data_i_filtered_run_soup,pt.size = 0.5, reduction = "umap",label = FALSE)

