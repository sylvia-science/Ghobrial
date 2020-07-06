library(Seurat)
library(ggplot2)
library(sctransform)
pbmc <- Read10X(data.dir = "/home/sujwary/Projects/PBMC/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc)

#pbmc = readRDS('/home/sujwary/Projects/PBMC/pbmc33k_identities.Rds', refhook = NULL)
#pbmc =ReadH5AD('/home/sujwary/Projects/PBMC/svensson_chromium_control.h5ad', layers = "data", verbose = TRUE)  
#pbmc = CreateSeuratObject(counts = pbmc@assays[["RNA"]]@counts, project = "BM", min.cells = 3, min.features = 1)

# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

####################
# run sctransform
pbmc_sctransform <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc_sctransform <- RunPCA(pbmc_sctransform, verbose = FALSE)
pbmc_sctransform <- RunUMAP(pbmc_sctransform, dims = 1:30, verbose = FALSE)

pbmc_sctransform <- FindNeighbors(pbmc_sctransform, dims = 1:30, verbose = FALSE)
pbmc_sctransform <- FindClusters(pbmc_sctransform, verbose = FALSE)
DimPlot(pbmc_sctransform, label = TRUE) + NoLegend()



rowMeans_list = rowMeans(pbmc_sctransform@assays[["SCT"]]@scale.data)
bin_list = as.numeric(cut_number(rowMeans_list,6))

colSums_list = colSums(pbmc_sctransform@assays[["RNA"]]@counts)

plot_list = vector('list', 6)
for (bin in 1:6){
  print(bin)
  y_list = pbmc_sctransform@assays[["SCT"]]@data[bin_list == bin,]
  df = as.data.frame(matrix(ncol = 2, nrow = length(colSums_list)))
  colnames(df) = c('colsum','exp')
  df$colsum = colSums_list # Unormalized Mean per cell
  df$exp = colMeans(y_list) # Norm Mean per cell
  
  
  plot = ggplot(df, aes(colsum, exp)) +
    geom_point(color="black") + 
    stat_smooth(aes(x=colsum,y=exp), method="loess", se=F, color="tomato2") +
    theme(text = element_text(size=20))
  
  plot_list[[bin]]= plot
  
  # plot mean counts on y
  
  
}


plot = plot_grid(
  plot_list[[1]], plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
  labels = "AUTO", ncol = 2)
print(plot)




################


pbmc_SeuratNorm = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_SeuratNorm = FindVariableFeatures(pbmc_SeuratNorm, selection.method = "vst", nfeatures = 2000)
pbmc_SeuratNorm = ScaleData(pbmc_SeuratNorm)
pbmc_SeuratNorm = RunPCA(pbmc_SeuratNorm,npcs = 30)
pbmc_SeuratNorm = FindNeighbors(pbmc_SeuratNorm, dims = 1:30)
pbmc_SeuratNorm = FindClusters(pbmc_SeuratNorm)
pbmc_SeuratNorm = RunUMAP(pbmc_SeuratNorm, dims = 1:30)

DimPlot(pbmc_SeuratNorm, label = TRUE) + NoLegend()



rowMeans_list = rowMeans(pbmc_SeuratNorm@assays[["RNA"]]@data)
bin_list = as.numeric(cut_number(rowMeans_list,6))

colSums_list = colSums(pbmc_SeuratNorm@assays[["RNA"]]@counts)

plot_list = vector('list', 6)
for (bin in 1:6){
  print(bin)
  y_list = pbmc_SeuratNorm@assays[["RNA"]]@data[bin_list == bin,]
  df = as.data.frame(matrix(ncol = 2, nrow = length(colSums_list)))
  colnames(df) = c('colsum','exp')
  df$colsum = colSums_list # Unormalized Mean per cell
  df$exp = colMeans(y_list) # Norm Mean per cell
  
  
  plot = ggplot(df, aes(colsum, exp)) +
    geom_point(color="black") + 
    stat_smooth(aes(x=colsum,y=exp), method="loess", se=F, color="tomato2") +
    theme(text = element_text(size=20))
  
  plot_list[[bin]]= plot
  
  # plot mean counts on y
  
  
}


plot = plot_grid(
  plot_list[[1]], plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
  labels = "AUTO", ncol = 2)
print(plot)
