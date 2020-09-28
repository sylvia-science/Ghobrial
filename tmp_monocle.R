# This script is designed to pipeline Seurat object output into Monocle3
# This script also contains all the necessary functionalities in Monocle3 as 
# This could also convert a 2D Seurat object and visualize/analyze it in 3D
# It will start from reading in a Seurat object with the count sparse matrix, UMAP coordinates, and cluster information
# This script is originally written for local machines but adaptations have also been included in annotations


### Require:: 'BiocManager', 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'reticulate', 'htmlwidgets'

# This is a required python package
#reticulate::py_install("louvain")

# This is installing the actual monocle3
#devtools::install_github('cole-trapnell-lab/monocle3')


### Installing the packages

library(Seurat)
library(monocle)
library(htmlwidgets)

### Building the necessary parts for a basic cds

# part one, gene annotations

seurat = data_run_subset
seurat = seurat[,!(Idents(seurat) %in% c(0,8,11,21,22))]
gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix

New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

cell_metadata <- new("AnnotatedDataFrame", data = cell_metadata)
gene_annotation <- new("AnnotatedDataFrame", data = gene_annotation)

seurat$Ident = Idents(seurat)
cds_raw = as.CellDataSet(seurat)
cds_run <- reduceDimension(cds_raw, max_components = 2,
                            method = 'tSNE')
cds_run <- orderCells(cds_run)
#cds_run <- orderCells(cds_run, root_state = 5 )

cds_run_input = cds_run
cds_run_input@reducedDimW=t(seurat@reductions[["pca"]]@cell.embeddings)
cds_run_input@reducedDimS=t(seurat@reductions[["umap"]]@cell.embeddings)
cds_run_input@reducedDimK=t(seurat@reductions[["umap"]]@cell.embeddings)

path = paste0(filepath_cluster,'Monocle/cds_run_subset','.Robj')
save(cds_run,file= path)
cds_run_input = loadRData(path)


#cds_run_input <- orderCells(cds_run_input)


cds_run_input = cds_run[,cds_run$Ident %in% c(0,1,2,3,4,5,6,7,8,9)]


path = paste0(filepath_cluster,'Monocle/cell_traj_tSNE','.png')
png(file=path,width=1000, height=500,res = 100)
plot = plot_cell_trajectory(cds_run_input, color_by = 'Ident', 
                            show_tree = TRUE, show_backbone = TRUE, cell_size	= 1,
                            show_state_number = F) #+facet_wrap(~State, nrow = 1)
print(plot)
dev.off()



###
cds_run_input = cds_raw
cds_run_input <- estimateSizeFactors(cds_run_input)
cds_run_input <- estimateDispersions(cds_run_input)
cds_run_input <- detectGenes(cds_run_input, min_expr = 0.1)
#print(head(fData(cds_run_input))
     
expressed_genes <- row.names(subset(fData(cds_run_input),
                                    num_cells_expressed >= 10))
      
# Log-transform each value in the expression matrix.
L <- log(exprs(cds_run_input[expressed_genes,]))
      
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
      
# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
        stat_function(fun = dnorm, size = 0.5, color = 'red') +
        xlab("Standardized log(FPKM)") +
        ylab("Density")

cds_run_input <- reduceDimension(cds_run_input, max_components = 2,
                            method = 'tSNE')
cds_run_input <- orderCells(cds_run_input)
plot_cell_trajectory(cds_run_input, color_by = "State")

plot_cell_trajectory(cds_run_input, color_by = "State") +
  facet_wrap(~State, nrow = 1)

path = paste0(folder_subcluster,'Monocle/cds_run_monocle_pipeline','.Robj')
save(cds_run,file= path)
cds_run_input = loadRData(path)

