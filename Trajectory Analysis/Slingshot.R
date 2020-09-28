library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

data_run_subset

seurat = data_run_subset
seurat = seurat[,!(Idents(seurat) %in% c(0,8,11,17,18,21,22))]
start_cluster = 6
sds <- slingshot(Embeddings(seurat, "umap"), clusterLabels = Idents(seurat), 
                 start.clus = start_cluster, stretch = 2)

output_base = paste0(filepath_cluster,'Slingshot/','Start',start_cluster,'/')
dir.create(output_base)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(Idents(seurat), brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(Idents(seurat), hue_pal())




path = paste0(output_base,'/lineage','.png')
png(file=path,width=1000, height=500,res = 100)
plot = plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
plot = plot + lines(sds, lwd = 2, type = 'lineages', col = 'black')
print(plot)
dev.off()

path = paste0(output_base,'/lineage_curve','.png')
png(file=path,width=1000, height=500,res = 100)
plot = plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
plot = plot +lines(sds, lwd = 2, col = 'black')
print(plot)
dev.off()



nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
#par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  #plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  #lines(sds, lwd = 2, col = 'black', type = 'lineages')
  
  path = paste0(output_base,'/lineage_',i,'.png')
  png(file=path,width=1000, height=500,res = 100)
  plot = plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  plot = plot + lines(sds, lwd = 2, col = 'black', type = 'lineages')
  print(plot)
  dev.off()
  
}

library(tradeSeq)
library(uwot)

DE_input = seurat

DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'

DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
                           DE_input$Best_Overall_Response, ' ', Idents(DE_input))
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


BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 12 # use 2 cores

# fit negative binomial GAM
matrix = seurat@assays[["RNA"]]@counts
matrix = matrix[seurat@assays[["RNA"]]@var.features,]
sds <- fitGAM(seurat@assays[["RNA"]]@counts ,sds = sds,  U = design,
              parallel=TRUE, BPPARAM = BPPARAM)

# test for dynamic expression
ATres <- associationTest(sds)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sim$slingPseudotime_1, na.last = NA)
heatdata <- assays(sim)$counts[topgenes, pst.ord]
heatclus <- sim$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


##

# Get top highly variable genes

top_hvg <- HVFInfo(seurat) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(variance)) %>% 
  top_n(300, variance) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(seurat, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

colnames(dat_use_df) = make.names(dat_use_df)



dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)

model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)

val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)
metrics(data = val_results, truth, estimate)

var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[1:6]
# Convert to gene symbol
top_gene_name <- gns$gene_name[match(top_genes, gns$gene)]


par(mfrow = c(3, 2))
for (i in seq_along(top_genes)) {
  colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
  plot(reducedDim(sds), col = colors, 
       pch = 16, cex = 0.5, main = top_gene_name[i])
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}