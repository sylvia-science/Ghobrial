############
## EdgeR
############
runEdgeR = function(DE_input,design, contrast,keep,folder_output, subfolder){
  #browser()
  
  base = paste0(folder_output, 'DE/','EdgeR/',subfolder,'/')
  dir.create(base,recursive = T)
  
  dge <- DGEList(DE_input@assays[["RNA"]]@counts, 
                 group = DE_input$DE_ident) 
  countpergene_mean <- rowMeans(dge$counts)
  
  isexpr = countpergene_mean > 0.02

  pathName <- paste0(base,'CountPerGeneHist','.png')
  png(file=pathName,width=500, height=500,res = 100)
  hist(log10(countpergene_mean), breaks=100, main="", col="grey",xlab=expression('log10 Counts per gene'))
  abline(v=log10(0.02), col="blue", lwd=2, lty=2)
  dev.off()

  sum(keep)
  sum(isexpr)
  
  dge = dge[keep,]

  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = design,robust=F)
  # Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
  # Conduct genewise statistical tests for a given coefficient or contrast.
  fit <- glmQLFit(dge, design = design, contrast = contrast,robust = T) 
  # glmQLFTest is similar to glmLRT except that it replaces likelihood ratio tests 
  # with empirical Bayes quasi-likelihood F-tests
  qlf <- glmQLFTest(fit,contrast = contrast) # Need to have contrast here
  #qlf_batch <- glmQLFTest(fit, coef=3:4)
  
  col = colnames(design)
  col[con !=0]
  
  # res_batch <- topTags(qlf_batch, n = Inf)
  # res_batch = res_batch$table
  # sum(res_batch$FDR < 0.05)
  # sum(res_batch$FDR > 0.05)
  
  #browser()
  res <- topTags(qlf, n = Inf, adjust.method = "BH")
  res = res$table
  res = res[order(-res$logFC),]
  res[1:20,]
  sum(res$FDR < 0.05)
  sum(res$FDR > 0.05)
  
  
  geneid = mapIds(org.Hs.eg.db, rownames(qlf), 'ENTREZID', 'SYMBOL')
  
  go <- goana(qlf,geneid = geneid, species="Hs")
  topG0 = topGO(go, sort="up",number = Inf)
  keg <- kegga(qlf,geneid = geneid, species="Hs")
  topKEGG = topKEGG(keg, sort="up",number = Inf)
  
  
  
  
  path = paste0(base,'qlf','.Robj')
  save(qlf,file= path)
  qlf = loadRData(path)
  
  path = paste0(base, 'DE_EdgeR ',ident1,' Vs ', ident2,'.csv')
  print(path)
  write.csv(res, file = path,row.names=TRUE)
  
  path = paste0(base, 'GO_EdgeR ',ident1,' Vs ', ident2,'.csv')
  write.csv(topG0, file = path,row.names=TRUE)
  
  path = paste0(base, 'KEGG_EdgeR ',ident1,' Vs ', ident2,'.csv')
  write.csv(topKEGG, file = path,row.names=TRUE)
  
  
  output <- list(dge, fit,qlf, res)
  return (output)
}



plotEdgeR = function(dge_edgeR, DE_ident,qlf,base, subfolder){
  #browser()
  
  DE_ident = factor(DE_ident) 
  base = paste0(base, subfolder,'/')
  pathName <- paste0(base,'BCV','.png')
  png(file=pathName,width=500, height=500,res = 100)
  plotBCV(dge_edgeR)
  dev.off()
  
  pathName <- paste0(base,'QLDDisp','.png')
  png(file=pathName,width=500, height=500,res = 100)
  plotQLDisp(qlf)
  dev.off()
  
  pathName <- paste0(base,'MDS','.png')
  png(file=pathName,width=500, height=500,res = 100)
  plot = limma:::plotMDS(dge_edgeR,top = 500, 
                         col=as.numeric(DE_ident) , pch = 16)
  print(plot)
  dev.off()
  
  pathName <- paste0(base,'Smear','.png')
  png(file=pathName,width=500, height=500,res = 100)
  plotSmear(qlf)
  dev.off()
}

############
## MAST
############

runMAST= function(DE_input,formula,ident_list , ident, keep, folder_output,subfolder){
  library(rsvd)
  library(data.table)
  library(GGally)
  base = paste0(folder_output, 'DE/','MAST/',subfolder,'/')
  dir.create(base,recursive = T)
  
  #browser()
  dge <- DGEList(DE_input@assays[["RNA"]]@counts , 
                 group = DE_input$DE_ident)
  
  countpergene <- rowMeans(dge$counts)
  isexpr = countpergene > 0.02
  
  pathName <- paste0(base,'CountPerGeneHist','.png')
  png(file=pathName,width=500, height=500,res = 100)
  hist(log10(countpergene), breaks=100, main="", col="grey",xlab=expression('log10 Counts per gene'))
  abline(v=log10(0.02), col="blue", lwd=2, lty=2)
  dev.off()
  

  
  dge = dge[keep,]
  
  
  
  dge <- calcNormFactors(dge)
  cpms <- cpm(dge)
  sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                    cData = data.frame(wellKey = names(ident_list), 
                                       grp = ident_list))
  plotPCAMast(sca, base)
  
  zlmdata <- zlm(formula, sca)

  
  #colnames(zlmdata@LMlike@modelMatrix)  <- gsub("ident", "", colnames(zlmdata@LMlike@modelMatrix))
  #colnames(zlmdata@LMlike@modelMatrix)
  #browser()
  summaryCond = summary(zlmdata, doLRT=paste0('DE_ident',ident))
  summaryDt <- summaryCond$datatable
  
  fcHurdle = summaryDt[contrast== paste0('DE_ident',ident),]
  fcHurdle = fcHurdle[component=='H',]
  fcHurdle$FDR = p.adjust(fcHurdle$`Pr(>Chisq)`, method = "BH")
  fcHurdle = fcHurdle[order(fcHurdle$FDR),]

  

  
  path = paste0(base,'zlmdata','.Robj')
  save(zlmdata,file= path)
  zlmdata = loadRData(path)
  
  path = paste0(base,'summaryCond','.Robj')
  save(summaryCond,file= path)
  #summaryCond = loadRData(path)
  
  path = paste0(base,'fcHurdle','.csv')
  write.csv(fcHurdle, file = path,row.names=TRUE)
  
  output <- list(dge,summaryCond, fcHurdle)
  return(output)
  
}
plotPCAMast <- function(sca_obj, base){

  projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  
  pathName <- paste0(base,'PCA','.png')
  png(file=pathName,width=500, height=500,res = 100)
  print(ggpairs(pca, columns=c("PC1","PC2","PC3"),
                mapping=aes(color=grp), upper=list(continuous='blank')))
  invisible(pca)
  dev.off()
}


############
## VoomLimma
############

runVoomLimma= function(DE_input,design,contrast, keep, folder_output, subfolder){
  #browser()
  base = paste0(folder_output, 'DE/','VoomLimma/',subfolder,'/')
  dir.create(base,recursive = T)
  
  dge <- DGEList(DE_input@assays[["RNA"]]@counts, 
                 group = Idents(DE_input)) 
  
  
  countpergene <- rowMeans(dge$counts)
  isexpr = countpergene > 0.02
  
  pathName <- paste0(base,'CountPerGeneHist','.png')
  png(file=pathName,width=500, height=500,res = 100)
  hist(log10(countpergene), breaks=100, main="", col="grey",xlab=expression('log10 Counts per gene'))
  abline(v=log10(0.02), col="blue", lwd=2, lty=2)
  dev.off()
  
  dge = dge[keep,]
  
  
  dge <- calcNormFactors(dge)
  
  vm <- voom(dge, design = design, plot = TRUE)
  fit <- lmFit(vm, design = design)
  
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  
  geneid = mapIds(org.Hs.eg.db, rownames(fit2), 'ENTREZID', 'SYMBOL')
  
  go <- goana(fit2,geneid = geneid, species="Hs")
  topG0 = topGO(go, sort="up",number = Inf)
  keg <- kegga(fit2,geneid = geneid, species="Hs")
  topKEGG = topKEGG(keg, sort="up",number = Inf)
  
  res <- limma::topTable(fit2, number = Inf,adjust.method = "BH", sort.by = "logFC")
  
  res = res[order(res$logFC),]
  res[1:10,]
  sum(res$adj.P.Val < 0.05)
  sum(res$adj.P.Val > 0.05)
  
  dir.create( base, recursive = TRUE)
  path = paste0(base, 'DE_Voom ',ident1,' Vs ', ident2,'.csv')
  print(path)
  write.csv(res, file = path,row.names=TRUE)
  
  path = paste0(base, 'GO_EdgeR ',ident1,' Vs ', ident2,'.csv')
  write.csv(topG0, file = path,row.names=TRUE)
  
  path = paste0(base, 'KEGG_EdgeR ',ident1,' Vs ', ident2,'.csv')
  write.csv(topKEGG, file = path,row.names=TRUE)
  
  output <- list(dge,fit2, res)
  return(output)
}


############
## DESeq2
############
  
runDESeq2= function(DE_input,design,contrast, keep, folder_output, subfolder){
  #browser()
  base = paste0(folder_output, 'DE/','DESeq2/',subfolder,'/')
  dir.create(base,recursive = T)
  
  counts = DE_input@assays[["RNA"]]@counts
  
  dds <- DESeqDataSetFromMatrix(trunc(counts + 1), 
                                colData = DE_input@meta.data, 
                                design = design)
  
  countpergene <- rowMeans(counts)
  isexpr = countpergene > 0.02
  
  
  pathName <- paste0(base,'CountPerGeneHist','.png')
  png(file=pathName,width=500, height=500,res = 100)
  hist(log10(countpergene), breaks=100, main="", col="grey",xlab=expression('log10 Counts per gene'))
  abline(v=log10(0.02), col="blue", lwd=2, lty=2)
  dev.off()
  
  dds = dds[keep,]
  
  run_rlog = F
  if (run_rlog){
    rlog <- DESeq2:::rlog(dds, blind=TRUE)
    DESeq2::plotPCA(rlog, intgroup = "DE_ident")
    
    # Extract the rlog matrix from the object and compute pairwise correlation values
    rlog_mat <- assay(rlog)
    rlog_cor <- cor(rlog)
    
    # Plot heatmap
    pheatmap(rlog_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
  }
  dds <- DESeq(dds)
  plotDispEsts(dds)
  
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05,  pAdjustMethod	 ='BH')
  
  res <- lfcShrink(dds,type ='ashr',contrast =  contrast,res=res)
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction

  res_table_thres = res_tbl[res_tbl$padj < 0.05,]
  ## Volcano plot
  #browser()
  pathName <- paste0(base,'Volcano','.png')
  png(file=pathName,width=500, height=500,res = 100)
  ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = 0.05)) +
    ggtitle("Volcano plot") +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  dev.off()
  
  path = paste0(base, 'DE_DESeq2 ',ident1,' Vs ', ident2,'.csv')
  print(path)
  write.csv(res, file = path,row.names=TRUE)
  
  output <- list(dds, res_tbl)
  return(output)
}
