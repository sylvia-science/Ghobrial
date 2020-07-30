GetDEVoomLimma = function(fit,design,nclust){
  counter <- 1
  
  for (i in seq_len(nclust)) {
    for (j in seq_len(i - 1)) {
      browser()
      con <- integer(ncol(design))
      con[i] <- 1
      con[j] <- -1
      # Compute data for set of contrasts
      fit2 <- limma::contrasts.fit(fit, con)
      fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
      
      res <- limma::topTable(fit2, number = Inf, sort.by = "none")
      all.results[[counter]] <- res
      all.pairs[[counter]] <- colnames(design)[c(i, j)]
      counter <- counter + 1
      
      ## Also filling the reverse comparison.
      res$logFC <- -res$logFC
      all.results[[counter]] <- res
      all.pairs[[counter]] <- colnames(design)[c(j, i)]
      counter <- counter + 1
    }
  }
  
  ## Combine results across all pairwise tests
  all.pairs <- do.call(rbind, all.pairs)
  combined <- scran::combineMarkers(all.results, all.pairs, 
                                    pval.field = "P.Value",
                                    pval.type = "any")
  return(combined)
}