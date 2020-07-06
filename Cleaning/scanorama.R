library(reticulate)
path_to_python <- "/home/sujwary/anaconda3/bin/python"
reticulate::use_python(path_to_python,required = TRUE)

system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sujwary/software/anaconda3/lib/')

scanorama <- import('scanorama')


# Integration.
integrated.data <- scanorama$integrate(datasets, genes_list)

# Batch correction.
corrected.data <- scanorama$correct(datasets, genes_list, return_dense=TRUE)

# Integration and batch correction.
integrated.corrected.data <- scanorama$correct(datasets, genes_list,
                                               return_dimred=TRUE, return_dense=TRUE)