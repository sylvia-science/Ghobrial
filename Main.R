
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')

# Run with Filter
filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'

metaData <- read_excel(filename_metaData)
sampleParam <- read_excel(filename_sampleParam)


sample_name <- 'GL1024BM' # Pre treatment
sample_name <- 'GL1290BM' # Post treatment
folder <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'
filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")

filter <- FALSE
data = run_pipeline(filename,folder,sample_name,sampleParam,filter)
filter <- TRUE
data = run_pipeline(filename,folder,sample_name,sampleParam,filter)
print(get_cellType(data,folder,sample_name,filter))

