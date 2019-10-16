
# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)

#library(biomaRt)

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')



sample_name <- 'GL1080BM' # 10 Pre treatment
sample_name <- 'GL1374BM' # 10 Post treatment

sample_name <- 'GL1024BM' # 5 Pre treatment
#sample_name <- 'GL1290BM' # 5 Post treatment

#sample_name <- metaData$Sample[6]

folder <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'
filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")

filter <- FALSE
data = run_pipeline(filename,folder,sample_name,sampleParam,filter)
                    
filter <- TRUE
data = run_pipeline(filename,folder,sample_name,sampleParam,filter)
print(get_cellType(data,folder,sample_name,filter))
label_cells(data,folder,sample_name,sampleParam,filter)



