# Libraries

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)


source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')

folder_base <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/'

filename_sampleParam <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

output = select(metaData, 'Sample','Dexa or not','Patient Number', 'Treatment', '10X kit')
output[,"Num Cells"] = NA

output[,"nFeature_RNA Mean"] = NA
output[,"nFeature_RNA Median"] = NA
output[,"nFeature_RNA STD"] = NA

output[,"nCount_RNA Mean"] = NA
output[,"nCount_RNA Median"] = NA
output[,"nCount_RNA STD"] = NA

output[,"percent.mt Mean"] = NA
output[,"percent.mt Median"] = NA
output[,"percent.mt STD"] = NA

output_nCount_RNA = vector()
output_nFeature_RNA = vector()
output_percent.mt = vector()

first_v2 = 1
first_v3 = 17
sample_name <- metaData$Sample[first_v2]
folder = makeFolders(folder_base,sample_name,FALSE,FALSE)
data_combined_v2 <- loadRData(paste0(folder,'data.Robj'))

sample_name <- metaData$Sample[first_v3]
folder = makeFolders(folder_base,sample_name,FALSE,FALSE)
data_combined_v3 <- loadRData(paste0(folder,'data.Robj'))

for(i in 1:32){
  sample_name <- metaData$Sample[i]
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  # Open up each unfiltered sample, get mean, median, std for all 3 variables (9 columns)
  
  print(paste0('Sample: ', sample_name))
  
  
  filter = FALSE
  regress_TF = FALSE
  
  folder = makeFolders(folder_base,sample_name,filter,regress_TF)
  print(paste0('folder: ', folder))
  
  data <- loadRData(paste0(folder,'data.Robj'))
  
  
  output[i,"Num Cells"] = length(data@meta.data[["nCount_RNA"]])
  
  output[i,"nFeature_RNA Mean"] = mean(data@meta.data[["nFeature_RNA"]])
  output[i,"nFeature_RNA Median"] = median(data@meta.data[["nFeature_RNA"]])
  output[i,"nFeature_RNA STD"] = sd(data@meta.data[["nFeature_RNA"]])
  
  output[i,"nCount_RNA Mean"] = mean(data@meta.data[["nCount_RNA"]])
  output[i,"nCount_RNA Median"] = median(data@meta.data[["nCount_RNA"]])
  output[i,"nCount_RNA STD"] = sd(data@meta.data[["nCount_RNA"]])
  
  output[i,"percent.mt Mean"] = mean(data@meta.data[["percent.mt"]])
  output[i,"percent.mt Median"] = median(data@meta.data[["percent.mt"]])
  output[i,"percent.mt STD"] = sd(data@meta.data[["percent.mt"]])
  
  output_nFeature_RNA = c(output_nFeature_RNA, data@meta.data[["nFeature_RNA"]])
  output_nCount_RNA = c(output_nCount_RNA, data@meta.data[["nCount_RNA"]])
  output_percent.mt = c(output_percent.mt, data@meta.data[["percent.mt"]])
  if (metaData$'10X kit'[i] == 'v2'){
    data_combined_v2 = merge(data_combined_v2, data)
  }else{
    data_combined_v3 = merge(data_combined_v3, data)
    }
}

# Subset output 
pre_v2 <- subset(output, Treatment =='Pre-treatment' & output$'10X kit' == 'v2')
post_v2 <- subset(output, Treatment =='Post-treatment' & output$'10X kit' == 'v2')

pre_v3 <- subset(output, Treatment =='Pre-treatment' & output$'10X kit' == 'v3')
post_v3 <- subset(output, Treatment =='Post-treatment' & output$'10X kit' == 'v3')

subset_list = list(pre_v2, post_v2, pre_v3, post_v3)
subset_name_list = list('pre_v2', 'post_v2', 'pre_v3', 'post_v3')

# Sumurize output by making mean of all variables in each subset
output_summary = data.frame(matrix(ncol = 10, nrow = 4))
rownames(output_summary) = subset_name_list
for (i in 1:4){ # Subsets
  subset_i = subset_list[[i]]
  col_names = colnames(subset_i)
  for (j in 6:15){ # Columns
    print(subset_name_list[i])
    print(j)
    col_name = paste0('mean(',col_names[j], ')' ) 
    colnames(output_summary)[j - 5] = col_name
    mean_val = mean(subset_i[[col_names[j] ]])
    output_summary[i,j-5] = mean_val
    
    print(  col_name )
    print(mean_val) 
    
  }
  
}

write.csv(output, 
          file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/QC_analysis.csv',row.names=TRUE)

write.csv(output_summary, 
          file = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/QC_analysis_summary.csv',row.names=TRUE)

# Boxplots for V2 and V3

pathName <- paste0(folder_base,'violin_v2.png')
png(file=pathName,width=600, height=600)
VlnPlot(data_combined_v2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
dev.off()

pathName <- paste0(folder_base,'violin_v3.png')
png(file=pathName,width=600, height=600)
VlnPlot(data_combined_v3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
dev.off()

pathName <- paste0(folder_base,'Box_v2.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(c(pre_v2$'nFeature_RNA Mean', post_v2$'nFeature_RNA Mean')
        , names= "nFeature_RNA", ylim = c(0, 1200), main="nFeature_RNA")

boxplot(c(pre_v2$'nCount_RNA Mean', post_v2$'nCount_RNA Mean')
        , names=c("nCount_RNA"), ylim = c(0, 4000), main="nCount_RNA")

boxplot(c(pre_v2$'percent.mt Mean', post_v2$'percent.mt Mean')
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")

dev.off()

pathName <- paste0(folder_base,'Box_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(c(pre_v3$'nFeature_RNA Mean', post_v3$'nFeature_RNA Mean')
        , names= "nFeature_RNA", ylim = c(0, 1200), main="nFeature_RNA")

boxplot(c(pre_v3$'nCount_RNA Mean', post_v3$'nCount_RNA Mean')
        , names=c("nCount_RNA"), ylim = c(0, 4000), main="nCount_RNA")

boxplot(c(pre_v3$'percent.mt Mean', post_v3$'percent.mt Mean')
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")
dev.off()

