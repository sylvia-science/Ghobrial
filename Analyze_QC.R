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

output_mono = select(metaData, 'Sample','Dexa or not','Patient Number', 'Treatment', '10X kit')
output_mono[,"Num Cells"] = NA

output_mono[,"nFeature_RNA Mean"] = NA
output_mono[,"nFeature_RNA STD"] = NA

output_mono[,"nCount_RNA Mean"] = NA
output_mono[,"nCount_RNA STD"] = NA

output_mono[,"percent.mt Mean"] = NA
output_mono[,"percent.mt STD"] = NA

output_NOTmono = select(metaData, 'Sample','Dexa or not','Patient Number', 'Treatment', '10X kit')
output_NOTmono[,"Num Cells"] = NA

output_NOTmono[,"nFeature_RNA Mean"] = NA
output_NOTmono[,"nCount_RNA Mean"] = NA
output_NOTmono[,"percent.mt Mean"] = NA


output_nCount_RNA = vector()
output_nFeature_RNA = vector()
output_percent.mt = vector()

# first_v2 = 1
# first_v3 = 17
# sample_name <- metaData$Sample[first_v2]
# folder = makeFolders(folder_base,sample_name,FALSE,FALSE, makeFolder_TF = FALSE)
# data_combined_v2 <- loadRData(paste0(folder,'data.Robj'))
# 
# sample_name <- metaData$Sample[first_v3]
# folder = makeFolders(folder_base,sample_name,FALSE,FALSE, makeFolder_TF = FALSE)
# data_combined_v3 <- loadRData(paste0(folder,'data.Robj'))
# 
# sample_name <- metaData$Sample[first_v2]
# folder = makeFolders(folder_base,sample_name,TRUE,TRUE, makeFolder_TF = FALSE)
# data_combined_label_v2 <- loadRData(paste0(folder,'data.Robj'))
# 
# sample_name <- metaData$Sample[first_v3]
# folder = makeFolders(folder_base,sample_name,TRUE,TRUE, makeFolder_TF = FALSE)
# data_combined_label_v3 <- loadRData(paste0(folder,'data.Robj'))
# 
# 

data_list = vector(length = 0)

for(i in 1:24){
  print(i)
  sample_name <- metaData$Sample[i]
  filename <- paste("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  # Open up each unfiltered sample, get mean, median, std for all 3 variables (9 columns)
  
  print(paste0('Sample: ', sample_name))
  
  
  filter = FALSE
  regress_TF = FALSE
  
  folder = makeFolders(folder_base,sample_name,filter,regress_TF,makeFolder_TF = FALSE)
  print(paste0('folder: ', folder))
  
  data <- loadRData(paste0(folder,'data.Robj'))
  
  data@meta.data$orig.ident = paste0('data_',metaData$Treatment[i])
  data@meta.data$data_pre_renamed = data@active.ident
  
  data@meta.data$sample_name = sample_name
  data@meta.data$dexa = metaData$'Dexa or not'[i]
  data@meta.data$'Patient Number' = metaData$'Patient Number'[i]
  data@meta.data$Response = metaData$'Response'[i]
  data@meta.data$'10X Kit' = metaData$'10X Kit'[i]
  
  
  
  
  output_nFeature_RNA = c(output_nFeature_RNA, data@meta.data[["nFeature_RNA"]])
  output_nCount_RNA = c(output_nCount_RNA, data@meta.data[["nCount_RNA"]])
  output_percent.mt = c(output_percent.mt, data@meta.data[["percent.mt"]])

  
  ##
  ## Monocytes
  ##
  folder = makeFolders(folder_base,sample_name,filter = TRUE,regress_TF = TRUE,makeFolder_TF = FALSE)
  print(paste0('folder: ', folder))
  
  PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  
  data <- loadRData(paste0(folder,'data.Robj'))
  data = getCluster (data,resolution_val, PCA_dim)
  
  cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
  data = label_cells(data,cluster_IDs)
  

  data_mono = SubsetData(object= data,  cells = data@active.ident == "CD14+ Mono")
  data_NOTmono = SubsetData(object= data,  cells = data@active.ident != "CD14+ Mono")
  
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
  
  output_mono[i,"Num Cells"] = length(data_mono@meta.data[["nCount_RNA"]])
  
  output_mono[i,"nFeature_RNA Mean"] = mean(data_mono@meta.data[["nFeature_RNA"]])
  output_mono[i,"nFeature_RNA STD"] = sd(data_mono@meta.data[["nFeature_RNA"]])
  
  output_mono[i,"nCount_RNA Mean"] = mean(data_mono@meta.data[["nCount_RNA"]])
  output_mono[i,"nCount_RNA STD"] = sd(data_mono@meta.data[["nCount_RNA"]])
  
  output_mono[i,"percent.mt Mean"] = mean(data_mono@meta.data[["percent.mt"]])
  output_mono[i,"percent.mt STD"] = sd(data_mono@meta.data[["percent.mt"]])
  
  ##
  output_NOTmono[i,"Num Cells"] = length(data_NOTmono@meta.data[["nCount_RNA"]])
  output_NOTmono[i,"nFeature_RNA Mean"] = mean(data_NOTmono@meta.data[["nFeature_RNA"]])
  output_NOTmono[i,"nCount_RNA Mean"] = mean(data_NOTmono@meta.data[["nCount_RNA"]])
  output_NOTmono[i,"percent.mt Mean"] = mean(data_NOTmono@meta.data[["percent.mt"]])

  
  stats = clusterStats(data@active.ident)
  print(stats)
  stats_Mono = stats[stats$Cluster == 'CD14+ Mono',]
  stats_TCell = stats[stats$Cluster == 'T Cell',]
  stats_NK = stats[stats$Cluster == 'NK',]
  stats_CD16Mono = stats[stats$Cluster == 'CD16+ Mono',]
  stats_DC = stats[stats$Cluster == 'DC',]
  if (nrow(stats_Mono) == 1){
    output[i,"Percent Mono"] = stats_Mono$Percent
    output[i,"Num Mono"] = stats_Mono$Num
  }
  if (nrow(stats_TCell) == 1){
    output[i,"Percent T"] = stats_TCell$Percent
    output[i,"Num T"] = stats_TCell$Num
  }
  if (nrow(stats_NK) == 1){
    output[i,"Percent NK"] = stats_NK$Percent
    output[i,"Num NK"] = stats_NK$Num
  }
  if (nrow(stats_CD16Mono) == 1){
    output[i,"Percent CD16Mono"] = stats_CD16Mono$Percent
    output[i,"Num CD16Mono"] = stats_CD16Mono$Num
  }
  if (nrow(stats_DC) == 1){
    output[i,"Percent DC"] = stats_DC$Percent
    output[i,"Num DC"] = stats_DC$Num
  }

  #browser()
}

# Subset output 
# data_merge = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
# for (i in 3:length(data_list)){
#   print(i)
#   data_merge = merge(x =  data_merge,y = data_list[[i]],merge.data = FALSE)
# }
# 
# data_mono = SubsetData(object = data_merge, cells = data_merge@active.ident == "CD14+ Mono")

v2 <- subset(output, output$'10X kit' == 'v2')
v3 <- subset(output, output$'10X kit' == 'v3')

v2_mono <- subset(output_mono, output$'10X kit' == 'v2')
v3_mono <- subset(output_mono, output$'10X kit' == 'v3')

v2_NOTmono <- subset(output_NOTmono, output$'10X kit' == 'v2')
v3_NOTmono <- subset(output_NOTmono, output$'10X kit' == 'v3')



wilcox_nfeature = wilcox.test(x = v2$`nFeature_RNA Mean`,y =  v3$`nFeature_RNA Mean`, paired = FALSE)
wilcox_ncount = wilcox.test(x = v2$`nCount_RNA Mean`,y =  v3$`nCount_RNA Mean`, paired = FALSE)
wilcox_MT = wilcox.test(x = v2$`percent.mt Mean`,y =  v3$`percent.mt Mean`, paired = FALSE)

wilcox_mono_nfeature = wilcox.test(x = v2_mono$`nFeature_RNA Mean`,y =  v3_mono$`nFeature_RNA Mean`, paired = FALSE)
wilcox_mono_ncount = wilcox.test(x = v2_mono$`nCount_RNA Mean`,y =  v3_mono$`nCount_RNA Mean`, paired = FALSE)
wilcox_monox_MT = wilcox.test(x = v2_mono$`percent.mt Mean`,y =  v3_mono$`percent.mt Mean`, paired = FALSE)

wilcox_NOTmono_nfeature = wilcox.test(x = v2_NOTmono$`nFeature_RNA Mean`,y =  v3_NOTmono$`nFeature_RNA Mean`, paired = FALSE)
wilcox_NOTmono_ncount = wilcox.test(x = v2_NOTmono$`nCount_RNA Mean`,y =  v3_NOTmono$`nCount_RNA Mean`, paired = FALSE)
wilcox_NOTmonox_MT = wilcox.test(x = v2_NOTmono$`percent.mt Mean`,y =  v3_NOTmono$`percent.mt Mean`, paired = FALSE)


wilcox_mono = wilcox.test(x = v2$`Percent Mono`,y =  v3$`Percent Mono`, paired = FALSE)
wilcox_NumMono = wilcox.test(x = v2$`Num Mono`,y =  v3$`Num Mono`, paired = FALSE)
wilcox_TCell = wilcox.test(x = v2$`Percent T`,y =  v3$`Percent T`, paired = FALSE)
wilcox_NumTCell = wilcox.test(x = v2$`Num T`,y =  v3$`Num T`, paired = FALSE)
wilcox_NK = wilcox.test(x = v2$`Percent NK`,y =  v3$`Percent NK`, paired = FALSE)
wilcox_NumNK = wilcox.test(x = v2$`Num NK`,y =  v3$`Num NK`, paired = FALSE)

wilcox_CD16Mono = wilcox.test(x = v2$`Percent CD16Mono`,y =  v3$`Percent CD16Mono`, paired = FALSE)
wilcox_NumCD16Mono = wilcox.test(x = v2$`Num CD16Mono`,y =  v3$`Num CD16Mono`, paired = FALSE)
wilcox_DC = wilcox.test(x = v2$`Percent DC`,y =  v3$`Percent DC`, paired = FALSE)
wilcox_NumDC = wilcox.test(x = v2$`Num DC`,y =  v3$`Num DC`, paired = FALSE)

pre_v3 <- subset(output, output$Treatment =='pre' & output$'10X kit' == 'v3')
post_v3 <- subset(output, output$Treatment =='post' & output$'10X kit' == 'v3')

pre_v2 <- subset(output, output$Treatment =='pre' & output$'10X kit' == 'v2')
post_v2 <- subset(output, output$Treatment =='post' & output$'10X kit' == 'v2')

pre_v3 <- subset(output, output$Treatment =='pre' & output$'10X kit' == 'v3')
post_v3 <- subset(output, output$Treatment =='post' & output$'10X kit' == 'v3')

subset_list = list(pre_v2, post_v2, pre_v3, post_v3)
subset_name_list = list('pre_v2', 'post_v2', 'pre_v3', 'post_v3')

# Sumurize output by making mean of all variables in each subset
output_summary = data.frame(matrix(ncol = 10, nrow = 4))
rownames(output_summary) = subset_name_list
for (i in 1:1){ # Subsets
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
###########################################################
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
boxplot(v2$'nFeature_RNA Mean'
        , names= "nFeature_RNA", ylim = c(0, 1200), main="nFeature_RNA")

boxplot(v2$'nCount_RNA Mean'
        , names=c("nCount_RNA"), ylim = c(0, 4000), main="nCount_RNA")

boxplot(v2$'percent.mt Mean'
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")

dev.off()

pathName <- paste0(folder_base,'Box_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(v3$'nFeature_RNA Mean'
        , names= "nFeature_RNA", ylim = c(0, 1200), main="nFeature_RNA")

boxplot(v3$'nCount_RNA Mean'
        , names=c("nCount_RNA"), ylim = c(0, 4000), main="nCount_RNA")

boxplot(v3$'percent.mt Mean'
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")
dev.off()
#####################
## Monocyte QC
#####################
pathName <- paste0(folder_base,'Box_mono_QC_v2.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(v2_mono$'nFeature_RNA Mean'
        , names= "nFeature_RNA", ylim = c(0, 1500), main="nFeature_RNA")

boxplot(v2_mono$'nCount_RNA Mean'
        , names=c("nCount_RNA"), ylim = c(0, 5000), main="nCount_RNA")

boxplot(v2_mono$'percent.mt Mean'
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")

dev.off()

pathName <- paste0(folder_base,'Box_mono_QC_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(v3_mono$'nFeature_RNA Mean'
        , names= "nFeature_RNA", ylim = c(0, 1500), main="nFeature_RNA")

boxplot(v3_mono$'nCount_RNA Mean'
        , names=c("nCount_RNA"), ylim = c(0, 5000), main="nCount_RNA")

boxplot(v3_mono$'percent.mt Mean'
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")
dev.off()


#####################
## NOT Monocyte QC
#####################
pathName <- paste0(folder_base,'Box_Notmono_QC_v2.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(v2_NOTmono$'nFeature_RNA Mean'
        , names= "nFeature_RNA", ylim = c(0, 1500), main="nFeature_RNA")

boxplot(v2_NOTmono$'nCount_RNA Mean'
        , names=c("nCount_RNA"), ylim = c(0, 5000), main="nCount_RNA")

boxplot(v2_NOTmono$'percent.mt Mean'
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")

dev.off()

pathName <- paste0(folder_base,'Box_Notmono_QC_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,3))
boxplot(v3_NOTmono$'nFeature_RNA Mean'
        , names= "nFeature_RNA", ylim = c(0, 1500), main="nFeature_RNA")

boxplot(v3_NOTmono$'nCount_RNA Mean'
        , names=c("nCount_RNA"), ylim = c(0, 5000), main="nCount_RNA")

boxplot(v3_NOTmono$'percent.mt Mean'
        , names= "percent.mt Mean", ylim = c(0, 40), main="percent.mt")
dev.off()


################################

pathName <- paste0(folder_base,'Box_PercentMono_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Percent Mono'
        , names= "% Monocytes V2", ylim = c(0, 100), main="% Monocytes V2")

boxplot(v3$'Percent Mono'
        , names= "% Monocytes V3", ylim = c(0, 100), main="% Monocytes V3")
dev.off()

pathName <- paste0(folder_base,'Box_NumMono_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Num Mono'
        , names= "Number Monocytes V2", ylim = c(0, 400), main="Number Monocytes V2")

boxplot(v3$'Num Mono'
        , names= "Number Monocytes V3", ylim = c(0, 400), main="Number Monocytes V3")
dev.off()


############################
pathName <- paste0(folder_base,'Box_PercentT_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Percent T'
        , names= "% T Cells V2", ylim = c(0, 100), main="% T Cells V2")

boxplot(v3$'Percent T'
        , names= "% T Cells V3", ylim = c(0, 100), main="% T Cells V3")
dev.off()

pathName <- paste0(folder_base,'Box_NumT_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Num T'
        , names= "Number T Cells V2", ylim = c(0, 1500), main="Number T Cells V2")

boxplot(v3$'Num T'
        , names= "Number T Cells V3", ylim = c(0, 1500), main="Number T Cells V3")
dev.off()

###########################
pathName <- paste0(folder_base,'Box_PercentNK_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Percent NK'
        , names= "% NK Cells V2", ylim = c(0, 100), main="% NK Cells V2")

boxplot(v3$'Percent NK'
        , names= "% NK Cells V3", ylim = c(0, 100), main="% NK Cells V3")
dev.off()

pathName <- paste0(folder_base,'Box_NumNK_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Num NK'
        , names= "Number NK Cells V2", ylim = c(0, 400), main="Number NK Cells V2")

boxplot(v3$'Num NK'
        , names= "Number NK Cells V3", ylim = c(0, 400), main="Number NK Cells V3")
dev.off()

#############################

pathName <- paste0(folder_base,'Box_PercentDC_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Percent DC'
        , names= "% DC V2", ylim = c(0, 20), main="% DC V2")

boxplot(v3$'Percent DC'
        , names= "% DC V3", ylim = c(0, 20), main="% DC V3")
dev.off()

pathName <- paste0(folder_base,'Box_NumDC_v2_v3.png')
png(file=pathName,width=600, height=600)
par(mfrow=c(1,2))
boxplot(v2$'Num DC'
        , names= "Number DC V2", ylim = c(0, 400), main="Number DC V2")

boxplot(v3$'Num DC'
        , names= "Number DC V3", ylim = c(0, 400), main="Number DC V3")
dev.off()

