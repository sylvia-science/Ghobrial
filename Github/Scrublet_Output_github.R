library(Matrix)
threshold =0.3
sample_name = metaData$Sample[4]
path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_scrublet',threshold,'.csv')
scrb = read.csv(path,header = T)


filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
data_i_filtered = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
data_i_filtered = CreateSeuratObject(counts = data_i_filtered, project = "BM", min.cells = 3, min.features = 1)

colSum_list = colSums(data_i_filtered ) # Needs to be from Matrix library
keep = colSum_list >= 100
data_i_filtered = data_i_filtered[,keep]



data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
data_i_filtered_run = ScaleData(data_i_filtered_run)
data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
data_i_filtered_run = FindClusters(data_i_filtered_run)
data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)

data_i_filtered_run$Scrublet_Boolean = as.logical(as.character(toupper(scrb$Scrublet_Boolean)))
data_i_filtered_run$doublet_scores = scrb$doublet_scores
data_i_filtered_run$predicted_doublet = as.logical(as.character(toupper(scrb$predicted_doublet)))


pathName = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_Scrublet_Boolean','.png')
png(file=pathName,width=1000, height=1000)
print(DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))
dev.off()


pathName = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_Scrublet_Boolean','.png')
png(file=pathName,width=1000, height=1000)
print(DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE,group.by  = 'Scrublet_Boolean'))
dev.off()

pathName = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'doublet_scores','.png')
png(file=pathName,width=1000, height=1000)
print(FeaturePlot(data_i_filtered_run,pt.size = 0.5, features = 'doublet_scores'))
dev.off()

pathName = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'predicted_doublet','.png')
png(file=pathName,width=1000, height=1000)
print(DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE,group.by  = 'predicted_doublet'))
dev.off()


data_i_filtered_run = data_i_filtered_run[,!(data_i_filtered_run$Scrublet_Boolean)]

path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_scrublet',threshold,'.Robj')
save(data_i_filtered_run,file= path)

