library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
library(sc)
library(scater)
library(dplyr)
library(scran)
library(stringr)
library(reshape2)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PipelineIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/PlotFunctionIntegrateAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/IntegrateAll_ClusterUmap.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')

source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Dmerge the samples, select highly variable genes, run PCA and then run it instead of findneighbors. esktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
sampleParam <- read_excel(filename_sampleParam)


downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  =downsample$x

i = 1

#folder = 'Intra-v3_1'
#folder = 'Inter-version'
folder = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")

sample_list = metaData$Sample
Umap_type = 'PCA_Umap'
Umap_type = 'SNN_Umap'
#folder_Scanorama = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/','Samples13/')

folder_Scanorama = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/',folder,'/',Umap_type,'/')

#folder_Scanorama = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/','Samples13_F2000/')

if (run){
  gene_list = read.csv(file = paste0(folder_Scanorama,'genes','.csv'),header = F)
  gene_list = gene_list$V2
  
  
  int = data.frame(matrix(ncol = 100, nrow = 0))
  data_list = vector(mode = "list",length = length(sample_list))
  

  for (i in 1:length(sample_list)){
    sample = sample_list[i]
    print(sample)
    #data_matrix_i_i = read.csv(file = paste0(folder_Scanorama,sample,'_corrected.csv'),header = F)
    #data_matrix_i_i = t(data_matrix_i_i)
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',sample,'/')
    #colnames(data_matrix_i_i) = readLines(paste0(folder,sample,'_colnames.txt'))
    #rownames(data_matrix_i_i) = gene_list
    
    #data_matrix = cbind(data_matrix,data_matrix_i_i)
    path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample,'/',sample,'.Robj')
    data_i = loadRData(path)
    
    #data_i = CreateSeuratObject(data_matrix_i_i, project = "BM", min.cells = 0)
    data_i$sample = sample
    
    
    path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample,'/cellIdents.csv')
    cellIdents1 = read.csv(path,sep = ',',row.names = 1)
    cellIdents1$x = paste0(cellIdents1$x, ' S',i)
    data_i$CellType = cellIdents1$x
    if (!is.na(downsample)){
      downsample = sub("_.*", "", downsample)
      cellnames = colnames(data_i)
      cellnames = sub("_.*", "", cellnames)
      data_i = data_i[,cellnames %in% downsample]
      #browser()
    }
    int_i =  read.csv(file = paste0(folder_Scanorama,sample,'_integrated.csv'),header = F)
    data_list[[i]] =data_i 
    int = rbind(int,int_i)
  }
  
  
  
  data_merge = merge(x =  data_list[[1]] ,y = data_list[2:length(data_list)], merge.data = T)

  data_matrix = matrix(data_merge@assays[["RNA"]]@data)
  
  colnames(int) <- paste0("PC_", 1:100)
  stdevs <- apply(int, MARGIN = 2, FUN = sd)
  

  rownames(int) = colnames(data_merge)
  
  data_merge$nCount_RNA = Matrix:::colSums(data_merge@assays[["RNA"]]@counts)
  #data_merge$percent.mt= PercentageFeatureSet(data_merge, pattern = "^MT-")
  
  
  #data_merge@assays[["RNA"]]@scale.data = (t((t(data_matrix) - colMedians(data_matrix))/colSds(data_matrix)) + 1)
  
  
  data_merge[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(int), stdev = stdevs, key = "PC_", assay = "RNA")
  
  
  
  data_run = ScaleData(data_merge)
  data_run = FindNeighbors(data_run, dims = 1:30) 
  data_run = FindClusters(data_run,resolution = 1.4)
  if (Umap_type == 'PCA_Umap'){
    data_run = RunUMAP(data_run, dims = 1:30)
  }else{
    data_run = RunUMAP(data_run,umap.method = 'umap-learn',graph= 'RNA_snn')
  }
  data_run = addMetaData(data_run, metaData)
  data_run = load_emptyDrops(data_run)
  data_run = load_Doublets(data_run)
  data_run = load_CellLabel(data_run)
  data_run$split_var = ''
  data_run$kit = data_run$`10X kit`
  
  
  folder_output = paste0(folder_Scanorama, '/')
  dir.create(folder_output)
  path = paste0(folder_output,'Data.Robj')
  save(data_run,file= path)
  
}else{

  folder_output = paste0(folder_Scanorama, '/')
  path = paste0(folder_output,'Data.Robj')
  data_run = loadRData(path)
}

folder_output = paste0(folder_Scanorama, '/')
data_run = addMetaData(data_run, metaData)
#data_run = load_emptyDrops(data_run)
#data_run = load_Doublets(data_run)
data_run = load_CellLabel(data_run)
data_run$split_var = ''
data_run$kit = data_run$`10X kit`


resolution_val = 1.4
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
groupBy_list = c('sample','Diagnosis','kit','Treatment','Batch','is_cell','LowCount','Doublet','GeneralCellType')
#groupBy_list = c('sample')
splitBy_list = NA
featurePlot_list = c('percentMT','nCount_RNA')
plotAll(data_run, folder = folder_output,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list, featurePlot_list = featurePlot_list,
        PCA_dim = 30,resolution_val = resolution_val)

filepath_cluster = paste0( folder_output, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',plotTogether = F)



fc_BCell = colorRampPalette(c("orange", "darkorange"))
fc_CD14 = colorRampPalette(c("purple", "purple4"))
fc_CD16 = colorRampPalette(c("lightgreen", "darkgreen"))
fc_CD8 = colorRampPalette(c("goldenrod1", "goldenrod4"))
fc_DC = colorRampPalette(c("yellow", "yellow4"))
fc_Eryth = colorRampPalette(c("red", "red4"))
fc_HSC = colorRampPalette(c("lightcyan", "darkcyan"))
fc_MK = 'darkseagreen'
fc_neutro = colorRampPalette(c("burlywood", "burlywood4"))
fc_NK = colorRampPalette(c("springgreen", "springgreen4"))
fc_pDC = colorRampPalette(c("lightpink", "pink4"))
fc_Plasma = 'plum'
fc_PreB = colorRampPalette(c("seashell", "seashell4"))
fc_T = colorRampPalette(c("lightblue", "darkblue"))




color_list = c(fc_BCell(9),fc_CD14(13),fc_CD16(5),
               fc_CD8(12),fc_DC(6),fc_Eryth(6),
               fc_HSC(10),fc_MK,fc_neutro(3),fc_NK(12),
               fc_pDC(2), fc_Plasma,fc_PreB(2),fc_T(13))


filepath_cluster = paste0( folder_Scanorama, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
cellList = sort(unique(data_run$CellType))
data_run$CellType = factor(data_run$CellType, levels = cellList)
group = 'CellType'
pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=5000, height=1500)
plot = DimPlot(data_run,pt.size = 1.5, reduction = "umap",label = FALSE,group.by  = group)

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
)



#plot = plot + scale_color_manual(values=color_list)

print(plot)

dev.off()

#########################
resolution_val = 1.4
filepath_cluster = paste0( folder_Scanorama, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )


k_num = 100
dim_num = 30
folder_output = paste0(filepath_cluster,'Entropy/')
dir.create(folder_output)
folder_output = paste0(folder_output,'Scanorama')

#data_run_downsample = data_run[,cellNames]
compute_entropy_Seurat(data_run, BBKNN = F,
                       corrected_assay = 'pca',folder = folder_output,
                       k_num = k_num, dim_num = dim_num)

file_entropy = paste0(folder_output,'_k',k_num,'_' ,'entropy.csv')
entropy = read.csv(file_entropy,sep = ',')

entropy$Method = 'Scanorama'

entropy$batch_entropy <- NULL
plotEntropy(entropy,folder_output)
   

#######################################################

# Baseline

patient1 = 20
patient2 = 30

folder_Scanorama = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/baseline/',
                'Patients',patient1,'_',patient2,'/')

if (run){
  gene_list = read.csv(file = paste0(folder_Scanorama,'genes','.csv'),header = F)
  gene_list = gene_list$V2
  samplename_baseline1 = metaData$Sample[metaData$`Patient Number` == patient1 & metaData$Treatment == 'baseline']
  data_matrix_i_i1 = read.csv(file = paste0(folder_Scanorama,samplename_baseline1,'_corrected.csv'),header = F)
  data_matrix_i_i1 = t(data_matrix_i_i1)
  int1 =  read.csv(file = paste0(folder_Scanorama,samplename_baseline1,'_integrated.csv'),header = F)
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',samplename_baseline1,'/')
  colnames(data_matrix_i_i1) = readLines(paste0(folder,samplename_baseline1,'_colnames.txt'))
  rownames(data_matrix_i_i1) = gene_list
  rownames(int1)= readLines(paste0(folder,samplename_baseline1,'_colnames.txt'))
  path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',samplename_baseline1,'/cellIdents.csv')
  cellIdents1 = read.csv(path,sep = ',',row.names = 1)
  cellIdents1$x = paste0(cellIdents1$x, ' S1')
  
  samplename_baseline2 = metaData$Sample[metaData$`Patient Number` == patient2 & metaData$Treatment == 'baseline']
  data_matrix_i_i2 = read.csv(file = paste0(folder_Scanorama,samplename_baseline2,'_corrected.csv'), header = F)
  data_matrix_i_i2 = t(data_matrix_i_i2)
  int2 =  read.csv(file = paste0(folder_Scanorama,samplename_baseline2,'_integrated.csv'),header = F)
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',samplename_baseline2,'/')
  colnames(data_matrix_i_i2) = readLines(paste0(folder,samplename_baseline2,'_colnames.txt'))
  rownames(data_matrix_i_i2) = gene_list
  rownames(int2) = readLines(paste0(folder,samplename_baseline2,'_colnames.txt'))
  path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',samplename_baseline2,'/cellIdents.csv')
  cellIdents2 = read.csv(path,sep = ',',row.names = 1)
  cellIdents2$x = paste0(cellIdents2$x, ' S2')
  
  int = cbind(t(int1),t(int2))
  rownames(int) <- paste0("PC_", 1:100)
  stdevs <- apply(int, MARGIN = 2, FUN = sd)
  
  data = cbind(data_matrix_i_i1,data_matrix_i_i2)
  data <- CreateSeuratObject(counts = data, assay = "RNA",  project = "BM",min.cells = 0, min.features = 0)
  
  rownames(data@meta.data) <- colnames(data)
  colnames(int) <- colnames(data)
  
  tmp = cbind(data_matrix_i_i1,data_matrix_i_i2)
  
  data@assays[["RNA"]]@scale.data = (t((t(tmp) - colMedians(tmp))/colSds(tmp)) + 1)
  
  
  data[["pca"]] <- CreateDimReducObject(embeddings = t(int), stdev = stdevs, key = "PC_", assay = "RNA")
  
  
  
  data$sample = ''
  data$CellType = ''
  cellnames_baseline1 = colnames(data_matrix_i_i1)
  #cellnames_baseline1 =  gsub("-","_",cellnames_baseline1,fixed = T)
  
  cellnames_baseline2 = colnames(data_matrix_i_i2)
  #cellnames_baseline2 =  gsub("-","_",cellnames_baseline2,fixed = T)
  
  data$sample[colnames(data) %in% cellnames_baseline1] = samplename_baseline1
  data$sample[colnames(data) %in% cellnames_baseline2] = samplename_baseline2
  
  data$CellType[colnames(data) %in% cellnames_baseline1] = cellIdents1$x
  data$CellType[colnames(data) %in% cellnames_baseline2] = cellIdents2$x
  
  #data_run = ScranNorm(data)
  #data_run = FindVariableFeatures(data_run, selection.method = "vst", nfeatures = 2000)
  #data_run = ScaleData(data_run)
  data_run = FindNeighbors(data, dims = 1:10) 
  data_run = FindClusters(data_run,resolution = 1.4)
  data_run = RunUMAP(data_run,umap.method = 'umap-learn',graph= 'RNA_snn')
  data_run$split_var = ''
  
  data_run[["percent.mt"]] <- PercentageFeatureSet(data_run, pattern = "^MT-")
  
  
  
  folder = paste0(folder_Scanorama,'SNN_Umap/')
  dir.create(folder, recursive = T)
  path = paste0(folder,'baseline','.Robj')
  save(data_run,file= path)
  

}else{
  folder = paste0(folder_Scanorama,'SNN_Umap/')
  path = paste0(folder,'/','baseline','.Robj')
  data_run = loadRData(path)
}
resolution_val = 1.4
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
groupBy_list = c('sample', 'CellType')
splitBy_list = c('sample')
plotAll(data_run, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, keepOldLabels = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = resolution_val)

filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '')

cellList = sort(unique(data_run$CellType))
data_run$CellType = factor(data_run$CellType, levels = cellList)
group = 'CellType'

pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=1000, height=1000)
plot = DimPlot(data_run,pt.size = 1.5, reduction = "umap",label = FALSE,group.by  = group)

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
  )


color_list = c("aquamarine", "aquamarine4", "orangered","orangered4","orchid1","orange1",
               "orange4", "snow4", "dodgerblue1","dodgerblue4","cyan","cyan4",
               "green","green4")



plot = plot + scale_color_manual(values=color_list)

print(plot)

dev.off()
