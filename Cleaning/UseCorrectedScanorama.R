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

source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

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

i = 1

run = T
patient_list = c(12, 16, 20)
for (i in 1:length(patient_list) ){
  if (run){
    patient = patient_list[i]
    folder_Scanorama = paste0('/home/sujwary/Desktop/scRNA/Output/Scanorama/patient',patient,'/')
    gene_list = read.csv(file = paste0(folder_Scanorama,'genes','.csv'),header = F)
    gene_list = gene_list$V2
    samplename_baseline = metaData$Sample[metaData$`Patient Number` == patient & metaData$Treatment == 'baseline']
    data_matrix_baseline = read.csv(file = paste0(folder_Scanorama,samplename_baseline,'.csv'),header = F)
    data_matrix_baseline = t(data_matrix_baseline)
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',samplename_baseline,'/')
    colnames(data_matrix_baseline) = readLines(paste0(folder,samplename_baseline,'_colnames.txt'))
    
    rownames(data_matrix_baseline) = gene_list
    
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Scanorama/patient',patient,'/')
    samplename_C9D1 = metaData$Sample[metaData$`Patient Number` == patient & metaData$Treatment == 'C9D1']
    data_matrix_C9D1 = read.csv(file = paste0(folder,samplename_C9D1,'.csv'), header = F)
    data_matrix_C9D1 = t(data_matrix_C9D1)
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',samplename_C9D1,'/')
    colnames(data_matrix_C9D1) = readLines(paste0(folder,samplename_C9D1,'_colnames.txt'))
    rownames(data_matrix_C9D1) = gene_list
    
    data = cbind(data_matrix_baseline,data_matrix_C9D1)
    
    data = ScranNorm(data)
    
    # Reproduce 6 panel plots. Divide each cell by size factor. That gives norm count without log. Do this for 3404
    colnames(data) <- sub('[.]', '_', make.names(colnames(data), unique=TRUE))
    
    data_i = as.Seurat(data)
    data_i$sample = ''
    cellnames_baseline = colnames(data_matrix_baseline)
    cellnames_baseline =  gsub("-","_",cellnames_baseline,fixed = T)
    
    cellnames_C9D1 = colnames(data_matrix_C9D1)
    cellnames_C9D1 =  gsub("-","_",cellnames_C9D1,fixed = T)
    
    #data_i$sample[colnames(data_i) %in% cellnames_baseline] = samplename_baseline
    #data_i$sample[colnames(data_i) %in% cellnames_C9D1] = samplename_C9D1
    
    data_i$sample[1:length(colnames(data_matrix_baseline))] = samplename_baseline
    data_i$sample[(length(colnames(data_matrix_baseline)) + 1):ncol(data_i)] = samplename_C9D1
    
    data_i$nCount_RNA = Matrix:::colSums(data_i@assays[["RNA"]]@counts)
    data_i$percent.mt= PercentageFeatureSet(data_i, pattern = "^MT-")
    
    
    data_i_filtered_run = FindVariableFeatures(data_i, selection.method = "vst", nfeatures = 2000)
    data_i_filtered_run = ScaleData(data_i_filtered_run)
    data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
    data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
    data_i_filtered_run = FindClusters(data_i_filtered_run,resolution = 1.4)
    data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
    
    data_i_filtered_run$split_var = ''
    
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Scanorama/patient',patient,'/')
    path = paste0(folder,'/',sample_name,'.Robj')
    save(data_i_filtered_run,file= path)
  }else{
    folder = paste0('/home/sujwary/Desktop/scRNA/Output/Scanorama/patient',patient,'/')
    path = paste0(folder,'/','Patient',patient,'.Robj')
    data_i_filtered_run = loadRData(path)
  }
  cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
  groupBy_list = c('sample')
  splitBy_list = c('sample')
  plotAll(data_i_filtered_run, folder = folder_Scanorama,
          sample_name,sampleParam = NA,
          cell_features = cell_features,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = T, keepOldLabels = F, 
          groupBy = groupBy_list, splitBy = splitBy_list,
          PCA_dim = 30,resolution_val = 1)
  
  filepath_cluster = paste0( folder_Scanorama, 'Cluster/', 'PCA',30,'/res',1,'/' )
  PlotKnownMarkers(data_i_filtered_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                   plotType ='FeaturePlotFix' , str = '')
}
##############

# Baseline
folder_Scanorama = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/Baseline/')

if (run){
  patient = patient_list[i]
  gene_list = read.csv(file = paste0(folder_Scanorama,'genes','.csv'),header = F)
  gene_list = gene_list$V2
  samplename_baseline1 = metaData$Sample[metaData$`Patient Number` == 12 & metaData$Treatment == 'baseline']
  data_matrix_baseline1 = read.csv(file = paste0(folder_Scanorama,samplename_baseline1,'_corrected.csv'),header = F)
  data_matrix_baseline1 = t(data_matrix_baseline1)
  int1 =  read.csv(file = paste0(folder_Scanorama,samplename_baseline1,'_integrated.csv'),header = F)
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',samplename_baseline1,'/')
  colnames(data_matrix_baseline1) = readLines(paste0(folder,samplename_baseline1,'_colnames.txt'))
  rownames(data_matrix_baseline1) = gene_list
  rownames(int1)= readLines(paste0(folder,samplename_baseline1,'_colnames.txt'))
  
  
  samplename_baseline2 = metaData$Sample[metaData$`Patient Number` == 20 & metaData$Treatment == 'baseline']
  data_matrix_baseline2 = read.csv(file = paste0(folder_Scanorama,samplename_baseline2,'_corrected.csv'), header = F)
  data_matrix_baseline2 = t(data_matrix_baseline2)
  int2 =  read.csv(file = paste0(folder_Scanorama,samplename_baseline2,'_integrated.csv'),header = F)
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup_MT_C100/',samplename_baseline2,'/')
  colnames(data_matrix_baseline2) = readLines(paste0(folder,samplename_baseline2,'_colnames.txt'))
  rownames(data_matrix_baseline2) = gene_list
  rownames(int2) = readLines(paste0(folder,samplename_baseline2,'_colnames.txt'))
  
  int = cbind(t(int1),t(int2))
  rownames(int) <- paste0("PC_", 1:100)
  stdevs <- apply(int, MARGIN = 2, FUN = sd)
  
  data = cbind(data_matrix_baseline1,data_matrix_baseline2)
  data <- CreateSeuratObject(counts = data, assay = "RNA",  project = "BM",min.cells = 0, min.features = 0)
  
  rownames(data@meta.data) <- colnames(data)
  colnames(int) <- colnames(data)
  
  tmp = cbind(data_matrix_baseline1,data_matrix_baseline2)
  
  data@assays[["RNA"]]@scale.data = log(t((t(tmp) - colMedians(tmp))/colSds(tmp)) + 1)
       
  
  data[["pca"]] <- CreateDimReducObject(embeddings = t(int), stdev = stdevs, key = "PC_", assay = "RNA")
  
  #data_run = ScranNorm(data)
  #data_run = FindVariableFeatures(data_run, selection.method = "vst", nfeatures = 2000)
  #data_run = ScaleData(data_run)
  data_run = FindNeighbors(data, dims = 1:10) 
  data_run = FindClusters(data_run,resolution = 1.4)
  data_run = RunUMAP(data_run, dims = 1:30)
  data_run$split_var = ''
  
  FeaturePlot(data_run, feature = 'SELL')
  
  data_run[["percent.mt"]] <- PercentageFeatureSet(data_run, pattern = "^MT-")
  
  data_run$sample = ''
  cellnames_baseline1 = colnames(data_matrix_baseline1)
  #cellnames_baseline1 =  gsub("-","_",cellnames_baseline1,fixed = T)
  
  cellnames_baseline2 = colnames(data_matrix_baseline2)
  #cellnames_baseline2 =  gsub("-","_",cellnames_baseline2,fixed = T)
  
  data_run$sample[colnames(data_run) %in% cellnames_baseline1] = samplename_baseline1
  data_run$sample[colnames(data_run) %in% cellnames_baseline2] = samplename_baseline2
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Scanorama/Baseline','/')
  path = paste0(folder,'/','baseline','.Robj')
  save(data_run,file= path)
  
  ###################################
  data_i = CreateSeuratObject(counts = data, project = "BM", min.cells = 0, min.features = 0)
  
  
  data_i = ScranNorm(data_i)
  #data_i$active.ident = 'RNA'
  #data_i@active.ident =  factor(data_i$active.ident)
  #DefaultAssay(data_i) <- "RNA"
  
  # Remove cells with Nan values
  tmp = data_i[[data_i@active.assay]]@data
  tmp = as.matrix(tmp)
  cell_name = names(which(is.na(colSums(tmp))))
  cell_val = colnames(data_i[[data_i@active.assay]]@data) != cell_name
  colname_list = colnames(data_i[[data_i@active.assay]]@data)
  # Stupid fix
  for ( i in 1:length(colname_list)){
    if(colname_list[i] %in% cell_name) {
      cell_val[i] = F
    }
    
  }
  #data_integrate_test = SubsetData(data, cells = cell_val)
  if (length(cell_val) !=0){
    data_i@assays[["RNA"]]@data = data_i@assays[["RNA"]]@data[,cell_val]
    data_i@meta.data = data_i@meta.data[cell_val,]
    #data_i@assays[["RNA"]]@scale.data = data_i@assays[["RNA"]]@scale.data[,cell_val]
  }
  #browser()
  
  
  
  
  # Reproduce 6 panel plots. Divide each cell by size factor. That gives norm count without log. Do this for 3404
  
  
  data_i$sample = ''
  cellnames_baseline1 = colnames(data_matrix_baseline1)
  cellnames_baseline1 =  gsub("-","_",cellnames_baseline1,fixed = T)
  
  cellnames_baseline2 = colnames(data_matrix_baseline2)
  cellnames_baseline2 =  gsub("-","_",cellnames_baseline2,fixed = T)
  
  data_i$sample[colnames(data_i) %in% cellnames_baseline1] = samplename_baseline1
  data_i$sample[colnames(data_i) %in% cellnames_baseline2] = samplename_baseline2
  
  #data_i$sample[1:length(colnames(data_matrix_baseline1))] = samplename_baseline1
  #data_i$sample[(length(colnames(data_matrix_baseline1)) + 1):ncol(data_i)] = samplename_baseline2
  
  data_i$nCount_RNA = Matrix:::colSums(data_i@assays[["RNA"]]@counts)
  data_i$percent.mt= PercentageFeatureSet(data_i, pattern = "^MT-")
  
  
  data_i_filtered_run = FindVariableFeatures(data_i, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run,features = rownames(data_i_filtered_run))
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run,resolution = 1.4)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  data_i_filtered_run$split_var = ''
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Scanorama/Baseline','/')
  path = paste0(folder,'/','baseline','.Robj')
  save(data_i_filtered_run,file= path)
}else{
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/Scanorama/Baseline','/')
  path = paste0(folder,'/','baseline','.Robj')
  data_run = loadRData(path)
}

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
groupBy_list = c('sample')
splitBy_list = c('sample')
plotAll(data_run, folder = folder_Scanorama,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
        clusterTF =F, markersTF = T, keepOldLabels = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = 30,resolution_val = 1)

filepath_cluster = paste0( folder_Scanorama, 'Cluster/', 'PCA',30,'/res',1,'/' )
PlotKnownMarkers(data_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '')

