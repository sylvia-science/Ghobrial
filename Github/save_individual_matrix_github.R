library(Matrix)
library(readxl)
library(Seurat)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sample_name = metaData$Sample[4]
filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")

path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_scrublet0.3.Robj')


data_i_scrublet = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
data_i_scrublet = CreateSeuratObject(counts = data_i_scrublet, project = "BM", min.cells = 3,min.features = 1)
colSum_list = colSums(data_i_scrublet)
keep = colSum_list >= 100

data_i_scrublet = data_i_scrublet[,keep]

data_matrix = data_i_scrublet@assays[["RNA"]]@counts

write(colnames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/RawMatrix/',sample_name,'_colnames.txt'))
write(rownames(data_matrix), file = paste0('/home/sujwary/Desktop/scRNA/Data/RawMatrix/',sample_name,'_rownames.txt'))
writeMM(data_matrix, file = paste0('/home/sujwary/Desktop/scRNA/Data/RawMatrix/',sample_name,'_matrix.txt'))



