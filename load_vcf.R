
library(Matrix)
library(VariantAnnotation)

matrix_dir = "/disk2/Projects/EloRD/Data/Bam/GL1024BM_demux_data_test/"
barcode.path = paste0(matrix_dir, "barcodes.tsv.gz")
features.path = paste0(matrix_dir, "features.tsv.gz")
matrix.path = paste0(matrix_dir, "matrix.mtx.gz")
alt = paste0(matrix_dir, "alt.mtx")
ref = paste0(matrix_dir, "ref.mtx")
souporcell_merged_sorted = paste0(matrix_dir,'souporcell_merged_sorted_vcf.vcf.gz')
cluster_genotypes = paste0(matrix_dir, 'cluster_genotypes.vcf')

alt = readMM(file = alt)
ref = readMM(file = ref)

cluster_genotypes = VariantAnnotation::readVcf(cluster_genotypes)
souporcell_merged_sorted = VariantAnnotation::readVcf(souporcell_merged_sorted)

mat = readMM(file = matrix.path)
features.names = read.delim(features.path, header = FALSE, 
                            stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = F, stringsAsFactors =  F)
colnames(mat) = barcode.names$V1
rownames(mat)= feature.names$V1

