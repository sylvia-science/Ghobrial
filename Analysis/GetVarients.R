
library(Rsamtools)
library(scruff)
library(Matrix)
library(bedr)
library(vcfR)
library(maftools)



base = '/disk2/Projects/EloRD/Data/Bam/'
sample = 'GL1024BM'

# barcodes
barcodes = read.csv(paste0(base, sample, '_out_cell_barcodes.csv'), header = F)
nrow(barcodes)

## MMTX

ref = readMM(paste0(base, sample,'_demux_data_test/ref.mtx'))
alt = readMM(paste0(base, sample,'_demux_data_test/alt.mtx'))
dim(ref)
dim(alt)

## VCF
souporcell_merged_sorted = read.vcf(paste0(base, sample,
                                           '_demux_data_test/',
                                           'souporcell_merged_sorted_vcf.vcf'))

nrow(souporcell_merged_sorted[["vcf"]])
#nrow(cluster_genotypes[["vcf"]])
#nrow(cluster_genotypes_vep[["vcf"]])

tmp = souporcell_merged_sorted[["vcf"]]


laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 

maf = read.maf(maf = paste0(base, sample,
                            '_demux_data_test/', sample,'.vep.maf'))

