BiocManager::install("Rsamtools")

library(Rsamtools)
library(scruff)
library(Matrix)
library(bedr)
library(vcfR)


base = '/disk2/Projects/EloRD/Data/Bam/'

#read in entire BAM file
bam <- scanBam("/disk2/Projects/EloRD/Data/Bam/GL1024BM.bam", 
               index = "/disk2/Projects/EloRD/Data/Bam/GL1024BM.bai")




.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam)
dim(bam_df)

# barcodes
barcodes = read.csv(paste0(base,'GL1024BM_out_cell_barcodes.csv'), header = F)
nrow(barcodes)

## MMTX

ref = readMM(paste0(base,'/GL1024BM_demux_data_test/ref.mtx'))
alt = readMM(paste0(base,'/GL1024BM_demux_data_test/alt.mtx'))
dim(ref)
dim(alt)

## VCF
souporcell_merged_sorted = read.vcf(paste0(base,
                                           '/GL1024BM_demux_data_test/',
                                           'souporcell_merged_sorted_vcf.vcf'))

cluster_genotypes = read.vcf(paste0(base,'/GL1024BM_demux_data_test/',
                'cluster_genotypes.vcf'))

cluster_genotypes_vep = read.vcf(paste0(base,'/GL1024BM_demux_data_test/',
                                    'cluster_genotypes.vep.vcf'))

nrow(souporcell_merged_sorted[["vcf"]])
nrow(cluster_genotypes[["vcf"]])
nrow(cluster_genotypes_vep[["vcf"]])

tmp = souporcell_merged_sorted[["vcf"]]

