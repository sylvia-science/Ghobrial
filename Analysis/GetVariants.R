
library(Rsamtools)
library(Matrix)
library(bedr)
library(vcfR)
library(maftools)
library(readxl)
library(numbers)

source('~/Desktop/scRNA/Code/Integrate All/LoadHarmonyData.R')


base = '/disk2/Projects/EloRD/Data/Bam/PBMC/'
sample = 'GL1024BM'

# Seurat
filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1 || metaData$`Sample Type` == 'PBMC',]
#metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)



data = LoadHarmonyData (metaData, sampleParam_combine)
data_label = data[[1]]
sample_list = unique(data_label$sample )

sample_list = metaData$Sample

for (i in 1:length(sample_list)){
  sample = sample_list[i]
  print(sample)
  data_label = data[[1]]
  filepath_cluster = data[[2]]
  
  #data_label = data_label[,data_label$sample == sample]
  #barcode_harmony = colnames(data_label)
  #barcode_harmony <-gsub("_.*","",barcode_harmony)
  
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample,"_raw_feature_bc_matrix.h5",sep = "")
  data_sample = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  barcode = colnames(data_sample)
  #barcode = unique(c(barcode,barcode_harmony))
  print(filename)
  file = paste0(base, sample,'_out_cell_barcodes.csv')
  write.table(barcode,file, sep = ',', row.names = F, col.names = F, quote = F)

}

## Load files and save maf data
str = '_harmony'
str = ''
for (i in 1:length(sample_list)){
  #sample_name = sampleParam$Sample[i]
  sample = sample_list[i]
  print(sample)

  data_label = data[[1]]
  filepath_cluster = data[[2]]
  if ( !file.exists(paste0(base, sample, '_out_cell_barcodes', str , '.csv'))){
    #next
  }
  #data_label = data_label[,data_label$sample == sample]
  #DimPlot(data_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  
  # barcodes
 
  barcodes = read.csv(paste0(base, sample, '_out_cell_barcodes',str,'.csv'), header = F)
  barcodes = as.character(barcodes$V1)
  length(barcodes)
  
  barcode_filter = read.csv(paste0(base,sample,"_out_cell_barcodes_cellranger.csv"), header = F)
  barcode_filter = as.character(barcode_filter$V1)
  #barcodes = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  #barcode_filter = colnames(data_sample)
  
  ## MMTX
  
  ref = readMM(paste0(base, sample,'_demux_data',str,'/ref.mtx'))
  alt = readMM(paste0(base, sample,'_demux_data',str,'/alt.mtx'))
  dim(ref)
  dim(alt)
  colnames(ref) = barcodes
  colnames(alt) = barcodes
  
  ref = ref[,barcode_filter]
  alt = alt[,barcode_filter]
  barcodes = colnames(ref)
  ## VCF

  file_vcf = paste0(base, sample,
                    '_demux_data', str,'/',
                    'souporcell_merged_sorted_vcf.vcf')
  souporcell_merged_sorted = read.vcf(file_vcf)
  
  nrow(souporcell_merged_sorted[["vcf"]])

  #nrow(cluster_genotypes[["vcf"]])
  #nrow(cluster_genotypes_vep[["vcf"]])
  
  vcf = souporcell_merged_sorted[["vcf"]]
  
  
  
  maf = read.maf(maf = paste0(base, sample,'_demux_data', str,'/', sample,'.vep.maf'),
                 useAll = TRUE)
  
  
  cluster_list = levels(unique(Idents(data_label)))
  

  #vcf_subset = vcf[rowSums(alt)!= 0 || rowSums(ref)!= 0,]
  
  #maf_subset = maf@data[maf@data[["Start_Position"]] %in% vcf_subset$POS]
  #maf_silent_subset = maf@maf.silent[maf@maf.silent[["Start_Position"]] %in% vcf_subset$POS]
  
  colnames = c('Barcode','Sample','CellType','Alt','Ref',colnames(maf@data))
  maf_output = setNames(data.frame(matrix(ncol = length(colnames), nrow=0)), colnames)
  
  ident_list = Idents(data_label)
  names(ident_list) <-gsub("_.*","",names(ident_list) )
  
  cnt = 1
  for(bc in barcodes){
    
    if (mod(cnt,100) == 0){
      print(cnt)
      #print(nrow(maf_output))
    }
    cnt = cnt + 1
    cellType = as.character(ident_list[names(ident_list) == bc])
    cellType = ''
    alt_list = alt[,bc]
    ref_list = ref[,bc]
    pos_list = vcf$POS
    
    maf_subset = maf@data
    maf_subset = maf_subset[maf_subset$Start_Position %in% pos_list]
    maf_pos = maf_subset$Start_Position
    pos_in_maf = pos_list %in% maf_pos
    alt_list = alt_list[pos_in_maf]
    ref_list = ref_list[pos_in_maf]
    pos_list = pos_list[pos_in_maf]
    
    match_pos_in_maf = match(pos_list,maf_pos)
    alt_list = alt_list[match_pos_in_maf]
    ref_list = ref_list[match_pos_in_maf]
    pos_list = pos_list[match_pos_in_maf]
    
    maf_subset$Barcode = bc
    maf_subset$Alt = alt_list
    maf_subset$Ref = ref_list
    maf_subset$CellType = cellType
    
    maf_subset = maf_subset[maf_subset$Alt>0 | maf_subset$Ref>0,]
    maf_subset$Sample = sample

    maf_output = rbind(maf_output,maf_subset)
      
  }
  
  mask = as.logical(colSums(is.na(maf_output)) != nrow(maf_output))
  col_list = colnames(maf_output)
  col_list = col_list[mask]
  maf_output = data.frame(maf_output)
  
  maf_output = maf_output[, col_list]
  
  maf_output = maf_output[maf_output$Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site'),]
  
  maf_output
  file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
  write.csv(maf_output, file = file)

  
  #file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
  #write.csv(output_maf, file = file)
}

########################

sample = sample_list[1]
print(sample)
file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
maf_output_1 = read.csv( file = file)
mask = as.logical(colSums(is.na(maf_output_1)) != nrow(maf_output_1))
col_list = colnames(maf_output_1)
col_list = col_list[mask]
maf_output_1 = data.frame(maf_output_1)

maf_output_1 = maf_output_1[, col_list]


colnames = c('Barcode','Sample','CellType','Alt','Ref',colnames(maf_output_1),"CLIN_SIG", "PUBMED")
maf_all = setNames(data.frame(matrix(ncol = length(colnames), nrow = 0)), colnames)

for (i in 1:length(sample_list)){
  sample = sample_list[i]
  print(sample)
  file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
  if (file.exists(file)){
    
    maf_output = read.csv( file = file)
    
    mask = as.logical(colSums(is.na(maf_output)) != nrow(maf_output))
    col_list = colnames(maf_output)
    col_list = col_list[mask]
    maf_output = data.frame(maf_output)
    

    Missing1 <- setdiff(colnames(maf_all), colnames(maf_output))  # Find names of missing columns
    maf_output[Missing1] <- NA
    
    Missing2 <- setdiff(colnames(maf_output), colnames(maf_all))  # Find names of missing columns
    maf_all[,Missing2] <- NA
    
    if (nrow(maf_output) > 0){
      maf_output$Sample = sample
      maf_all = rbind(maf_all,maf_output)
    }
  }
}

maf_all = maf_all[!(maf_all$CellType %in% 'Remove'),]
maf_all$VAF = maf_all$Alt/(maf_all$Alt + maf_all$Ref)

#patient_list = unique(metaData$`Patient Number`[metaData$`Sample Type` == 'PBMC'])

#metaData_subset = metaData[metaData$`Patient Number` %in% patient_list,]

maf_all = maf_all[maf_all$VAF > 0.5,]
maf_all = merge(maf_all, metaData, by="Sample")
maf_all$CellType = as.character(maf_all$CellType)
maf_all$CellType[maf_all$'Sample Type' == 'PBMC'] = 'PBMC'

sample_list = unique(maf_all$Sample)
patient_list = unique(maf_all$`Patient Number`)

pos_list = unique(maf_all$Start_Position)
HGVSc_list = unique(maf_all$HGVSc)

#tmp = apply(maf_all, 2, function(x) length(unique(x)))


maf_list <- vector("list", length = length(patient_list))
cnt = 1

colnames = c('Patient','CellType','Effect','PatientCellTypeEffect')
mutation_df = setNames(data.frame(matrix(ncol = length(colnames), nrow = 0)), colnames)


for (patient in patient_list){
  print(patient)
  maf_patient = maf_all[ maf_all$`Patient Number` == patient,]
  cellType_list = unique(maf_patient$CellType)
  cellType_list = cellType_list[!(cellType_list %in% c('Remove'))]
  
  effect_list = unique(maf_patient$all_effects)
  maf_matrix = (data.frame(matrix(ncol = length(cellType_list), nrow = length(effect_list))))
  colnames(maf_matrix) = cellType_list
  rownames(maf_matrix) = effect_list
  for (cellType in cellType_list){
    maf_sample_cellType = maf_patient[maf_patient$CellType == cellType,]
    
    effect_cellType = maf_sample_cellType$all_effects
    maf_matrix[,cellType] = effect_list %in% effect_cellType
    
    colnames = c('Patient','CellType','Effect','PatientCellTypeEffect')
    mutation_df_patient = setNames(data.frame(matrix(ncol = length(colnames), nrow = length(unique(maf_sample_cellType$all_effects)))), colnames)
    mutation_df_patient$Patient = patient
    mutation_df_patient$CellType = cellType
    mutation_df_patient$Effect = unique(maf_sample_cellType$all_effects)
    
    mutation_df_patient$PatientCellTypeEffect = paste(mutation_df_patient$Patient,
                                                      mutation_df_patient$CellType, 
                                                      mutation_df_patient$Effect)
    
    mutation_df = rbind(mutation_df,mutation_df_patient)
  }
  
  effect_per_cell = rowSums(maf_matrix)
  effect_keep = effect_list[effect_per_cell< 4]
  maf_matrix_keep = maf_matrix[effect_per_cell< 4,]
  
  
  
  
  maf_list[[cnt]] = maf_matrix_keep
  cnt = cnt + 1
  
  file = paste0('/disk2/Projects/EloRD/Output/MafResults/Filter/Patient',patient,'.csv')
  write.csv(maf_matrix_keep, file = file)
}

nrow(mutation_df)
length(unique(mutation_df$PatientCellTypeEffect))

mutation_df$keep = F


for (i in 1:length(maf_list)){
  print(i)
  patient = patient_list[i]
  maf_matrix_keep = maf_list[[i]]
  
  for (cellType in colnames(maf_matrix_keep)){
    for (effect in rownames(maf_matrix_keep)){
      val = maf_matrix_keep[[cellType]][rownames(maf_matrix_keep) == effect]
  
      if (val){
        mutation_df$keep[mutation_df$Patient == patient & 
                           mutation_df$CellType == cellType & 
                           mutation_df$Effect == effect ] = T
      }
    }
  }
}
metaData$Patient = metaData$`Patient Number`
mutation_df = merge(mutation_df, metaData, by="Patient")

mutation_df = mutation_df[mutation_df$keep,]
mutation_df_SMM = mutation_df[mutation_df$`Diagnosis` == 'High Risk SMM',]
mutation_df_NBM = mutation_df[mutation_df$`Diagnosis` == 'NBM',]
  
variant_SMM = mutation_df_SMM[ !(mutation_df_SMM$Effect %in% mutation_df_NBM$Effect), ]

tmp = unique(variant_SMM$Effect)

