

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
library(readxl)

filename_upreg <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa Genes/UpReg.xlsx'
filename_downreg <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa Genes/DownReg.xlsx'

genes_upreg <- read_excel(filename_upreg)
genes_downreg <- read_excel(filename_downreg)

filename_sample_Integrate_pairs <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Integrate_pairs.xlsx'
sample_Integrate_pairs <- read_excel(filename_sample_Integrate_pairs)

filename_metaData <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData <- read_excel(filename_metaData)

Dexa_MinResponse = metaData[metaData$Treatment == 'Pre-treatment' 
                                    & metaData$`Dexa or not` == 'Yes' 
                                    & metaData$Response == 	'Minimal Response'
                                    & metaData$'Sample Type' == 'Bone Marrow',]

Dexa_MRTP = metaData[metaData$Treatment == 'Pre-treatment' 
                            & metaData$`Dexa or not` == 'Yes' 
                            & metaData$Response == 	'Minimal Response then Progression'
                            & metaData$'Sample Type' == 'Bone Marrow',]

Dexa_VGPR = metaData[metaData$Treatment == 'Pre-treatment' 
                     & metaData$`Dexa or not` == 'Yes' 
                     & metaData$Response == 	'VGPR'
                     & metaData$'Sample Type' == 'Bone Marrow',]

NoDexa_MinResponse = metaData[metaData$Treatment == 'Pre-treatment' 
                            & metaData$`Dexa or not` == 'No' 
                            & metaData$Response == 	'Minimal Response'
                            & metaData$'Sample Type' == 'Bone Marrow',]

NoDexa_MRTP = metaData[metaData$Treatment == 'Pre-treatment' 
                     & metaData$`Dexa or not` == 'No' 
                     & metaData$Response == 	'Minimal Response then Progression'
                     & metaData$'Sample Type' == 'Bone Marrow',]

NoDexa_VGPR = metaData[metaData$Treatment == 'Pre-treatment' 
                     & metaData$`Dexa or not` == 'No' 
                     & metaData$Response == 	'VGPR'
                     & metaData$'Sample Type' == 'Bone Marrow',]

folder_base_output <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/'

sample_type = 'BM'
#sample_type = 'PB'
patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all


for(patient in patient_list){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  
  file = paste0(folder_pre,'diffExprGene_Patient',patient,'.csv')
  print(file)
  diffExprGene <- read.csv(file)
  combine_upreg <- merge(diffExprGene,genes_upreg,by="Gene")
  combine_downreg <- merge(diffExprGene,genes_downreg,by="Gene")
  combine_upreg = combine_upreg[order(combine_upreg$'Gene'),]
  combine_downreg = combine_downreg[order(combine_upreg$'Gene'),]
  
  print(combine_upreg[,1:8])
  write.csv(combine_upreg, file = paste0(folder_pre,'diffExprGene_upreg_patient',patient,'.csv'),row.names=FALSE)
  write.csv(combine_downreg, file = paste0(folder_pre,'diffExprGene_downreg_patient',patient,'.csv'),row.names=FALSE)
  
}
#######################
for(patient in Dexa_VGPR$`Patient Number`){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  
  file = paste0(folder_pre,'diffExprGene_Patient',patient,'.csv')
  print(file)
  diffExprGene <- read.csv(file)
  combine_upreg <- merge(diffExprGene,genes_upreg,by="Gene")
  combine_downreg <- merge(diffExprGene,genes_downreg,by="Gene")
  combine_upreg = combine_upreg[order(combine_upreg$'Gene'),]
  combine_downreg = combine_downreg[order(combine_upreg$'Gene'),]
  
  print(combine_upreg[,1:8])
  
}
