
# source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/AnalyzePairDE.R')

source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Functions.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Plot_func.R')
source('C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Integration/FunctionsIntegrate.R')

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(stringr)

filename_sampleParam_integrate <- 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_integrate_parameters.xlsx'
sampleParam_integrate <- read_excel(filename_sampleParam_integrate)

base_path_gene = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa Genes/'

filename_list = c('AD2009_general','AD2009_specific1','AD2009_specific2'
                  ,'AG2018','AM2012','HJ2012'
                  ,'JG2001_1','JJ2019_1','JJ2019_2','KK2017'
                  ,'KM2009', 'MB2012', 'MG2018_1'
                  ,'MG2018_2','MG2018_3')
Log2FC_str = c('Log2FC', 'Log2FC_vivo','Log2FC'
               ,NA,'Log2FC1','Log2FC'
               ,'log2FC_ave', NA, NA, 'log2FC_ave'
               ,'Log2FC', 'Log2FC', 'Log2FC'
               ,'log2FC_FP','Log2FC_FP')

folder_base_output = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Code/Output/Integrate Pair/'
sample_type = 'BM'

filename_sample_Integrate_pairs = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/sample_Integrate_pairs.xlsx'
sample_Integrate_pairs = read_excel(filename_sample_Integrate_pairs)

filename_metaData = 'C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx'
metaData = read_excel(filename_metaData)


###########################################################
## Combine all paper gene data together
###########################################################

data_paper = data.frame(matrix(ncol = 4, nrow = 0))
names(data_paper) = c('Gene','Log2FC','paper','paper_cell') 

for (i in 1:length(filename_list)){
  filename = paste0(base_path_gene,filename_list[i],'.xlsx')
  data_orig = read_excel(filename)
  print(filename_list[i])
  print(colnames(data_orig))
  print('')
  
  data = subset(data_orig, select = c(Gene))
  data[,'Log2FC'] = NA
  data[,'paper'] = NA
  data[,'paper_cell'] = NA
  if (is.na(Log2FC_str[i]) ){
    data$Log2FC = NA
  }else{
    col_name = Log2FC_str[i]
    print(col_name)
    log2FC  = data_orig[,col_name]
    log2FC = as.data.frame(log2FC)
    log2FC = log2FC[,col_name]
    log2FC = gsub(',','.',log2FC)
    data$Log2FC = log2FC
    
  }
  
  if (filename_list[i] =='JJ2019_1' || filename_list[i] == 'JJ2019_2'){
    data$Gene = gsub("\\|.*","",data$Gene) # Get rid of things after "|'
  }
  
  data$paper =  filename_list[i]
  data$paper_cell = data_orig$Cell
  #browser()
  data_paper = rbind(data_paper, data)
}

data_paper = unique(data_paper)


###############################################################################
## Combine data from DE and from data_paper
###############################################################################

#combine_all = data.frame(matrix(ncol = ncol(combine), nrow = 0))
#names(combine_all) = names(combine)
patient_list = c(10, 5, 20, 12, 34, 28, 21, 31, 16, 51, 6, 40) # all

combine_all = c()

for(patient in patient_list){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  
  file = paste0(folder_pre,'DE/DE_Patient',patient,'.csv')
  print(file)
  diffExprGene = read.csv(file)
  combine = merge(diffExprGene,data_paper,by="Gene") # Take only genes that appear in literature
  #combine = combine[order(combine$'Gene'),]
  combine = combine[combine$p_val_adj < 0.1,]
  combine = combine[order(combine$Cell),]

  combine  = combine[order(combine$Gene),]
  combine  = levels(unique(combine$Gene))
  combine_all = c(combine_all,combine)
}
combine_all_count = as.data.frame(table(combine_all))

combine_all_count_filter = combine_all_count[combine_all_count$Freq >2,] # Genes that appear in DE, lit, and at least twice

############################################################################
## Save only genes that appear at least twice in all lists
############################################################################



for(patient in patient_list){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  
  file = paste0(folder_pre,'DE/DE_Patient',patient,'.csv')
  print(file)
  diffExprGene = read.csv(file)
  
  combine = merge(diffExprGene,data_paper,by="Gene")
  combine = combine[combine$p_val_adj < 0.1,]
  combine = combine[order(combine$Cell),]
  gene_intersect = intersect(combine$Gene,combine_all_count_filter$combine_all) 
  
  # Take only genes that are in gene_intersect
  combine_filter = combine
  for (i in 1:length(combine$Gene)){
    if (!any(gene_intersect == combine$Gene[i])){ # if a combine gene is not in the intersect, remove it
      #print('not found')
      combine_filter = combine_filter[combine_filter$Gene != combine$Gene[i],]
    }else{
      #print('found')
    }
  }
  #combine_filter = combine[combine$Gene == gene_intersect,]
  #print(combine[,1:8])
  write.csv(combine, file = paste0(folder_pre,'DE/DE_Paper_patient',patient,'.csv'),row.names=FALSE)
  write.csv(combine_filter, file = paste0(folder_pre,'DE/DE_Paper_filter_patient',patient,'.csv'),row.names=FALSE)
}





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




#sample_type = 'PB'
output = getStatSummary (Dexa_MinResponse$`Patient Number`, sample_Integrate_pairs,folder_base_output)
  

for(patient in patient_list){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  
  file = paste0(folder_pre,'DE/DE_Patient',patient,'.csv')
  print(file)
  diffExprGene = read.csv(file)
  
  combine = merge(diffExprGene,data_paper,by="Gene")
  combine = combine[combine$p_val_adj < 0.1,]
  combine = combine[order(combine$Cell),]
  gene_intersect = intersect(combine$Gene,combine_all_count_filter$combine_all) 
  
  # Take only genes that are in gene_intersect
  combine_filter = combine
  for (i in 1:length(combine$Gene)){
    if (!any(gene_intersect == combine$Gene[i])){ # if a combine gene is not in the intersect, remove it
      print('not found')
      combine_filter = combine_filter[combine_filter$Gene != combine$Gene[i],]
    }else{
      print('found')
    }
  }
  #combine_filter = combine[combine$Gene == gene_intersect,]
  #print(combine[,1:8])
  write.csv(combine, file = paste0(folder_pre,'DE/DE_Paper_patient',patient,'.csv'),row.names=FALSE)
  write.csv(combine_filter, file = paste0(folder_pre,'DE/DE_Paper_filter_patient',patient,'.csv'),row.names=FALSE)
}

# Get cluster DE for pre vs post
# Make folder DE/FeaturePlot
# Save MakeFeaturePlot FeaturePlots

for(patient in 40){ # Patient numbers 
  print('')
  print(patient)
  pair_list =  sample_Integrate_pairs[ sample_Integrate_pairs$'Patient Number' == patient, ]
  
  sample_name_pre = pair_list[[paste0('Sample Pre ', sample_type)]]
  folder_pre = makeFolders(folder_base_output,sample_name_pre,filter = TRUE,regress_TF = TRUE,FALSE)
  folder_output = paste0(folder_pre,'DE/FeaturePlot/')
  dir.create( folder_output, recursive = TRUE)
  

  data = loadRData(paste0(folder_pre,'data.Robj'))
  cluster_IDs = sampleParam_integrate[['Cluster_IDs']][sampleParam_integrate['Patient Number'] == patient]
  data = label_cells(data,cluster_IDs)
  
  
  file = paste0(folder_pre,'DE/DE_Patient',patient,'.csv')
  print(file)
  diffExprGene = read.csv(file)

  for (i in 1:(nrow(diffExprGene))){
    feature = as.character(diffExprGene$Gene[i])
    cell_str = as.character(diffExprGene$Cell[i])
    print(i)
    print(feature)
    browser()
    
    FeaturePlotFix(data, feature,folder_output,cell_str, split = TRUE)
  }
  # cnt = 1
  # for (i in 1:(length(diffExprGene)/3)){
  #   cnt = i*3
  #   print(cnt)
  #   if (cnt < length(diffExprGene)){
  #     gene_idx = (cnt-2):cnt
  #     gene_list = diffExprGene[gene_idx]
  #   }else{
  #     gene_idx = (cnt-2):length(diffExprGene)
  #     gene_list = diffExprGene[gene_idx]
  #   }
  #   gene_list = paste( unlist(gene_list), collapse=', ')
  #   
  #   cell_features = data.frame(matrix(ncol = 2, nrow = 1))
  #   colnames(cell_features) = c("Cell",'Markers')
  #   cell_features$Cell =  paste0('Gene',gene_idx[1],'_',gene_idx[length(gene_idx)])
  #   cell_features$Markers =  gene_list
  #   
  #   MakeFeaturePlot(data,data,folder_output,cell_features, split = TRUE)
  # }
  
  

}
