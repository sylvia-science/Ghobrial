
library(SingleCellSignalR)


library(readxl)
library(Seurat)

source('/home/sujwary/Desktop/scRNA/Code/Functions.R')


filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)


sample_type = 'Harmony_AllSamples_Sample_Kit'
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 

folder_name = 'AllSamples'
#folder_name = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample

folder = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                '/Batch_Sample_Kit/','/')
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )

setwd(filepath_cluster)

folder_name = 'AllSamples'
folder_main = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/',folder_name,
                     '/Batch_Sample_Kit/','/')
path = paste0(folder,'data_run','.Robj')
data_merge_run = loadRData(path)

tmp = data_merge_run@meta.data[paste0('RNA_snn_res.', resolution_val)]
tmp = tmp[,1]
tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
Idents(data_merge_run) = tmp

data_main_label = label_cells(data_merge_run,cluster_IDs)


celltype = 'Mono_DC'
cell_list = c('CD14+ Mono','CD16+ Mono','DC')
resolution_val_subset = 1.6
cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]

folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
path = paste0(folder_subcluster,'data_run','.Robj')
data_run_subset = loadRData(path)
tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
tmp = tmp[,1]
tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
Idents(data_run_subset) = tmp

data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)


Ident_main = colnames(data_main_label)
Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]

Ident_subset = colnames(data_run_subset_label)
Ident_subset_match = match(Ident_subset, Ident_main)
cell_subset = Idents(data_run_subset_label)
Ident_subset[1]
Ident_main[1]

newIdents = as.character(Idents(data_main_label))
newIdents2= newIdents
newIdents2[colnames(data_main_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
Idents(data_main_label) = newIdents2

plot = DimPlot(data_main_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)

print(plot)
##



patient_list = unique(data_main_label$`Patient Number`[data_main_label$Treatment %in% c('baseline','NBM')])

for (patient in patient_list){
  
  data_input = data_main_label[,data_main_label$Treatment %in% c('baseline','NBM')]
  print(patient)
  data_input = data_input[,data_input$`Patient Number` == patient]
  Treatment = unique(data_input$Treatment)
  data_input$Best_Overall_Response[data_input$Best_Overall_Response == 'MR' ] = 'PR'
  data_input$Best_Overall_Response[data_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  data_input$input_ident = paste0('P',data_input$'Patient Number', ' ',
                                  data_input$Treatment, ' ', Idents(data_input))
  
  #data_input$input_ident = paste0(data_input$Treatment, ' ', Idents(data_input))
  data_input$input_ident = gsub("NBM NA", "NBM", data_input$input_ident)
  
  #data_input$input_ident =  Idents(data_input)
  
  
  #data_input = data_input[,data_input$Best_Overall_Response == 'PR']
  remove_list = c('Remove','DC/T-Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL')
  data_input = data_input[,!(Idents(data_input) %in% remove_list )]
  
  
  data_matrix = as.matrix(data_input@assays[["RNA"]]@data)
  data_matrix = data_matrix[Matrix::rowSums(data_matrix)!=0,]
  umap = data_input@reductions[["umap"]]@cell.embeddings
  labels = factor(data_input$input_ident)
  
  
  # class = cell_classifier(data=as.data.frame(data_matrix), genes=rownames(data_matrix), 
  #                         markers = markers(c("immune")), 
  #                         tsne=umap,
  #                         plot.details=TRUE,write = FALSE)
  
  
  try(signal_auto <- cell_signaling(data = as.data.frame(data_matrix), 
                                genes = rownames(data_matrix), 
                                cluster = as.numeric(labels), 
                                c.names = levels(labels),
                                int.type = c( "autocrine"), write = T))
  
  file.rename(paste0(filepath_cluster, 'cell-signaling/'), 
              paste0(filepath_cluster, 'cell-signaling autocrine ', Treatment,' ',patient,'/'))
  
  try(signal_para <- cell_signaling(data = as.data.frame(data_matrix), 
                                genes = rownames(data_matrix), 
                                cluster = as.numeric(labels), 
                                c.names = levels(labels),
                                int.type = c("paracrine"), write = T))
  
  file.rename(paste0(filepath_cluster, 'cell-signaling/'), paste0(filepath_cluster, 'cell-signaling paracrine ', Treatment,' ',patient,'/'))
  
}




patients_SMM = unique(data_main_label$`Patient Number`[data_main_label$Treatment == 'baseline'])
patients_NBM  = unique(data_main_label$`Patient Number`[data_main_label$Treatment == 'NBM'])

unique(Idents(data_main_label) )

source('~/Desktop/scRNA/Code/singleCellRHelper.R')

cell_list = unique(Idents(data_main_label) )
cell_list = cell_list[!(cell_list %in% c('Remove','DC/T-Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL','14','42')) ]

colnames = c("Cell1","Cell2","Cell1_gene", "Cell2_gene", 
             "interaction.type", "LRscore", "Patient", "Treatment", 
             "interactionName", "frac_high_baseline", "frac_high_NBM" ,"frac_diff",
             "sig_high_baseline", "sig_high_NBM","p.value", "FDR"     )
data_interaction_all = setNames(data.frame(matrix(ncol = 16, nrow = 0)), colnames)

for (cell in cell_list){
  ident1='CD8+ T Cell'
  ident2 = cell
  data_interaction = singleCellRHelper(ident1, ident2,data_main_label, metaData, filepath_cluster)
  #data_interaction
  if (nrow(data_interaction) > 0){
    data_interaction$Cell1 = ident1
    data_interaction$Cell2 = ident2
    colnames(data_interaction)[1]= 'Cell1_gene'
    colnames(data_interaction)[2] = 'Cell2_gene'
    data_interaction_all = rbind(data_interaction_all,data_interaction)
  }
  if (sum(data_interaction$sig_high_baseline) + sum(data_interaction$sig_high_NBM) > 0) {
    data= data_interaction[c('interactionName','frac_high_baseline', 'frac_high_NBM','FDR', 'sig_high_baseline','sig_high_NBM')]
    data = data[order(data$FDR),]
    print('High Baseline')
    print(unique(data[data$sig_high_baseline == T,]))
    print('High NBM')
    print(unique(data[data$sig_high_NBM == T,]))
    
    
    }
  print(' ')
}

 


library(dplyr)
col_idx <-  which(names(data_interaction_all) == 'Cell2' )
data_interaction_all <- data_interaction_all[, c(col_idx, (1:ncol(data_interaction_all))[-col_idx])]
col_idx <-  which(names(data_interaction_all) == 'Cell1' )
data_interaction_all <- data_interaction_all[, c(col_idx, (1:ncol(data_interaction_all))[-col_idx])]

data= unique(data_interaction_all[c('Cell1','Cell2','interactionName','interaction.type','p.value')])
data$FDR_all = p.adjust(data$p.value, method="BH") 

plot = ggplot(data, aes(x = p.value))  + geom_histogram(bins = 30) + scale_y_log10()


print(plot)

for (i in 1:nrow(data_interaction_all )){
  row = data_interaction_all[i,]
  data_interaction_all$FDR_all[i] = data$FDR_all[data$Cell1 == row$Cell1 & data$Cell2 == row$Cell2 & data$interactionName == row$interactionName]
}


data_interaction_all$sig_high_baseline = (data_interaction_all$frac_diff > 0 & data_interaction_all$p.value < 0.05)
data_interaction_all$sig_high_NBM = (data_interaction_all$frac_diff < 0 & data_interaction_all$p.value < 0.05)

data_interaction_all_sig = data_interaction_all[data_interaction_all$p.value < 0.05,]

data_interaction_all_sig = data_interaction_all_sig[, c("Cell1" ,"Cell2", "Cell1_gene", "Cell2_gene", 
                                                        "interaction.type" ,"interactionName",
                                                        "frac_high_baseline","frac_high_NBM"  ,
                                                        "sig_high_baseline", "sig_high_NBM", "p.value" ,"FDR" ,'FDR_all')]
data_interaction_all_sig = unique(data_interaction_all_sig)

data_interaction_all_sig = data_interaction_all_sig[order(data_interaction_all_sig$Cell2, 
                                                          data_interaction_all_sig$sig_high_baseline,
                                                          data_interaction_all_sig$p.value),]

#write.csv(data_interaction_all_sig,paste0(filepath_cluster,
#                                          '/SingleCellSignalR/Baseline/','All',' Vs ', ident2,' Sig.csv'))

write.csv(data_interaction_all_sig,paste0(filepath_cluster,
                                          '/SingleCellSignalR/Baseline/',ident1,' Vs All Sig.csv'))



ident2= 'CD16+ Mono'
ident1 ='NK'

data_interaction = singleCellRHelper(ident1, ident2,data_main_label, metaData, filepath_cluster)

unique(data_interaction$interactionName[data_interaction$sig_high_baseline == T])
unique(data_interaction$interactionName[data_interaction$sig_high_NBM == T])

colnames = c(ident1,ident2, "interaction.type", "LRscore", 'Patient','Treatment')
data_interaction = setNames(data.frame(matrix(ncol = 6, nrow = 0)), colnames)




visualize_interactions(signal)

ident1='T Cell'
ident2 ='CD16+ Mono'

visualize_interactions(signal, show.in=paste0(ident1,'-',ident2))


signal[[paste0(ident1,'-',ident2)]]

for (i in 1:600){
  signal_i = signal_auto[[i]]
  cols = colnames(signal_i)
  if (nrow(signal_i) >20){
    print(paste0(cols[1],' ', cols[2]))
    print(nrow(signal_i))
    #print(signal[[i]])
  }
  
}

ident1='CD16+ Mono'
ident2 ='T Cell'
ident1='CD14+ Mono'
ident2 ='CD8+ T Cell'

visualize_interactions(signal_auto, show.in=paste0(ident1,'-',ident2))

signal_auto[[paste0(ident1,'-',ident2)]]


for (i in 1:1000){
  signal_i = signal_para[[i]]
  cols = colnames(signal_i)
  if ('IFNG' %in% signal_i[,1] || 'IFNG' %in% signal_i[,2]){
    print(paste0(cols[1],' ', cols[2]))
    #print(nrow(signal_i))
    #print(signal[[i]])
  }
  
}

ident1='CD16+ Mono'
ident2 ='T Cell'

#visualize_interactions(signal, show.in=paste0(ident1,'-',ident2))
signal[[paste0(ident1,'-',ident2)]]

ident1='NK'
ident2 ='Plasma Cell'
ident1='CD16+ Mono'
ident2 = 'NK'
signal[[paste0(ident1,'-',ident2)]]
