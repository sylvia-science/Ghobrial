
singleCellRHelper = function(ident1, ident2,data,metaData, filepath_cluster){
  #browser()
  print(paste0(ident1,' ', ident2))
  colnames = c(ident1,ident2, "interaction.type", "LRscore", 'Patient','Treatment')
  data_interaction = setNames(data.frame(matrix(ncol = 6, nrow = 0)), colnames)
  patient_list = unique(data$`Patient Number`[data$Treatment %in% c('baseline','NBM')])
  
  
  for (patient in patient_list){
    #print(patient)
    Treatment = metaData$Treatment[metaData$`Patient Number` == patient]
    Treatment = Treatment[Treatment %in% c('baseline','NBM')]
    path = paste0(filepath_cluster, 'SingleCellSignalR/Baseline/', 'cell-signaling paracrine ', Treatment,' ',patient,'/')
    filename = paste0('LR_interactions_','P',patient,' ', Treatment,' ', ident1
                      ,'-','P',patient,' ', Treatment, ' ', ident2)
    if (file.exists(paste0(path,filename,'-paracrine.txt'))){
      data_para = read.delim(paste0(path,filename,'-paracrine.txt'))
      data_para$Patient = patient
      data_para$Treatment = Treatment
      colnames(data_para) = c(ident1,ident2, "interaction.type", "LRscore", 'Patient', 'Treatment')
      data_interaction = rbind(data_interaction, data_para)
      
    }else{
      #print("Auto Doesn't Exist")
    }
    
    if (file.exists(paste0(path,filename,'-autocrine.txt'))){
      data_auto = read.delim(paste0(path,filename,'-autocrine.txt'))
      data_auto$Patient = patient
      data_auto$Treatment = Treatment
      colnames(data_auto) = c(ident1,ident2, "interaction.type", "LRscore", 'Patient', 'Treatment')
      data_interaction = rbind(data_interaction, data_auto)
    }else{
      #print("Para Doesn't Exist")
    }
    
  }
  
  if (nrow(data_interaction) == 0){
    return( data_interaction)
  }
  
  data_interaction = unique(data_interaction)
  data_interaction$interactionName = paste0(data_interaction[,1],' ',data_interaction[,2])
  data_interaction$frac_high_baseline = NA
  data_interaction$frac_high_NBM = NA
  data_interaction$frac_diff = NA
  data_interaction$sig_high_baseline = F
  data_interaction$sig_high_NBM = F
  data_interaction$p.value = 1
  
  all_baseline = length(unique(data_interaction$Patient[data_interaction$Treatment =='baseline']))
  all_NBM = length(unique(data_interaction$Patient[data_interaction$Treatment =='NBM']))
  #browser()
  for (interaction in unique(data_interaction$interactionName) ){
    #browser()
    
    #browser()
    tmp = data_interaction[data_interaction$interactionName == interaction &  
                             data_interaction$Treatment =='baseline',]
    n_occur <- data.frame(table(tmp$Patient))
    if (any(n_occur$Freq > 1)){
      print(tmp)
      browser()
    }
    
    
    num_baseline = length(unique(data_interaction$Patient[data_interaction$interactionName == interaction &  
                                                     data_interaction$Treatment =='baseline']))
    num_NBM = length(unique(data_interaction$Patient[data_interaction$interactionName == interaction &  
                                                       data_interaction$Treatment =='NBM']))
    
    
    frac_baseline = num_baseline/all_baseline
    frac_NBM = num_NBM/all_NBM
    data_interaction$frac_high_baseline[data_interaction$interactionName == interaction]  = frac_baseline
    data_interaction$frac_high_NBM[data_interaction$interactionName == interaction]  = frac_NBM
    data_interaction$frac_diff[data_interaction$interactionName == interaction] = frac_baseline - frac_NBM
    #browser()
    #print(all_NBM)
    #print(all_baseline)
    #print(num_NBM)
    #print(num_baseline)

    table_input = rbind(c(all_NBM- num_NBM,num_NBM), 
                        c(all_baseline - num_baseline,num_baseline))
    colnames(table_input) = c('NotUp','Up')
    rownames(table_input) = c('NBM','SMM')
    #browser()
    library(rcompanion)
    
    test = fisher.test(table_input)
    test2 = pairwiseNominalIndependence(table_input,
                                fisher = TRUE,
                                gtest  = FALSE,
                                chisq  = T,
                                digits = 3)
    #browser()
    data_interaction$p.value[data_interaction$interactionName == interaction]  = test$p.value
    #data_interaction$chisq[data_interaction$interactionName == interaction]  =test2$p.Chisq
    
    
    mean_frac = mean(c(frac_baseline,frac_NBM))
  
    if ( test$p.value < 0.05 && (frac_baseline - frac_NBM)  > 0){
      #browser()
      #print(interaction)
      data_interaction$sig_high_baseline[data_interaction$interactionName == interaction] = T
      
      
    }else if (test$p.value < 0.05 && (frac_baseline - frac_NBM ) < 0){
      #browser()
      #print(interaction)
      data_interaction$sig_high_NBM[data_interaction$interactionName == interaction] = T
    }
    
  }
  
  #browser()
  
  data_interaction$p.value = as.numeric(as.character(data_interaction$p.value))
  
  data= unique(data_interaction[c('interactionName','interaction.type','p.value')])
  data$FDR = p.adjust(as.numeric(data$p.value), method="BH", n= length(data$p.value)) 
  
  data_interaction$FDR = 1
  for (i in 1:nrow(data_interaction )){
    row = data_interaction[i,]
    data_interaction$FDR[i] = data$FDR[data$interactionName == row$interactionName & 
                                         data$interaction.type == row$interaction.type ]
  }
  
  #print(hist(data$p.value ), breaks = 12)
  #print(nrow(data))
  #print(sum(data$FDR < 0.05))
  #browser()
  #data_interaction$FDR = p.adjust(data_interaction$p.value, method="BH")                            
  
  
  #browser()
  
  return(data_interaction)
  
  #unique(data_interaction$interactionName[data_interaction$sig_high_baseline == T])
  #unique(data_interaction$interactionName[data_interaction$sig_high_NBM == T])

}