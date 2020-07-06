#Import meta table
cell = 'CD14+ Mono'

stats_summary = rbind(stats_category[[1]][, c(cell,'category')],stats_category[[2]][, c(cell,'category')],stats_category[[3]][, c(cell,'category')])
stats_summary$category = factor(stats_summary$category , levels = c('NBM','Pre','Post'))


meta <- read_excel("C:/Users/Sylvia/Dropbox (Partners HealthCare)/Sylvia_Romanos/scRNASeq/Data/Dexa_meta.xlsx")
stats_summary$Sample <- rownames(stats_summary)
meta2 <- merge(stats_summary, meta, by= "Sample")

meta2$Category <-paste0(meta2$category, "_",meta2$`10X kit`)
meta2$Category[meta2$Category == "NBM_NBM"] <- "NBM"
meta2$Category <- factor(meta2$Category, levels = c("NBM", "Pre_v2", "Post_v2", "Pre_v3", "Post_v3"))

wilcox.test(meta2$`CD14+ Mono`[meta2$Category== "Pre_v3"], meta2$`CD14+ Mono`[meta2$Category== "Post_v3"])


x_name = 'category'
plot = ggplot(meta2, aes(x = Category, y = !!ensym(cell))) + 
  geom_boxplot()+ geom_jitter()+
  coord_cartesian(ylim = c(0, 50))+
  ggtitle(cell)+
  xlab("") + ylab("%")

plot = plot + theme(
  plot.title = element_text(color="black", size=24, face="bold.italic"),
  axis.title.x = element_text(color="black", size=24, face="bold"),
  axis.title.y = element_text(color="black", size=24, face="bold"),
  axis.text=element_text(size=24)
)



plot = plot + theme(plot.title = element_text(hjust = 0.5))