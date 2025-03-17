library(dplyr)
library(SeuratDisk)
library(hdf5r)
library(Seurat) 
library(ggplot2)
library(RColorBrewer)


cs9_embryo_somite <- readRDS(file="cs9_embryo_somite.rds")
N70<-cs9_embryo_somite[,which(cs9_embryo_somite@meta.data$sample_num == "N70")]
meta_data_N70 <- N70@meta.data


##spatial plot for the location of NMP-N and NMP-M in slice 70
ggplot(data =meta_data_N70,aes(x=spatial_X,y=spatial_Y,colour=as.factor(somie_celltype)))+geom_point(size=2)+scale_color_manual(values=c("NMP>meso"="#3365aa","NMP>neural"="#fbeb93"))+guides(color = guide_legend(override.aes = list(size = 5)))+
theme_bw()+theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90,size=7),panel.grid.major=element_line(colour=NA),panel.grid.minor = element_blank(),legend.title = element_text(size = 16),
                 legend.text = element_text(size = 12)) +
                 guides(color = guide_legend(override.aes = list(size = 5)))


##density curve chart for the distribution of these two NMP along the D-V axis
meta_data_N70_sub <- meta_data_N70_sub[which(meta_data_N70_sub$somie_celltype %in% c("NMP>meso","NMP>neural"))]
ggplotmeta_data_N70_sub, aes(x=spatial_Y,fill=somie_celltype,colour =somie_celltype))+ geom_density(alpha=0.5,adjust=2,lwd = 1)+theme_bw()+
  scale_fill_manual(values = c("NMP>meso"="#3365aa","NMP>neural"="#fbeb93"))+
  scale_color_manual(values = c("NMP>meso"="#3365aa","NMP>neural"="#fbeb93"))+
  theme(panel.grid =element_blank()) + 
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(size=0.75, colour = "black")) 


##2D density plot
library(Nebulosa)

marker_gene<-c("WNT3A","WNT5B","FGF8","FGF17","FGF3","FGFR1")
for (i in 1:length(marker_gene)){
  pdf(file="NMP_pathway_spatial.pdf",width=3,height=3)
  plot_density(N70,features = marker_gene[1], joint =F,pal = "viridis")
  plot_density(N70,features = marker_gene[2], joint =F,pal = "viridis")
  plot_density(N70,features = marker_gene[3], joint =F,pal = "viridis")
  plot_density(N70,features = marker_gene[4], joint =F,pal = "viridis")
  plot_density(N70,features = marker_gene[5], joint =F,pal = "viridis")
  plot_density(N70,features = marker_gene[6], joint =F,pal = "viridis")
dev.off()