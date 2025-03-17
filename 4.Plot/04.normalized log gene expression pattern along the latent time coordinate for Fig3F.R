library(dplyr)
library(SeuratDisk)
library(hdf5r)
library(Seurat) 
library(ggplot2)
library(RColorBrewer)

cs9_embryo_somite <- readRDS(file="cs9_embryo_somite.rds")
final_maker<-c("ALDH1A2",'SOX2','NKX1-2','CYP26A1','TBXT','HES7','TBX6','MSGN1','FOXC1')

##extract expression matrix
gene_matrix<- cs9_embryo_somite@assays$RNA@data[unique(final_maker),]
gene_matrix<- t(as.matrix(gene_matrix))

##normalize log gene expression
for(i in 1:dim(gene_matrix)[2])
{
  temp<- (gene_matrix[,i]-min(gene_matrix[,i]))/(max(gene_matrix[,i])-min(gene_matrix[,i]))
  
  gene_matrix[,i]<-temp
  
}
gene_matrix<- cbind(cs9_embryo_somite@meta.data,gene_matrix)
gene_matrix<-data.frame(gene_matrix)


##plot
my_umap_Plot <- function(i){
  ggplot(gene_matrix, aes(latent_time_new,gene_matrix[,i + dim(cs9_embryo_somite@meta.data)[2]]))+
    geom_point(aes(color=somie_celltype),size=1)+scale_color_manual(values=somite)+
    geom_smooth(method = "gam",se=FALSE,color="black")+ theme_classic()+theme(axis.line.x=element_line(linetype="dotted",color="black",size=0.8),axis.line.y=element_line(linetype="dotted",color="black",size=0.8))+ggtitle(check_marker[i])
}

umap_list <- lapply(1:9, my_umap_Plot)


library(cowplot)
p<-plot_grid(plotlist = umap_list, align = "h", ncol = 1)
print(p)
