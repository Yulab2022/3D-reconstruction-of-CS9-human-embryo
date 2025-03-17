library(dplyr)
library(SeuratDisk)
library(hdf5r)
library(Seurat) 
library(ggplot2)
library(RColorBrewer)
library(ggplot2) 
library(cowplot)

cs9<- readRDS("cs9_embryo.rds")
meta_data_embryo<-read.table(file="meta_data_cs9_embryo.csv",head=T,sep=",")
meta_data_embryo_nei<-meta_data_embryo[which(meta_data_embryo$clusters_1.2_adj %in% c(10,21,2,15,18,8,12,5,17,14,11,16,13,0,16)),]
cs9_new_embryo<-cs9[,meta_data_embryo_nei$cells]


cs9_embryo_new <- NormalizeData(cs9_embryo_new)
cs9_embryo_new <- FindVariableFeatures(cs9_embryo_new,selection.method = "vst", nfeatures = 5000)
##VariableFeaturePlot(cs9_embryo_new_new)
all.genes <- rownames(cs9_embryo_new)
cs9_embryo_new <- ScaleData(cs9_embryo_new, features = all.genes)
cs9_embryo_new <- RunPCA(cs9_embryo_new, features = VariableFeatures(object = cs9_embryo_new))

cs9_embryo_new <- FindNeighbors(cs9_embryo_new, dims = 1:20)
cs9_embryo_new <- FindClusters(cs9_embryo_new, resolution = 1.9,algorithm = 2)


cs9_embryo_new <- RunUMAP(cs9_embryo_new, dims = 1:20, min.dist = 0.3,spread = 1.2)
pdf(file="cs9_cluster_embryo_umap_dim20mindist=0.3spread=1.2.pdf",width=15,height=10)
DimPlot(cs9_embryo_new, reduction = "umap",group.by=c("RNA_snn_res.1.9"),label = TRUE,label.size=5)
dev.off()

saveRDS(cs9_embryo_new,file="cs9_embryo_new.rds")