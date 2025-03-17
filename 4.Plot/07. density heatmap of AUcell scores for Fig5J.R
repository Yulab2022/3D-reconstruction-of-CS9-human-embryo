library(irGSEA)
library(dplyr)
library(SeuratDisk)
library(hdf5r)
library(Seurat) 
library(ggplot2)
library(RColorBrewer)

cs9_embryo_LPM<-readRDS(file="cs9_embryo_LPM.rds")

sce <- irGSEA.score(object = cs9_embryo_LPM,assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "C5",  
                             subcategory = "CP:KEGG", geneid = "symbol",
                             method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')


##Integrate differential gene set
result.dge_customclassif <- irGSEA.integrate(object = sce_new,
                               group.by = "LPM_celltype",
                               method = c("AUCell","UCell","singscore","ssgsea"))


##density heatmap
densityheatmap <- irGSEA.densityheatmap(object = sce_new,
                                        method = "AUCell",
                                        show.geneset = "KEGG-TGF-BETA-SIGNALING-PATHWAY")
densityheatmap

densityheatmap <- irGSEA.densityheatmap(object = sce_new,
                                        method = "AUCell",
                                        show.geneset = "KEGG-ECM-RECEPTOR-INTERACTION")
densityheatmap