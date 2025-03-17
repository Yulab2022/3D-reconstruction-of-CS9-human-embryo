library(Seurat)

####import data
cs9_embryo_somite <- readRDS(file="cs9_embryo_somite.rds")
NMP_monkey <- readRDS(file="monkey_NMP.rds")
NMP_mouse <- readRDS(file="mouse_NMP.rds")
NMP_fish <- readRDS(file="fish_NMP.rds")

####prepare orthologous genes across species
library(biomart)
library(biomaRt)


umart = useMart('ensembl')
datalist = listDatasets(umart)
searchDatasets(mart = umart, pattern = "Human")


human_BM<-getBM(attributes = c("ensembl_gene_id","external_gene_name","mmusculus_homolog_ensembl_gene","mfascicularis_homolog_ensembl_gene"),mart = human,uniqueRows = T)
macaque_BM<-getBM(attributes = c("ensembl_gene_id","external_gene_name"),mart = macaque,uniqueRows = T)
mouse_BM<-getBM(attributes = c("ensembl_gene_id","external_gene_name"),mart = mouse,uniqueRows = T)
zeberafish_BM<-getBM(attributes = c("ensembl_gene_id","external_gene_name"),mart = zeberafish,uniqueRows = T)


macaque <- useMart("ensembl", dataset = "mfascicularis_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
zeberafish <- useMart("ensembl", dataset = "drerio_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 


gene.monkey2human <- getLDS(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name"),filters = "external_gene_name",values = unique(macaque_BM$external_gene_name),mart = macaque,attributesL = c("ensembl_gene_id","external_gene_name","chromosome_name"),martL = human,uniqueRows = T)
gene.mouse2human <- getLDS(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name"),filters = "external_gene_name",values = unique(mouse_BM$external_gene_name),mart = mouse,attributesL = c("ensembl_gene_id","external_gene_name","chromosome_name"),martL = human,uniqueRows = T)
gene.zeberafish2huamn <- getLDS(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name"),filters = "external_gene_name",values = unique(zeberafish_BM$external_gene_name),mart = zeberafish,attributesL = c("ensembl_gene_id","external_gene_name","chromosome_name"),martL = human,uniqueRows = T)


save.image(file="biomart_human_mouse_monkey.RData")


####one to one orthologous genes
load(file="biomart_human_mouse_monkey.RData")

monkey2human <- unique(gene.monkey2human[,c(2,5)])
monkey2human <- subset(monkey2human, !duplicated(mouse2human$Gene.name))
monkey2human <- subset(monkey2human, !duplicated(mouse2human$Gene.name.1))
monkey2human <- monkey2human[which(monkey2human$Gene.name.1 != ""),]
rownames(monkey2human) <- monkey2human$Gene.name.1

mouse2human <- unique(gene.mouse2human[,c(2,5)])
mouse2human <- subset(mouse2human, !duplicated(mouse2human$Gene.name))
mouse2human <- subset(mouse2human, !duplicated(mouse2human$Gene.name.1))
mouse2human <- mouse2human[which(mouse2human$Gene.name.1 != ""),]
rownames(mouse2human) <- mouse2human$Gene.name.1


zeberafish2human <- unique(gene.zeberafish2huamn[,c(2,5)])
zeberafish2human <- subset(zeberafish2human, !duplicated(zeberafish2human$Gene.name))
zeberafish2human <- subset(zeberafish2human, !duplicated(zeberafish2human$Gene.name.1))
zeberafish2human <- zeberafish2human[which(zeberafish2human$Gene.name.1 != ""),]
rownames(zeberafish2human) <- zeberafish2human$Gene.name.1



common_gene<-Reduce(intersect,list(mouse2human$Gene.name.1,zeberafish2human$Gene.name.1,monkey2human$Gene.name.1,rownames(cs9_embryo_somite)),accumulate =TRUE)


####human matrix###########
gene_cell_exp_human <- AverageExpression(cs9_embryo_somite,
                                         features = common_gene,
                                         group.by = 'somite_celltype',
                                         slot = 'data')
gene_cell_exp_human <- as.data.frame(gene_cell_exp_human$RNA)
gene_cell_exp_human <- gene_cell_exp_human[,c("NG","NMP","PSM","aPSM","Somite")]
gene_cell_exp_human <- t(scale(t(gene_cell_exp_human),scale = TRUE,center = TRUE))
rownames(gene_cell_exp_human) <- paste("human_",rownames(gene_cell_exp_human),sep="")



####monkey matrix####################
gene_cell_exp_monkey <- AverageExpression(NMP_monkey,
                                          features = monkey2human[common_gene_1,1],
                                          group.by = 'predicted.celltype',
                                          slot = 'data')
gene_cell_exp_monkey <- as.data.frame(gene_cell_exp_monkey$RNA)
gene_cell_exp_monkey <- gene_cell_exp_monkey[,c("NG","NMP","PSM","aPSM","Somite")]
gene_cell_exp_monkey <- t(scale(t(gene_cell_exp_monkey),scale = TRUE,center = TRUE))
rownames(gene_cell_exp_monkey) <- c(monkey2human[which(monkey2human$Gene.name %in% rownames(gene_cell_exp_monkey)),2])
rownames(gene_cell_exp_monkey) <- paste("monkey_",rownames(gene_cell_exp_monkey),sep="")


####mouse matrix###########
gene_cell_exp_mouse <- AverageExpression(NMP_mouse,
                                         features = mouse2human[common_gene_1,1],
                                         group.by = 'predicted.celltype',
                                         slot = 'data')
gene_cell_exp_mouse <- as.data.frame(gene_cell_exp_mouse$RNA)
gene_cell_exp_mouse <- gene_cell_exp_mouse[,c("NG","NMP","PSM","aPSM","Somite")]
gene_cell_exp_mouse <- t(scale(t(gene_cell_exp_mouse),scale = TRUE,center = TRUE))
rownames(gene_cell_exp_mouse) <- c(mouse2human[which(mouse2human$Gene.name %in% rownames(gene_cell_exp_mouse)),2])
rownames(gene_cell_exp_mouse) <- paste("mouse_",rownames(gene_cell_exp_mouse),sep="")


####zeberafish matrix###########################

gene_cell_exp_fish <- AverageExpression(NMP_fish,
                                        features = zeberafish2human[common_gene_1,2],
                                        group.by = 'predicted.celltype',
                                        slot = 'data')
gene_cell_exp_fish <- as.data.frame(gene_cell_exp_fish$RNA)
gene_cell_exp_fish <- gene_cell_exp_fish[,c("NG","NMP","PSM","aPSM","Somite")]
gene_cell_exp_fish <- t(scale(t(gene_cell_exp_fish),scale = TRUE,center = TRUE))
rownames(gene_cell_exp_fish) <- c(zeberafish2human[which(zeberafish2human$Gene.name %in% rownames(gene_cell_exp_fish)),2])
rownames(gene_cell_exp_fish) <- paste("fish_",rownames(gene_cell_exp_fish),sep="")


final<-rbind(gene_cell_exp_human,gene_cell_exp_monkey,gene_cell_exp_mouse,gene_cell_exp_fish)



####Mfuzz ################

library(Mfuzz)

#create object
mfuzz_class<- ExpressionSet(assayData = as.matrix(final))

#nreprocess missing values or outliers
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

#normalize
mfuzz_class <- standardise(mfuzz_class)


set.seed(123)
cluster_num <- 5
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

#plot
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(1, 5), time.labels = colnames(mfuzz_class))
