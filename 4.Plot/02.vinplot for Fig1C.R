library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(reshape2)
library(Seurat)
library(cowplot)
library(tidyr)
library(readxl)


cs9_embryo_new <- readRDS("cs9_embryo_new.rds")

mylevels<-c("FB","MB","HB","Sur.Ecto","SC-a","SC-p","Noto","NMP","Hea.Meso","Pha.Meso","Par.Meso","Somite","aPSM","pPSM","Car.Meso","MC","ECT","Spl.LPM","Som.LPM","Int.Meso","mLPM","FG","MG","HG","AT","Stalk","AM")
gene_list <- c("GABRP","TFAP2A","VTCN1","WNT6","SOX2","SFRP1","HESX1","IRX2","CRABP1","DLX5","PAX6","SHH","FOXA2","NKX1-2","TBXT","HES7","FGF17","TBX6","RIPPLY1","TCF15","MEOX1","SIX1","TCF21","FOXC2","PITX2","GATA4","GATA6","MYL7","KDR","CDH5","HAND1","HAND2","MSX1","OSR1","NID2","PAX1","APOA2","LCN15","DCC")
vln.dat=FetchData(cs9_embryo_new,c(gene_list,"cell_type"),slot="data")
vln.dat$Cell <- rownames(vln.dat)

#convert matrix format from wide to long
vln.dat.melt <- reshape2::melt(vln.dat, id.vars = c("Cell","cell_type"), 
                               measure.vars = gene_list,
                               variable.name = "gene",
                               value.name = "Expr")


vln.dat.melt$cell_type <-factor(vln.dat.melt$cell_type, levels= mylevels)



p<- ggplot(vln.dat.melt, aes(cell_type,Expr, fill = cell_type)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE,draw_quantiles = 0.5) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  scale_fill_manual(values = my_col) + 
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =element_blank()
        ##axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  ) 
p
