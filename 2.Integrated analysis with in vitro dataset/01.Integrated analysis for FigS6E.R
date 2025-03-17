# Load required R packages
library(scater)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(tidyverse)

## Load cs9 data
cs9 = readRDS("cs9_embryo_somite.rds")
segmentoid = readRDS("3d_segmentoid.rds")

my_palette <- c("#8665ae","#6e8db8","#305c99","#f6e5a0","#fdbf70","#fba5a3","#c999ca")

colnames(cs9@meta.data)

# Rename cs9 metadata columns
cs9@meta.data$orig.ident <- "CS9"

cs9@meta.data$celltype <- paste(cs9@meta.data$orig.ident, cs9@meta.data$somite_celltype, sep = "_")
unique(cs9@meta.data$celltype)

cs9@meta.data$celltype <- as.factor(cs9@meta.data$celltype)
levels(cs9@meta.data$celltype)

# Extract specific cell subpopulations
Idents(cs9) <- "celltype"
interesting_groups <- c("CS9_aPSM", "CS9_Early somite", "CS9_Later somite", "CS9_NMP", "CS9_pPSM")
cs9 <- subset(cs9, idents = interesting_groups)
unique(cs9@meta.data$celltype)

# Create UMAP plot
plot <- DimPlot(cs9, reduction = "umap", group.by = "celltype", label = TRUE) +
  labs(x = "umap_1", y = "umap_2", title = "UMAP", color = "Version") +
  scale_color_manual(values = my_palette) +  # Set custom color palette
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("new_cluster_planA")
plot

## CCA-based integration method
yu.list <- list(
  segmentoid = segmentoid,
  cs9 = cs9
)

yu.list <- lapply(X = yu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

features <- SelectIntegrationFeatures(object.list = yu.list)

# Perform data integration
anchors <- FindIntegrationAnchors(object.list = yu.list, anchor.features = features)
yu.integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(yu.integrated) <- "integrated"

# Process integrated dataset
all.genes <- rownames(yu.integrated)
yu.integrated <- ScaleData(yu.integrated, features = all.genes)


yu.integrated <- RunPCA(yu.integrated, npcs = 20, verbose = FALSE)
yu.integrated <- RunUMAP(yu.integrated, reduction = "pca", dims = 1:20, min.dist = 0.7, spread = 1.2)
yu.integrated <- FindNeighbors(yu.integrated, reduction = "pca", dims = 1:20)
yu.integrated <- FindClusters(yu.integrated, resolution = 0.4)

DimPlot(yu.integrated, reduction = "umap", group.by = "celltype", label = TRUE)
saveRDS(yu.integrated, file = "3d_CCA.rds")

# Load 3D dataset
segmentoid3d = readRDS("3d_CCA.rds")

DefaultAssay(segmentoid3d) <- "RNA"
Idents(segmentoid3d) <- "celltype"
unique(segmentoid3d@meta.data$celltype)

# Define color scheme
somitecolor <- c("CS9_Late.somite"="#305c99",
                 "CS9_Early.somite"="#6e8db8",
                 "CS9_aPSM"="#8665ae",
                 "CS9_pPSM"="#c999ca",
                 "CS9_NMP-N"="#fb8073",
                 "CS9_NMP-M"="#fba5a3",
                 "CS9_SC-1"="#fdbf70",
                 "CS9_SC-2"="#f6e5a0",
                 "iPSC"="#1a9e76",
                 "NMP"="#b85816",
                 "PSM Post"="#e94749",
                 "PSM Ant"="#ff8000",
                 "Somite"="#8dd3c7",
                 "Neural"="#f681bf")

# Retrieve factor levels
factor_levels <- levels(segmentoid3d@meta.data$celltype)

# Reorder color vector based on factor levels
somitecolor_reordered <- somitecolor[factor_levels]

# Create DimPlot
plot <- DimPlot(segmentoid3d, reduction = "umap", group.by = "celltype") +
  labs(x = "umap_1", y = "umap_2", title = "UMAP", color = "Version") +
  scale_color_manual(values = somitecolor_reordered) +  # Apply custom color scheme
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("celltype")
plot

# Save as PDF
ggsave("umap_plot_3d.pdf", plot=plot, width=6, height=4)





