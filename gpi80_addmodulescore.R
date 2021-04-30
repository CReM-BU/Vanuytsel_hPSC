
library(Seurat)
library(ggplot2)
library(gridExtra)
library(data.table)

cd34 <- readRDS("../Data/sc_CD34P_Bulk.Rds")
# gpi80 <- readRDS("../Data/sc_GPI80P.Rds")
sc <- readRDS("../Data/sc.Rds")

seurat.markers <- FindMarkers(sc,
    test.use = "MAST",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.10,
    verbose = TRUE,
    group.by = "orig.ident",
    ident.1 = "GM-KV_CD34P-GPI80P-mRNA",
    ident.2 = "GM-KV_CD34P-Bulk-mRNA")
seurat.markers <- seurat.markers[order(seurat.markers[,"avg_logFC"], decreasing = TRUE), ]
write.csv(seurat.markers, file = "deg_gpi80_vs_cd34_pos_only.csv")

seurat.markers2 <- FindMarkers(sc,
    test.use = "MAST",
    only.pos = FALSE,
    min.pct = 0.10,
    logfc.threshold = 0.10,
    verbose = TRUE,
    group.by = "orig.ident",
    ident.1 = "GM-KV_CD34P-GPI80P-mRNA",
    ident.2 = "GM-KV_CD34P-Bulk-mRNA")
seurat.markers2 <- seurat.markers2[order(seurat.markers2[,"avg_logFC"],
    decreasing = TRUE), ]
write.csv(seurat.markers2, file = "deg_gpi80_vs_cd34_all.csv")

gpi80deg <- fread("deg_gpi80_vs_cd34_pos.csv")
gpi80deg <- gpi80deg[order(avg_logFC, decreasing = TRUE), ]
gpi80pro <- gpi80deg[seq(30), V1]

cd34rnaams <- AddModuleScore(cd34, features = list(gpi80pro),
    name = "gpi80_30deg", assay = "SCT")
ams30 <- cd34rnaams@meta.data[, "gpi80_30deg1", drop = FALSE]
fwrite(ams30, file = "gpi80_addmodulescore_in_cd34.csv",
    row.names = TRUE)
