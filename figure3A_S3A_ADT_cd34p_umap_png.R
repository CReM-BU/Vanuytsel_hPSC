
library(Seurat)
library(data.table)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

cd34 <- readRDS("./Data/sc_CD34P_Bulk.Rds")

cd34pcohk <- fread("./Data/cd34p_8735cells_ADT_CLR_subtracted.csv")
cd34hsc <- fread("./Data/CD34_HSC_cell_list_6838.csv")
cd34hsc <- cd34hsc$CD34_HSC

cd34pcoh <- as.matrix(cd34pcohk[, -1])
rownames(cd34pcoh) <- cd34pcohk$rn

cd34pcohk <- cd34pcoh
cd34pcohk <- cd34pcohk[, tstrsplit(cd34hsc, "_")[[1]]]
colnames(cd34pcohk) <- paste0(colnames(cd34pcohk), "_3")

# plot rna
dtkeep <- cd34pcohk

mrnakeep <- cd34@assays$RNA@counts
markertable <- data.table(marker = rownames(dtkeep),
  gene = c("CD34", "CD38", "THY1", "PTPRC", "ITGA6",
    "PROM1", "PROCR", "CDH5", "TEK", "FLT4",
    "CLEC1B", "GYPA", "CXCR4", "HLA-DR-DP-DQ", "Pan-HLA-I",
    "THBD", "ENG", "CD164", "PECAM1"))

mtfound <- markertable[which(markertable$gene %in% rownames(mrnakeep)), ]


cd34all <- readRDS(paste0("./Data/GM-KV_CD34P-Bulk-mRNA_filter",
    "_3_CCregressed/processed.seurat.object/sc.Rds"))
DefaultAssay(cd34all) <- "RNA"

cd34final <- cd34all[, tstrsplit(cd34hsc, "_")[[1]]]
cd34final <- NormalizeData(cd34final)

ce <- cd34@reductions$umap@cell.embeddings
rownames(ce) <- tstrsplit(rownames(ce), "_")[[1]]
cd34final@reductions$umap@cell.embeddings <- ce

#DimPlot(cd34final)

# mRNA
for (i in seq(nrow(mtfound))) {
    p <- FeaturePlot(cd34final, features = mtfound$gene[i],
        cols = c("#f0f0f0", brewer.pal(9,"OrRd")), order = TRUE)
    p <- p & NoAxes()
    p2 <- p & NoLegend() + theme(plot.title = element_blank())
    ggsave(plot = p,
        device = "png",
        dpi = 1200,
        filename = paste0("./cd34p/mRNA/cd34p_umap_rna_",
            mtfound$gene[i], ".png"),
        height = 4,
        width = 4)
    ggsave(plot = p2,
        device = "png",
        dpi = 1200,
        filename = paste0("./cd34p/mRNA/cd34p_umap_rna_",
            mtfound$gene[i], "_noleg.png"),
        height = 4,
        width = 4)
}

cd34[["ADT"]] <- CreateAssayObject(data = cd34pcohk)
DefaultAssay(cd34) <- "ADT"

# ADT
for (i in seq(nrow(markertable))) {
    p1 <- FeaturePlot(cd34, features = markertable$marker[i],
        cols = c("#f0f0f0", brewer.pal(9,"OrRd")), order = TRUE)
    p2 <- FeaturePlot(cd34, features = markertable$marker[i],
        max.cutoff = "q99",
        cols = c("#f0f0f0", brewer.pal(9,"OrRd")), order = TRUE)
    p1 <- p1 & NoAxes()
    p1.2 <- p1 & NoLegend() + theme(plot.title = element_blank())
    p2 <- p2 & NoAxes()
    p2.2 <- p2 & NoLegend() + theme(plot.title = element_blank())

    ggsave(plot = p1,
        device = "png",
        dpi = 1200,
        filename = paste0("./cd34p/ADT/cd34p_umap_adt_",
            markertable$marker[i], ".png"),
        height = 4,
        width = 4)
    ggsave(plot = p1.2,
        device = "png",
        dpi = 1200,
        filename = paste0("./cd34p/ADT/cd34p_umap_adt_",
            markertable$marker[i], "_noleg.png"),
        height = 4,
        width = 4)

    ggsave(plot = p2,
        device = "png",
        dpi = 1200,
        filename = paste0("./cd34p/ADT_q99/cd34p_umap_adt_q99_",
            markertable$marker[i], ".png"),
        height = 4,
        width = 4)
    ggsave(plot = p2.2,
        device = "png",
        dpi = 1200,
        filename = paste0("./cd34p/ADT_q99/cd34p_umap_adt_q99_",
            markertable$marker[i], "_noleg.png"),
        height = 4,
        width = 4)
}
