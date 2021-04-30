
library(data.table)
library(Seurat)
library(ggplot2)
library(patchwork)

.themePublication <- function(base_size = 12,
    base_family = "sans") {
    (ggthemes::theme_foundation(base_size = base_size,
        base_family = base_family) +
            ggplot2::theme(plot.title = ggplot2::element_text(
                face = "bold",
                size = ggplot2::rel(1),
                hjust = 0.5),
                text = ggplot2::element_text(),
                panel.background = ggplot2::element_rect(color = NA),
                plot.background = ggplot2::element_rect(color = NA),
                panel.border = ggplot2::element_rect(color = NA),
                axis.title = ggplot2::element_text(
                    face = "bold",
                    size = ggplot2::rel(1)),
                axis.title.y = ggplot2::element_text(angle = 90,
                    vjust = 2),
                axis.title.x = ggplot2::element_text(vjust = -0.2),
                axis.text = ggplot2::element_text(),
                axis.line = ggplot2::element_line(color = "black"),
                axis.ticks = ggplot2::element_line(),
                panel.grid.major = ggplot2::element_line(color = "#f0f0f0"),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(color = NA),
                legend.position = "right",
                legend.direction = "vertical",
                # legend.key.size = ggplot2::unit(0.2, "cm"),
                # legend.margin = ggplot2::margin(0),
                # legend.title = ggplot2::element_text(face = "bold"),
                plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                strip.background = ggplot2::element_rect(
                    color = "#f0f0f0", fill = "#f0f0f0"),
                strip.text = ggplot2::element_text(face = "bold")
            ))
}

cd34rna <- readRDS("./Data/sc_CD34P_Bulk.Rds")

classical270 <- fread("./Data/pg_CD90pCD49fp_270cells.csv")
sumide3 <- fread("./Data/pg_CD133_1_gpi80_top3_207cells.csv")
epcr184 <- fread("./Data/pg_CD201p_1_184cells.csv")

idents <- Idents(cd34rna)

cd34rnarn <- tstrsplit(colnames(cd34rna), "_")[[1]]

classical270id <- rep(NA, ncol(cd34rna))
classical270id[cd34rnarn %in% classical270$rn] <- 1
classical270id[!cd34rnarn %in% classical270$rn] <- 0

idsumide3 <- rep(NA, ncol(cd34rna))
idsumide3[cd34rnarn %in% sumide3$rn] <- 1
idsumide3[!cd34rnarn %in% sumide3$rn] <- 0

idepcr184 <- rep(NA, ncol(cd34rna))
idepcr184[cd34rnarn %in% epcr184$rn] <- 1
idepcr184[!cd34rnarn %in% epcr184$rn] <- 0

cd34rna@meta.data$classical270 <- classical270id
cd34rna@meta.data$idsumide3 <- idsumide3
cd34rna@meta.data$epcr184 <- idepcr184


dtc <- data.table(rn = cd34rnarn)

dtc[rn %in% classical270$rn, classical270 := 1]
dtc[is.na(classical270), classical270 := 0]

dtc[rn %in% sumide3$rn, sumide3 := 1]
dtc[is.na(sumide3), sumide3 := 0]

dtc[rn %in% epcr184$rn, epcr184 := 1]
dtc[is.na(epcr184), epcr184 := 0]

size <- 0.7

Idents(cd34rna) <- classical270id
hl2 <- WhichCells(cd34rna, idents = 1)
g2 <- DimPlot(cd34rna,
    reduction = "umap",
    group.by = "classical270",
    pt.size = size,
    order = TRUE,
    cells.highlight = hl2) +
    scale_color_manual(labels = c("Other", "Classical"),
        values = c("lightgrey", "#187AF0")) +
    ggtitle("CD90+CD49f+ (top 270 cells)")

Idents(cd34rna) <- idsumide3
hl3 <- WhichCells(cd34rna, idents = 1)
g3 <- DimPlot(cd34rna,
    reduction = "umap",
    group.by = "idsumide3",
    pt.size = size,
    order = TRUE,
    cells.highlight = hl3) +
    scale_color_manual(labels = c("Other", "Sumide"),
        values = c("lightgrey", "#FF12FA")) +
    ggtitle("Sumide top3% (n = 207)")

Idents(cd34rna) <- idepcr184
hl5 <- WhichCells(cd34rna, idents = 1)
g5 <- DimPlot(cd34rna,
    reduction = "umap",
    group.by = "epcr184",
    pt.size = size,
    order = TRUE,
    cells.highlight = hl5) +
    scale_color_manual(labels = c("Other", "EPCR+"),
        values = c("lightgrey", "springgreen4")) +
    ggtitle("CD201+CD34+ (n = 184)")

pdf("umap_overlay_w10.pdf", width = 10)
print(g2)
print(g3)
print(g5)
dev.off()

pdf("umap_overlay2_noleg.pdf")
print(g2 + NoLegend())
print(g3 + NoLegend())
print(g5 + NoLegend())
dev.off()
