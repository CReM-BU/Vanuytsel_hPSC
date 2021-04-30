
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

cd34rna <- readRDS(paste0("./Data/GM-KV_CD34P-Bulk-mRNA_filter",
    "_3_CCregressed/processed.seurat.object/sc.Rds"))


rns <- fread("./Data/cd34_colnames.csv")

classical1 <- fread("./Data/pg_CD90_1_CD49f_1.csv")
classical270 <- fread("./Data/pg_CD90pCD49fp_270cells.csv")
sumide3 <- fread("./Data/pg_CD133_1_gpi80_top3_207cells.csv")
epcr119 <- fread("./Data/pg_CD201_1_119cells.csv")
epcr184 <- fread("./Data/pg_CD201p_1_184cells.csv")

cd34rna <- cd34rna[, rns$rn]
cd34rnarn <- colnames(cd34rna)

idents <- Idents(cd34rna)

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

levels(cd34rna)
a <- Idents(cd34rna)

# DE
# classical270
id2 <- factor(rep("other", ncol(cd34rna)),
    levels = c("classical270", "other"))
id2[cd34rnarn %in% classical270$rn] <- "classical270"
Idents(cd34rna) <- id2
mks2 <- FindMarkers(cd34rna,
    assay = "RNA",
    test.use = "MAST",
    min.pct = 0.10,
    logfc.threshold = 0.10,
    ident.1 = "classical270",
    only.pos = FALSE)

# Sumide3
id3 <- factor(rep("other", ncol(cd34rna)),
    levels = c("sumide3", "other"))
id3[cd34rnarn %in% sumide3$rn] <- "sumide3"
Idents(cd34rna) <- id3
mks3 <- FindMarkers(cd34rna,
    assay = "RNA",
    test.use = "MAST",
    min.pct = 0.10,
    logfc.threshold = 0.10,
    ident.1 = "sumide3",
    only.pos = FALSE)

# EPCR 184 cells
id5 <- factor(rep("other", ncol(cd34rna)),
    levels = c("idepcr184", "other"))
id5[cd34rnarn %in% epcr184$rn] <- "idepcr184"
Idents(cd34rna) <- id5
mks5 <- FindMarkers(cd34rna,
    assay = "RNA",
    test.use = "MAST",
    min.pct = 0.10,
    logfc.threshold = 0.10,
    ident.1 = "idepcr184",
    only.pos = FALSE)

fwrite(mks2, file = "deg_classical270.csv", row.names = TRUE)
fwrite(mks3, file = "deg_sumide3.csv", row.names = TRUE)
fwrite(mks5, file = "deg_epcr184.csv", row.names = TRUE)

