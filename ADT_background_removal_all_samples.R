
library(Seurat)
library(data.table)
library(Matrix)
library(ggplot2)
library(gridExtra)


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
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(color = NA),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.size = ggplot2::unit(0.2, "cm"),
                legend.margin = ggplot2::margin(0),
                legend.title = ggplot2::element_text(face = "bold"),
                plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                strip.background = ggplot2::element_rect(
                    color = "#f0f0f0", fill = "#f0f0f0"),
                strip.text = ggplot2::element_text(face = "bold")
            ))
}

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


celllist <- readRDS("./Data/new_cell_list/final.cell.list.Rds")
cd235ankeep <- celllist[tstrsplit(celllist, "_")[[2]] == 1]
cd235apkeep <- celllist[tstrsplit(celllist, "_")[[2]] == 2]
cd34pkeep <- celllist[tstrsplit(celllist, "_")[[2]] == 3]
gpi80pkeep <- celllist[tstrsplit(celllist, "_")[[2]] == 4]


#ADT cd235ankeep
hscdt1 <- readMM("./Data/GM-KV_CD34N-CD235AN-mRNA/raw.data/antibody/matrix.mtx.gz")
rn1 <- fread("./Data/GM-KV_CD34N-CD235AN-mRNA/raw.data/antibody/features.tsv.gz", header = FALSE)
cn1 <- fread("./Data/GM-KV_CD34N-CD235AN-mRNA/raw.data/antibody/barcodes.tsv.gz", header = FALSE)
colnames(hscdt1) <- paste0(cn1[[1]], "_1")
rownames(hscdt1) <- rn1[[1]]
hscdtm1 <- as.matrix(hscdt1)
hscdtm1 <- hscdtm1[-20, ]

cd235anmouse <- readRDS("./Data/GM-KV_CD34N-CD235AN-mRNA_mouse.cells.Rds")
cd235anmouse <- paste0(cd235anmouse, "_1")

#ADT cd235apkeep
hscdt2 <- readMM("./Data/GM-KV_CD34N-CD235AP-mRNA/raw.data/antibody/matrix.mtx.gz")
rn2 <- fread("./Data/GM-KV_CD34N-CD235AP-mRNA/raw.data/antibody/features.tsv.gz", header = FALSE)
cn2 <- fread("./Data/GM-KV_CD34N-CD235AP-mRNA/raw.data/antibody/barcodes.tsv.gz", header = FALSE)
colnames(hscdt2) <- paste0(cn2[[1]], "_2")
rownames(hscdt2) <- rn2[[1]]
hscdtm2 <- as.matrix(hscdt2)
hscdtm2 <- hscdtm2[-20, ]

cd235apmouse <- readRDS("./Data/GM-KV_CD34N-CD235AP-mRNA_mouse.cells.Rds")
cd235apmouse <- paste0(cd235apmouse, "_2")


#ADT cd34
hscdt3 <- readMM("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/matrix.mtx.gz")
rn3 <- fread("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/features.tsv.gz", header = FALSE)
cn3 <- fread("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/barcodes.tsv.gz", header = FALSE)
colnames(hscdt3) <- paste0(cn3[[1]], "_3")
rownames(hscdt3) <- rn3[[1]]
hscdtm3 <- as.matrix(hscdt3)
hscdtm3 <- hscdtm3[-20, ]

cd34pmouse <- readRDS("./Data/GM-KV_CD34P-Bulk-mRNA_mouse.cells.Rds")
cd34pmouse <- paste0(cd34pmouse, "_3")


#ADT gpi80
hscdt4 <- readMM("./Data/GM-KV_CD34P-GPI80P-mRNA/raw.data/antibody/matrix.mtx.gz")
rn4 <- fread("./Data/GM-KV_CD34P-GPI80P-mRNA/raw.data/antibody/features.tsv.gz", header = FALSE)
cn4 <- fread("./Data/GM-KV_CD34P-GPI80P-mRNA/raw.data/antibody/barcodes.tsv.gz", header = FALSE)
colnames(hscdt4) <- paste0(cn4[[1]], "_4")
rownames(hscdt4) <- rn4[[1]]
hscdtm4 <- as.matrix(hscdt4)
hscdtm4 <- hscdtm4[-20, ]

gpi80pmouse <- readRDS("./Data/GM-KV_CD34P-GPI80P-mRNA_mouse.cells.Rds")
gpi80pmouse <- paste0(gpi80pmouse, "_4")


mats <- cbind(hscdtm1, hscdtm2, hscdtm3, hscdtm4)

rntemp <- rownames(mats)
rnfinal <- tstrsplit(rntemp, "-")[[1]]
rnfinal[14] <- "HLA-DR-DP-DQ"
rnfinal[15] <- "Pan-HLA-I"
rownames(mats) <- rnfinal

hscdtmn <- NormalizeData(mats, normalization.method = "CLR", margin = 1)

cut <- data.table(ADT = factor(rownames(hscdtmn)),
    cutoff = apply(hscdtmn[,
        c(cd235anmouse, cd235apmouse, cd34pmouse, gpi80pmouse)],
        1, function(x) {mean(x) + sd(x)}))

gdt <- as.data.table(t(hscdtmn), keep.rownames = TRUE)
gdt[rn %in% c(cd235anmouse, cd235apmouse, cd34pmouse, gpi80pmouse), Species := "Mouse"]
gdt[rn %in% c(cd235ankeep, cd235apkeep, cd34pkeep, gpi80pkeep),
    Species := "Human"]
gdt[, Species := factor(Species, levels = c("Mouse", "Human"))]

gdtm <- melt(gdt, id.vars = c("rn", "Species"),
    variable.name = "ADT")

gdtm <- gdtm[order(Species, decreasing = TRUE), ]

colors <- gg_color_hue(2)
bw <- 0.15

g1 <- ggplot() + geom_density(data = gdtm, aes(x = value, color = Species,
    fill = Species), alpha = 0.1, size = 0.5, bw = bw) +
    geom_vline(data = cut,
        mapping = aes(xintercept = cutoff),
        color = "blue", linetype = "dashed", size = 0.5) +
    facet_wrap(~ADT, scales = "free", shrink = FALSE) +
    xlab("CLR transformed ADT counts") +
    xlim(-0.5, NA) +
    ylim(0, 3) +
    .themePublication() +
    scale_color_manual(limits = c("Human", "Mouse"),
        values = c(colors[1], colors[2])) +
    scale_fill_manual(limits = c("Human", "Mouse"),
        values = c(colors[1], colors[2])) + ylab("Density")

pdf("all_samples_adt_background_cutoff_byrow_ylim3_bw015.pdf")
print(g1)
dev.off()

pdf("all_samples_adt_background_cutoff_byrow_ylim3_bw015_w10.pdf",
    width = 10)
print(g1)
dev.off()


# apply cutoff
sub <- copy(cut)
sub[ADT == "CD49f", cutoff := 0]

cd34pco <- hscdtmn - sub[[2]]
cd34pco[cd34pco < 0] <- 0
cd34pcoh <- cd34pco[, c(cd235ankeep, cd235apkeep, cd34pkeep, gpi80pkeep)]

fwrite(as.data.table(cd34pcoh, keep.rownames = TRUE),
    file = "all_samples_19x26407_ADT_CLR_subtracted.csv")

# cluster labels
clusts <- fread("./Data/sample_clust_ref.csv")

normadtdt <- as.data.table(t(cd34pcoh), keep.rownames = TRUE)

normadtdt <- merge(normadtdt, clusts,
    by.x = "rn", by.y = "cn", all = TRUE, sort = FALSE)
colnames(normadtdt)[1] <- "cellID"
setcolorder(normadtdt, c(1, seq(ncol(normadtdt) - 2, ncol(normadtdt)),
    seq(2, ncol(normadtdt) - 3)))

fwrite(normadtdt, file = "all_samples_norm_adt_26407x23.csv")
