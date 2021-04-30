
library(data.table)
library(Matrix)
library(ggplot2)
library(Seurat)


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
                # legend.key.size = ggplot2::unit(0.2, "cm"),
                # legend.margin = ggplot2::margin(0),
                # legend.title = ggplot2::element_text(face = "bold"),
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

n = 2
cols = gg_color_hue(n)


# CD34+
cd34pkeep <- fread("./Data/cd34_cell_list_8735.csv")
cd34pkeep <- cd34pkeep[[1]]

hscdt <- readMM("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/matrix.mtx.gz")
rn <- fread("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/features.tsv.gz", header = FALSE)
cn <- fread("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/barcodes.tsv.gz", header = FALSE)
colnames(hscdt) <- cn[[1]]
rownames(hscdt) <- rn[[1]]
hscdtm <- as.matrix(hscdt)
hscdtm <- hscdtm[-20, ]

cd34pmouse <- readRDS("./Data/GM-KV_CD34P-Bulk-mRNA_mouse.cells.Rds")
cd34pmouse <- paste0(cd34pmouse, "_3")

rntemp <- rownames(hscdtm)
rnfinal <- tstrsplit(rntemp, "-")[[1]]
rnfinal[14] <- "HLA-DR-DP-DQ"
rnfinal[15] <- "Pan-HLA-I"
rownames(hscdtm) <- rnfinal

colnames(hscdtm) <- paste0(colnames(hscdtm), "_3")

# gpi80
gpi80keep <- fread("../gpi80_cell_list_7235.csv")
gpi80keep <- gpi80keep[[1]]

hscdt2 <- readMM("./Data/GM-KV_CD34P-GPI80P-mRNA/raw.data/antibody/matrix.mtx.gz")
rn2 <- fread("./Data/GM-KV_CD34P-GPI80P-mRNA/raw.data/antibody/features.tsv.gz", header = FALSE)
cn2 <- fread("./Data/GM-KV_CD34P-GPI80P-mRNA/raw.data/antibody/barcodes.tsv.gz", header = FALSE)
colnames(hscdt2) <- cn2[[1]]
rownames(hscdt2) <- rn2[[1]]
hscdtm2 <- as.matrix(hscdt2)
hscdtm2 <- hscdtm2[-20, ]

gpi80pmouse <- readRDS("./Data/GM-KV_CD34P-GPI80P-mRNA_mouse.cells.Rds")
gpi80pmouse <- paste0(gpi80pmouse, "_4")

rntemp2 <- rownames(hscdtm2)
rnfinal2 <- tstrsplit(rntemp2, "-")[[1]]
rnfinal2[14] <- "HLA-DR-DP-DQ"
rnfinal2[15] <- "Pan-HLA-I"
rownames(hscdtm2) <- rnfinal2

colnames(hscdtm2) <- paste0(colnames(hscdtm2), "_4")

mat <- cbind(hscdtm[, c(cd34pkeep, cd34pmouse)],
    hscdtm2[, c(gpi80keep, gpi80pmouse)])

hscdtmn <- NormalizeData(mat, normalization.method = "CLR", margin = 1)

cut <- data.table(ADT = factor(rownames(hscdtmn)),
    cutoff = apply(hscdtmn[, c(cd34pmouse, gpi80pmouse)],
        1, function(x) {mean(x) + sd(x)}))

gdt <- as.data.table(t(hscdtmn), keep.rownames = TRUE)
gdt[rn %in% c(cd34pmouse, gpi80pmouse), Species := "Mouse"]
gdt[rn %in% c(cd34pkeep, gpi80keep), Species := "Human"]
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

sub <- copy(cut)
sub[ADT == "CD49f", cutoff := 0]

hscdtmnc <- hscdtmn - sub[[2]]
hscdtmnc[hscdtmnc < 0] <- 0
hscdtmnc <- hscdtmnc[, c(cd34pkeep, gpi80keep)]


hscdtmf <- hscdtmnc
dtr <- data.table(t(hscdtmf), keep.rownames = TRUE)
dtr[rn %in% colnames(hscdtm), sample := "CD34+"]
dtr[rn %in% colnames(hscdtm2), sample := "GPI-80+"]

cd34cluster <- fread("./Data/selected_cluster_cd34.csv")
gpi80cluster <- fread("./Data/selected_cluster_gpi80.csv")
clstsr <- sort(unique(cd34cluster$selected_cluster))


for (i in seq(length(clstsr))) {
    dtr[rn %in% cd34cluster[cd34cluster$selected_cluster ==
            clstsr[i], cn], cluster := clstsr[i]]
    dtr[rn %in% gpi80cluster[gpi80cluster$selected_cluster ==
            clstsr[i], cn], cluster := clstsr[i]]
}

dtfr <- dtr[!is.na(cluster), ]

clsts <- clstsr
# dot plot
avgexp <- vector("numeric", length = nrow(hscdtmn) * length(clsts))

hscdtmfcd34 <- hscdtmf[, cd34cluster$cn]

for (i in seq(length(clsts))) {
    aexp <- apply(hscdtmfcd34[,
        cd34cluster$selected_cluster == clsts[i]], 1, mean)
    avgexp[((i - 1) * nrow(hscdtmfcd34)) + seq(nrow(hscdtmfcd34))] <- aexp
}

nexp <- vector("numeric", length = nrow(hscdtmfcd34) * length(clsts))
for (i in seq(length(clsts))) {
    nexp[((i - 1) * nrow(hscdtmfcd34)) + seq(nrow(hscdtmfcd34))] <-
        apply(hscdtmfcd34[, cd34cluster$selected_cluster == clsts[i]], 1,
            function(x) {sum(x > 0)})
}

medexp <- vector("numeric", length = nrow(hscdtmfcd34) * length(clsts))
for (i in seq(length(clsts))) {
    mexp <- apply(hscdtmfcd34[,
        cd34cluster$selected_cluster == clsts[i]], 1, median)
    medexp[((i - 1) * nrow(hscdtmfcd34)) + seq(nrow(hscdtmfcd34))] <- mexp
}

dt2 <- data.table(ADT = rep(rownames(hscdtmfcd34), length(clsts)),
    cluster = rep(clsts, each = nrow(hscdtmfcd34)),
    number_of_cells = rep(table(cd34cluster$selected_cluster),
        each = nrow(hscdtmfcd34)),
    number_of_expressed = nexp,
    percent_expressed = nexp * 100 /
        rep(table(cd34cluster$selected_cluster),
            each = nrow(hscdtmfcd34)),
    avgexp = avgexp,
    medexp = medexp)


# gpi80
#cd34pcohscale <- ScaleData(cd34pcoh)
clsts2 <- sort(unique(gpi80cluster$selected_cluster))
hscdtmfgpi80 <- hscdtmf[, gpi80cluster$cn]
avgexp2 <- vector("numeric", length = nrow(hscdtmfgpi80) * length(clsts2))

for (i in seq(length(clsts2))) {
    aexp2 <- apply(hscdtmfgpi80[,
        gpi80cluster$selected_cluster == clsts2[i]], 1, mean)
    avgexp2[((i - 1) * nrow(hscdtmfgpi80)) + seq(nrow(hscdtmfgpi80))] <- aexp2
}


nexp2 <- vector("numeric", length = nrow(hscdtmfgpi80) * length(clsts2))
for (i in seq(length(clsts2))) {
    nexp2[((i - 1) * nrow(hscdtmfgpi80)) + seq(nrow(hscdtmfgpi80))] <-
        apply(hscdtmfgpi80[, gpi80cluster$selected_cluster == clsts2[i]], 1,
            function(x) {sum(x > 0)})
}

medexp2 <- vector("numeric", length = nrow(hscdtmfgpi80) * length(clsts2))
for (i in seq(length(clsts2))) {
    mexp2 <- apply(hscdtmfgpi80[,
        gpi80cluster$selected_cluster == clsts2[i]], 1, median)
    medexp2[((i - 1) * nrow(hscdtmfgpi80)) + seq(nrow(hscdtmfgpi80))] <- mexp2
}

dt3 <- data.table(ADT = rep(rownames(hscdtmfgpi80), length(clsts2)),
    cluster = rep(clsts2, each = nrow(hscdtmfgpi80)),
    number_of_cells = rep(table(gpi80cluster$selected_cluster),
        each = nrow(hscdtmfgpi80)),
    number_of_expressed = nexp2,
    percent_expressed = nexp2 * 100 /
        rep(table(gpi80cluster$selected_cluster),
            each = nrow(hscdtmfgpi80)),
    avgexp = avgexp2,
    medexp = medexp2)

dt2[, sample := "CD34+"]
dt3[, sample := "GPI-80+"]

dtf <- rbind(dt2, dt3)

cols = c("#f0f0f0", RColorBrewer::brewer.pal(9,"OrRd"))

#cut(dtf$avgexp, 10)

# avgexp scale
for (i in seq(nrow(hscdtmfcd34))) {
    mclsexp <- mean(dtf[ADT == rownames(hscdtmfcd34)[i], avgexp])
    sd <- sd(dtf[ADT == rownames(hscdtmfcd34)[i], avgexp])
    dtf[ADT == rownames(hscdtmfcd34)[i],
        avgexpscale := scales::rescale(avgexp, c(0, 10))]
}

dtf[, ADT := factor(ADT, levels = unique(ADT))]


mks <- c("CD164", "CD133", "CD34", "CD90", "CD49f",
    "CD201", "CD45", "ENG", "HLA-DR-DP-DQ", "Pan-HLA-I")
dts <- dtf[ADT %in% mks, ]
dts[, ADT := factor(ADT, levels = mks)]

g4 <- ggplot(dts, aes(ADT, cluster)) +
    geom_point(aes(color = avgexpscale, size = percent_expressed)) +
    scale_color_gradientn(colors = cols,
        guide = guide_colorbar(order = 1)) +
    # scale_colour_gradient2(low = "blue4",
    #     high = "firebrick1",
    #     mid = "floralwhite",
    #     midpoint = 5,
    #     breaks = seq(10, 0, -2.5),
    #     limits = c(0, 10),
    #     guide = guide_colorbar(order = 1)) +
    # scale_color_continuous(low = "blue", high = "red",
    #     breaks = seq(-2, 2, 0.5)) +
    scale_size(breaks = seq(100, 0, -25),
        range = c(0, 6),
        limits = c(0, 100),
        guide = guide_legend(order = 2)) +
    facet_wrap(~sample, ncol = 1) +
    .themePublication() +
    theme(axis.text.x = element_text(angle = 45,
        hjust = 1, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        panel.spacing = unit(1, "lines")) +
    labs(color = expression("Rescaled average\nexpression"),
        size = "Percent expressed")

g5 <- ggplot(dts, aes(ADT, cluster)) +
    geom_point(aes(color = avgexpscale, size = percent_expressed)) +
    scale_colour_gradient2(low = "blue4",
        high = "firebrick1",
        mid = "floralwhite",
        midpoint = 5,
        breaks = seq(10, 0, -2.5),
        limits = c(0, 10),
        guide = guide_colorbar(order = 1)) +
    # scale_color_continuous(low = "blue", high = "red",
    #     breaks = seq(-2, 2, 0.5)) +
    scale_size(breaks = seq(100, 0, -25),
        range = c(0, 6),
        limits = c(0, 100),
        guide = guide_legend(order = 2)) +
    facet_wrap(~sample, ncol = 1) +
    .themePublication() +
    theme(axis.text.x = element_text(angle = 45,
        hjust = 1, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        panel.spacing = unit(1, "lines")) +
    labs(color = expression("Rescaled average\nexpression"),
        size = expression("Percent expressed"))

pdf("dotplot_combine_w10.pdf", width = 10)
print(g4)
print(g5)
dev.off()

g5 <- g4 + theme(legend.position = "none")

pdf("dotplot_combine_nolegend_w10.pdf", width = 10)
print(g5)
dev.off()





