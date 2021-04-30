
library(data.table)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggthemes)


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


get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

plothex <- function(x, y, matrix, bins = 200) {
    m <- matrix[c(x, y), ]
    tm <- data.table(t(m), keep.rownames = TRUE)
    g <- ggplot(tm, aes_string(as.name(x), as.name(y))) +
        geom_bin2d(bins = bins) +
        scale_color_continuous(type = "viridis") +
        .themePublication() +
        ggtitle(paste0(x, " vs. ", y))
    return(g)
}


plothex2 <- function(x, y, matrix, n = 100, ...) {
    m <- matrix[c(x, y), ]
    tm <- data.table(t(m), keep.rownames = TRUE)
    tm$density <- get_density(tm[[x]], tm[[y]], n = n, ...)

    g <- ggplot(tm, aes_string(as.name(x), as.name(y), color = "density")) +
        geom_point(size = 1) +
        scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")),
            trans = "sqrt") +
        .themePublication() +
        ggtitle(paste0(x, " vs. ", y))
    return(g)
}

cd34pcohscale <- fread("./Data/cd34p_8735cells_ADT_CLR_subtracted.csv")
cd34m <- as.matrix(cd34pcohscale[, -1])
rownames(cd34m) <- cd34pcohscale$rn
cd34pcohscale <- cd34m

# cd38(rna) vs cd34(adt)
cd34rna <- readRDS(paste0("./Data/GM-KV_CD34P-Bulk-mRNA_filter",
    "_3_CCregressed/processed.seurat.object/sc.Rds"))

rnascale <- cd34rna@assays$RNA@data
rnascale <- rnascale[, colnames(cd34pcohscale)]

cd38 <- t(rnascale["CD38", ])
rownames(cd38) <- "cd38_rna_log"
colnames(cd38) <- tstrsplit(colnames(cd38), "_")[[1]]

cd34pcohscalebind3 <- rbind(cd34pcohscale[, colnames(cd38)], cd38)

cd34cell <- fread("./Data/CD34_HSC_cell_list_6838.csv")
cd34cell <- tstrsplit(cd34cell$CD34_HSC, "_")[[1]]

cd34pcohscalebind3 <- cd34pcohscalebind3[, cd34cell]


x <- "cd38_rna_log"
y <- rownames(cd34pcohscalebind3)
tcd34 <- t(cd34pcohscalebind3)
tcd34 <- data.table(tcd34, keep.rownames = TRUE)
tcd34t <- tcd34[CD34 > 0.5 & CD34 < 2 & cd38_rna_log == 0, ]

#0.5 < cd34 < 2
g1 <- plothex2(x, y[1], cd34pcohscalebind3, n = 500, h = 0.5)
g1 <- g1 + geom_segment(aes(x = 0.05, y = 0.5, xend = 0.05, yend = 2),
    color = "black", size = 1) +
    geom_segment(aes(x = -0.05, y = 0.5, xend = 0.05, yend = 0.5),
        color = "black", size = 1) +
    geom_segment(aes(x = -0.05, y = 2, xend = 0.05, yend = 2),
        color = "black", size = 1) +
    geom_segment(aes(x = -0.05, y = 0.5, xend = -0.05, yend = 2),
        color = "black", size = 1)
g1 <- g1 + ggtitle(paste0("0.5 < CD34 < 2 (n = ",
    nrow(tcd34t),
    ", N = ", ncol(cd34pcohscalebind3), ")")) + xlab("CD38_RNA")


ams30 <- fread("./Data/gpi80_addmodulescore_in_cd34.csv")
ams30$V1 <- tstrsplit(ams30$V1, "_")[[1]]

cd34pcohscalebind <- rbind(cd34pcohscale[, ams30$V1],
    gpi80_30deg1 = ams30$gpi80_30deg1)
cd34pcohscalebind <- cd34pcohscalebind[, tcd34t$rn]


xs <- c("CD49f", "gpi80_30deg1", "CD201")
ys <- c("CD90", "CD133", "CD34")

# CD49f x CD90

dt <- as.data.table(t(cd34pcohscalebind), keep.rownames = TRUE)


dt2 <- dt[CD90 >= 1 & CD49f >= 1, ]
g2 <- plothex2(xs[1], ys[1], cd34pcohscalebind, n = 500)

# 273 double positive cells
min(dt2[, CD90])
min(dt2[, CD49f])

slope1 <- (3.9 - 1) / (3.1 - 1)
dt3 <- dt[CD49f >= (1 + 0.42) & CD90 >= (1 + 0.42 * slope1), ]

xseg3 <- 1.43
yseg3 <- 1.6


g3 <- plothex2(xs[1], ys[1], cd34pcohscalebind, n = 500)
g3 <- g3 + geom_segment(aes(x = xseg3, y = yseg3, xend = 3.1, yend = yseg3),
    color = "black", size = 1) +
    geom_segment(aes(x = xseg3, y = yseg3, xend = xseg3, yend = 3.9),
        color = "black", size = 1) +
    geom_segment(aes(x = xseg3, y = 3.9, xend = 3.1, yend = 3.9),
        color = "black", size = 1) +
    geom_segment(aes(x = 3.1, y = yseg3, xend = 3.1, yend = 3.9),
        color = "black", size = 1)
g3 <- g3 + ggtitle(paste0("CD49f+CD90+ (n = ",
    nrow(dt3), ")"))



# CD133+ GPI-80+
# CD133 > 1 top 3% GPI-80+ 207 cells

dt4 <- dt[CD133 >= 1,
    .(rn, CD133, gpi80_30deg1)][order(gpi80_30deg1, decreasing = T)]

# GPI-80 > 0.516
dt4 <- dt4[gpi80_30deg1 > 0.516]

xseg4 <- 0.516
yseg4 <- 1

g4 <- plothex2(xs[2], ys[2], cd34pcohscalebind, n = 500)
g4 <- g4 + geom_segment(aes(x = xseg4, y = yseg4, xend = 1, yend = yseg4),
    color = "black", size = 1) +
    geom_segment(aes(x = xseg4, y = yseg4, xend = xseg4, yend = 2.5),
        color = "black", size = 1) +
    geom_segment(aes(x = xseg4, y = 2.5, xend = 1, yend = 2.5),
        color = "black", size = 1) +
    geom_segment(aes(x = 1, y = yseg4, xend = 1, yend = 2.5),
        color = "black", size = 1)
g4 <- g4 + ggtitle(paste0("Top 3% GPI-80+ & CD133 > 1 (n = ",
    nrow(dt4), ")"))


# CD34+ CD201+
dt6 <- dt[order(CD201, decreasing = T), .(rn, CD201, CD34)]

# 26 unique values
table(dt6$CD201)

dt61 <- dt6[CD201 > 0.93, ] # 184 cells


xseg6 <- 0.8
yseg6 <- 0.5

g6 <- plothex2(xs[3], ys[3], cd34pcohscalebind, n = 500) + ylim(0, NA)
g6 <- g6 + geom_segment(aes(x = xseg6, y = yseg6, xend = 2.8, yend = yseg6),
    color = "black", size = 1) +
    geom_segment(aes(x = xseg6, y = yseg6, xend = xseg6, yend = 2),
        color = "black", size = 1) +
    geom_segment(aes(x = xseg6, y = 2, xend = 2.8, yend = 2),
        color = "black", size = 1) +
    geom_segment(aes(x = 2.8, y = yseg6, xend = 2.8, yend = 2),
        color = "black", size = 1)
g6 <- g6 + ggtitle(paste0("CD201+CD34+ (n = ",
    nrow(dt61), ")"))


# save plots
pdf("psort_gates_w10.pdf", width = 10)
print(g1)
print(g3)
print(g4)
print(g6)
dev.off()

pdf("psort_gates_noleg.pdf")
print(g1 + theme(legend.position = "none"))
print(g3 + theme(legend.position = "none"))
print(g4 + theme(legend.position = "none"))
print(g6 + theme(legend.position = "none"))
dev.off()

fwrite(dt3, file = "pg_CD90pCD49fp_270cells.csv")
fwrite(dt4, file = "pg_CD133_1_gpi80_top3_207cells.csv")
fwrite(dt61, file = "pg_CD201p_1_184cells.csv")
