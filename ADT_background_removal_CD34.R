
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

cd34pkeep <- fread("./Data/cd34_cell_list_8735.csv")
cd34pkeep <- tstrsplit(cd34pkeep[[1]], "_")[[1]]

#ADT
hscdt <- readMM("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/matrix.mtx.gz")
rn <- fread("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/features.tsv.gz",
    header = FALSE)
cn <- fread("./Data/GM-KV_CD34P-Bulk-mRNA/raw.data/antibody/barcodes.tsv.gz",
    header = FALSE)
colnames(hscdt) <- cn[[1]]
rownames(hscdt) <- rn[[1]]
hscdtm <- as.matrix(hscdt)
hscdtm <- hscdtm[-20, ]

cd34pmouse <- readRDS("./Data/GM-KV_CD34P-Bulk-mRNA_mouse.cells.Rds")

rntemp <- rownames(hscdtm)
rnfinal <- tstrsplit(rntemp, "-")[[1]]
rnfinal[14] <- "HLA-DR-DP-DQ"
rnfinal[15] <- "Pan-HLA-I"
rownames(hscdtm) <- rnfinal

# by row
hscdtmn2 <- NormalizeData(hscdtm[, c(cd34pkeep, cd34pmouse)],
    normalization.method = "CLR", margin = 1)

cut2 <- data.table(ADT = factor(rownames(hscdtmn2)),
    cutoff = apply(hscdtmn2[, cd34pmouse], 1, function(x) {mean(x) + sd(x)}))

gdt <- as.data.table(t(hscdtmn2), keep.rownames = TRUE)
gdt[rn %in% cd34pmouse, Species := "Mouse"]
gdt[rn %in% cd34pkeep, Species := "Human"]
gdt[, Species := factor(Species, levels = c("Mouse", "Human"))]


gdtm <- melt(gdt, id.vars = c("rn", "Species"),
    variable.name = "ADT")

colors <- gg_color_hue(2)
bw <- 0.15

g0.5 <- ggplot() + geom_density(data = gdtm,
    aes(x = value, color = Species,
    fill = Species), alpha = 0.1, size = 0.5, bw = bw) +
    geom_vline(data = cut2,
        mapping = aes(xintercept = cutoff),
        color = "blue", linetype = "dashed", size = 0.5) +
    facet_wrap(~ADT, scales = "free", shrink = FALSE) +
    xlab("CLR transformed ADT counts") +
    xlim(-0.5, NA) +
    #ylim(0, 4) +
    .themePublication() +
    scale_color_manual(limits = c("Human", "Mouse"),
        values = c(colors[1], colors[2])) +
    scale_fill_manual(limits = c("Human", "Mouse"),
        values = c(colors[1], colors[2])) + ylab("Density")


g1 <- ggplot() + geom_density(data = gdtm, aes(x = value, color = Species,
    fill = Species), alpha = 0.1, size = 0.5, bw = bw) +
    geom_vline(data = cut2,
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

pdf("cd34p_adt_background_cutoff_byrow_bw015.pdf")
print(g0.5)
dev.off()

pdf("cd34p_adt_background_cutoff_byrow_bw015_w10.pdf", width = 10)
print(g0.5)
dev.off()

pdf("cd34p_adt_background_cutoff_byrow_ylim3_bw015.pdf")
print(g1)
dev.off()

pdf("cd34p_adt_background_cutoff_byrow_ylim3_bw015_w10.pdf",
    width = 10)
print(g1)
dev.off()


# Discussed ADTs

dadt <- c("CD34", "CD38", "CD90", "CD45", "CD49f", "CD133", "CD201",
    "ENG", "CD164")

gdtm2 <- melt(gdt[, c("rn", dadt, "Species"),
    with = FALSE], id.vars = c("rn", "Species"),
    variable.name = "ADT")

g2 <- ggplot() + geom_density(data = gdtm2,
    aes(x = value, color = Species,
        fill = Species), alpha = 0.1, size = 0.5, bw = bw) +
    geom_vline(data = cut2[ADT %in% dadt, ],
        mapping = aes(xintercept = cutoff),
        color = "blue", linetype = "dashed", size = 0.5) +
    facet_wrap(~ADT, scales = "free", shrink = FALSE) +
    xlab("CLR transformed ADT counts") +
    xlim(-0.5, NA) +
    #ylim(0, 4) +
    .themePublication() +
    scale_color_manual(limits = c("Human", "Mouse"),
        values = c(colors[1], colors[2])) +
    scale_fill_manual(limits = c("Human", "Mouse"),
        values = c(colors[1], colors[2])) + ylab("Density")


g3 <- ggplot() + geom_density(data = gdtm2, aes(x = value, color = Species,
    fill = Species), alpha = 0.1, size = 0.5, bw = bw) +
    geom_vline(data = cut2[ADT %in% dadt, ],
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

pdf("cd34p_adt_background_cutoff_byrow_bw015_selected.pdf")
print(g2)
dev.off()

pdf("cd34p_adt_background_cutoff_byrow_bw015_selected_w10.pdf", width = 10)
print(g2)
dev.off()

pdf("cd34p_adt_background_cutoff_byrow_ylim3_bw015_selected.pdf")
print(g3)
dev.off()

pdf("cd34p_adt_background_cutoff_byrow_ylim3_bw015_selected_w10.pdf",
    width = 10)
print(g3)
dev.off()






# apply cutoff

sub <- copy(cut2)
sub[ADT == "CD49f", cutoff := 0]

cd34pco <- hscdtmn2 - sub[[2]]
cd34pco[cd34pco < 0] <- 0
cd34pcoh <- cd34pco[, cd34pkeep]

fwrite(as.data.table(cd34pcoh, keep.rownames = TRUE),
   file = "cd34p_8735cells_ADT_CLR_subtracted.csv")










