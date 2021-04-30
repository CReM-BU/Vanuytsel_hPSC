library(VennDiagram)
library(data.table)
library(ggplot2)
library(grDevices)


.themePublication <- function(base_size = 20,
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


mks2 <- fread("./Data/20200930_deg_classical270.csv")
mks3 <- fread("./Data/20200930_deg_sumide3.csv")
mks5 <- fread("./Data/20200930_deg_epcr184.csv")

nrow(mks2) # 551
nrow(mks2[avg_logFC > 0, ]) # 336
nrow(mks2[avg_logFC < 0, ]) # 215

nrow(mks3) # 1468
nrow(mks3[avg_logFC > 0, ]) # 743
nrow(mks3[avg_logFC < 0, ]) # 725

nrow(mks5) # 644
nrow(mks5[avg_logFC > 0, ]) # 398
nrow(mks5[avg_logFC < 0, ]) # 246

length(intersect(mks2[avg_logFC > 0, V1], mks3[avg_logFC > 0, V1]))
# [1] 272
length(intersect(mks2[avg_logFC > 0, V1], mks5[avg_logFC > 0, V1]))
# [1] 198
length(intersect(mks3[avg_logFC > 0, V1], mks5[avg_logFC > 0, V1]))
# [1] 316

cols1 <- c("#187AF0", "#FF12FA", "#58D03D")

# venn 4
# Classical (n = 270) vs Sumide 3% (n = 207) vs EPCR+ (n = 184)

# up-regulated
venn.diagram(x = list(mks2[avg_logFC > 0, V1],
    mks3[avg_logFC > 0, V1],
    mks5[avg_logFC > 0, V1]),
    category.names = c(expression(bold(atop("Classical vs. other", "(336)"))),
        expression(bold(atop("Sumide 3% vs. other", "(743)"))),
        expression(bold(atop("EPCR+ vs. other", "(398)")))),
    filename = 'venn4_Classical270_sumide_epcr184_up.tiff',
    fill = cols1,
    alpha = 0.5,
    output = TRUE,
    # direct.area = TRUE,
    # area.vector = c(8, 70, 505, 0, 107, 34, 2),
    scaled = FALSE,
    euler.d = FALSE,
    lty = 'blank',
    cex = 1.8,
    cat.pos = c(-28, 28, 180),
    cat.dist = c(0.07, 0.07, 0.07),
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans")

# down-regulated
venn.diagram(x = list(mks2[avg_logFC < 0, V1],
    mks3[avg_logFC < 0, V1],
    mks5[avg_logFC < 0, V1]),
    category.names = c(expression(bold(atop("Classical vs. other", "(215)"))),
        expression(bold(atop("Sumide 3% vs. other", "(725)"))),
        expression(bold(atop("EPCR+ vs. other", "(246)")))),
    filename = 'venn4_Classical270_sumide_epcr184_down.tiff',
    fill = cols1,
    alpha = 0.5,
    output = TRUE,
    lty = 'blank',
    cex = 1.8,
    cat.pos = c(-28, 28, 180),
    cat.dist = c(0.07, 0.07, 0.07),
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans")

# all
venn.diagram(x = list(mks2[, V1],
    mks3[, V1],
    mks5[, V1]),
    category.names = c(expression(bold(atop("Classical vs. other", "(551)"))),
        expression(bold(atop("Sumide 3% vs. other", "(1468)"))),
        expression(bold(atop("EPCR+ vs. other", "(644)")))),
    filename = 'venn4_Classical270_sumide_epcr184_all.tiff',
    fill = cols1,
    alpha = 0.5,
    output = TRUE,
    lty = 'blank',
    cex = 1.8,
    cat.pos = c(-28, 28, 180),
    cat.dist = c(0.07, 0.07, 0.07),
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans")


