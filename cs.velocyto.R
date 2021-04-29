module load R/3.6.0
module load python3/3.6.5 # UMAP
module load pandoc/2.5 # knitr

# for velocyto:
module load boost
module load openmpi
R

params <- list()
params$prefix <- "GM-KV_CD34P-GPI80P-mRNA"
params$resDEG <- "RNA_snn_res.0.5"
params$percent.mito <- 0.25

library(future)
plan("multiprocess", workers = 14)
options(future.globals.maxSize = 8000 * 1024^2)

ifelse(dir.exists("/mnt/scc/project_workspace/"),
       cremwd <- "/mnt/scc/",
       cremwd <- "/restricted/projectnb/crem-bioinfo/")
master_path <- file.path(cremwd, "/project_workspace/19_06_24_kim_citeseq/calculations/analysis/")
prefix <- params$prefix
calculations_path <- paste0(master_path, prefix, "/")
plots_path <- file.path(calculations_path, "plots/")
rds_path <- file.path(calculations_path, "rds/")
cache_path <- file.path(calculations_path, "cache/") # it must end in a slash

library(Seurat)
# devtools::install_github("velocyto-team/velocyto.R")
library(velocyto.R)
# devtools::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

# If you don't have velocyto's example mouse bone marrow dataset, download with the CURL command
# curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile
# = '~/Downloads/SCG71.loom')
ldat <- ReadVelocity(file = paste0(cremwd, "project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-KV_CD34P-GPI80P-mRNA/velocyto/GM-KV_CD34P-GPI80P-mRNA.loom"), engine = "hdf5r")
sc <- as.Seurat(x = ldat)
final.cells.to.keep <- readRDS(paste0(master_path, prefix, "/rds/final.cells.to.keep.Rds"))
sc <- sc[, paste0("GM-KV_CD34P-GPI80P-mRNA:",final.cells.to.keep,"x")]

sc <- SCTransform(object = sc, assay = "spliced")
sc <- RunPCA(object = sc, verbose = FALSE)
sc <- FindNeighbors(object = sc, dims = 1:20)
sc <- FindClusters(object = sc)
sc <- RunUMAP(object = sc, dims = 1:20)
sc <- RunVelocity(object = sc, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sc)))
names(x = ident.colors) <- levels(x = sc)
cell.colors <- ident.colors[Idents(object = sc)]
names(x = cell.colors) <- colnames(x = sc)
show.velocity.on.embedding.cor(emb = Embeddings(object = sc, reduction = "umap"), vel = Tool(object = sc, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1, n.cores = 14)
# Sys.setenv(DISPLAY="10.48.225.54:114.0")
