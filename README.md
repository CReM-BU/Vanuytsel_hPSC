
# 19_06_24_kim_citeseq

## Genome reference

Mixture of human and mouse genomes.
No reporters.

cite-seq: ADT + 10X 3'

ADT (Antibody Derived Tags) libraries have different library indexes than cDNA libraries
ADT use 6mer indexes
10X 3' uses standard 8mer indexes

# Design
https://media.nature.com/original/nature-assets/nmeth/journal/v14/n9/extref/nmeth.4380-S3.pdf
How much coverage for ADT libraries?
Sequencing CITE-seq libraries: We estimate that an average of 100 molecules per ADT per cell is sufficient to achieveuseful information. The numberof readsrequired to obtain 100 molecules depends on the complexity of the sequencing library (e.g. duplication rate).ADT and cDNA sequencing libraries can be pooled at desired proportions. To obtain sufficient read coverage for both libraries we typically sequence ADT libraries in 10% of a lane and cDNA library fraction at 90% of a lane (HiSeq Rapid Run Lane).

https://citeseq.files.wordpress.com/2019/02/cite-seq_190213.pdf
ADT and  cDNA  sequencing  libraries  can  be  pooled  at  desired  proportions  and  sequenced  in   parallel.  We typically  allocate  ~50K  reads  per  cell  for  the  cDNA  library  as  recommended  by  10x  Genomics.  We  typically sequence the ADT library at ~2-5K reads per cell depending on the size of the antibody panel.

# Analysis

Get reference that combines mouse and human
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-and-mm10-3.0.0.tar.gz

all: standard mkfastq demultiplexing
cDNA -> align, count, etc
ADT -> https://hoohm.github.io/CITE-seq-Count/Running-the-script/

Option 1 for demultiplexing:
use mkfastq with the simple csv sheet or (preferred) IEM option (mkfastq supports oligo sequences explicitly as input). Caveat: it will expect a 8mer index, so add AT, which is the following bp in biolegend total seq A (from CITE-seq faq).

Option 2:
use bcl2fasq directly, with a IEM (Illumina Experiment Manager) sample sheet specifying a customized base mask for bcl2fastq

Customizing bcl2fastq parameters:
https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/bcl2fastq-direct
For example, omiting extra bases from reads:
bcl2fastq --use-bases-mask=Y150,I8,Y150

cd /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/ADT/HHW2LBGXB/outs/fastq_path
rename "_S1_" "_S0_" GM*/*gz
rename "_S2_" "_S0_" GM*/*gz
rename "_S3_" "_S0_" GM*/*gz
rename "_S4_" "_S0_" GM*/*gz

# Process

1. mkfastq mRNAs
2. mkfastq ADTs (project == sample name)
3. count mRNAs
4. cat ADT fastq lanes
5. CITE-seq count ADTs
6. multimodal analysis 

cd /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/ADT/H77VWBGXC/outs/fastq_path
mkdir citeseq
cp GM*/*R1* citeseq
cp GM*/*R2* citeseq
cd citeseq
# put sample names in file (to iterate)  
ls -1 *R1*.gz | awk -F '_' '{print $1 "_" $2}' | sort | uniq > ID  
# to cat ADT lanes
for i in `cat ./ID`; do echo $i*_R1_*; done 
for i in `cat ./ID`; do echo $i*_R2_*; done 
for i in `cat ./ID`; do cat $i*_R1_* > ${i}_R1.fastq.gz; done
for i in `cat ./ID`; do cat $i*_R2_* > ${i}_R2.fastq.gz; done 
# to cat ADT runs
cd /restricted/projectnbls/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/ADT/aggr
for i in `cat ./ID`; do echo ../H*/outs/fas*/citeseq/${i}_R1.fastq.gz ; done
for i in `cat ./ID`; do ll ../H*/outs/fas*/citeseq/${i}_R1.fastq.gz ; done
for i in `cat ./ID`; do cat ../H*/outs/fas*/citeseq/${i}_R1.fastq.gz > ${i}_R1.fastq.gz; done
for i in `cat ./ID`; do cat ../H*/outs/fas*/citeseq/${i}_R2.fastq.gz > ${i}_R2.fastq.gz; done

# get wl files
for f in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/* ; do
pre=`basename $f`
zcat /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/$pre/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | cut -f 1 -d "-"  > /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/CITE-count-seq/${pre}.wl
done

# Questions


# Pending:
FeaturePlot: formatted RNA and ADT levels of genes that were targetted in ADTs (list will be supplied by Kim) in combined analysis
COMET results
Add ADT levels to SPRING annotations
Distribution of ADTs re-formatted: histograms w/ limit of x-axis at 1000
Re-draw feature plots with color scale limits at  max.cutoff = "q80"

Re-analysis w/ mouse data:
    what's the background noise of ADTs in mouse?
    How similar is endothelial cluster
POSTPONE: In GPI80P sample, regressing out any cell cycle differences -> keep previous cluster annotation, to be able to visualize them

# For COMET (selection of marker combinations):
module load python3
# pip install COMETSC --user 
/usr3/bustaff/cvmar/.local/bin/Comet -h

/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples
/usr3/bustaff/cvmar/.local/bin/Comet -C 16 comet.expr.txt comet.coord.txt comet.cluster.txt /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples/comet

/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples
/usr3/bustaff/cvmar/.local/bin/Comet -C 16 comet.expr.txt comet.coord.txt comet.cluster.txt /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples/comet_

Hypergeometric marker detection. Finds markers identifying a cluster.
Documentation available at https://hgmd.readthedocs.io/en/latest/index.html

positional arguments:
  marker              Marker file input
  vis                 vis file input
  cluster             Cluster file input
  output_path         the output directory where output files should go

optional arguments:
  -h, --help          show this help message and exit
  -g [G]              Optional Gene list
  -C [C]              Num of cores avail for parallelization
  -X [X]              X argument for XL-mHG
  -L [L]              L argument for XL-mHG
  -Abbrev [ABBREV]    Choose between abbreviated or full 3-gene computation
  -K [K]              K-gene combinations to include
  -Down [DOWN]        Downsample
  -Trim [TRIM]        Trim output files
  -Count [COUNT]      Set to True when count data is being used, for
                      visualizations.
  -tenx [TENX]        Set to True when count data is being used, for
                      visualizations.
  -online [ONLINE]    Set to True for online version.
  -skipvis [SKIPVIS]  Set to True to skip visualizations.

CD34N-CD235AP
GM-KV_CD34N-CD235AN/run_report.yaml:Reads processed: 7523211
GM-KV_CD34N-CD235AN/run_report.yaml:Percentage mapped: 84
GM-KV_CD34N-CD235AN/run_report.yaml:    Expected cells: 8236
GM-KV_CD34N-CD235AP/run_report.yaml:Reads processed: 8337770
GM-KV_CD34N-CD235AP/run_report.yaml:Percentage mapped: 79
GM-KV_CD34N-CD235AP/run_report.yaml:    Expected cells: 5281
GM-KV_CD34P-Bulk/run_report.yaml:Reads processed: 10646298
GM-KV_CD34P-Bulk/run_report.yaml:Percentage mapped: 90
GM-KV_CD34P-Bulk/run_report.yaml:       Expected cells: 9576
GM-KV_CD34P-GPI80P/run_report.yaml:Reads processed: 10664297
GM-KV_CD34P-GPI80P/run_report.yaml:Percentage mapped: 93
GM-KV_CD34P-GPI80P/run_report.yaml:     Expected cells: 8129

cell bc  1-16     UMI 17-26
GTCCACTCATTGCCTC TTAGAGTTTT AA


TCATGCCGTAGCTTGT GGGTGCTGGC CT


GTTTCCTTGACCA NNTAAAAAANNNAAAAAAANNNAAAAAANNNAACAAAANNNAAAAAAANNNAAAAAGTNTAGGGGGNNNNGGGGGGGN

## SPRING plots
cvmar
└── 19_06_24_kim_citeseq
    ├── GM-KV_CD34N-CD235AN-mRNA
    │   ├── cellranger.mRNA.summary
    │   │   └── web_summary.html
    │   ├── processed.seurat.ADT.mRNA.object
    │   │   └── sc.Rds
    │   └── raw.data
    │       ├── antibody
    │       │   ├── barcodes.tsv.gz
    │       │   ├── features.tsv.gz
    │       │   └── matrix.mtx.gz
    │       └── mRNA
    │           └── filtered_feature_bc_matrix.h5
    ├── GM-KV_CD34N-CD235AP-mRNA
    │   ├── cellranger.mRNA.summary
    │   │   └── web_summary.html
    │   ├── processed.seurat.ADT.mRNA.object
    │   │   └── sc.Rds
    │   └── raw.data
    │       ├── antibody
    │       │   ├── barcodes.tsv.gz
    │       │   ├── features.tsv.gz
    │       │   └── matrix.mtx.gz
    │       └── mRNA
    │           └── filtered_feature_bc_matrix.h5
    ├── GM-KV_CD34P-Bulk-mRNA
    │   ├── cellranger.mRNA.summary
    │   │   └── web_summary.html
    │   ├── processed.seurat.ADT.mRNA.object
    │   │   └── sc.Rds
    │   └── raw.data
    │       ├── antibody
    │       │   ├── barcodes.tsv.gz
    │       │   ├── features.tsv.gz
    │       │   └── matrix.mtx.gz
    │       └── mRNA
    │           └── filtered_feature_bc_matrix.h5
    └── GM-KV_CD34P-GPI80P-mRNA
        ├── cellranger.mRNA.summary
        │   └── web_summary.html
        ├── processed.seurat.ADT.mRNA.object
        │   └── sc.Rds
        └── raw.data
            ├── antibody
            │   ├── barcodes.tsv.gz
            │   ├── features.tsv.gz
            │   └── matrix.mtx.gz
            └── mRNA
                └── filtered_feature_bc_matrix.h5


# RunALRA() seurat imputation

# Velocity:
http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/velocity.html
1. get looms for each individual library (velocyto.sh)
2. run velocyto.R for each combination of loom and pre-analyzed individual Seurat objects (velocyto.Rmd and parameters)
3. run velocyto.R for pre-analyzed combined analysis and all individual looms (velocyto.sh on combined data after merge.sctransform...)


# Clean dataset
Sample  Cell_number
CD34N-CD235AN 6899
CD34N-CD235AP 3645
CD34P-Bulk  8739
CD34P-GPI80P  7237


# requests
1. new annot column grouping some clusters as per email in certain samples
2. for analysis of 2samples: import CCregress_res.0.15 column from individual CCregress analyses, then DEG between CD34+bulk (@0.15 res w/o #2,#3,#4,#5; keep only #0 and #1) and GPI-80 (@0.15 res w/o #3,#4)
3. new VlnPlots of  ITGA6, PROCR, PROM1, CD164, MLLT3, HLF, HMGA2, LMO2 for bulk vs gpi80, and also for the subsetted bulk vs gpi80 (see above)

10.48.225.54:152.0


# UMAP.combined.updated.plots.Rmd
    Purpose: 
        + To recreate the old UMAP clustered by sample identity with the new CCRegressed UMAP embeddings and overlay interesting genes onto this UMAP
        + Requested for ISSCR poster presentation (reason for ISSCR in plot names)
    Samples: 
        + All Samples (Combined)
    What is done: 
        + Generate UMAP of combined CCregressed analysis clustered by sample identity with specified colors
        + Overlay HSC.MPP signal over above UMAP
        + Overlay interesting genes over the UMAP
    Path to Plots: /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.all.after.final.filt_CCregressed/plots/
        + UMAP.colored.by.orig.ident.ISSCR.poster.pdf
        + UMAP.colored.by.HSC.MPP.ISSCR.poster.pdf
        + UMAP.specific.feats.all.samples.pdf
    Google Drive Location: 
        + UMAP.colored.by.orig.ident.ISSCR.poster.pdf: https://drive.google.com/drive/folders/1ncVDr400-KbLyR0jxe8WwWEnNgkQWcVa
        + UMAP.colored.by.HSC.MPP.ISSCR.poster.pdf: https://drive.google.com/drive/folders/1ncVDr400-KbLyR0jxe8WwWEnNgkQWcVa
        + UMAP.specific.feats.all.samples.pdf: https://drive.google.com/drive/folders/1GU2FryIKkshxx4Qu6HqlbeBfpcEQynbx 

# ISSCR_updated_plots.Rmd
    Purpose: 
        + Generate plots and analysis requested from Kim for the ISSCR poster presentation*
        + Apply cell lineage labels from the combined analysis to the individual samples and generate a number of different plots for exploration and presentation
        + Also generate visualizations and analysis with the cluster labels created by Kim for the individual samples CD34P_Bulk and GPI80P
        + Visualize data without cells from the high mitochondria cluster
        * Some of the plots have been made obsolete by new requests and different approaches taken in later analysis
    Samples:
        + All Samples (Individually)
    What is done:
        + Get cell lineage labels from combined analysis and apply to the seurat object from the individual analysis
        + Remove cells from high.mito cluster for visualizations
        + Plot UMAPs for each of the individual samples using the cell lineage labels from the combined analysis
        + Plot UMAPS for CD34P_Bulk sample and GPI80P sample using Kim's cluster labels (based on CCregress_res.0.5)
        + Generate annotation slot in seurat objects for CD34P_Bulk and GPI80P samples called annot_20200701 containing Kims cluster 
        + Find DEG for CD34P_Bulk and GPI80P samples with Kim's clustering labels
        + Generate Dotplots with important features for CD34P_Bulk and GPI80P samples
        + Generate Module Scores using top 100 and top 30 DEG from the combined 2 sample analysis (merged.sctransform.2samples.after.final.filt_CCregressed)
            - Plot Heatmaps of CD34P_Bulk and GPI80P samples
        + Select subsets of the CD34P_Bulk and GPI80P samples using the selections already made for bulk.cd34pos.cells and gpi80.cd34pos.cells
                + bulk.cd34pos <- readRDS("/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-KV_CD34P-Bulk-mRNA_filter_3_CCregressed/rds/bulk.cd34pos.Rds")
                + gpi80.cd34pos <- readRDS("/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/                GM-KV_CD34P-GPI80P-mRNA_filter_3_CCregressed/rds/gpi80.cd34pos.Rds")
            - Using the above two populations and the whole of CD34P_Bulk population generate two stacked barplots
                + The first stacked barplot shows the percentage of cells at different cell phases among the three populations
                + The second stacked barplot shows the percentage of cells in each cell lineage (labels from Kim) among the three populations
        + Plot UMAP of the CD34P_Bulk and GPI80P samples pre Cell Cycle regression with cell Phase overlayed

    Path to Plots: /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/
        + GM-KV_CD34N-CD235AN-mRNA_filter_3_CCregressed/plots/
            - UMAP.color.by.cell_lineage_from_combined_sc_CD34N_CD235AN.pdf
        + GM-KV_CD34N-CD235AP-mRNA_filter_3_CCregressed/plots/
            - UMAP.color.by.cell_lineage_from_combined_sc_CD34N_CD235AP.pdf
        + GM-KV_CD34P-Bulk-mRNA_filter_3_CCregressed/plots/
            - UMAP.color.by.cell_lineage_from_combined_sc_CD34P_Bulk.pdf
            - UMAP.color.by.cell_lineage_sc_CD34P_Bulk.pdf
            - Dot.annot.sc_CD34P_Bulk.cell_lineage.pdf
            - heatmap_enrichment_scores_CD34P_Bulk_ISSCR_poster.pdf
            - UMAP.color.by.pre.CCRegress.Phase.Bulk.pdf
        + GM-KV_CD34P-GPI80P-mRNA_filter_3_CCregressed/plots/
            - UMAP.color.by.cell_lineage_from_combined_sc_CD34P_GPI80P.pdf
            - UMAP.color.by.cell_lineage_sc_CD34P_GPI80P.pdf
            - Dot.annot.sc_CD34P_GPI80P.cell_lineage.pdf
            - heatmap_enrichment_scores_CD34P_GPI80P_ISSCR_poster.pdf
            - UMAP.color.by.pre.CCRegress.Phase.GPI80P.pdf
        + merged.sctransform.2samples.after.final.filt_CCregressed/plots/
            - barplot_pct_phase.pdf
            - barplot_pct_lineage.pdf
    Google Drive Location: https://drive.google.com/drive/folders/1ncVDr400-KbLyR0jxe8WwWEnNgkQWcVa 

# progenitors.from.all.Rmd
    Purpose: 
        + Analyze and Visualize the CD34P_Bulk and GPI80P samples from the HSC/progenitor cluster in the combined analysis
    Samples: 
        + CD34P_Bulk and GPI80P
    What is done: 
        + Generate UMAP of the combined analysis with all the progenitor clusters (prog.0 - prog.5) labeled HSC/progenitors
        + From the combined analysis extract the cells that belong to the HSC/progenitor clusters from the combined analysis and either the CD34P_Bulk or GPI80P sample
            - Clusters included are HSC/progenitors (prog.#) and pDC.DC.precursor and neutrophil.myeloid.progenitor
        + Create three seurat objects*
            - sc: contains all of the cells from the HSC/progenitor + pDC.DC.precursor + neutrophil.myeloid.progenitor clusters
            - sc_CD34P_Bulk: contains the CD34P_Bulk cells from the same clusters
            - sc_GPI80P: contains the GPI80P cells from the same clusters
        + Normalize each of the objects using the standard 10X analysis procedure (SCTransform)
        + Reduce dimensions, plot heatmap of top 9 PC
        + Find Clusters and run UMAP 
            - plot UMAPs with orig.ident, cell lineage labels from the combined analysis, and cell phase
        + Run cell-cycle stage classification using the standard 10X analysis procedure
        + Generate molecular signature scores for genes sets in link below:
            - /restricted/projectnb/crem-bioinfo/reference_data/gene_sets/kim.citeseq.20190826.csv
        + Plot UMAP with clustering at all resolutions
        + Test how many cells from each cell lineage label from the combined analysis are in each cluster
        + Generate a UMAP with preselected features highlighted
        + Generate a UMAP of other interesting features to highlight
        + Plot clusters from the combined analysis onto the new UMAPS
        + Perform DEG analysis on all of the clustering resolutions
        + Create annotation with recommendations from Kim (SCT_snn_res.0.5 with cluster 3 and 4 combined into one cluster)
        + Generate stacked barplots:
            - For each object plot the percentage of each cells in each cluster one a shared plot
            - Create similar barplot for cell phase percentage in each object
            - Create barplots of cell phase percentage in each cluster for sc_CD34P_Bulk and sc_GPI80P
        + Generate Dot plots with specific features highlighted
        + Look at LMNA expression in cluster 0 of sc_GPI80P and generate correlation heatmaps for highly variable genes
        + Generate table of top 20 DEG for the selected clusters
        + Generate Heatmaps for the top 20 DEG for the selected clusters
        * The analysis steps for each of the 3 seurat objects is the same unless specifically stated
    Path to Plots: /restricted/projectnb/crem-bioinfo//project_workspace/19_06_24_kim_citeseq/calculations/analysis/progenitors.from.all/plots/
        + UMAP.full.color.by.cell.lineage.pdf (1)
        + elbow.pdf (1)
        + elbow.CD34P_Bulk.pdf (1)
        + elbow.GPI80.pdf (1)
        + heatmap_top_9_PCA_sc.pdf (2)
        + heatmap_top_9_PCA_sc_GPI80P.pdf (2)
        + heatmap_top_9_PCA_sc_CD34P_Bulk.pdf (2)
        + UMAP.color.by.CD34_Bulk.pdf (1)
        + UMAP.color.by.GPI80P.pdf (1)
        + UMAP.color.by.orig.ident.pdf (1)
        + UMAP.color.by.cell_lineage.from.full.pdf (1)
        + UMAP.color.by.Phase.from.full.pdf (2)
        + UMAP.color.by.orig.ident.GPI80P.pdf (1)
        + UMAP.color.by.cell_lineage.from.full.GPI80P.pdf (1)
        + UMAP.color.by.Phase.from.full.GPI80P.pdf (2)
        + UMAP.color.by.orig.ident.CD34P_Bulk.pdf (1)
        + UMAP.color.by.cell_lineage.from.full.CD34P_Bulk.temp.pdf (1)
        + UMAP.color.by.Phase.from.full.CD34P_Bulk.temp.pdf (2)
        + UMAP.clust.pdf (1)
        + UMAP.preselected.feats.pdf (1)
        + UMAP.preselected.feats.GPI80.pdf (1)
        + UMAP.preselected.feats.CD34P_Bulk.pdf (1)
        + UMAP.specific.feats.sc.pdf (2)
        + UMAP.specific.feats.sc_GPI80.pdf (2)
        + UMAP.specific.feats.sc_CD34P.pdf (2)
        + UMAP.clust.from.combined.GPI80P.pdf (1)
        + UMAP.clust.from.combined.CD34P_Bulk.pdf (1)
        + UMAP.color.by.selected_clusters.combined.pdf (1)
        + UMAP.color.by.selected_clusters.CD34P_Bulk.pdf (1)
        + UMAP.color.by.selected_clusters.GPI80P.pdf (1)
        + table_lineage_pct.pdf (2)
        + barplot_pct_lineage.pdf (2)
        + barplot_pct_phase.pdf (2)
        + barplot.pct.phase.by.cluster.sc.pdf (2)
        + barplot.pct.phase.by.cluster.sc_CD34P_Bulk.pdf (2)
        + barplot.pct.phase.by.cluster.sc_GPI80P.pdf (2)
        + Dot.annot.sc_CD34P_Bulk.selected_clusters.pdf (2)
        + Dot.annot.sc_CD34P_GPI80P.selected_clusters.pdf (2)
        + corrplot_GPI80P_c0_LMNA.pdf (2)
        + corrplot_GPI80P_c0_top1000.pdf (2)
        + heatmap_enrichment_scores.pdf (2)
	+ heatmap_enrichment_scores_30_and_100deg.pdf (2)
	+ heatmap_top_20_DEG_selected_clusters_sc.pdf (2)
        + heatmap_top_20_DEG_selected_clusters_sc_GPI80.pdf (2)
        + heatmap_top_20_DEG_selected_clusters_sc_CD34P_Bulk.pdf (2)
    Google Drive Location: 
        + (1): https://drive.google.com/drive/folders/19WT3E9jkWIFBHjOkwI5ZRAdQxXHCmuy5
        + (2): https://drive.google.com/drive/folders/1GU2FryIKkshxx4Qu6HqlbeBfpcEQynbx
        + Other objects: https://drive.google.com/drive/folders/1ZlXIMAH-04R5RV7Sa_HJMwvKczQ0l9na 

# violin_and_barplots.Rmd
    Purpose: 
        + Create Violin and barplots of expression of specific genes for CD34P_Bulk and GPI80P samples to compare expression via violin plots and percent of cells that express each gene via barplots. 
    Samples: 
        + CD34P_Bulk and GPI80P
    What is done: 
        + Remove cells from high mitochondria cluster from combined analysis
        + Calculate the percentage of cells in each sample that expresses the genes of interest
        + Plot violin plots and barplots side by side in a figure for each of the genes of interest
    Path to Plots: /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples.after.final.filt_CCregressed//plots/ISSCR_poster/
        + All plots are in the same folder with no other plots there
    Google Drive Location: https://drive.google.com/drive/folders/1z4nxKGNqHIZ22HcF8mITPX8HhgXgVWLM

# monocle3.combined.Rmd
    Purpose: 
        + Perform trajectory analysis on the combined analysis using Monocle 3
    Samples: 
        + All Samples (Combined)
    What is done: 
        + Load combined seurat object from merged.sctransform.all.after.final.filt_CCregressed analysis
            - Label clusters according to Kims cell lineage labels
            - Remove high mitochondrial cells
        + Generate a monocle object (cds) to be filled with embeddings and information from loaded seurat object
        + Fill monocle object with UMAP embeddings, clustering information, and partition information
            - Create partition based on UMAP select connected clusters for each partition
            - Plot UMAPs confirming clustering and lineage labels transfered properly to monocle object
        + Learn graph: Calculate pseudotime and order the cells in a trajectory based on a root node
            - root node chosen with unbiased function from monocle from within provided category of cells (GM-KV_CD34P-GPI80P-mRNA)
            - Generate UMAPs
        + Create subset monocle objects for the 4 arms of the analysis: Erythroid, Lymphoid B, Lymphoid T/NK, and Myeloid
            - Follow same steps as when creating the full monocle object from the combined analysis: generate object, fill with seurat embeddings and clusters, assign partitions, generate pseudotime and order cells.
            - Some cells needed to be excluded because they were disconnected from their groups of cells and disrupted the trajectory. umap_cell_lineage_arm.pdf vs umap_cell_lineage_arm_post_filter.pdf show the differences.
        + Plot UMAP of branches with pseudotime over the combined analysis UMAP
            - Cells in arm follow pseudotime color gradient, other cells in grey
        + Find the top markers expressed in each cluster of the full monocle object and plot DotPlots of results
            - Also create dotplot of specific set of genes of interest
        + Setup a number of functions that are no longer active in Monocle3, but that are essential for generating the pseudotime heatmaps
        + Find DEG over pseudotime for each of the subset monocle objects and the full monocle object
            - Two methods: 
                + fit_models(cds_subset, model_formula_str = "~Pseudotime", cores = 4 )
                + graph_test(merged.cds, neighbor_graph="principal_graph", cores=4)
            - Find gene modules from DEG and generate heatmaps of DEG modules
                + Heatmaps with columns of pseudotime bins and cell lineages
            - Plot pseudotime heatmaps of top 50 genes (and later of specific genes of interest)
        + For full object plot pseudotime heatmap of genes from the top 10 Principal Components
        + Generate pseudotime heatmaps with specific genes for each of the arms (subset monocle objects)
    Path to Plots: /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.all.after.final.filt_CCregressed/monocle/plots/
        + umap_cell_lineage.pdf (1)
        + umap_cluster.pdf (1)
        + umap_cell_lineage_no_labels.pdf/umap_cell_lineage_no_labels.png (1)
        + umap_partition.pdf (1)
        + umap_trajectory_plot_pseudotime.pdf (1)
        + umap_cell_lineage_erythroid_arm.pdf (1)
        + umap_cell_lineage_erythroid_arm_post_filter.pdf (1)
        + umap_trajectory_plot_pseudotime_erythroid_arm.pdf (1)
        + umap_cell_lineage_lymphoid_B_arm.pdf (1)
        + umap_cell_lineage_lymphoid_B_arm_post_filter.pdf (1)
        + umap_trajectory_plot_pseudotime_lymphoid_B_arm.pdf (1)
        + umap_cell_lineage_lymphoid_T_NK_arm.pdf (1)
        + umap_cell_lineage_lymphoid_T_NK_arm_post_filter.pdf (1)
        + umap_trajectory_plot_pseudotime_lymphoid_T_NK_arm.pdf (1)
        + umap_cell_lineage_myeloid_arm.pdf (1)
        + umap_cell_lineage_myeloid_arm_post_filter.pdf (1)
        + umap_trajectory_plot_pseudotime_myeloid_arm.pdf (1)
        + umap_trajectory_plot_pseudotime_erythroid_over_combined.pdf (1,2)
        + umap_trajectory_plot_pseudotime_lymphoid_B_over_combined.pdf (1,2)
        + umap_trajectory_plot_pseudotime_lymphoid_T_NK_over_combined.pdf (1,2)
        + umap_trajectory_plot_pseudotime_myeloid_over_combined.pdf (1,2)
        + plot_genes_by_group_gene_list.pdf (1)
        + plot_genes_by_group_top_markers_cluster.pdf (1)
        + plot_genes_by_group_top_markers_cluster_top10.pdf (1)
        + plot_genes_by_group_gene_list_all.pdf (1)
        + heatmap_deg_pseudotime_modules_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_bin_fit_models.pdf (1)
        + umap_modules_highlighted_fit_models.pdf (1)
        + heatmap_DEG_pseudotime_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_bin_graph_test.pdf (1)
        + umap_modules_highlighted_graph_test.pdf (1)
        + heatmap_DEG_pseudotime_graph_test.pdf (1)
        + heatmap_topPCA_genes_pseudotime.pdf (1)
        + heatmap_deg_pseudotime_modules_erythroid_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_erythroid_bin_fit_models.pdf (1)
        + heatmap_DEG_pseudotime_erythroid_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_erythroid_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_erythroid_bin_graph_test.pdf (1)
        + heatmap_DEG_pseudotime_erythroid_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_B_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_B_bin_fit_models.pdf (1)
        + heatmap_DEG_pseudotime_lymphoid_B_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_B_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_B_bin_graph_test.pdf (1)
        + heatmap_DEG_pseudotime_lymphoid_B_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_T_NK_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_T_NK_bin_fit_models.pdf (1)
        + heatmap_DEG_pseudotime_lymphoid_T_NK_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_T_NK_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_lymphoid_T_NK_bin_graph_test.pdf (1)
        + heatmap_DEG_pseudotime_lymphoid_T_NK_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_myeloid_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_myeloid_bin_fit_models.pdf (1)
        + heatmap_DEG_pseudotime_myeloid_fit_models.pdf (1)
        + heatmap_deg_pseudotime_modules_myeloid_graph_test.pdf (1)
        + heatmap_deg_pseudotime_modules_myeloid_bin_graph_test.pdf (1)
        + heatmap_DEG_pseudotime_myeloid_graph_test.pdf (1)
        + heatmap_pseudotime_erythroid_genes.pdf (1,2)
        + heatmap_pseudotime_lymphoid_B_genes.pdf (1,2)
        + heatmap_pseudotime_lymphoid_T_NK_genes.pdf (1,2)
        + heatmap_pseudotime_myeloid_genes.pdf (1,2)
    Google Drive Location: 
        + (1): https://drive.google.com/drive/folders/1Gp2q4DVTyyu_6YcvjHs2YKihBIIwzXaI 
        + (2): https://drive.google.com/drive/folders/1ncVDr400-KbLyR0jxe8WwWEnNgkQWcVa 
        + Other Objects: https://drive.google.com/drive/folders/1mPF8XyVD4jLIfBoE6w8yWxc7yq5-9_fh 



