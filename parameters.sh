module load R/3.6.0
# module load python3/3.6.5 # UMAP
module load pandoc/2.5 # knitr
module load hdf5

# for velocyto:
module load boost
module load openmpi
R
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis
# for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA  GM-KV_CD34P-GPI80P-mRNA; do
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.5', percent.mito= 0.25 ))"  &
sleep 20
done 


for pref in GM-KV_CD34P-GPI80P-mRNA GM-KV_CD34P-Bulk-mRNA; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.5', percent.mito= 0.25 ))"  &
done 

params <- list()
params$prefix <- "GM-KV_CD34P-GPI80P-mRNA"
params$prefix <- "GM-KV_CD34P-Bulk-mRNA"
params$resDEG <- "RNA_snn_res.0.5"
params$percent.mito <- 0.25
params$regress.batch <- TRUE
params$sc.transform <- TRUE

# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d` ; done
# cp htmls
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp $d/../`basename $d`.html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp -R plots
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp -R $d/plots /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp -R $d/plots/feature.plot.all.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/plots/. ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp  $d/plots/*colScaled.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/plots/. ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp  $d/plots/rdim.res.0.5.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/plots/. ; done
# cp xlsx
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp -R $d/*xlsx /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp web_summary and loupe
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp -R $d/outs/*html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp -R $d/outs/*loupe /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp SPRING files
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*/; do cp -R $d/spring.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp final cells to keep
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*mRNA/; do cp -R $d/rds/final.cells.to.keep.Rds /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`_final.cells.to.keep.Rds ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*mRNA/; do cp -R $d/rds/mouse.cells.Rds /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`_mouse.cells.Rds ; done

# merge

Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.html',\
params= list( prefix= 'merged', resDEG= 'RNA_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='FALSE'))" &

Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.regress.html',\
params= list( prefix= 'merged.regress', resDEG= 'RNA_snn_res.0.25', percent.mito= 0.25, regress.batch='TRUE', sc.transform ='FALSE'))" &


# combine sctransform

Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.regress.sctransform.html',\
params= list( prefix= 'merged.regress.sctransform', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='TRUE', sc.transform ='TRUE'))" &


Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.html',\
params= list( prefix= 'merged.sctransform', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

# combine sctransform 2 samples

Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.2samples.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.regress.sctransform.2samples.html',\
params= list( prefix= 'merged.regress.sctransform.2samples', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='TRUE', sc.transform ='TRUE'))" &

# final:
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.2samples.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples.html',\
params= list( prefix= 'merged.sctransform.2samples', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

params <- list()
params$prefix <- "merged.sctransform.2samples"
params$resDEG <- "SCT_snn_res.0.5"
params$percent.mito <- 0.25
params$regress.batch <- FALSE
params$sc.transform <- TRUE

# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples*/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d` ; done
# cp htmls
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples*/; do cp $d/../merged.sctransform.2samples.html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp -R plots
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples*/; do cp -R $d/plots /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp xlsx
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples*/; do cp -R $d/*xlsx /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp SPRING files
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples/; do cp -R $d/spring.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp COMET files
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples/; do cp -R $d/comet.*txt /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp COMET results
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.2samples/; do cp -R $d/COMET_web_version /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done

https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?client_datasets/cite_seq_2sample_sctranform/cite_seq_2sample_sctranform

# exact same cell numbers as individual libraries
https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?client_datasets/cite_seq_2sample_sctranform_cell_num/cite_seq_2sample_sctranform_cell_num 



# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/raw.data ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/raw.data/antibody ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/raw.data/mRNA ; done
# cp h5 mRNA
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/`basename $d`/outs/filtered_feature_bc_matrix.h5 /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/raw.data/mRNA/. ; done
# cp ADT
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/CITE-count-seq//${$(basename $d)%-mRNA}/read_count/* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/raw.data/antibody/. ; done


# for Josh Campbell

for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/raw.data/antibody ; done
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/raw.data/mRNA ; done
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/processed.seurat.ADT.mRNA.object ; done
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/cellranger.mRNA.summary ; done
# cp h5 mRNA
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/`basename $d`/outs/filtered_feature_bc_matrix.h5 /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/raw.data/mRNA/. ; done
# cp ADT
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/CITE-count-seq//${$(basename $d)%-mRNA}/read_count/* /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/raw.data/antibody/. ; done
# cp seurat object
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/`basename $d`/rds/sc.Rds /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/processed.seurat.ADT.mRNA.object/. ; done
# cp cellranger summaries
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/`basename $d`/outs/web_summary.html /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq/`basename $d`/cellranger.mRNA.summary/. ; done

# new updated seurat objects for Josh / Zhe
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq_2020-02-21/`basename $d`_filter_3/processed.seurat.object ; done
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do mkdir -p /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq_2020-02-21/`basename $d`_filter_3_CCregressed/processed.seurat.object ; done
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/`basename $d`_filter_3/rds/sc.Rds /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq_2020-02-21/`basename $d`_filter_3/processed.seurat.object/. ; done
for d in /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/GM-*/; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/`basename $d`_filter_3_CCregressed/rds/sc.Rds /restricted/projectnb/camplab/home/cvmar/19_06_24_kim_citeseq_2020-02-21/`basename $d`_filter_3_CCregressed/processed.seurat.object/. ; done



# filter top 8 in total ADT
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x_filter.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.5', percent.mito= 0.25 ))"  &
sleep 20
done 

# filter top 3 in total ADT and top 5 in ADT CD235a
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x_filter_2.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter_2.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.5', percent.mito= 0.25 ))"  &
sleep 20
done 

# mouse <- Not good. PCA and UMAPs are nonsensical
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x_mouse.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_mouse.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.5', percent.mito= 0.25 ))"  &
sleep 20
done 

# mouse 2 (dim red with human var features + extra scatter plot)
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x_mouse_2.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_mouse_2.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.5', percent.mito= 0.25 ))"  &
sleep 20
done 

# filter_3 = QC filt, doublet filt, mouse filt, outlier filt. Also changed normalization method to SCT
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter_3
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x_filter_3.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter_3/${pref}_filter_3.html',\
params= list( prefix= '$pref', resDEG= 'SCT_snn_res.0.5', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'  ))"  &
sleep 20
done 

params <- list()
params$prefix <- "GM-KV_CD34P-GPI80P-mRNA"
params$resDEG <- "SCT_snn_res.0.5"
params$percent.mito <- 25 # range 0 to 100
params$regress.batch <- FALSE
params$sc.transform <- TRUE

# filter_3_CCregressed = does Cell Cycle regression on the final filtered dataset (filter_3)
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter_3_CCregressed
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x_filter_3_CCregressed.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter_3_CCregressed/${pref}_filter_3_CCregressed.html',\
params= list( prefix= '$pref', resDEG= 'CCregress_res.0.5', percent.mito= 25, regress.batch='FALSE', sc.transform ='TRUE'  ))"  &
sleep 20
done 

params <- list()
params$prefix <- "GM-KV_CD34P-GPI80P-mRNA"
params$resDEG <- "CCregress_res.0.5"
params$percent.mito <- 25 # range 0 to 100
params$regress.batch <- FALSE
params$sc.transform <- TRUE

# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*ssed/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d` ; done
# cp htmls
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*ssed/; do cp $d/`basename $d`.html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp -R plots
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*ssed/; do cp -R $d/plots /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp xlsx
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*ssed/; do cp -R $d/*xlsx /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp plots/velocyto
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*/; do cp -R $d/plots/*velo* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/plots/. ; done
# cp SPRING files
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*ssed/; do cp -R $d/spring.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_filter_3*ssed/; do cp -R $d/spring.custom.color.tracks.csv /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done


# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_mouse_2/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d` ; done
# cp htmls
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_mouse_2/; do cp $d/../`basename $d`.html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp -R plots
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_mouse_2/; do cp -R $d/plots /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp xlsx
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/GM-*_mouse_2/; do cp -R $d/*xlsx /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done



# final:
pref="merged.sctransform.all.after.final.filt"
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.all.after.final.filt.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref/$pref.html',\
params= list( prefix= 'merged.sctransform.all.after.final.filt', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

params <- list()
params$prefix <- "merged.sctransform.all.after.final.filt"  # version of all samples together, with signatures from popescu scored, with final filtered list
params$resDEG <- "SCT_snn_res.0.25"





# final of combination of 2 samples:
pref="merged.sctransform.2samples.after.final.filt"
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.2samples.after.final.filt.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref/$pref.html',\
params= list( prefix= 'merged.sctransform.2samples.after.final.filt', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

params <- list()
params$prefix <- "merged.sctransform.2samples.after.final.filt"  # version of all samples together, with signatures from popescu scored, with final filtered list
params$resDEG <- "SCT_snn_res.0.25"

# final of combination of 2 samples CC regressed:
pref="merged.sctransform.2samples.after.final.filt_CCregressed"
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.2samples.after.final.filt_CCregressed.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref/$pref.html',\
params= list( prefix= '$pref', resDEG= 'CCregress_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

params <- list()
params$prefix <- "merged.sctransform.2samples.after.final.filt_CCregressed"  # version of all samples together, with signatures from popescu scored, with final filtered list
params$resDEG <- "CCregress_res.0.25"

# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*final.filt/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d` ; done
# cp htmls
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*final.filt/; do cp $d/`basename $d`.html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp -R plots
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*final.filt/; do cp -R $d/plots /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp xlsx
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*final.filt/; do cp -R $d/*xlsx /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp SPRING files
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*final.filt/; do cp -R $d/spring.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp plots/velocyto
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*final.filt/; do cp -R $d/plots/*velo* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/plots/. ; done

module load R/3.6.0
module load python3/3.6.5 # UMAP
module load pandoc/2.5 # knitr
# for velocyto:
module load boost
module load openmpi

# run velocyto
for pref in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/velocyto.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${pref}_filter_3_CCregressed/velocyto.html',\
params= list( prefix= '$pref', resDEG= 'SCT_snn_res.0.5'))" &
sleep 20
done

params <- list()
params$prefix <- "GM-KV_CD34N-CD235AN-mRNA"
params$resDEG <- "SCT_snn_res.0.5"



# final of combination of all samples CC regressed:
pref="merged.sctransform.all.after.final.filt_CCregressed"
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_06_24_kim_citeseq/cs.10x.merge.sctransform.all.after.final.filt_CCregressed.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/$pref/$pref.html',\
params= list( prefix= '$pref', resDEG= 'CCregress_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

params <- list()
params$prefix <- "merged.sctransform.all.after.final.filt_CCregressed"  # version of all samples together, with signatures from popescu scored, with final filtered list
params$resDEG <- "CCregress_res.0.25"

# mkdirs
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do mkdir -p /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d` ; done
# cp htmls
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do cp $d/`basename $d`.html /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp -R plots
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do cp -R $d/plots /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp xlsx
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do cp -R $d/*xlsx /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp SPRING files
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do cp -R $d/spring.* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do cp -R $d/spring.groupings.csv /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/. ; done
# cp plots/velocyto
for d in /mnt/scc/project_workspace/19_06_24_kim_citeseq/calculations/analysis/*all.after.final.filt_CCregressed/; do cp -R $d/plots/*velo* /Volumes/GoogleDrive/My\ Drive/CReM/Murphy_Lab/19_06_24_kim_citeseq/`basename $d`/plots/. ; done

# transfer for A.Belkina
for d in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do mkdir -p /projectnb/crem-transfer/19_06_24_kim_citeseq/$d; done
for d in GM-KV_CD34N-CD235AN  GM-KV_CD34N-CD235AP  GM-KV_CD34P-Bulk GM-KV_CD34P-GPI80P; do mkdir -p /projectnb/crem-transfer/19_06_24_kim_citeseq/${d}-ADT; done
for d in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/cellranger_count_all3runs/$d/outs/filtered_feature_bc_matrix.h5 /projectnb/crem-transfer/19_06_24_kim_citeseq/$d/. ; done
for d in GM-KV_CD34N-CD235AN-mRNA  GM-KV_CD34N-CD235AP-mRNA  GM-KV_CD34P-Bulk-mRNA GM-KV_CD34P-GPI80P-mRNA; do cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/${d}_filter_3_CCregressed/rds/sc.Rds  /projectnb/crem-transfer/19_06_24_kim_citeseq/$d/seurat.Rds; done
for d in GM-KV_CD34N-CD235AN  GM-KV_CD34N-CD235AP  GM-KV_CD34P-Bulk GM-KV_CD34P-GPI80P; do cp -R /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/CITE-count-seq/$d/umi_count  /projectnb/crem-transfer/19_06_24_kim_citeseq/${d}-ADT/ADT_umi_count; done

/rprojectnb/camplab/home/cvmar/progenitors.Seurat
cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/progenitors.from.all/rds/sc.Rds .
cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/progenitors.from.all/rds/sc_CD34P_Bulk.Rds .
cp /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/progenitors.from.all/rds/sc_GPI80P.Rds .

# browser
module load python3
pip install --user cellbrowser
export PATH=$PATH:~/.local/bin
cd /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations
cbImportSeurat -i /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/analysis/merged.sctransform.all.after.final.filt_CCregressed/rds/sc.definitive.ucsc.Rds -o /restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/CReM-BU.github.io --htmlDir=/restricted/projectnb/crem-bioinfo/project_workspace/19_06_24_kim_citeseq/calculations/CReM-BU.github.io
cbBuild -o .
cbBuild -h

module load ncftp
ncftp
set passive on
set so-bufsize 33554432
open ftp://geoftp:rebUzyi1@ftp-private.ncbi.nlm.nih.gov
cd uploads/cvillamar_2ngVqg0y
put -R 19_06_24_kim_citeseq