#!/bin/bash -l

#SBATCH -A snic2021-5-436
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 36:00:00
#SBATCH -J mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.umap.louvain.2500.2000.50.10.velocity
#SBATCH -o mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.umap.louvain.2500.2000.50.10.velocity.log

module load bioinfo-tools
module load java

##########################
##########################
##########################
####               Variables                 #####
##########################
##########################
##########################

#Inputs for preprocessing
#Master DB file with columns: source, bam_root_fldr, and out_root_fldr
mstrDBfls='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.mstrDB.tsv'
#Repeat masker GTF file. The file "mm10_rmsk.gtf" was obatined from the UCSC browser, repeat masker track for mm10.
rmskFl='mm10_rmsk.gtf'
#Annotation GTF file.
anntGTFFl='mm10.genecodeV18Comp.ERCCYFPiCre.cfflinks.gnNms.biotyp.gtf'
#Output loom merge
outLoomMrgd='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.mrgd.loom'
#Prefix for the project
scPrfx='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.umap.louvain.2500.2000.50.10.velocity'
#Sufix for the file name
sfxNm=''
#Minimum shared counts
min_shared_counts=2500
#Number of top genes to build velocity
n_top_genes=2000
#Number of genes to plot
n_genes_plot=100
#Gene expression space or reduced PCA
n_pcs=50
#kNN neighbors
n_neighbors=10
#Reference genes file
file_genes_ref='frlnTst.gnsRef.mm10.txt'
#Interest genes file
file_gene_interest='frlnTst.gnsIntrst.final.mm10.txt'
#Species
species='Mouse'
#Input counts h5ad file
cntsFl='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk_cnts.h5ad'
#Input embedding h5ad file
emdngFl='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk_scnpy.depNmtchndrlRgrs_ccRgrssn.pca.data.louvain_n_leiden.h5ad'
#Output prefix h5ad file with adata to plot
outFlPrfx='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk_cnts.scnpy.depNmtchndrlRgrs_ccRgrssn.pca.data'
#Output folder to write adata merged file
outFldr='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.d'
#Input embedding h5ad file or #'umap'
embdgNm='umap'
#Cluster name in the embedding file, or #louvain 
clstrNm='louvain'
#Ouput h5ad from counts and embedding files
h5adEmbddngFl='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk_cnts.scnpy.depNmtchndrlRgrs_ccRgrssn.pca.data.umap.louvain.h5ad'
#Ancestral gene
anctrlGn='Erbb3'
######################
#First methods to run#
######################
# If True will run preprocessing
preprocs=FALSE
# If True will run scVelo
scvelo=TRUE
#Build a new h5ad from counts and embedding files
mkNwH5adFrm2h5ads=FALSE


######################
#      Execute       #
######################

python3 li_2022_velocytoAnalysis.py -m=$mstrDBfls -r=$rmskFl -g=$anntGTFFl -l=$outLoomMrgd -p=$preprocs -v=$scvelo -i=$outFlPrfx -s=$scPrfx -M=$n_top_genes -S=$min_shared_counts -P=$n_genes_plot -R=$file_genes_ref -I=$file_gene_interest -e=$species -a=$anctrlGn -n=$mkNwH5adFrm2h5ads -C=$cntsFl -T=$emdngFl -f=$outFldr -H=$h5adEmbddngFl -b=$embdgNm -N=$clstrNm
