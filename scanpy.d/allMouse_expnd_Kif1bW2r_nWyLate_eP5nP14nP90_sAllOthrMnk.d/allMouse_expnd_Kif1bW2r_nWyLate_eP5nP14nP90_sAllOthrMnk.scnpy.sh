#!/bin/bash -l

#SBATCH -A snic2021-5-436
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 72:00:00
#SBATCH -J s.allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk.scnpy
#SBATCH -o allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk.scnpy.log
#SBATCH -M rackham

module load bioinfo-tools
module load java

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################

#Input h5ad file from that is going to be processed, include raw counts and previous PAGODA embedding for all cells in allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90 samples
inH5ad='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_cnts.h5ad'
#Samples to exclude from inH5ad as they are not in "sAllOthrMnk" dataset
excldSmpls='SS2_19_194|SS2_19_318|SS2_19_320|SS2_19_190|SS2_19_196|SS2_20_009|SS2_20_011|SS2_20_049|SS2_20_055|SS2_20_057|SS2_20_043|SS2_20_051|SS2_20_053'
#output h5ad to write information for a subset of samples excluding $excldSmpls
outH5ad='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk_cnts.h5ad'
#Output prefix
depot_prefix='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk'
#Files to make dictionaries to plot genes of interests and cell types of references
flTypGnsToPlt='mm10_gnIntrstToPlt.txt'
#File of two columns with Tirosh's genes for cell cycle
flTrshCllCcl='melanomaCllCycl_mm10.tsv'
#max_genes is the threshold  maximum number of genes with expression allowed
max_genes='10000'
#max_prcntg_mito is the maximum percentage of mitochondrial genes allowed
max_prcntg_mito='0.1'

######################
#First methods to run#
######################
# If True will run Scanpy's
run_scanpy=TRUE
#If True will write a h5ad file with the basics counts and PAGODA embedding
w_bscs=TRUE

#Excute
python3 li_2022_scnpyAnalysis.py -b=$w_bscs -i=$inH5ad -o=$outH5ad -d=$depot_prefix -p=$flTypGnsToPlt -t=$flTrshCllCcl -S=$run_scanpy -M=$max_prcntg_mito -m=$max_genes -X=$excldSmpls
