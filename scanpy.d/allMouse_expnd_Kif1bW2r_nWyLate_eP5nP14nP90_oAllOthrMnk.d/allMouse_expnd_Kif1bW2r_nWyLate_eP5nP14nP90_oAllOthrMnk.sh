#!/bin/bash -l

#SBATCH -A snic2021-5-436
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 72:00:00
#SBATCH -J s.allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.scnpy
#SBATCH -o allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk.scnpy.log
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
#Samples to exclude from inH5ad as they are not in "oAllOthrMnk" dataset
excldSmpls='SS2_17_343|SS2_17_344|SS2_17_345|SS2_17_346|SS2_17_347|SS2_17_348|SS2_17_368|SS2_17_369|SS2_17_369_2|SS2_17_370|SS2_17_371|SS2_17_371_2|SS2_18_111|SS2_18_112|SS2_18_114|SS2_18_115|SS2_18_116|SS2_20_077|SS2_20_101|SS2_20_103|SS2_20_107|SS2_20_109|SS2_20_111|SS2_20_113|SS2_20_115|SS2_20_117|SS2_20_119|SS2_20_149|SS2_20_155'
#output h5ad to write information for a subset of samples excluding $excldSmpls
outH5ad='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk_cnts.h5ad'
#Output prefix
depot_prefix='allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk'
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

#Execute
python3 li_2022_scnpyAnalysis.py -b=$w_bscs -i=$inH5ad -o=$outH5ad -d=$depot_prefix -p=$flTypGnsToPlt -t=$flTrshCllCcl -S=$run_scanpy -M=$max_prcntg_mito -m=$max_genes -X=$excldSmpls
