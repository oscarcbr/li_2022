#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  li_2022_scnpyAnalysis.py
#  
#  Copyright 2019 oscar <oscar@oscar-J53kOiSr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


##################################
# Work on scanpy-based analysis  #                
##################################

import argparse,os
import scanpy as sc

from numpy import array,inf,float32
from singlecell.li_2022_scnpyAnalysis_mthds import mk_adata_frm_csv,rn_scnpy,addPAGODAtSNE

#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3

#################################
#         Parse inputs          #                
#################################
def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def str2float(v):
	value = float32(v)
	return value
	
def str2int(v):
	value = int(v)
	return value
	


parser = argparse.ArgumentParser()
 
#input parameters  
parser.add_argument('-i', '--inH5ad', default=None,help='input h5ad file from that is going to be processed')
parser.add_argument('-o', '--outH5ad', default=None,help='output h5ad to write information for a subset of samples or cells')
parser.add_argument('-d', '--depot_prefix', default=None,help='is an output preffix to write results')
#Files to make dictionaries
parser.add_argument('-p', '--flTypGnsToPlt', default=None,help='file of two columns with gene and cluster names to plot')
parser.add_argument('-t', '--flTrshCllCcl', default=None,help="file of two columns with Tirosh's genes for cell cycle")
#Filtering parameters
parser.add_argument('-m', '--max_genes', default=inf,type=str2float,help='max_genes is the threshold  maximum number of genes with expression allowed')
parser.add_argument('-M', '--max_prcntg_mito', default=inf,type=str2float,help='max_prcntg_mito is the maximum percentage of mitochondrial genes allowed')
parser.add_argument('-X','--excldSmpls',default=None,help='samples to exclude from the input h5ad file (i.e. inH5ad). Joined multiple samples with "|"  (optional)')
parser.add_argument('-C','--incldCllNms',default=None,help='cells to include from the input h5ad file (i.e. inH5ad). Joined multiple cells with "|" (optional)')
#Write output file with subset cells/samples
parser.add_argument('-b','--w_bscs',default=False, type=str2bool,help='If True will write a h5ad file with the basics counts and embedding for a subset of cells/samples (i.e. outH5ad)')
#Methods to run
parser.add_argument('-S', '--scanpy', default=False, type=str2bool,help='If True will run the scanpy optimized pipelines')
#Test differences between groups using PAGODA embedding
parser.add_argument('-e', '--tstPAGODA_clstr',default=True, type=str2bool,help='Test differences between groups using PAGODA embedding')
parser.add_argument('-f', '--anlyzOnlyDpRgrssdMtchndrl',default=False, type=str2bool,help='If True will also analyze also non-cell cycle corrected results')

args = parser.parse_args()

if args:
	#input paramenters
	inH5ad = args.inH5ad
	outH5ad = args.outH5ad
	depot_prefix = args.depot_prefix
	anlyzOnlyDpRgrssdMtchndrl = args.anlyzOnlyDpRgrssdMtchndrl
	#Files to make dictionaries
	flTypGnsToPlt = args.flTypGnsToPlt
	flTrshCllCcl = args.flTrshCllCcl
	excldSmpls = args.excldSmpls
	incldCllNms = args.incldCllNms
	#Interactive parameters
	max_genes = args.max_genes#inf#Extreme high value, change after exploring figure
	max_prcntg_mito = args.max_prcntg_mito#inf#After exploring figure, change
	#Optional PAGODA embeding results
	w_bscs = args.w_bscs
	#Methods to run
	scnpy = args.scanpy
	#Test differences between groups using PAGODA embedding
	tstPAGODA_clstr = args.tstPAGODA_clstr
	

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################
#Other parameters
s_genes,g2m_genes = [],[]
for l in open(flTrshCllCcl).read().splitlines():
	if l.strip() and l[0]!='#':
		gnNm,cllTyp=l.split('\t')
		if cllTyp=='G1/S':
			s_genes.append(gnNm)
		elif cllTyp=='G2/M':
			g2m_genes.append(gnNm)


lGnsToPlt = [l.split()[0] for l in open(flTypGnsToPlt).read().splitlines() \
if l.strip()]

#Get list of samples to exclude
assert (incldCllNms is not None or excldSmpls is not None)
if excldSmpls is not None:
	excldSmpls = sorted(set([v for v in excldSmpls.split('|') if v.strip()]))
elif incldCllNms is not None:
	incldCllNms = sorted(set([v for v in incldCllNms.split('|') if v.strip()]))


###########################
###########################
###########################
######   Switches   #######
###########################
###########################
###########################
#Write basic info filke
wrtBscs=0
if w_bscs:
	wrtBscs = 1

#Scanpy's 
run_scnpy=0
assert scnpy
if scnpy:
	run_scnpy = 1 #If True will run Scanpy's (2015 == Optimized seurat) pipeline on adrenal gland (i.e. nuc-seq)

##########################
##########################
##########################
######   Execute   #######
##########################
##########################
##########################

####################
####################
#Make adata set
####################
####################
if run_scnpy or wrtBscs:
	#Starting with h5ad file that have raw counts and previously computed tSNE/UMAP embedding
	adata = sc.read(inH5ad)
	#The method "mk_adata_frm_csv" can be used to build "inH5ad" from raw counts.
	#The method "addPAGODAtSNE" can be used to add tSNE in a csv file to "inH5ad".
	if excldSmpls is not None:
		sFllSmpls = set(adata.obs['samples'].unique())
		assert len(excldSmpls)==len(sFllSmpls.intersection(set(excldSmpls)))
		adata = adata[~adata.obs['samples'].isin(excldSmpls)]
	elif incldCllNms is not None:
		sFllCllNms = set(adata.obs_names.unique())
		assert len(incldCllNms)==len(sFllCllNms.intersection(set(incldCllNms)))
		adata = adata[adata.obs_names.isin(incldCllNms)]
	#Write basic shape
	if wrtBscs:
		adata.write(outH5ad)


####################
####################
#Scanpy's 
####################
####################
#Run Scanpy's (2015 == Optimized seurat) pipeline
if run_scnpy:
	dt_nmPrfx_scnpy = '%s_scnpy'%depot_prefix
	print(max_genes,max_prcntg_mito)
	rn_scnpy(adata,dt_nmPrfx_scnpy,max_genes=max_genes, \
	max_prcntg_mito=max_prcntg_mito,geneRefs=lGnsToPlt+['PAGODA_hc'], \
	s_genes=s_genes,g2m_genes=g2m_genes,tstPAGODA_clstr= \
	tstPAGODA_clstr,anlyzOnlyDpRgrssdMtchndrl=anlyzOnlyDpRgrssdMtchndrl)
	
####################

