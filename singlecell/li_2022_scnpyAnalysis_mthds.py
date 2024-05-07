#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  li_2022_scnpyAnalysis_mthds.py
#  
#  Copyright 2019 Oscar C. Bedoya Reina <oscar.bedoya.reina@ki.se>
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

"""
Group of methods to work with Scanpy. 

Methods developed to run the single-cell analysis included in the section "Cell clustering and expression analysis" 
in Li W et al. manuscript entitled "Chromaffin to neuroblast cell state transitions drive tumor plasticity in NF1 and 
KIF1Bβ deficient neuroblastoma, pheochromocytoma and composite tumors." included in Li W. 2022 thesis 
"Exploring Inter-and Intraheterogeneity in Childhood Neuroblastoma and Pheochromocytoma".
"""

#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3


import matplotlib
matplotlib.use('Agg')
matplotlib.rc('font', family='sans-serif')
matplotlib.rcParams.update({'font.size': 22})
import palantir
import rpy2.robjects.numpy2ri
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

from numpy import array,inf,isinf,nan,where,zeros,float32,isnan
from rpy2 import robjects
from rpy2.robjects.packages import importr


sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
#Other parameters
sc.settings.set_figure_params(dpi=80)

###################
###################
###################
####    Preprocessing   ####
####       methods        ####
###################
###################
###################

def rnStjScnpy(adata,dt_nmPrfx,min_genes=200,min_cells=3,min_mean=0.0125, \
	max_mean=3,min_disp=0.5,clip_max_value=10,max_genes=inf, \
	max_prcntg_mito=inf,copy=False,regrssByNmbrCnts=True, \
	regrssByPrcntgMitchdrl=True):
	"""
	Method to run Satija's (Seurat) pipeline with the improvements of
	scanpy.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add prefix 
	to all output files. min_genes is the threshold minimum number of genes to filter 
	out cells. min_cells is the threshold minimum number of cells to filter out genes. 
	min_mean and max_mean are the minimum and maximum thresholds for mean 
	gene expression to filter out highly variable genes. min_disp is the threshold 
	minimum dispersion of normalized expression to filter genes variable genes. 
	clip_max_value is a threshold to filter genes exceeding a standard deviation of 10. 
	max_genes is the threshold  maximum number of genes with high expression 
	allowed, by scanpy receipt 2500, by Satija's "inf". max_prcntg_mito is the 
	maximum percentage of mitochondrial genes allowed, by scanpy receipt 0.05, by
	Satija's "inf" (Remove cells that have too many mitochondrial genes expressed or
	too many total counts). copy is a switch if so will return a copy of the adata 
	structure with Satija's processing. regrssByNmbrCnts is a switch to regress out 
	effects of total counts  per cell and regrssByPrcntgMitchdrl by the percentage of
	mitochondrial genes expressed (both False in Satija's and True in scanpy).
	Output: adata with updated fields.
	"""
	if copy: adata = adata.copy()
	#Plot top 20 highest expressed genes
	sc.pl.highest_expr_genes(adata,n_top=20,save= \
	'.%s.top20HghstExprssdGns.svg'%dt_nmPrfx)
	#Filter cell and genes with low numbers (basic)
	sc.pp.filter_cells(adata, min_genes=min_genes)
	sc.pp.filter_genes(adata, min_cells=min_cells)
	#Process mitochondrial info
	mito_genes = adata.var_names.str.startswith('MT-')
	if not mito_genes.sum():
		print('Trying to obtain mitochondrial-encoded genes for mouse')
		mito_genes = adata.var_names.str.startswith('mt-')
		if not mito_genes.sum():
			print('Warning: No mitochondrial-encoded genes found')
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, \
	axis=1) / np.sum(adata.X, axis=1)
	adata.obs['n_counts'] = adata.X.sum(axis=1)
	#Plot data info preprocessing
	#1
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], \
	jitter=0.4,multi_panel=True,save='.%s.nmbrGnsNcntsNprnctgMitchndrl.svg' \
	%dt_nmPrfx)
	#2
	sc.pl.scatter(adata, x='n_counts', y='percent_mito', \
	save='.%s.prnctMitchndrlExprssd.svg'%dt_nmPrfx)
	#3
	sc.pl.scatter(adata, x='n_counts', y='n_genes',save='.%s.gnExprssd.svg'% \
	dt_nmPrfx)
	#Filter data for maximum number of genes and mitochondrial
	adata = adata[adata.obs['n_genes'] < max_genes, :]
	adata = adata[adata.obs['percent_mito'] < max_prcntg_mito, :]
	#Total-count normalize (library-size correct) the data matrix X to 
	#10,000 reads per cell, so that counts become comparable among cells.
	sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
	#Logarithmize the data.
	sc.pp.log1p(adata)
	#Set this as raw data
	adata.raw = adata
	#Identify highly-variable genes and plot
	sc.pp.highly_variable_genes(adata, min_mean=min_mean, \
	max_mean=max_mean, min_disp=min_disp)
	sc.pl.highly_variable_genes(adata,save='.%s.hghlyVrblGns.svg'%dt_nmPrfx)
	#Select highly variable genes
	adata = adata[:, adata.var['highly_variable']]
	################################
	#The order of the following two steps might be debatable...
	#This order follow Scanpy's tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
	#Yet might make more sense to first scale and then regress...
	#See: https://stackoverflow.com/questions/77185897/scanpy-single-cell-rna-sequencing-pre-processing-of-individual-vs-integrate
	################################
	#Regress out effects of total counts per cell and the percentage of 
	#mitochondrial genes expressed.
	if regrssByNmbrCnts and regrssByPrcntgMitchdrl:
		sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
	elif regrssByNmbrCnts:
		sc.pp.regress_out(adata, ['n_counts'])
	elif regrssByPrcntgMitchdrl:
		sc.pp.regress_out(adata, ['percent_mito'])
	#Scale each gene to unit variance. Clip values exceeding standard 
	#deviation "clip_max_value".
	sc.pp.scale(adata, max_value=clip_max_value)
	return adata if copy else None


###################
###################
###################
####   Variation and       ####
####  denoise methods  ####
###################
###################
###################
def rn_vrntn_estmtn(adata,dt_nmPrfx,geneRefs=['TH','PNMT','MYCN','PRRX1', \
	'CLDN11'],copy=False,cSeed=1234,keysToPlot=['louvain','leiden', \
	'PAGODA_hc','stage','outcome','samples','S_score','G2M_score', \
	'phase']):
	"""
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data 
	to add prefix to all output files. geneRef is a list of genes to 
	plot in pca plots. copy is a switch if so will return a copy of the 
	adata structure. cSeed is a seed to replicate results.  keysToPlot is
	a list of variables of 'obs' in adata to plot in the embedding.
	Output: Plots with variance estimation and a PCA objtect 
	calculation.
	"""
	if copy: adata = adata.copy()
	#Reduce the dimensionality of the data by running PCA, which reveals
	#the main axes of variation and denoises the data
	sc.tl.pca(adata, svd_solver='arpack',random_state=cSeed)
	#List of factors to plot
	lPltInsrts=[]
	for key in keysToPlot:
		if key in adata.obs and len(adata.obs[key].unique(). \
		tolist())>1:
			lPltInsrts.append(key)
	lPltInsrts.extend([g for g in geneRefs if g in set(adata.var_names)])
	#Scatter plot in the PCA coordinates
	sc.pl.pca(adata, color=lPltInsrts,save='.%s.pca_2Cmpnts_gnRef.svg'% \
	dt_nmPrfx)
	#Inspect the contribution of single PCs to the total variance in the
	#data. This gives information about how many PCs should be
	#considered in order to compute the neighborhood relations of cells,
	#e.g. used in the clustering function sc.tl.louvain() or tSNE 
	#sc.tl.tsne(). In author's experience, often, a rough estimate of 
	#the number of PCs does fine.
	sc.pl.pca_variance_ratio(adata,log=True,save='.%s.pca_cmpnts_logVrnc.svg'% \
	dt_nmPrfx)
	return adata if copy else None


###################
###################
###################
####   Embedding the   ####
### neighborhood graph  ###
###################
###################
###################
def cmpt_nghbrnhd_grph(adata,dt_nmPrfx,mthd='pca',n_pcs=40,copy=False, \
	cSeed=1234):
	"""
	Method to compute neighbor graph.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data  to add prefix 
	to all output files. mthd is the method to compute the neighborhood graph, the 
	options are 'pca' and 'diffmap'. n_pcs is the number of components to calculate 
	the 'pca' graph. copy is a switch if so will return a copy of the adata structure. 
	cSeed is a seed to replicate results.
	Output: It returns the graph computed as part of the adata structure.
	"""
	assert mthd in {'pca','diffmap'}
	if copy: adata = adata.copy()
	#Compute the neighborhood graph of cells using the PCA 
	#representation of the data matrix. You might simply use default 
	#values here. For the sake of reproducing Seurat’s results, let’s 
	#take the following values.
	if mthd=='pca':
		sc.pp.neighbors(adata, n_neighbors=10,n_pcs=n_pcs,random_state= \
		cSeed)
	#Denoise using few of the first spectral components.
	else:
		sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs,random_state= \
		cSeed)
		sc.tl.diffmap(adata)
		sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap', \
		random_state=cSeed)
		sc.tl.draw_graph(adata,save='.%s.dffmap.svg'%dt_nmPrfx, \
		random_state=cSeed)
	return adata if copy else None

def cmpt_embedding(adata,dt_nmPrfx,geneRefs=['TH','PNMT','MYCN','PRRX1', \
	'CLDN11'],crtrEmbdng='both',copy=False,pltRdctns=False, \
	cSeed=1234,keysToPlot=['louvain','leiden','PAGODA_hc','stage', \
	'outcome','samples','S_score','G2M_score','phase']):
	"""
	Method to compute embeddings.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add prefix 
	to all output files.  geneRef is a list of genes to plot in intermediate embedding 
	plots. crtrEmbdng is the algorithm to make the embedding, the options are 
	umap (recommended), tsne, and both. copy is a switch if so will return a copy of
	the adata  structure. pltRdctns is a switch, if True will save the plots. cSeed is a 
	seed to replicate results.  keysToPlot is a list of variables of 'obs' in adata to plot
	in the embedding.
	Output: Includes the embedding and plots them.
	"""
	assert crtrEmbdng in {'umap','tsne','both'}
	if copy: adata = adata.copy()
	#Calculate embedding. The .raw attribute of adata, the previous 
	#plots showed the “raw” (normalized, logarithmized, but uncorrected) 
	#gene expression. 
	lPltInsrts=[]
	for key in keysToPlot:
		if key in adata.obs and len(adata.obs[key].unique(). \
		tolist())>1:
			lPltInsrts.append(key)
	lPltInsrts.extend([g for g in geneRefs if g in set(adata.var_names)])
	#
	if crtrEmbdng=='tsne' or crtrEmbdng=='both':
		sc.tl.tsne(adata,random_state=cSeed)
		#1. Plot normalized,logarithmized, & uncorrected gene expression
		sc.pl.tsne(adata,color=lPltInsrts)#,save='.%s.tsne.raw.svg'%dt_nmPrfx)
		if pltRdctns:
			#2. Plot corrected gene expression
			sc.pl.tsne(adata,color=lPltInsrts,use_raw=False,save= \
			'.%s.tsne.svg'%dt_nmPrfx)
	if crtrEmbdng=='umap' or crtrEmbdng=='both':
		sc.tl.umap(adata)
		#1. Plot normalized,logarithmized, & uncorrected gene expression
		sc.pl.umap(adata,color=lPltInsrts)#,save='.%s.umap.raw.svg'%dt_nmPrfx)
		if pltRdctns:
			#2. Plot corrected gene expression
			sc.pl.umap(adata, color=lPltInsrts,use_raw=False,save= \
			'.%s.umap.svg'%dt_nmPrfx)
	return adata if copy else None

def paga_embdg(adata,dt_nmPrfx,geneRefs=['TH','PNMT','MYCN','PRRX1', \
	'CLDN11'],copy=False,cSeed=1234,mthdsToPlot=['louvain','leiden', \
	'PAGODA_hc'],keysToPlot=['louvain','leiden','PAGODA_hc','stage', \
	'outcome','samples','S_score','G2M_score','phase']):
	"""
	Method to compute PAGA embedding.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add 
	prefix to all output files.  geneRef is a list of genes to plot in intermediate 
	embedding plots.  copy is a switch if so will return a copy of the adata 
	structure. cSeed is a seed to replicate results.  keysToPlot is a list of 
	variables of 'obs' in adata to plot in the embedding. mthdsToPlot is a list 
	of embedding methods to run paga.
	Output: Includes the embedding and plots them.
	"""
	if copy: adata = adata.copy()
	# In some ocassions, it might still observe disconnected clusters 
	#and similar connectivity violations. They can usually be remedied 
	#by running:
	#1 Make keys to run PAGA
	lPltInsrts=[]
	for key in keysToPlot:
		if key in adata.obs and len(adata.obs[key].unique(). \
		tolist())>1:
			lPltInsrts.append(key)
	lPltInsrts.extend([g for g in geneRefs if g in set(adata.var_names)])
	#Calculate PAGA for each method
	for mthd in mthdsToPlot:
		if mthd in adata.obs and len(adata.obs[mthd].unique().tolist())>1:
			sc.tl.paga(adata,groups=mthd)
			sc.pl.paga(adata, plot=True,save='.%s.%s.dffmap.svg'% \
			(dt_nmPrfx,mthd))
			# remove `plot=False` if you want to see the coarse-grained graph
			# as next...
			#2
			#~ plt.figure()
			#~ sc.pl.paga(adata)
			#~ plt.savefig('%s.dffmap.coarseGrn.svg'%dt_nmPrfx)
			#~ plt.close()
			#3
			sc.tl.draw_graph(adata,init_pos='paga',random_state=cSeed)
			sc.pl.draw_graph(adata,color=lPltInsrts,legend_loc='on data', \
			save='.%s.%s.dffmap.rcmptd.svg'%(dt_nmPrfx,mthd))
	return adata if copy else None

def prcss_palantir(adata,geneRefs,dt_nmPrfx_rgrs,strtMrkrs=['NANOG', \
	'PROM1','SOX10','LGR5','POU5F1','SOX2','NOTCH1','NOTCH3','WNT1'], \
	plt_tsne=False,plt_tsne_by_cell_sizes=False,plt_gene_expression=True, \
	plt_diffusion_components=False,keysToPlot=['louvain','leiden', \
	'PAGODA_hc','stage','outcome','samples','S_score','G2M_score','phase']):
	"""
	Wrapper to run palantir, all the way. adata is an AnnData object. geneRef 
	is a list of genes to plot in intermediate embedding plots.  dt_nmPrfx 
	is a name of the data to add prefix to all output files. strtMrkrs is list of
	genes expected to be expressed in the most "ancestral" cell. plt_tsne is
	a switch to plot palantir tSNE. plt_tsne_by_cell_sizes is a switch to plot 
	palantir tSNE by cell sizes. plt_gene_expression is a switch to plot  gene
	expression. plt_diffusion_components is a switch to plot the difussion
	components. keysToPlot is a list of variables of 'obs' in adata to plot in the
	embedding. 
	Output: Includes the embedding and plots them.
	"""
	gnNms = set(adata.var_names).intersection(set(strtMrkrs))
	#Variables to plot
	lPltInsrts=[]
	for key in keysToPlot:
		if key in adata.obs and len(adata.obs[key].unique(). \
		tolist())>1:
			lPltInsrts.append(key)
	lPltInsrts.extend([g for g in geneRefs if g in set(adata.var_names)])
	#
	while gnNms:
		strtMrkr = gnNms.pop()
		sce.tl.palantir(adata)
		#select starting cell
		allVls = adata.X[:,where(adata.var_names==strtMrkr)[0]]
		strtngCll = adata.obs_names[where(allVls==max(allVls))[0]][0]
		#Run Palantir
		if plt_tsne:
			fig, ax = palantir.plot.plot_tsne(adata.uns['palantir_tsne'])
			plt.savefig('figures/%s.palantir.tsne.%s.svg'%(dt_nmPrfx_rgrs, \
			strtMrkr))
			plt.close()
		if plt_tsne_by_cell_sizes:
			fig, ax = palantir.plot.plot_tsne_by_cell_sizes(adata.uns \
			['palantir_imp_df'],adata.uns['palantir_tsne'])
			plt.savefig('figures/%s.palantir.tsneCllSzs.%s.svg'%(dt_nmPrfx_rgrs, \
			strtMrkr))
			plt.close()
		#Gene expression on tSNE
		if plt_gene_expression:
			palantir.plot.plot_gene_expression(adata.uns['palantir_imp_df'], \
			adata.uns['palantir_tsne'],lPltInsrts)
			plt.savefig('figures/%s.palantir.gnExprssn.tsne.%s.svg'% \
			(dt_nmPrfx_rgrs,strtMrkr))
			plt.close()
		#Diffusion maps
		if plt_diffusion_components:
			palantir.plot.plot_diffusion_components(adata.uns \
			['palantir_tsne'],adata.uns['palantir_diff_maps'])
			plt.savefig('figures/%s.palantir.dffsnMap.%s.svg'% \
			(dt_nmPrfx_rgrs,strtMrkr))
			plt.close()
		#Visualize results
		print('Starting cell %s'%strtngCll)
		try:
			pr_res = palantir.core.run_palantir(adata.uns['palantir_ms_data'], \
			strtngCll, num_waypoints=500)
			palantir.plot.plot_palantir_results(pr_res, \
			adata.uns['palantir_tsne'])
			plt.savefig('figures/%s.palantir.dffsnMap.rslt.%s.svg'% \
			(dt_nmPrfx_rgrs,strtMrkr))
			plt.close()
		except:
			pass
	return 0


##################
##################
##################
####      Wrappers     ####
##################
##################
##################
def wpr_dnstrm1(adata,dt_nmPrfx_rgrs,geneRefs=['TH','PNMT','MYCN','PRRX1', \
	'CLDN11'],adata_pca=None,adata_diffmap=None,adata_phate=None, \
	adata_palantir=None,cSeed=1234,mkPaga=True,fllAnalysis=True, \
	keysToPlot=['louvain','leiden','PAGODA_hc','stage','outcome','samples', \
	'S_score','G2M_score','phase'],tstPAGODA_clstr=True,lFtrsAnls=['outcome', \
	'samples','stage']):
	"""
	Wrapper to run the core of the analysis.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add prefix
	to all output files. geneRef is a list of genes to plot in intermediate embedding 
	plots. adata_pca is an adata computed using pca. adata_diffmap is an adata 
	computed using diffmap. adata_phate is an adata compute with phate. 
	adata_palantir is an adata computed with palantir. cSeed is a seed to replicate 
	results.  If mkPaga the PAGA embedding will be plot. fllAnalysis is a switch to run
	all the statistical comparisons for clusters, for every annotation of interest, and 
	include additional analysis. keysToPlot is a list of variables of 'obs' in adata to 
	plot in the embedding.  tstPAGODA_clstr will run statistical tests in the 
	"PAGODA_hc" clusters included in the ".obs" object. lFtrsAnls is a list of features 
	of interests for which the statistical differences will be computed as well.
	Output: Includes the embedding and plots them, and output folders with statistics
	for gene expression for different groups of cells in embeddings and lFtrsAnls.
	"""
	#Embedding
	if adata_pca is None:
		adata_pca = cmpt_nghbrnhd_grph(adata,dt_nmPrfx_rgrs,mthd='pca', \
		n_pcs=40,copy=True,cSeed=cSeed)
	cmpt_embedding(adata_pca,dt_nmPrfx_rgrs,geneRefs=geneRefs, \
	crtrEmbdng='both',copy=False,cSeed=cSeed)
	if adata_diffmap is None:
		adata_diffmap = cmpt_nghbrnhd_grph(adata,dt_nmPrfx_rgrs,mthd= \
		'diffmap',n_pcs=40,copy=True,cSeed=cSeed)
	cmpt_embedding(adata_diffmap,dt_nmPrfx_rgrs,geneRefs=geneRefs, \
	crtrEmbdng='both',copy=False,cSeed=cSeed)
	if adata_phate is None:
		sce.tl.phate(adata,random_state=cSeed)
		sce.pl.phate(adata,color=[g for g in geneRefs if g in \
		set(adata.var_names).union(set(['PAGODA_hc','stage','outcome', \
		'samples']))],save='.%s.phate.svg'%dt_nmPrfx_rgrs)
	if adata_palantir is None:
		prcss_palantir(adata,geneRefs,dt_nmPrfx_rgrs)
	#
	for adata_embd,addPfrx_embd in ((adata_pca,'pca'),(adata_diffmap, \
	'diffmap')):
		dt_nmPrfx_rgrs_embd = '%s.%s'%(dt_nmPrfx_rgrs,addPfrx_embd)
		#Clustering
		sc.tl.louvain(adata_embd,random_state=cSeed)
		sc.tl.leiden(adata_embd,random_state=cSeed)
		#List of things to plot
		lPltInsrts=[]
		for key in keysToPlot:
			if key in adata_embd.obs and len(adata_embd.obs[key].unique(). \
			tolist())>1:
				lPltInsrts.append(key)
		lPltInsrts.extend([g for g in geneRefs if g in set(adata_embd. \
		var_names)])
		# ~ [g for g in geneRefs if g in set(adata.var_names).union(set(['PAGODA_hc', \
		# ~ 'stage','outcome','samples']))]
		#Plot Raw
		if adata_embd.raw is not None:
			sc.pl.umap(adata_embd,color=lPltInsrts,use_raw=True,save= \
			'.%s.umap.clusterng.raw.svg'%dt_nmPrfx_rgrs_embd)
			sc.pl.tsne(adata_embd,color=lPltInsrts,use_raw=True,save= \
			'.%s.tsne.clusterng.raw.svg'%dt_nmPrfx_rgrs_embd)
		#Plot transformed
		sc.pl.umap(adata_embd,color=lPltInsrts,use_raw=False,save= \
		'.%s.umap.clusterng.svg'%dt_nmPrfx_rgrs_embd)
		sc.pl.tsne(adata_embd,color=lPltInsrts,use_raw=False,save= \
		'.%s.tsne.clusterng.svg'%dt_nmPrfx_rgrs_embd)
		#Run the full analysis
		if fllAnalysis:
			#Write data
			if 'palantir_pca_results' in adata_embd.uns:
				toDel=adata_embd.uns.pop('palantir_pca_results')
				toDel=adata_embd.uns.pop('palantir_diff_maps')
				toDel=adata_embd.uns.pop('palantir_ms_data')
				toDel=adata_embd.uns.pop('palantir_imp_df')
				print('Warning: part of "palantir_pca_results" and "palantir_diff_maps" were removed from file %s.data.h5ad and would have to be built from scratch'% \
				dt_nmPrfx_rgrs_embd)
			#write data in case is not written later by paga
			# ~ if not mkPaga:
				# ~ adata_embd.write('%s.data.h5ad'%dt_nmPrfx_rgrs_embd)
			#
			for mthd,addPrfx_exprssd,statsc in (('wilcoxon','wlcnx','pvals_adj'), \
			('t-test_overestim_var','ttst_ovar','pvals_adj'),('logreg', \
			'logreg','pvals'),('t-test','welch','pvals_adj')):#updated 23/09/2021
				dt_nmPrfx_rgrs_embd_exprssd = '%s.%s'%(dt_nmPrfx_rgrs_embd, \
				addPrfx_exprssd)
				#Finding marker genes
				#1
				adata_embd = sc.tl.rank_genes_groups(adata_embd, \
				'louvain',method=mthd,copy=True,n_genes=4000000)
				pltNwrtGnGrps(adata_embd,mthd,'%s.%s'% \
				(dt_nmPrfx_rgrs_embd_exprssd,'louvain'))
				#cluster vs. all others
				all_grps=list(adata_embd.uns['rank_genes_groups']['names']. \
				dtype.names)
				for ref in all_grps:
					sc.tl.rank_genes_groups(adata_embd,'louvain',groups= \
					all_grps,reference=ref,method=mthd,n_genes=4000000)
					pltNwrtGnGrpsVEchOthr(adata_embd,mthd, \
					dt_nmPrfx_rgrs_embd_exprssd,ref,'louvain')
				#2
				adata_embd = sc.tl.rank_genes_groups(adata_embd, \
				'leiden',method=mthd,copy=True,n_genes=4000000)
				pltNwrtGnGrps(adata_embd,mthd,'%s.%s'% \
				(dt_nmPrfx_rgrs_embd_exprssd,'leiden'))
				#cluster vs. all others
				all_grps=list(adata_embd.uns['rank_genes_groups']['names']. \
				dtype.names)
				for ref in all_grps:
					sc.tl.rank_genes_groups(adata_embd,'leiden',groups= \
					all_grps,reference=ref,method=mthd,n_genes=4000000)
					pltNwrtGnGrpsVEchOthr(adata_embd,mthd, \
					dt_nmPrfx_rgrs_embd_exprssd,ref,'leiden')
				#3
				if 'PAGODA_hc' in set(adata_embd.obs_keys()) and tstPAGODA_clstr:
					##############
					#next line added to confirm there is enough samples in each group
					clstrINPagoda = list(adata_embd.obs['PAGODA_hc'].cat.categories)
					clstrINPagodaMoreThn1Smpl = [c for c in clstrINPagoda \
					if len(adata_embd.obs['PAGODA_hc'][adata_embd.obs['PAGODA_hc']==c])>1]
					##############
					if len(clstrINPagodaMoreThn1Smpl)==len(clstrINPagoda) \
					and len(clstrINPagodaMoreThn1Smpl)>1:
						adata_embd = sc.tl.rank_genes_groups(adata_embd, \
						'PAGODA_hc',method=mthd,copy=True,n_genes=4000000)
						pltNwrtGnGrps(adata_embd,mthd, \
						'%s.%s'%(dt_nmPrfx_rgrs_embd_exprssd,'PAGODA_hc'))
						#cluster vs. all others
						all_grps=list(adata_embd.uns['rank_genes_groups']['names']. \
						dtype.names)
						if len(all_grps)>1:
							for ref in all_grps:
								sc.tl.rank_genes_groups(adata_embd,'PAGODA_hc', \
								groups=all_grps,reference=ref,method=mthd, \
								n_genes=4000000)
								pltNwrtGnGrpsVEchOthr(adata_embd,mthd, \
								dt_nmPrfx_rgrs_embd_exprssd,ref,'PAGODA_hc')
				##################
				#Added to simplify and analyze features of
				#interest
				##################
				for ftrIntrst in lFtrsAnls:
					if ftrIntrst in set(adata_embd.obs_keys()):
						if len(adata_embd.obs[ftrIntrst].unique())>1:
							if [len(adata_embd[adata_embd.obs[ftrIntrst]==vrbl]) \
							for vrbl in adata_embd.obs[ftrIntrst].unique()]:
								print('Feature %s has less than 2 samples for some cases'%ftrIntrst)
							else:
								adata_embd = sc.tl.rank_genes_groups(adata_embd, \
								ftrIntrst,method=mthd,copy=True,n_genes=4000000)
								pltNwrtGnGrps(adata_embd,mthd, \
								'%s.%s'%(dt_nmPrfx_rgrs_embd_exprssd,ftrIntrst))
								#cluster vs. all others
								all_grps=list(adata_embd.uns['rank_genes_groups']['names']. \
								dtype.names)
								if len(all_grps)>1:
									for ref in all_grps:
										sc.tl.rank_genes_groups(adata_embd,ftrIntrst, \
										groups=all_grps,reference=ref,method=mthd, \
										n_genes=4000000)
										pltNwrtGnGrpsVEchOthr(adata_embd,mthd, \
										dt_nmPrfx_rgrs_embd_exprssd,ref,ftrIntrst)
				##################
				if mkPaga:#The following will only write for wilcoxon
					mkPaga = False
					#Run one standard each vs. average to save 100 top genes
					sc.tl.rank_genes_groups(adata_embd, \
					'louvain',method=mthd,copy=False,n_genes=100)
					#PAGA embedding
					paga_embdg(adata_embd,'%s.louvain_n_leiden'% \
					dt_nmPrfx_rgrs_embd,geneRefs=geneRefs,cSeed=cSeed)
				#Write data
				adata_embd.write('%s.data.louvain_n_leiden.h5ad'% \
				dt_nmPrfx_rgrs_embd)
	return 0


def wrtGnTblsWlcnx(adata,dt_nmPrfx,thrshld=0.01,statsc='pvals_adj'):
	"""
	Method to write gene tables, and the significance of their expression in cell
	of one cluster compared to cells in all the other clusters, using Wilcoxon.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add 
	prefix to all output files. thrshld is the p-adjusted threshold to call significance.
	statsc is the key of the statistic to test the significance.
	Output: A folder with files included the statistics is going to be written.
	"""
	result = adata.uns['rank_genes_groups']
	groups = result['names'].dtype.names
	allDta = pd.DataFrame({group + '_' + key[:1]: result[key][group] for \
	group in groups for key in ['names', statsc,'logfoldchanges']}).to_numpy()
	nGns,nClmns = allDta.shape
	nClstrs = nClmns/3
	outFldr='%s.avrg.d'%'.'.join(dt_nmPrfx.split('.')[:-2])
	if not os.path.exists(outFldr):
		os.mkdir(outFldr)
	for cntClstr in range(int(nClstrs)):
		clstrPos=cntClstr*3
		crrntClstrGns = allDta[:,array([clstrPos+1,clstrPos,clstrPos+2])]
		# Significant
		crrntClstrGns_sgfcnt = crrntClstrGns[where(crrntClstrGns[:,0]< \
		thrshld)]
		crrntClstrGns_sgfcnt = '\n'.join(['%s\t%s\t%s'%(l[1],l[0],l[2]) for l in \
		sorted(crrntClstrGns_sgfcnt.tolist())])
		#write
		outFl = open(os.path.join(outFldr,'%sclstr_%s_avrgSgnfcnt.tsv'% (dt_nmPrfx, \
		cntClstr)),'w')
		outFl.write('%s\n'%'\t'.join(['#Name','FDR(BH)','logFold_change']))
		outFl.write(crrntClstrGns_sgfcnt)
		outFl.close()
		# Greater
		crrntClstrGns_grtr = crrntClstrGns[where((crrntClstrGns[:,0]< \
		thrshld) & (crrntClstrGns[:,2]>0))]
		crrntClstrGns_grtr = '\n'.join(['%s\t%s\t%s'%(l[1],l[0],l[2]) for l in \
		sorted(crrntClstrGns_grtr.tolist())])
		outFl = open(os.path.join(outFldr,'%sclstr_%s_avrgSgnfcntGrtr.tsv'% (dt_nmPrfx, \
		cntClstr)),'w')
		outFl.write('%s\n'%'\t'.join(['#Name','FDR(BH)','logFold_change']))
		outFl.write(crrntClstrGns_grtr)
		outFl.close()		
		# Lower
		crrntClstrGns_lwr = crrntClstrGns[where((crrntClstrGns[:,0]< \
		thrshld) & (crrntClstrGns[:,2]<0))]
		crrntClstrGns_lwr = '\n'.join(['%s\t%s\t%s'%(l[1],l[0],l[2]) for l in \
		sorted(crrntClstrGns_lwr.tolist())])
		outFl = open(os.path.join(outFldr,'%sclstr_%s_avrgSgnfcntLwr.tsv'% (dt_nmPrfx, \
		cntClstr)),'w')
		outFl.write('%s\n'%'\t'.join(['#Name','FDR(BH)','logFold_change']))
		outFl.write(crrntClstrGns_lwr)
		outFl.close()
	return 0

def wrtGnTblsWlcnxVEchOthr(adata,dt_nmPrfx,groupIn,thrshld=0.01,statsc='pvals_adj'):
	"""
	Method to write gene tables, and the significance of their expression in cell
	of one cluster compared to cells in each of all the other clusters, using Wilcoxon.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add 
	prefix to all output files. thrshld is the p-adjusted threshold to call significance.
	statsc is the key of the statistic to test the significance.
	Output: A folder with files included the statistics is going to be written.
	"""
	result = adata.uns['rank_genes_groups']
	groups = result['names'].dtype.names	
	allDtaFldNms = pd.DataFrame({group + '_' + key[:1]: result[key][group] for \
	group in groups for key in ['names']}).to_numpy()
	srtdGnNms = sorted(allDtaFldNms[:,0])
	#
	outFldr='%s.spfc.d'%'.'.join(dt_nmPrfx.split('.')[:-2]). \
	replace('.wilcoxon','').replace('.t-test_overestim_var',''). \
	replace('.t-test','')
	if not os.path.exists(outFldr):
		os.mkdir(outFldr)
	#
	toSrt_nms = []
	for lNms in allDtaFldNms.T:
		dGnNmsPos = dict([(g,p) for p,g in enumerate(lNms)])
		toSrt_nms.append([dGnNmsPos[nm] for nm in srtdGnNms])
	toSrt_nms = array(toSrt_nms).T	
	#
	allDtaStstcs = pd.DataFrame({group + '_' + key[:1]: result[key][group] for \
	group in groups for key in [statsc]}).to_numpy()
	allDtaFldChng = pd.DataFrame({group + '_' + key[:1]: result[key][group] for \
	group in groups for key in ['logfoldchanges']}).to_numpy()
	#sort
	for pos in range(len(groups)):
		allDtaStstcs[:,pos] = allDtaStstcs[:,pos][toSrt_nms[:,pos]]
		allDtaFldChng[:,pos] = allDtaFldChng[:,pos][toSrt_nms[:,pos]]
		allDtaFldNms[:,pos] = allDtaFldNms[:,pos][toSrt_nms[:,pos]]
	#Select upregulated
	#From the scanpy code: foldchanges = (expm1_func(mean_group) + 1e-9) / (expm1_func(mean_rest) + 1e-9) # add small value to remove 0's
	ar_IntrstngPos = where((np.max(allDtaStstcs,axis=1)<thrshld) & \
	(np.min(allDtaFldChng,axis=1)>0))#significant and the lower change is > 0 (e.g. greater in the rest than in the group) confirmed by looking into the markers
	#Get all	
	vsOtherFldChng=array(['|'.join(['%s:%s'%(groups[p],g) for p,g in \
	enumerate(fldChng)]) for fldChng in allDtaFldChng[ar_IntrstngPos]])	
	vsOtherFldFDRs=array(['|'.join(['%s:%s'%(groups[p],g) for p,g in \
	enumerate(vsFDRs)]) for vsFDRs in allDtaStstcs[ar_IntrstngPos]])
	allDta = array([np.mean(allDtaStstcs,axis=1)[ar_IntrstngPos], \
	allDtaFldNms[:,0][ar_IntrstngPos],vsOtherFldFDRs,vsOtherFldChng]).T
	crrntClstrGns = '\n'.join(['%s\t%s\t%s\t%s'%(l[1],l[0],l[2],l[3]) for \
	l in sorted(allDta.tolist())])
	#Write
	outFl = open(os.path.join(outFldr,'%sclstr_vsEchOthrGrtr.dwnrgltd.tsv'% \
	dt_nmPrfx),'w')
	outFl.write('%s\n'%'\t'.join(['#Name','avrg_FDR(BH)','FDR(BH)_vs_echOthr', \
	'logFoldChng_vs_echOthr']))
	outFl.write(crrntClstrGns)
	outFl.close()
	#Select upregulated
	ar_IntrstngPos = where((np.max(allDtaStstcs,axis=1)<thrshld) & \
	(np.max(allDtaFldChng,axis=1)<0))#significant and the higher change is < 0 (e.g. greater in the group than in the rest) confirmed by looking into the markers
	#Get all	
	vsOtherFldChng=array(['|'.join(['%s:%s'%(groups[p],g) for p,g in \
	enumerate(fldChng)]) for fldChng in allDtaFldChng[ar_IntrstngPos]])	
	vsOtherFldFDRs=array(['|'.join(['%s:%s'%(groups[p],g) for p,g in \
	enumerate(vsFDRs)]) for vsFDRs in allDtaStstcs[ar_IntrstngPos]])
	allDta = array([np.mean(allDtaStstcs,axis=1)[ar_IntrstngPos], \
	allDtaFldNms[:,0][ar_IntrstngPos],vsOtherFldFDRs,vsOtherFldChng]).T
	crrntClstrGns = '\n'.join(['%s\t%s\t%s\t%s'%(l[1],l[0],l[2],l[3]) for \
	l in sorted(allDta.tolist())])
	#Write
	outFl = open(os.path.join(outFldr,'%sclstr_vsEchOthrGrtr.uprgltd.tsv'% \
	dt_nmPrfx),'w')
	outFl.write('%s\n'%'\t'.join(['#Name','avrg_FDR(BH)','FDR(BH)_vs_echOthr', \
	'logFoldChng_vs_echOthr']))
	outFl.write(crrntClstrGns)
	outFl.close()
	#Test
	ar_IntrstngPos = where((np.max(allDtaStstcs,axis=1)<thrshld))
	vsOtherFldChng=array(['|'.join(['%s:%s'%(groups[p],g) for p,g in \
	enumerate(fldChng)]) for fldChng in allDtaFldChng[ar_IntrstngPos]])	
	vsOtherFldFDRs=array(['|'.join(['%s:%s'%(groups[p],g) for p,g in \
	enumerate(vsFDRs)]) for vsFDRs in allDtaStstcs[ar_IntrstngPos]])
	allDta = array([np.mean(allDtaStstcs,axis=1)[ar_IntrstngPos], \
	allDtaFldNms[:,0][ar_IntrstngPos],vsOtherFldFDRs,vsOtherFldChng]).T
	return 0

def wrtGnTblsLogreg(adata,dt_nmPrfx,statsc='scores'):
	"""
	Method to write gene tables, and the expression in cell of one cluster 
	compared to cells in each of all the other clusters, using log regrssion.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add 
	prefix to all output files. statsc is the key of the statistic to return results.
	Output: A folder with files included the statistics is going to be written.
	"""
	result = adata.uns['rank_genes_groups']
	groups = result['names'].dtype.names
	allDta = pd.DataFrame({group + '_' + key[:1]: result[key][group] for \
	group in groups for key in ['names', statsc]}).to_numpy()
	nGns,nClmns = allDta.shape
	nClstrs = nClmns/2
	#
	outFldr='%s.logReg.d'%'.'.join(dt_nmPrfx.split('.')[:-2])
	if not os.path.exists(outFldr):
		os.mkdir(outFldr)
	for cntClstr in range(int(nClstrs)):
		clstrPos=cntClstr*2
		crrntClstrGns = allDta[:,array([clstrPos+1,clstrPos])]
		# Greater
		crrntClstrGns_grtr = crrntClstrGns[where(crrntClstrGns[:,0]>0)]
		crrntClstrGns_grtr = sorted(crrntClstrGns_grtr.tolist())
		crrntClstrGns_grtr.reverse()
		crrntClstrGns_grtr = '\n'.join(['%s\t%s'%(l[1],l[0]) for l in \
		crrntClstrGns_grtr])
		outFl = open(os.path.join(outFldr,'%sclstr_%s_avrgSgnfcntGrtr.tsv'% \
		(dt_nmPrfx,cntClstr)),'w')
		outFl.write('%s\n'%'\t'.join(['#Name','score']))
		outFl.write(crrntClstrGns_grtr)
		outFl.close()		
		# Lower
		crrntClstrGns_lwr = crrntClstrGns[where(crrntClstrGns[:,0]<0)]
		crrntClstrGns_lwr = '\n'.join(['%s\t%s'%(l[1],l[0]) for l in \
		sorted(crrntClstrGns_lwr.tolist())])
		outFl = open(os.path.join(outFldr,'%sclstr_%s_avrgSgnfcntLwr.tsv'% \
		(dt_nmPrfx,cntClstr)),'w')
		outFl.write('%s\n'%'\t'.join(['#Name','score']))
		outFl.write(crrntClstrGns_lwr)
		outFl.close()
	return 0

def pltNwrtGnGrps(adata,mthd,nmPrfx):
	"""
	Method to plot and write gene tables, and the expression in cell of one 
	cluster compared to cells of all the other clusters.
	Input: adata is an AnnData object. mthd is the approach to compute the 
	signfiicance, including "Wilcoxon", "T-test with variance correction",  
	"Welch's t-test", and "logreg". nmPrfx is a name of the data to add prefix 
	to all output files.
	For more info about the tests look in:
	https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
	Output: A folder with files including the statistics and plots is going to be 
	written.
	"""
	sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False, \
	save='.top25DffExprssdGns.%s.%s.svg'% \
	(mthd,nmPrfx))
	if mthd in {'wilcoxon','t-test_overestim_var','t-test'}:
		wrtGnTblsWlcnx(adata,'%s.%sClstrs.'% \
		(nmPrfx,mthd))
	else:
		wrtGnTblsLogreg(adata,'%s.%sClstrs.'% \
		(nmPrfx,mthd))
	return None

def pltNwrtGnGrpsVEchOthr(adata,mthd,nmPrfx,group,embdg,plot=False):
	"""
	Method to plot and write gene tables, and the expression in cell of one 
	cluster compared to cells of each of all of the other clusters.
	Input: adata is an AnnData object. mthd is the approach to compute the 
	signfiicance, including "Wilcoxon", "T-test with variance correction",  
	"Welch's t-test", and "logreg". nmPrfx is a name of the data to add prefix 
	to all output files.
	For more info about the tests look in:
	https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
	Output: A folder with files including the statistics and plots is going to be 
	written.
	"""
	if plot:
		sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False, \
		save='.top25DffExprssdGns.%s.%s.%s.%s.svg'%(embdg,nmPrfx,mthd,group))
	if mthd in {'wilcoxon','t-test_overestim_var','t-test'}:
		wrtGnTblsWlcnxVEchOthr(adata,'%s.%s.%s.%sClstrs.'% \
		(nmPrfx,embdg,mthd,group),groupIn=group)
	else:
		pass
	return None
	

###########################
###########################
###########################
####    Cell-cycle     ####
###########################
###########################
###########################
def cc_rgrsn(adata,dt_nmPrfx,s_genes,g2m_genes,copy=False,cSeed=1234, \
	rgrss=True):
	"""
	Method to regress or score the cell-cycle associated genes.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add 
	prefix to all output files. s_genes is the list of genes in S phase, g2m_genes
	is the list of genes in G2/M phase. If copy will return a copy of adata will
	the newly added objects. cSeed is a seed to replicate results. If rgrss is 
	True it will regress the data following the cc scores, otherwise will just
	score the cell cycle.
	Output: adata and plots with cell-cyle regressed or scored.
	"""
	if copy: adata = adata.copy()
	cell_cycle_genes = sorted(set([x for x in s_genes+g2m_genes if x in \
	adata.var_names]))
	if cell_cycle_genes:
		#score cell cycle
		sc.tl.score_genes_cell_cycle(adata,s_genes=s_genes, \
		g2m_genes=g2m_genes)
		adata_cc_genes = adata[:, cell_cycle_genes]
		if rgrss:
			sc.tl.pca(adata_cc_genes,svd_solver='arpack',random_state= \
			cSeed)
			sc.pl.pca_scatter(adata_cc_genes, color='phase',save= \
			'%s.cell_cycle_phase.svg'%(dt_nmPrfx))
			#regress out
			sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
			sc.pp.scale(adata)
			#project again
			adata_cc_genes = adata[:, cell_cycle_genes]
			sc.tl.pca(adata_cc_genes,svd_solver='arpack',random_state= \
			cSeed)
			sc.pl.pca_scatter(adata_cc_genes, color='phase',save= \
			'.%s.cell_cycle_phase.rgrssd.svg'%(dt_nmPrfx))
		return adata if copy else None
	else:
		return None

def cc_rmv(adata,dt_nmPrfx,s_genes,g2m_genes,copy=False):
	"""
	Method to cell cycle correction by removing cell-cycle genes.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add 
	prefix to all output files. s_genes is the list of genes in S phase, g2m_genes
	is the list of genes in G2/M phase. If copy will return a copy of adata will
	the newly added objects. If copy the method will return a copy of the adata
	without the cell cycle genes.
	Output: adata excluding the cell cycle genes.
	"""
	if copy: adata = adata.copy()
	not_cell_cycle_genes = [x for x in adata.var_names if x not in \
	set(s_genes+g2m_genes)]
	#remove cell cycle
	adata_not_cc_genes = adata[:, not_cell_cycle_genes]
	return adata_not_cc_genes if copy else None

################################
################################
################################
####   Pipeline wrappers    ####
################################
################################
################################

############
#Scanpy's. Iterate to define suitable values for "max_genes" and "max_prcntg_mito"
############
def rn_scnpy(adata,dt_nmPrfx,max_genes=inf,max_prcntg_mito=inf, \
	geneRefs=['TH','PNMT','MYCN','PRRX1','CLDN11'],crrctFrBtch=True, \
	s_genes=None,g2m_genes=None,cSeed=1234,tstPAGODA_clstr=True, \
	anlyzOnlyDpRgrssdMtchndrl=False):
	"""
	Wrapper to run Satija's (Seurat) pipeline with the improvements of
	scanpy.
	Input: adata is an AnnData object. dt_nmPrfx is a name of the data to add prefix 
	to all output files. min_genes is the threshold minimum number of genes to filter 
	out cells. min_cells is the threshold minimum number of cells to filter out genes. 
	min_mean and max_mean are the minimum and maximum thresholds for mean 
	gene expression to filter out highly variable genes. min_disp is the threshold 
	minimum dispersion of normalized expression to filter genes variable genes. 
	clip_max_value is a threshold to filter genes exceeding a standard deviation of 10. 
	max_genes is the threshold  maximum number of genes with high expression 
	allowed, by scanpy receipt 2500, by Satija's "inf". max_prcntg_mito is the 
	maximum percentage of mitochondrial genes allowed, by scanpy receipt 0.05, by
	Satija's "inf" (Remove cells that have too many mitochondrial genes expressed or
	too many total counts). copy is a switch if so will return a copy of the adata 
	structure with Satija's processing. regrssByNmbrCnts is a switch to regress out 
	effects of total counts  per cell and regrssByPrcntgMitchdrl by the percentage of
	mitochondrial genes expressed (both False in Satija's and True in scanpy). 
	anlyzOnlyDpRgrssdMtchndrl is a switch, if If True will also analyze also non-cell 
	cycle corrected results.
	Output: adata with updated fields.
	
	"""
	#Preprocessing
	adata_rgrssdCDepNMtchndrl=rnStjScnpy(adata,dt_nmPrfx,max_genes= \
	max_genes,max_prcntg_mito=max_prcntg_mito,copy=True)
	if isinf(max_genes) or isinf(max_prcntg_mito):
		raise(Exception('Please update the "max_genes" and "max_prcntg_mito" values from the *Mitchndrl*.svg files...'))
	##########
	#added analyze only deep corrected results
	##########
	if anlyzOnlyDpRgrssdMtchndrl:
		dt_nmPrfx_rgrs1 = '%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs')
		wpr_dnstrm1(adata_rgrssdCDepNMtchndrl,dt_nmPrfx_rgrs1,geneRefs, \
		tstPAGODA_clstr=tstPAGODA_clstr)
		#depRgrs
		adata_rgrssdCDep=rnStjScnpy(adata,dt_nmPrfx,max_genes= \
		max_genes,max_prcntg_mito=max_prcntg_mito,copy=True, \
		regrssByPrcntgMitchdrl=False)
		dt_nmPrfx_rgrs1 = '%s.%s'%(dt_nmPrfx,'depRgrs')
		wpr_dnstrm1(adata_rgrssdCDep,dt_nmPrfx_rgrs1,geneRefs, \
		tstPAGODA_clstr=tstPAGODA_clstr)
	##########
	##########
	#Cell cycle score
	##########
	adata_rgrssdCDepNMtchndrl_ccScr = cc_rgrsn(adata_rgrssdCDepNMtchndrl, \
	'%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs'),s_genes,g2m_genes,copy=True, \
	cSeed=cSeed,rgrss=False)
	#Cell regression
	adata_rgrssdCDepNMtchndrl_ccRgrssn = cc_rgrsn(adata_rgrssdCDepNMtchndrl, \
	'%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs'),s_genes,g2m_genes,copy=True, \
	cSeed=cSeed)
	#Cell removal
	adata_rgrssdCDepNMtchndrl_ccRmv = cc_rmv(adata_rgrssdCDepNMtchndrl, \
	'%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs_ccRmv'),s_genes,g2m_genes,copy=True)
	#Downstream pipeline 1
	#Only plot cell cycle scores
	for adataP,dt_nmPrfx_rgrs1 in [(adata_rgrssdCDepNMtchndrl_ccScr, \
	'%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs_ccScr'))]:
		if adataP is not None:
			wpr_dnstrm1(adataP,dt_nmPrfx_rgrs1,geneRefs,fllAnalysis=False, \
			tstPAGODA_clstr=tstPAGODA_clstr)
	#Work on cell cycle regression
	for adataP,dt_nmPrfx_rgrs1 in ((adata_rgrssdCDepNMtchndrl_ccRgrssn, \
	'%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs_ccRgrssn')),(adata_rgrssdCDepNMtchndrl_ccRmv, \
	'%s.%s'%(dt_nmPrfx,'depNmtchndrlRgrs_ccRmv'))):
		if adataP is not None:
			wpr_dnstrm1(adataP,dt_nmPrfx_rgrs1,geneRefs,tstPAGODA_clstr= \
			tstPAGODA_clstr)
	##########


#######################################################
#######################################################
################                                  		              ###################
################      Methods to format data          ###################
################      in h5da data structure            ###################
################                                                      ###################
#######################################################
#######################################################


###################
###################
###################
####   Preprocessing   ####
####       methods       ####
###################
###################
###################
def mk_adata_frm_csv(infl):
	"""
	Method to make adata from csv infile.
	Input: infl is a csv file with the raw genec counts.
	Output: adata is an h5ad data structure with raw counts excluding
	ERCC prefixed gene names (i.e. spike-ins).
	"""
	adata = sc.read_csv(infl).T
	gnNms = adata.var_names
	noERCCpos = np.array([p for p in range(len(gnNms)) if \
	gnNms[p].find('ERCC')==-1])#exclude spike ins
	adata = adata[:,noERCCpos]
	return adata
	
def addPAGODAtSNE(adata,pagoda_tSNEFl,flippingTsne=False,addSfx=None, \
	flippingHrzntlTsne=False):
	"""
	Short method to add a tSNE from pagoda into the scanpy file 
	Input: adata is an AnnData object. pagoda_tSNEFl is a csv file with the 2-D 
	coordinates of the tSNE. flippingTsne is a switch, if True it will flip over the tSNE. 
	addSfx is a sufix to add to the embedding name in the AnnData object. 
	flippingHrzntlTsne is a switch, if True it will flip horizontally the tSNE. 
	Output: adata with the tSNE/embedding included as a "obsm" object.
	"""
	#Introduce the tSNE from PAGODA
	tsne_pagoda = mk_adata_frm_csv(pagoda_tSNEFl)
	tsne_pagoda_cllNms = array(tsne_pagoda.var_names)
	tsne_pagoda = array(tsne_pagoda.X).T
	if flippingTsne:
		tsne_pagoda[:,1]=-tsne_pagoda[:,1]
	if flippingHrzntlTsne:
		tsne_pagoda[:,0]=-tsne_pagoda[:,0]
	tsne_pagoda = list(tsne_pagoda)
	#############
	#Consider cells not present in the tsne
	#############
	tsne_pagoda.append([nan,nan])
	#############
	dCllNmPos = dict([(cllNm,p) for p,cllNm in enumerate(tsne_pagoda_cllNms)])
	#############
	#Consider different cells
	#############
	lstPos = len(tsne_pagoda_cllNms)
	#############
	arryCllPos = array([dCllNmPos.get(cllNm,lstPos) for cllNm in adata.obs_names])#this can be a source of errors and was designed to be sure all cell are included under the expectation that the number of cells is not greatly different
	if addSfx is None:
		adata.obsm['X_tsne'] = array(tsne_pagoda)[arryCllPos,:]
	else:
		adata.obsm['X_tsne_%s'%addSfx] = array(tsne_pagoda)[arryCllPos,:]
	return 0
