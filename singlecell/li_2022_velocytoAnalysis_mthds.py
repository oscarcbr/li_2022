#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  li_2022_velocytoAnalysis_mthds.py
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

"""
Group of methods to work with velocyto

Methods developed to run the velocity analysis included in the section "Cell clustering and expression analysis" 
in Li W et al. manuscript entitled "Chromaffin to neuroblast cell state transitions drive tumor plasticity in NF1 and 
KIF1Bβ deficient neuroblastoma, pheochromocytoma and composite tumors." included in Li W. 2022 thesis 
"Exploring Inter-and Intraheterogeneity in Childhood Neuroblastoma and Pheochromocytoma".
"""

#Assert is running on python 3
import os,sys,shutil
assert sys.version_info.major>2#Assert is running on python 3

import matplotlib
import pandas as pd
matplotlib.use('Agg')
matplotlib.rc('font', family='sans-serif')
matplotlib.rcParams.update({'font.size': 22})

import loompy
import scvelo as scv
import scanpy.external as sce
import palantir
import matplotlib.pyplot as plt

from velocyto.commands._run import _run
from pandas import Index
from numpy import array,isnan,where

scv.settings.set_figure_params('scvelo')


###################
###################
###################
####   Preprocessing    ####
####       methods        ####
###################
###################
###################
def wrprRun_smartseq2(output,rmskFl,prfxNm,anntGTFFl,lBAMFls,dtype="uint32", \
	dump='0',verbose=1,additional_ca={}):
	"""
	Wrapper to perform the read counting for UMI-less, not stranded, full-length 
	techniques such as SmartSeq2.
	Input: output is the output folder, if it does not exist it will be created. rmskFl is a 
	repeat mask .gtf file containing intervals to mask. prfxNm is the sample name that 
	will be used as the filename of the output. anntGTFFl is the annotation gtf files 
	for genes in  the genome. lBAMFls is the list of BAM files. dtype is the type of the 
	loom file layers - if more than 6000 molecules/reads per gene per cell are expected 
	set uint32 to avoid truncation (default in run_smartseq2: uint32)". dump is for 
	debugging purposes only: it will dump a molecular mapping report to hdf5. --dump 
	N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle 
	report is printed (default: 0). verbose sets the verbosity level: -v (only warnings) 
	-vv (warnings and info) -vvv (warnings, info and debug).
	Output: output is the output folder with the counts for spliced and  unspliced values.
	The main result file is a 4-layered loom file: prfxNm.loom. 
	"""
	bamfiles = tuple(lBAMFls)
	#Execute
	_run(bamfile=bamfiles,gtffile=anntGTFFl,bcfile=None, outputfolder=output, \
	sampleid=prfxNm,metadatatable=None,repmask=rmskFl,onefilepercell=True, \
	logic="SmartSeq2",without_umi=True,umi_extension="no",multimap=False, \
	test=False,samtools_threads=1,samtools_memory=1,dump=dump, \
	loom_numeric_dtype=dtype,verbose=verbose,additional_ca= additional_ca)
	return 0


def preprcss_mltplBAMFldrs(lBAMFldrs,lPrfxNm,lOutFldr,rmskFl,anntGTFFl, \
	outLoomMrgd,ovrWrtRslts=False):
	"""
	Wrapper to preprocess multiple BAM files for velocyto.
	Input: lBAMFldrs is a list of BAM folders for different samples.  lPrfxNm is a list of 
	prefix names (e.g. samples for the BAM folders) that must be in the same order as 
	the folders in lBAMFldrs. lOutFldr is a list of output folders in the same order. rmskFl 
	is a repeat mask .gtf file containing intervals to mask. anntGTFFl is the annotation 
	gtf files for genes in the genome. outLoomMrgd is the output loom with merged files.
	ovrWrtRslts is a switch if True will overwrite the previous files.
	Output: outLoomMrgd is the output loom with merged files.
	"""
	#Test that lBAMFldrs and lPrfxNm have the same samples
	assert len(lBAMFldrs)==len(lPrfxNm)
	#
	if not os.path.exists(outLoomMrgd) or ovrWrtRslts:
		lLoomFls=[]
		for p in range(len(lBAMFldrs)):
			prfxNm = lPrfxNm[p]
			output = lOutFldr[p]
			loomFl=os.path.join(output,'%s.loom'%prfxNm)
			if not os.path.exists(loomFl) or ovrWrtRslts:
				BAMFldr = lBAMFldrs[p]
				lBAMFls = (os.path.join(BAMFldr,f) for f in \
				os.listdir(BAMFldr) if f.find('.bam')>-1)
				wrprRun_smartseq2(output,rmskFl,prfxNm,anntGTFFl,lBAMFls)
			#append output files
			lLoomFls.append(loomFl)
		#Merge the samples
		shutil.copyfile(lLoomFls[0],outLoomMrgd)
		ds = loompy.connect(outLoomMrgd)
		#~ print(lLoomFls[0])
		for fn in lLoomFls[1:]:
			ds.add_loom(fn)
		ds.close()
	return 0

def wrpr_scvelo(adata,outLoomMrgd,dt_nmPrfx,ovrWrtRslts=False,min_shared_counts= \
	2000,n_top_genes=3000,n_pcs=30,n_neighbors=30,keysToPlot=['louvain','leiden', \
	'PAGODA_hc','stage','outcome','samples','S_score','G2M_score','phase'], \
	tstGnPrgrssn=True,geneRefs=['MYCN','PRRX1','CLDN11','EPAS1','MBOAT1','NTRK2','TH', \
	'MUC3A', 'PNMT', 'SOX6', 'BRCA1', 'SOX5', 'PIAS2','NTRK1'],lGnsIntrst=['MYCN','PRRX1', \
	'CLDN11','EPAS1','MBOAT1','NTRK2','TH','MUC3A', 'PNMT', 'SOX6', 'BRCA1', 'SOX5',  \
	'PIAS2','NTRK1'],frmCnts=True,n_genes_plot=100,flTrshCllCcl='melanomaCllCycl.tsv', \
	mostAcntrlGn='RTTN',retain_genes=None,clstr_basis='X_tsne_PAGODA_hc',clr_basis= \
	'PAGODA_hc'):
	"""
	Wrapper to run scVelo
	Input: adata is a anndata dataset annotated. outLoomMrgd is the output loom with merged 
	files. dt_nmPrfx is the prefix for all the files. ovrWrtRslts is a switch if True will overwrite the 
	previous files. min_shared_counts is the minimum number of cells expressed required to pass
	filtering (unspliced). n_top_genes is the number of genes to keep. n_pcs is the number of 
	principal components to use to compute the first- and second-order moments. n_neighbors 
	is the number of neighbors to use to compute moments. keysToPlot is a list of variables of 
	'obs' in adata to plot in the embedding. tstGnPrgrssn is a switch, if True will test gene 
	progressions with different models. geneRefs is the genes of references to plot in embeddings
	and plots. lGnsIntrst is the list of genes to plot for progression. frmCnts is a switch, if True will
	force to run log1p. n_genes_plot is the number of genes to plot for progression.  flTrshCllCcl is 
	the file of Tirosh's genes for gene expression. mostAcntrlGn is the gene that is expected to be
	expected by the most ancestral cell. retain_genes is a list of gene names to be retained 
	independent of thresholds. clstr_basis is the embedding name to compute moments and plot
	this and other results. clr_basis is the color basis to plot embedding, as included in the input 
	ha5d file.
	Ouput: results of the analysis for scvelo, including loom files and objects storing moments,
	and spliced/unspliced quantification.
	"""
	ouFldrDflt = os.path.join(os.path.split(outLoomMrgd)[0],'figures')
	if not os.path.exists(ouFldrDflt):
		os.mkdir(ouFldrDflt)
	#preprocessed adata object merge with the spliced/unspliced counts	
	ldata = scv.read(outLoomMrgd)#, cache=True)
	#correct the cell names accordingly	
	#~ ldata.obs.index=Index([c.replace('%s.'%'SK_N_SH','').replace('.bam',''). \
	#~ replace(':','.') for c in ldata.obs.index])	
	ldata.obs.index=Index(['.'.join(c.split('.')[1:-1]).replace(':','.') \
	for c in ldata.obs.index])
	#In case of nan in tnse, remove points
	if isnan(adata.obsm[clstr_basis].max()):
		adata = adata[adata.obs_names[array([p for p,v in \
		enumerate(adata.obsm[clstr_basis][:,0]) if not isnan(v)])]]
	adataVelo = scv.utils.merge(adata, ldata)
	#Preprocessing that is necessary consists of : - gene selection by 
	#detection (detected with a minimum number of counts) and high 
	#variability (dispersion). - normalizing every cell by its initial 
	#size and logarithmizing X.
	#Filtering and normalization is applied in the same vein to 
	#spliced/unspliced counts and X. Logarithmizing is only applied to 
	#X. If X is already preprocessed from former analysis, it won’t 
	#touch it.
	adataVelo_0 = adataVelo.copy()
	#All of this is summarized in a single function pp.filter_and_normalize:
	scv.pp.filter_and_normalize(adataVelo, min_shared_counts= \
	min_shared_counts,log=False,n_top_genes=n_top_genes,retain_genes= \
	retain_genes)
	if frmCnts:
		scv.pp.log1p(adataVelo)
	print(set(adataVelo.var_names).intersection(set(lGnsIntrst)))
	print(set(lGnsIntrst).difference(set(adataVelo.var_names)))
	print(len(lGnsIntrst))
	print(len(set(adataVelo.var_names).intersection(set(lGnsIntrst))))
	#After basic preprocessing (gene selection and normalization is 
	#sufficient), we compute the first- and second-order moments 
	#(basically means and variances) for velocity estimation:
	###############
	#Other parameters
	if flTrshCllCcl is not None:
		s_genes,g2m_genes = [],[]
		for l in open(flTrshCllCcl).read().splitlines():
			if l.strip() and l[0]!='#':
				gnNm,cllTyp=l.split('\t')
				if cllTyp=='G1/S':
					s_genes.append(gnNm)
				elif cllTyp=='G2/M':
					g2m_genes.append(gnNm)
	###############
	#List of factors to plot
	lPltInsrts=[]
	for key in keysToPlot:
		if key in adataVelo.obs and len(adataVelo.obs[key].unique(). \
		tolist())>1:
			lPltInsrts.append(key)
	lPltInsrts.extend([g for g in geneRefs if g in set(adataVelo.var_names)])	
	###############
	scv.pp.moments(adataVelo,n_pcs=n_pcs,n_neighbors=n_neighbors, \
	use_rep=clstr_basis)
	##############
	##############
	#Palantir
	sce.tl.palantir(adataVelo)
	fig = plt.figure()
	palantir.plot.plot_gene_expression(adataVelo.uns['palantir_imp_df'], \
	adataVelo.uns['palantir_tsne'],lGnsIntrst)
	plt.savefig(os.path.join(ouFldrDflt,'%s.palantir.gnIntrst.svg'%dt_nmPrfx))
	plt.close()
	#Pick up the start cell with the highest expression of gene mostAcntrlGn
	start_cell=adataVelo.obs_names[where(adataVelo.X[:,adataVelo.var_names== \
	mostAcntrlGn].T==adataVelo.X[:,adataVelo.var_names==mostAcntrlGn].T.max())[1]]
	pr_res = palantir.core.run_palantir(adataVelo.uns['palantir_ms_data'], \
	start_cell)
	fig = plt.figure()
	palantir.plot.plot_palantir_results(pr_res, adataVelo.uns['palantir_tsne'])
	plt.savefig(os.path.join(ouFldrDflt,'%s.palantir.embdngRslts.svg'%dt_nmPrfx))
	plt.close()
	assert (adataVelo.obs_names==pr_res.pseudotime.keys()).all() and \
	(adataVelo.obs_names==pr_res.entropy.keys()).all()
	adataVelo.obs['palantir_psdtime']=pr_res.pseudotime.values
	adataVelo.obs['palantir_entropy']=pr_res.entropy.values			
	##############
	##############
	#Copy data with moments calculated for dynamic model
	adataVelo_ori = adataVelo.copy()
	#The core of the software is the efficient and robust estimation of 
	#velocities, obtained with:
	scv.tl.velocity(adataVelo,basis=clstr_basis,groupby=clr_basis)
	#The probabilities of one cell transitioning into another cell are 
	#computed using cosine correlation and are stored in a matrix 
	#denoted as velocity graph:
	scv.tl.velocity_graph(adataVelo)
	#Finally the velocities can be projected and visualized in any 
	#embedding (e.g. UMAP) on single cell level, grid level, or as 
	#streamplot:
	scv.tl.velocity_confidence(adataVelo)
	print('\n#########Velocity confidence mean: %s, median: %s#########\n'% \
	(adataVelo.obs['velocity_confidence'].mean(), \
	adataVelo.obs['velocity_confidence'].median()))
	scv.pl.scatter(adataVelo,basis=clstr_basis,color= \
	'velocity_confidence', perc=[2, 98],save='.%s.scvelo.confidence.svg'% \
	dt_nmPrfx)
	#Plot embeddings
	for basis in [embd for embd in adataVelo.obsm.keys()]:
		basis=basis.split('X_')[1]
		scv.pl.velocity_embedding(adataVelo, basis=basis,color=lPltInsrts, \
		save='.%s.scvelo.%s.svg'%(dt_nmPrfx,basis),arrow_size=3,arrow_length=3)
		scv.pl.velocity_embedding_grid(adataVelo,basis=basis,color=lPltInsrts, \
		save='.%s.scvelo.grid.%s.svg'%(dt_nmPrfx,basis))
		scv.pl.velocity_embedding_stream(adataVelo,basis=basis,color=lPltInsrts, \
		save='.%s.scvelo.strm.%s.svg'%(dt_nmPrfx,basis))
		#Next two for palantir
		scv.pl.scatter(adataVelo,basis=basis,color='palantir_psdtime', \
		save='.%s.palantir.pseudotime.%s.svg'%(dt_nmPrfx,basis), \
		perc=[2, 98], colorbar=True,figsize=(35,35),size=2000,color_map='gnuplot', \
		rescale_color=[0,1])
		scv.pl.scatter(adataVelo,basis=basis,color='palantir_entropy', \
		save='.%s.palantir.entropy.%s.svg'%(dt_nmPrfx,basis), \
		perc=[2, 98], colorbar=True,figsize=(35,35),size=2000,color_map='gnuplot', \
		rescale_color=[0,1])
	#
	#For every tool module there is a plotting counterpart, which allows
	# you to examine your results in detail, e.g.:
	#scv.pl.velocity(adataVelo, var_names=['gene_A', 'gene_B'], **params)
	#scv.pl.velocity_graph(adataVelo, **params)	
	if tstGnPrgrssn:
		#Recover general stats
		adataVelo = adataVelo_ori.copy()#Reset the velocity
		expndStst(adataVelo,dt_nmPrfx,lGnsIntrst,lPltInsrts,flTrshCllCcl, \
		g2m_genes,s_genes,n_genes_plot,smplX='gnrl',clstr_basis=clstr_basis, \
		clr_basis=clr_basis)
		#With differentially kinetics test
		adataVelo = adataVelo_ori.copy()
		expndStst(adataVelo,'%s.fit_diff_kinetics'%dt_nmPrfx,lGnsIntrst,lPltInsrts,flTrshCllCcl, \
		g2m_genes,s_genes,n_genes_plot,smplX='gnrl',diff_kintcs=True, \
		clstr_basis=clstr_basis,clr_basis=clr_basis)
		#Calculate dynamical model
		adataVelo = adataVelo_ori.copy()
		expndStst(adataVelo,'%s.dynamical'%dt_nmPrfx,lGnsIntrst,lPltInsrts,flTrshCllCcl, \
		g2m_genes,s_genes,n_genes_plot,smplX='gnrl',diff_kintcs=True, \
		mode='dynamical',clstr_basis=clstr_basis,clr_basis=clr_basis)
		##############
		##############
		##############
		##############
	return adataVelo
	

def expndStst(adataVeloTmp,dt_nmPrfx,lGnsIntrst,lPltInsrts,flTrshCllCcl, \
	g2m_genes,s_genes,n_genes_plot,smplX='',diff_kintcs=False,mode='stochastic', \
	mostAcntrlGn=None,clstr_basis='X_tsne_PAGODA_hc',clr_basis='PAGODA_hc'):
	"""
	Compute extended statistics for scvelo.
	Input: adataVeloTmp is a anndata dataset annotated. dt_nmPrfx is the prefix for all the files. 
	lGnsIntrst is the list of genes to plot for progression. lPltInsrts is the list of genes/variables to 
	plot in embeddigns. flTrshCllCcl is the file of Tirosh's genes for gene expression. g2m_genes
	is the list of G2/M genes, s_genes is the list of S genes. n_genes_plot is the number of genes
	to plot for progression. smplX is a sufix for the case study to save the results. diff_kintcs is a 
	switch, if True will test different kinetics. mode is the mode to run the estimation using the 
	steady-state/deterministic, stochastic or dynamical model of transcriptional dynamics.
	mostAcntrlGn is the gene that is expected to be expected by the most ancestral cell.
	clstr_basis is the embedding name to compute moments and plot this and other results.
	clr_basis is the color basis to plot embedding, as included in the input.
	Ouput: results from velocity and Palantir.
	"""
	scv.tl.recover_dynamics(adataVeloTmp)
	#Test different kinetics
	if diff_kintcs:
		top_genes = adataVeloTmp.var['fit_likelihood'].sort_values(ascending=False).index[:100]
		scv.tl.differential_kinetic_test(adataVeloTmp, var_names=top_genes, groupby=clr_basis)
		kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
		scv.pl.scatter(adataVeloTmp, basis=top_genes[:20], ncols=5,add_outline='fit_diff_kinetics', \
		save='.%s.%s.scvelo.fit_likelihood.svg'%(dt_nmPrfx,smplX), **kwargs)
	#Copy data with moments calculated for model allowing different dynamics
	scv.tl.velocity(adataVeloTmp,basis=clstr_basis,diff_kinetics=diff_kintcs,groupby=clr_basis, \
	mode=mode)
	scv.tl.velocity_graph(adataVeloTmp)
	#############
	#Gather stats
	scv.tl.velocity_confidence(adataVeloTmp)
	print('\n#########Velocity confidence expndStst mean: %s, median: %s#########\n'% \
	(adataVeloTmp.obs['velocity_confidence'].mean(),adataVeloTmp.obs['velocity_confidence']. \
	median()))
	scv.pl.scatter(adataVeloTmp,basis=clstr_basis,color=('velocity_length','velocity_confidence'), \
	perc=[2, 98], \
	save='.%s.scvelo.confidence.%s.svg'%(dt_nmPrfx,smplX),cmap='coolwarm')
	#Plot embeddings
	for basis in [embd for embd in adataVeloTmp.obsm.keys() if embd.find('X_')>-1]:
		basis=basis.split('X_')[1]
		scv.pl.velocity_embedding(adataVeloTmp,basis=basis,color=lPltInsrts,save= \
		'.%s.scvelo.%s.%s.svg'%(dt_nmPrfx, \
		basis,smplX),arrow_size=3,arrow_length=3)
		scv.pl.velocity_embedding_grid(adataVeloTmp,basis=basis,color=lPltInsrts,save= \
		'.%s.scvelo.grid.%s.%s.svg'%(dt_nmPrfx,basis,smplX))
		scv.pl.velocity_embedding_stream(adataVeloTmp,basis=basis,color=lPltInsrts,save= \
		'.%s.scvelo.strm.%s.%s.svg'%(dt_nmPrfx,basis,smplX))
	#Plot genes explaining transitions in each clsuter the best
	scv.tl.rank_velocity_genes(adataVeloTmp, groupby=clr_basis)
	df = pd.DataFrame(adataVeloTmp.uns['rank_velocity_genes']['names'])
	df.head()
	#Plot these genes
	for clstrX in adataVeloTmp.obs[clr_basis].unique():
		scv.pl.scatter(adataVeloTmp,basis=df[clstrX][:5].tolist(),ylabel=clstrX,frameon=False,save= \
		'.%s.scvelo.%s.drctnlty.%s.%s.%s.svg'%(dt_nmPrfx,clr_basis,'top5rnkdGns',smplX,clstrX), \
		color=clr_basis)
	#Plot cell cycle
	if flTrshCllCcl is not None:
		cell_cycle_genes = sorted(set([x for x in s_genes+g2m_genes if x in adataVeloTmp. \
		var_names]))
		if cell_cycle_genes:
			#score cell cycle
			scv.tl.score_genes_cell_cycle(adataVeloTmp,s_genes=s_genes,g2m_genes=g2m_genes)
			s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adataVeloTmp)
			s_genes = scv.get_df(adataVeloTmp[:, s_genes], 'spearmans_score', sort_values=True).index
			g2m_genes = scv.get_df(adataVeloTmp[:, g2m_genes], 'spearmans_score', sort_values= \
			True).index
			kwargs = dict(frameon=False, ylabel='cell cycle genes')
			scv.pl.scatter(adataVeloTmp,  list(s_genes[:2]) + list(g2m_genes[:3]),	**kwargs,save= \
			'.%s.scvelo.%s.%s.%s.svg'%(dt_nmPrfx,clr_basis,'cllCycl_genes',smplX),color=clr_basis)
	############
	############
	#Calculate time dynamics
	adataVeloLtnt = adataVeloTmp.copy()
	adataVeloPsd = adataVeloTmp.copy()
	#select start cell
	if mostAcntrlGn:
		start_cellLtnt = adataVeloLtnt.obs_names[where(adataVeloLtnt.X[:,adataVeloLtnt.var_names== \
		mostAcntrlGn].T==adataVeloLtnt.X[:,adataVeloLtnt.var_names==mostAcntrlGn].T.max())[1]]
		start_cellPsd = adataVeloPsd.obs_names[where(adataVeloPsd.X[:,adataVeloPsd.var_names== \
		mostAcntrlGn].T==adataVeloPsd.X[:,adataVeloPsd.var_names==mostAcntrlGn].T.max())[1]]
		scv.tl.recover_latent_time(adataVeloLtnt,root_key=start_cellLtnt)
		scv.tl.velocity_pseudotime(adataVeloPsd,root_key=start_cellPsd)
	#
	else:
		scv.tl.recover_latent_time(adataVeloLtnt)
		scv.tl.velocity_pseudotime(adataVeloPsd)
	for adataVeloN,tkey in ((adataVeloPsd,'velocity_pseudotime'),(adataVeloLtnt,'latent_time')):
		#Plot gene dynamics progression over all
		scv.pl.scatter(adataVeloN,basis=clstr_basis, color=tkey,fontsize=24, size=100,color_map= \
		'gnuplot', perc=[2, 98], colorbar=True,rescale_color=[0,1],save='.%s.scvelo.%s.%s.svg'% \
		(dt_nmPrfx,tkey,smplX))	
		adataVeloN.uns['neighbors']['distances'] = adataVeloN.obsp['distances']
		adataVeloN.uns['neighbors']['connectivities'] = adataVeloN.obsp['connectivities']			
		#Calculate PAGA
		scv.tl.paga(adataVeloN, groups=clr_basis)
		df = scv.get_df(adataVeloN, 'paga/transitions_confidence', precision=2).T
		df.style.background_gradient(cmap='Blues').format('{:.2g}')
		scv.pl.paga(adataVeloN,basis='tsne_PAGODA_hc',size=100,perc=[2, 98],alpha=.1, \
		min_edge_width=2,node_size_scale=1.5,save='.%s.scvelo.paga.%s.%s.svg'%(dt_nmPrfx, \
		tkey,smplX),color=clr_basis)
		#Plot progression of genes
		scv.pl.heatmap(adataVeloN, var_names=sorted(set(lGnsIntrst)),tkey='velocity_pseudotime', \
		n_convolve=10,col_color=clr_basis,save='.%s.scvelo.gnPrgssn.%s.lIntrst.%s.%s.%s.svg'% \
		(dt_nmPrfx,clr_basis,tkey,'lGnsIntrst',smplX),figsize=(30,100))
		##############		
		top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:300]
		if dt_nmPrfx.split('.')[0]!='11NBs_noMitchnd_KEql12':
			for cluster in lPltInsrts:
				if cluster in ['louvain','leiden',clr_basis]:
					scv.pl.heatmap(adataVeloN, var_names=top_genes,tkey=tkey,	n_convolve=100, \
					col_color=cluster,save='.%s.scvelo.gnPrgssn.%s.%s.%s.svg'%(dt_nmPrfx,tkey,cluster, \
					smplX))
		##############
		#Plot Susanne's plus top genes
		top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:n_genes_plot]. \
		tolist()
		top_genes.extend(lGnsIntrst)
		scv.pl.heatmap(adataVeloN, var_names=sorted(set(top_genes)),tkey='velocity_pseudotime', \
		n_convolve=10,col_color=clr_basis, \
		save='.%s.scvelo.gnPrgssn.%s.topGnsNlIntrst.%s.%s.%s.svg'%(dt_nmPrfx,	clr_basis,tkey, \
		n_genes_plot,smplX),figsize=(30,30))
		#Expand to 300
		top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:300].tolist()
		top_genes.extend(lGnsIntrst)
		scv.pl.heatmap(adataVeloN, var_names=sorted(set(top_genes)),tkey='velocity_pseudotime', \
		n_convolve=10,col_color=clr_basis,save='.%s.scvelo.gnPrgssn.%s.top300GnsNlIntrst.%s.%s.%s.svg'% \
		(dt_nmPrfx,clr_basis,tkey,300,smplX),figsize=(30,100))			
		#Plot dynamics
		slctd_gns = sorted(set(lGnsIntrst).intersection(set(adataVeloN.var_names)))[:20]
		scv.pl.scatter(adataVeloN,basis=slctd_gns,ncols=5,frameon=False,save= \
		'.%s.scvelo.gnPrgssn.%s.lIntrst20Dynamics.%s.%s.%s.svg'%(dt_nmPrfx,clr_basis,tkey,'dynamics', \
		smplX),color=clr_basis)	
		##############
		##############
		##############
		#Special cases for velocities reconstruction in Li et al. 2022
		###############
		###############
		###############
		# Wenyu paper
		if dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn1')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['3','4'])]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
		#
		elif dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn2')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['5'])]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
		#
		elif dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn3')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['7','10'])]
			#Select only precursos in proximity to chromaffin	
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
		#
		elif dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn4')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['7','1'])]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
		#
		elif dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn5')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['5','8'])]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
		#
		elif dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn6')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['5','19'])]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
		#
		elif dt_nmPrfx.find('allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk')>-1 and dt_nmPrfx. \
			find('trnstn7')>-1:
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[adataVeloN.obs[clr_basis].isin(['5','17'])]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis,clr_basis= \
			clr_basis)
			##############
		elif dt_nmPrfx.split('.')[0]=='mm14_allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk':
			assert lGnsIntrst is not None
			#next special case for chromaffin an precursors
			adataVeloN = adataVeloN[(adataVeloN.obs[clr_basis]=='5') | (adataVeloN.obs[clr_basis]=='8') | \
			(adataVeloN.obs[clr_basis]=='15') | (adataVeloN.obs[clr_basis]=='17') | (adataVeloN.obs[clr_basis]=='18') | \
			(adataVeloN.obs[clr_basis]=='19')]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloN.obs_names])
			adataVeloN = adataVeloN[adataVeloTmp2Pos]
			analyseMouseAGIncldPlntr(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis=clstr_basis, \
			clr_basis=clr_basis)
			####
			#next special case for excluding c18
			adataVeloE = adataVeloN[(adataVeloN.obs[clr_basis]=='5') | (adataVeloN.obs[clr_basis]=='8') | \
			(adataVeloN.obs[clr_basis]=='15') | \
			(adataVeloN.obs[clr_basis]=='17') | (adataVeloN.obs[clr_basis]=='19')]
			#Select only precursos in proximity to chromaffin
			adataVeloTmp2Pos=array([nm for nm in adataVeloE.obs_names])
			adataVeloE = adataVeloE[adataVeloTmp2Pos]
			analyseMouseAGIncldPlntr(adataVeloE,lGnsIntrst,tkey,'%s.excldC18'%dt_nmPrfx,n_genes_plot,clstr_basis= \
			clstr_basis,clr_basis=clr_basis)
			##############
			####
			#next special case for excluding WTs
			adataVeloE = adataVeloN[(adataVeloN.obs[clr_basis]=='5') | (adataVeloN.obs[clr_basis]=='8') | \
			(adataVeloN.obs[clr_basis]=='15') | (adataVeloN.obs[clr_basis]=='17') | (adataVeloN.obs[clr_basis]=='18') | \
			(adataVeloN.obs[clr_basis]=='19')]
			adataVeloE = adataVeloE[adataVeloE.obs['outcome'].isin(['DKO','KIF_KO','KIF_KO_AG_3m_wyL','NF1_KO', \
			'NF1_KO_AG_3m_wyL'])]
			#Select only precursos in proximity to chromaffin
			analyseMouseAGIncldPlntr(adataVeloE,lGnsIntrst,tkey,'%s.excldWT'%dt_nmPrfx,n_genes_plot,clstr_basis= \
			clstr_basis,clr_basis=clr_basis)
			##############
			##############
			####
			#special case for each genotype
			for addPrfx,stSmpls in (('onlyWT',['SS2_17_343','SS2_17_344','SS2_18_111','SS2_18_114','SS2_18_116', \
			'SS2_20_101','SS2_20_103','SS2_20_107']),('onlyKIF_KO',['SS2_17_345','SS2_17_346','SS2_18_112', \
			'SS2_18_115','SS2_20_113','SS2_20_149','SS2_20_155']),('onlyNF1_KO',['SS2_17_368','SS2_17_369', \
			'SS2_17_369_2','SS2_20_111','SS2_20_115','SS2_20_119']),('onlyDKO',['SS2_17_347','SS2_17_348', \
			'SS2_17_370','SS2_17_371','SS2_17_371_2','SS2_20_077','SS2_20_109','SS2_20_117'])):
				adataVeloE = adataVeloN[(adataVeloN.obs[clr_basis]=='5') | \
				(adataVeloN.obs[clr_basis]=='8') | (adataVeloN.obs[clr_basis]=='15') | (adataVeloN.obs[clr_basis]=='17') | \
				(adataVeloN.obs[clr_basis]=='18') | (adataVeloN.obs[clr_basis]=='19')]
				adataVeloE = adataVeloE[adataVeloE.obs['samples'].isin(stSmpls)]
				analyseMouseAGIncldPlntr(adataVeloE,lGnsIntrst,tkey,'%s.%s'%(dt_nmPrfx,addPrfx),n_genes_plot, \
				clstr_basis=clstr_basis,clr_basis=clr_basis)
	############
	############
	return 0


def analyseMouseAG(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot,clstr_basis='X_tsne_PAGODA_hc', \
	clr_basis='PAGODA_hc'):
	"""
	Analyze adrenal glands data.
	Input: adataVeloN is a anndata dataset annotated. lGnsIntrst is the list of genes to plot for progression. 
	dt_nmPrfx  is the prefix for all the files. n_genes_plot is the number of genes to plot for progression. 
	clstr_basis is the embedding name to compute moments and plot this and other results. clr_basis is the 
	color basis to plot embedding, as included in the input.
	Ouput: results for gene progression
	"""
	##############
	##############
	scv.pl.scatter(adataVeloN,basis=clstr_basis,color=tkey,fontsize=24, size=100,color_map='gnuplot', \
	perc=[2, 98], colorbar=True,rescale_color=[0,1],save='.%s.scvelo.%s.2clstrs.%s.svg'%(dt_nmPrfx,tkey,'full_ori'))
	#List interest Susanne
	scv.pl.heatmap(adataVeloN, var_names=sorted(set(lGnsIntrst)),tkey=tkey,n_convolve=10,col_color=clr_basis, \
	save='.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey,'lGnsIntrst'),figsize=(30,30))
	#Plot Susanne's plus top genes
	top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:n_genes_plot].tolist()
	top_genes.extend(lGnsIntrst)
	scv.pl.heatmap(adataVeloN, var_names=sorted(set(top_genes)),tkey=tkey,n_convolve=10,col_color=clr_basis, \
	save='.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey,n_genes_plot),figsize=(30,30))
	#Expand to 300
	top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:300].tolist()
	top_genes.extend(lGnsIntrst)
	scv.pl.heatmap(adataVeloN, var_names=sorted(set(top_genes)),tkey=tkey,n_convolve=10,col_color=clr_basis, \
	save='.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey,300),figsize=(30,100))			
	#Plot dynamics
	slctd_gns = sorted(set(lGnsIntrst).intersection(set(adataVeloN.var_names)))[:20]
	scv.pl.scatter(adataVeloN,basis=slctd_gns,ncols=5,frameon=False,save= \
	'.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey,'dynamics'),color=clr_basis)
	##############
	##############
	

def analyseMouseAGIncldPlntr(adataVeloN,lGnsIntrst,tkey,dt_nmPrfx,n_genes_plot, \
	clstr_basis='X_tsne_PAGODA_hc',clr_basis='PAGODA_hc'):
	"""
	Analyze adrenal glands data, including palantir analysis.
	Input: adataVeloN is a anndata dataset annotated. lGnsIntrst is the list of genes to plot for progression. 
	dt_nmPrfx  is the prefix for all the files. n_genes_plot is the number of genes to plot for progression. 
	clstr_basis is the embedding name to compute moments and plot this and other results. clr_basis is the 
	color basis to plot embedding, as included in the input.
	Ouput: results for gene progression and palatir analysis
	"""
	##############
	##############
	scv.pl.scatter(adataVeloN,basis=clstr_basis,color=tkey,fontsize=24,color_map='gnuplot', \
	perc=[2, 98], colorbar=True,rescale_color=[0,1],save='.%s.scvelo.%s.2clstrs.%s.svg'% \
	(dt_nmPrfx,tkey,'full_ori'), \
	figsize=(75,30),size=12000)
	#List interest Susanne
	scv.pl.heatmap(adataVeloN, var_names=sorted(set(lGnsIntrst)),tkey=tkey,n_convolve=10, \
	col_color=clr_basis,save='.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey, \
	'lGnsIntrst'),figsize=(30,30))
	#Plot Susanne's plus top genes
	top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:n_genes_plot].tolist()
	top_genes.extend(lGnsIntrst)
	scv.pl.heatmap(adataVeloN, var_names=sorted(set(top_genes)),tkey=tkey,n_convolve=10, \
	col_color=clr_basis,save='.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey, \
	n_genes_plot),figsize=(30,30))
	#Expand to 300
	top_genes = adataVeloN.var_names[adataVeloN.var.fit_likelihood.argsort()[::-1]][:300].tolist()
	top_genes.extend(lGnsIntrst)
	scv.pl.heatmap(adataVeloN, var_names=sorted(set(top_genes)),tkey=tkey,n_convolve=10, \
	col_color=clr_basis,save='.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey, \
	300),figsize=(30,100))			
	#Plot dynamics
	slctd_gns = sorted(set(lGnsIntrst).intersection(set(adataVeloN.var_names)))[:20]
	scv.pl.scatter(adataVeloN,basis=slctd_gns,ncols=5,frameon=False,save= \
	'.%s.scvelo.gnPrgssn.%s.%s.2clstrs.%s.svg'%(dt_nmPrfx,clr_basis,tkey,'dynamics'),color= \
	clr_basis)
	##############
	##############
	#Palantir
	scv.pl.scatter(adataVeloN,basis=clstr_basis,color='palantir_psdtime',save= \
	'.%s.palantir.pseudotime.%s.svg'%(dt_nmPrfx,'2clstrs'),perc=[2, 98], colorbar=True, \
	color_map='gnuplot',rescale_color=[0,1],figsize=(75,30),size=12000)
	scv.pl.scatter(adataVeloN,basis=clstr_basis,color='palantir_entropy', \
	save='.%s.palantir.entropy.%s.svg'%(dt_nmPrfx,'2clstrs'),perc=[2, 98], colorbar=True,figsize=(75, \
	30),size=12000,rescale_color=[0,1])
