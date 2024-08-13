## AuxFunc.py

import anndata
import argparse
import cellphonedb
import datetime
import gzip
import inspect
import itertools
import math
import matplotlib
import numpy
import openpyxl
import os
import pandas
import re
import scanpy
import scipy
import seaborn
import string
import struct
import sys
import time

import Common


def annotation(adata, dirOut, keyAdd, keyClu, clus2cell):
	try:
		## Annotate the cells.
		adata.obs[keyAdd] = adata.obs[keyClu].map(clus2cell).astype('category')
		## UMAP
		titlePlot = ('%s: %s cell types' % (keyClu, adata.obs[keyAdd].nunique()))
		um = scanpy.pl.umap(adata, color=keyAdd, use_raw=True, legend_loc='on data', legend_fontsize=6, \
							legend_fontoutline=1, title=titlePlot, return_fig=True)
		dirFigClu = os.path.join(dirOut, 'figures', 'clusters')
		Common.dirExist(dirFigClu)
		fileUMAP = os.path.join(dirFigClu, 'umap_%s.pdf' % (keyAdd))
		um.savefig(fileUMAP, bbox_inches='tight')
		Common.fileExist(fileUMAP)
		## UMAP for each cell type
		celltypes = adata.obs[keyAdd].unique()
		numCelltypes = len(celltypes)
		numCols = min(10, numCelltypes)
		numRows = math.ceil(numCelltypes/numCols)
		fig, axes = matplotlib.pyplot.subplots(numRows, numCols, figsize=(3.65 * numCols, 3.95 * numRows))
		axes = axes.flatten()
		for ax, celltype in zip(axes, celltypes):
			scanpy.pl.umap(adata, color=keyAdd, use_raw=True, groups=[celltype], legend_loc=None, na_color='white', na_in_legend=False, \
						   title=celltype, show=False, ax=ax)
		for i in range(numCelltypes, numRows * numCols):
			fig.delaxes(axes[i])
		matplotlib.pyplot.tight_layout()
		fileUMAP = os.path.join(dirFigClu, 'umap_split%s.png' % (keyAdd))
		matplotlib.pyplot.savefig(fileUMAP, dpi=400)
		return adata
	except Exception as e:
		Common.exception(e, 'could not do the cell annotation')


def clustering(adata, resol, dirOut, splitClusters=False, markers=Common.markersKnownMouse):
	try:
		print('\nUsing resolution %s ...' % (resol))
		## Compute clusters using the leiden method and store the results with the name 'leiden_res.{resol:3.1f}'.
		keyAdd = f'leiden_res.{resol:3.1f}'
		scanpy.tl.leiden(adata, resolution=resol, key_added=keyAdd, directed=False, n_iterations=2, flavor='igraph')
		## UMAP
		titlePlot = ('%s: %s clusters' % (keyAdd, adata.obs[keyAdd].nunique()))
		um = scanpy.pl.umap(adata, color=keyAdd, use_raw=True, legend_loc='on data', legend_fontsize=6, \
							legend_fontoutline=1, title=titlePlot, return_fig=True)
		dirFigClu = os.path.join(dirOut, 'figures', 'clusters')
		Common.dirExist(dirFigClu)
		fileUMAP = os.path.join(dirFigClu, 'umap_%s.png' % (keyAdd))
		um.savefig(fileUMAP, bbox_inches='tight')
		Common.fileExist(fileUMAP)
		## dotplot
		#for k in markers.keys():
		#	markers[k] = [g for g in markers[k] if g in adata.var_names]
		#dp = scanpy.pl.dotplot(adata, var_names=markers, groupby=keyAdd, use_raw=False, dendrogram=False, \
		#					   title=titlePlot, return_fig=True)
		#fileDotPlot = os.path.join(dirFigClu, 'dotplot_%s.pdf' % (keyAdd))
		#dp.style(cmap='viridis_r', dot_edge_color=None, grid=True).savefig(fileDotPlot, bbox_inches='tight')
		#Common.fileExist(fileDotPlot)

		## stacked-violin plot
		#sv = scanpy.pl.stacked_violin(adata, var_names=markers, groupby=keyAdd, use_raw=True, dendrogram=False, \
		#							  title=titlePlot, return_fig=True)
		#fileStkViol = os.path.join(dirFigClu, 'stacked_violin_%s.pdf' % (keyAdd))
		#sv.style(cmap='viridis_r', linewidth=0).savefig(fileStkViol, bbox_inches='tight')
		#Common.fileExist(fileStkViol)
		if (splitClusters):
			## UMAP for each cluster
			clusters = adata.obs[keyAdd].unique()
			numClusters = len(clusters)
			numCols = min(10, numClusters)
			numRows = math.ceil(numClusters/numCols)
			fig, axes = matplotlib.pyplot.subplots(numRows, numCols, figsize=(3.65 * numCols, 3.95 * numRows))
			axes = axes.flatten()
			for ax, cluster in zip(axes, sorted(clusters, key=int)):
				scanpy.pl.umap(adata, color=keyAdd, use_raw=True, groups=[cluster], legend_loc=None, na_color='white', na_in_legend=False, \
							   title=f'Cluster {cluster}', show=False, ax=ax)
			for i in range(numClusters, numRows * numCols):
				fig.delaxes(axes[i])
			matplotlib.pyplot.tight_layout()
			fileUMAP = os.path.join(dirFigClu, 'umap_splitclusters_%s.png' % (keyAdd))
			matplotlib.pyplot.savefig(fileUMAP, dpi=400)
		return adata
	except Exception as e:
		Common.exception(e, 'could not do the clustering')


def degs(adata, dirOut, keyClu, keyGroup, groupsUse, refUse, nameExcel, methodTest='wilcoxon'):
	try:
		keys = ['names', 'pvals', 'logfoldchanges', 'pvals_adj', 'scores']
		group_dfs = {}
		cluAll = adata.obs[keyClu].values.categories
		for cluCurr in cluAll:
			print('\nProcessing \'%s\' ...' % (cluCurr))
			adataClu = adata[adata.obs[keyClu] == cluCurr].copy()
			scanpy.tl.rank_genes_groups(adataClu, groupby=keyGroup, use_raw=True, groups=groupsUse, reference=refUse, \
										method=methodTest, corr_method='bonferroni', pts=True)
			result = adataClu.uns['rank_genes_groups']
			groups = result['names'].dtype.names
			for group in groups:
				groupAdd = '%s in %s' % (group, cluCurr)
				print('Processing \'%s\' ...' % (groupAdd))
				group_data = {key: result[key][group] for key in keys}
				group_df = pandas.DataFrame(group_data)
				pts = result['pts'][group]
				pts_rest = None
				if ('pts_rest' in result):
					pts_rest = result['pts_rest'][group]
				else:
					pts_rest = result['pts'][refUse]
				names = group_df['names']
				pts = pts[names]
				pts_rest = pts_rest[names]
				group_df['pts'] = list(pts)
				group_df['pts_rest'] = list(pts_rest)
				group_dfs[groupAdd] = group_df
		print(group_dfs)
		dirResDEG = os.path.join(dirOut, 'results', 'DEG')
		Common.dirExist(dirResDEG)
		fileExcel = os.path.join(dirResDEG, nameExcel)
		if os.path.isfile(fileExcel):
			os.remove(fileExcel)
		with pandas.ExcelWriter(fileExcel, mode='w') as writer:
			for group, df in group_dfs.items():
				df.to_excel(writer, sheet_name=group, header=True, index=True)
		return 0
	except Exception as e:
		Common.exception(e, 'could not calculate DEGs')


def markers(adata, dirOut, keyGroup, keyAdd, nameExcel, groupsUse='all', refUse='rest', methodTest='wilcoxon', preSheet=''):
	try:
		## Compute a ranking for the highly differential genes in each cluster.
		## Use 'bonferroni' as the p-value correction method, which is what Seurat::FindMarkers is using.
		scanpy.tl.rank_genes_groups(adata, groupby=keyGroup, use_raw=True, groups=groupsUse, reference=refUse, \
									method=methodTest, corr_method='bonferroni', pts=True, key_added=keyAdd)
		## Save the DEGs of each cluster.
		result = adata.uns[keyAdd]
		groups = result['names'].dtype.names
		print(groups)
		keys = ['names', 'pvals', 'logfoldchanges', 'pvals_adj', 'scores']
		group_dfs = {}
		for group in groups:
			## Extract data for the current group
			group_data = {key: result[key][group] for key in keys}
			group_df = pandas.DataFrame(group_data)
			pts = result['pts'][group]
			pts_rest = result['pts_rest'][group]
			names = group_df['names']
			pts = pts[names]
			pts_rest = pts_rest[names]
			group_df['pts'] = list(pts)
			group_df['pts_rest'] = list(pts_rest)
			group_dfs[group] = group_df
		print(group_dfs)
		dirResDEG = os.path.join(dirOut, 'results', 'DEG')
		Common.dirExist(dirResDEG)
		fileExcel = os.path.join(dirResDEG, nameExcel)
		if os.path.isfile(fileExcel):
			os.remove(fileExcel)
		with pandas.ExcelWriter(fileExcel, mode='w') as writer:
			## Write each group's DataFrame to a different sheet
			for group, df in group_dfs.items():
				df.to_excel(writer, sheet_name=preSheet+group, header=True, index=True)
		return adata
	except Exception as e:
		Common.exception(e, 'could not calculate marker genes')


def harmony(adata, keyBatch, keyPCA='X_pca', keyHmy='X_harmony'):
	try:
		## Integrate data using Harmony
		scanpy.external.pp.harmony_integrate(adata, key=keyBatch, basis=keyPCA, adjusted_basis=keyHmy, max_iter_harmony=100)
		return adata
	except Exception as e:
		Common.exception(e, 'could not run Harmony')


def neighbor(adata, dirOut, repUse, npcs=51, keyBatch='orig.ident'):
	try:
		print('\nUsing \'%s\' representation ...' % (repUse))
		print('Using %s PCs ...' % (npcs))
		## Compute the neighborhood graph
		scanpy.pp.neighbors(adata, n_neighbors=30, n_pcs=npcs, use_rep=repUse)
		## Embed the graph in 2 dimensions using UMAP
		scanpy.tl.umap(adata)
		## UMAP for each batch
		batches = adata.obs[keyBatch].unique()
		numBatches = len(batches)
		numCols = min(4, numBatches)
		numRows = math.ceil(numBatches/numCols)
		fig, axes = matplotlib.pyplot.subplots(numRows, numCols, figsize=(3.65 * numCols, 3.95 * numRows))
		axes = axes.flatten()
		#sizePoint = 120000/adata.n_obs
		for ax, batch in zip(axes, batches):
			#scanpy.pl.umap(adata[adata.obs[keyBatch] == batch], color=keyBatch, use_raw=True, legend_loc=None, size=sizePoint, \
			#			   na_in_legend=False, title=batch, show=False, ax=ax)
			scanpy.pl.umap(adata, color=keyBatch, use_raw=True, groups=[batch], legend_loc=None, na_color='white', na_in_legend=False, \
						   title=batch, show=False, ax=ax)
		for i in range(numBatches, numRows * numCols):
			fig.delaxes(axes[i])
		matplotlib.pyplot.tight_layout()
		dirFigClu = os.path.join(dirOut, 'figures', 'clusters')
		Common.dirExist(dirFigClu)
		fileUMAP = os.path.join(dirFigClu, 'umap_splitbatches.png')
		matplotlib.pyplot.savefig(fileUMAP, dpi=400)
		return adata
	except Exception as e:
		Common.exception(e, 'could not compute the neighborhood graph')


def pca(adata, npcs=51):
	try:
		print('\nUsing %s PCs ...' % (npcs))
		## Scale each gene to unit variance. 
		## Clip values exceeding standard deviation 10. Setting this can help reduce 
		## the effects of features that are only expressed in a very small number of cells.
		scanpy.pp.scale(adata, max_value=10)
		## Reduce the dimensionality of the data by running principal component analysis (PCA), 
		## which reveals the main axes of variation and denoises the data.
		scanpy.pp.pca(adata, n_comps=npcs)
		return adata
	except Exception as e:
		Common.exception(e, 'could not run PCA')


def preprocess(adata, keyBatch):
	try:
		## Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, 
		## so that counts become comparable among cells. 
		scanpy.pp.normalize_total(adata, target_sum=1e4)
		## Logarithmize the data.
		scanpy.pp.log1p(adata)
		## Identify highly-variable genes.
		## Highly-variable genes are selected within each batch separately and merged.
		## Avoid the selection of batch-specific genes and acts as a lightweight batch correction method.
		scanpy.pp.highly_variable_genes(adata, flavor='seurat', subset=False, batch_key=keyBatch)
		## Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene expression 
		## for later use in differential testing and visualizations of gene expression. 
		## This simply freezes the state of the AnnData object.
		adata.raw = adata
		return adata
	except Exception as e:
		Common.exception(e, 'could not do the preprocessing')

