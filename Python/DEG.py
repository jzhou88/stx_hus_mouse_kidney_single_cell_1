## DEG.py
print('\nRunning \'DEG.py\' ...')

calcClusterMarkers = False
calcCelltypeDegs = True
calcCelltypeMarkers = False

if (calcClusterMarkers or calcCelltypeMarkers):
	resol = 1.1
	nameInput = f'after_clustering_res{resol:3.1f}.h5ad'
	keyGroup = f'leiden_res.{resol:3.1f}'
	keyAdd = f'rank_genes_leiden_res.{resol:3.1f}'
	groupsUse = 'all'
	refUse = 'rest'
	methodTest = 'wilcoxon'
	preSheet = 'Cluster '
	nameExcel = f'cluster_{methodTest}_markers_leiden_res.{resol:3.1f}.xlsx'
	nameOutput = f'after_DETesting_clusters_res{resol:3.1f}.h5ad'

	if (calcCelltypeMarkers):
		resol = 1.1
		nameInput = f'after_annotation_res{resol:3.1f}.h5ad'
		keyGroup = f'celltypes_leiden_res.{resol:3.1f}'
		keyAdd = f'rank_genes_celltypes_leiden_res.{resol:3.1f}'
		groupsUse = 'all'
		refUse = 'rest'
		methodTest = 'wilcoxon'
		preSheet = ''
		nameExcel = f'celltype_{methodTest}_markers_leiden_res.{resol:3.1f}.xlsx'
		nameOutput = f'after_DETesting_celltypes_res{resol:3.1f}.h5ad'

	fileInput = os.path.join(dirRDS, nameInput)
	Common.fileExist(fileInput)
	print('\nReading \'%s\' ...' % (fileInput))
	adata = scanpy.read(fileInput, cache=True)
	print(adata)
	print(adata.obs)
	print(adata.var)

	adata = AuxFunc.markers(adata, dirOut, keyGroup=keyGroup, keyAdd=keyAdd, nameExcel=nameExcel, \
							groupsUse=groupsUse, refUse=refUse, methodTest=methodTest, preSheet=preSheet)
	print(adata)
	print(adata.obs)
	print(adata.var)
	fileOutput = os.path.join(dirRDS, nameOutput)
	adata.write(fileOutput)
	Common.fileExist(fileOutput)

if calcCelltypeDegs:
	resol = 1.1
	nameInput = f'after_DETesting_celltypes_res{resol:3.1f}.h5ad'
	keyClu = f'celltypes_leiden_res.{resol:3.1f}'
	keyGroup = 'model'
	groupsUse = ['36h']
	refUse = '0h'
	methodTest = 'wilcoxon'
	nameExcel = f'celltype_{keyGroup}VS{refUse}_{methodTest}_DEGs_leiden_res.{resol:3.1f}.xlsx'

	fileInput = os.path.join(dirRDS, nameInput)
	Common.fileExist(fileInput)
	print('\nReading \'%s\' ...' % (fileInput))
	adata = scanpy.read(fileInput, cache=True)
	print(adata)
	print(adata.obs)
	print(adata.var)

	AuxFunc.degs(adata, dirOut, keyClu=keyClu, keyGroup=keyGroup, groupsUse=groupsUse, refUse=refUse, nameExcel=nameExcel, methodTest=methodTest)

print('\n\'DEG.py\' completed successfully!')

