## Annotate.py
print('\nRunning \'Annotate.py\' ...')

resol = 1.1
keyAdd = f'celltypes_leiden_res.{resol:3.1f}'
keyClu = f'leiden_res.{resol:3.1f}'
clus2cell = {
	'0': 'PT', \
	'1': 'Endo', \
	'2': 'PT', \
	'3': 'T lymph', \
	'4': 'PT', \
	'5': 'CD IC', \
	'6': 'CD IC', \
	'7': 'DCT', \
	'8': 'Macro', \
	'9': 'PT', \
	'10': 'Neutro', \
	'11': 'Endo', \
	'12': 'Endo', \
	'13': 'SMC', \
	'14': 'LOH', \
	'15': 'Uro', \
	'16': 'Macro', \
	'17': 'JG', \
	'18': 'TEC', \
	'19': 'Mono', \
	'20': 'B lymph', \
	'21': 'T lymph', \
	'22': 'T lymph', \
	'23': 'NK', \
	'24': 'T lymph', \
	'25': 'PEC', \
	'26': 'Endo', \
	'27': 'Fib', \
	'28': 'Podo', \
	'29': 'CD PC'
}

fileInput = os.path.join(dirRDS, 'after_DETesting_clusters_res%s.h5ad' % (resol))
Common.fileExist(fileInput)
print('\nReading \'%s\' ...' % (fileInput))
adata = scanpy.read(fileInput, cache=True)
print(adata)
print(adata.obs)
print(adata.var)

adata = AuxFunc.annotation(adata, dirOut, keyAdd=keyAdd, keyClu=keyClu, clus2cell=clus2cell)
print(adata)
print(adata.obs)
print(adata.var)
fileOutput = os.path.join(dirRDS, 'after_annotation_res%s.h5ad' % (resol))
adata.write(fileOutput)
Common.fileExist(fileOutput)

print('\n\'Annotate.py\' completed successfully!')

