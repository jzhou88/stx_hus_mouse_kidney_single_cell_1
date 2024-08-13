## Cluster.py
print('\nRunning \'Cluster.py\' ...')

saveClusters = True
splitClusters = True

fileInput = os.path.join(dirRDS, 'after_neighbor.h5ad')
Common.fileExist(fileInput)
print('\nReading \'%s\' ...' % (fileInput))
adata = scanpy.read(fileInput, cache=True)
print(adata)
print(adata.obs)
print(adata.var)

adata = AuxFunc.clustering(adata, resolCurr, dirOut, splitClusters)
print(adata)
print(adata.obs)
print(adata.var)
if (saveClusters):
	fileOutput = os.path.join(dirRDS, 'after_clustering_res%s.h5ad' % (resolCurr))
	adata.write(fileOutput)
	Common.fileExist(fileOutput)

print('\n\'Cluster.py\' completed successfully!')

