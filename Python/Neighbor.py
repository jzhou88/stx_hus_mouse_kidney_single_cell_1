## Neighbor.py
print('\nRunning \'Neighbor.py\' ...')

repUse = 'X_harmony'
npcs = 51 # 9,503 cells in total
keyBatch = 'orig.ident'

fileInput = os.path.join(dirRDS, 'after_harmony.h5ad')
Common.fileExist(fileInput)
print('\nReading \'%s\' ...' % (fileInput))
adata = scanpy.read(fileInput, cache=True)
print(adata)
print(adata.obs)
print(adata.var)

adata = AuxFunc.neighbor(adata, dirOut, repUse=repUse, npcs=npcs, keyBatch=keyBatch)
print(adata)
print(adata.obs)
print(adata.var)
fileOutput = os.path.join(dirRDS, 'after_neighbor.h5ad')
adata.write(fileOutput)
Common.fileExist(fileOutput)

print('\n\'Neighbor.py\' completed successfully!')

