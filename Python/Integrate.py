## Integrate.py
print('\nRunning \'Integrate.py\' ...')

keyBatch = 'orig.ident'
npcs = 51 # 9,503 cells in total

fileInput = '/home/zlab/projects/kidney/2024-05-14/MouseStxKidneyCellRnaSeurat3.1/RDS/after_QC.h5ad'
Common.fileExist(fileInput)
print('\nReading \'%s\' ...' % (fileInput))
adata = scanpy.read(fileInput, cache=True)
print(adata)
print(adata.obs)
print(adata.var)

adata = AuxFunc.preprocess(adata, keyBatch)
print(adata)
print(adata.obs)
print(adata.var)
#fileOutput = os.path.join(dirRDS, 'after_HVG.h5ad')
#adata.write(fileOutput)
#Common.fileExist(fileOutput)

adata = AuxFunc.pca(adata, npcs)
print(adata)
print(adata.obs)
print(adata.var)
#fileOutput = os.path.join(dirRDS, 'after_pca.h5ad')
#adata.write(fileOutput)
#Common.fileExist(fileOutput)

adata = AuxFunc.harmony(adata, keyBatch=keyBatch)
print(adata)
print(adata.obs)
print(adata.var)
fileOutput = os.path.join(dirRDS, 'after_harmony.h5ad')
adata.write(fileOutput)
Common.fileExist(fileOutput)

print('\n\'Integrate.py\' completed successfully!')

