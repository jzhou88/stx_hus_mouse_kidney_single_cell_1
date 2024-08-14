# QC.R
base::cat("\nRunning \"QC.R\" ...\n")

min.cell <- 10
mt.pattern <- "^mt-"
max.mito <- 0.5
min.umi <- 500
max.umi <- 15000
min.gene <- 200
max.gene <- 4000
min.complex <- 0.25

cellranger.ver <- "cellranger-8.0.0"
cellranger.dir <- base::file.path(datasets.dir, "matrix", cellranger.ver)

split.sobj <- base::list()
for (curr.mouse in curr.mice) {
  base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.mouse))
  curr.outs.dir <- base::file.path(cellranger.dir, curr.mouse, "outs")
  # 1_RnaQC.R
  base::cat("\nRunning \"1_RnaQC.R\" ...\n")
  
  ### Run SoupX to remove ambient mRNA ------
  adj.counts <- NULL
  adj.count.file <- base::file.path(rds.dir, base::sprintf(fmt = "soupx_corrected_counts_%s.RDS", curr.mouse))
  if (base::file.exists(adj.count.file)) {
    adj.counts <- base::readRDS(file = adj.count.file)
  } else {
    sc <- SoupX::load10X(dataDir = curr.outs.dir)
    grDevices::jpeg(file = base::file.path(fig.qc.dir, base::sprintf("SoupX_estimates_density_%s.jpg", curr.mouse)), 
                    width = 5, height = 5, units = "in", pointsize = 12, quality = 100, res = 300)
    sc <- SoupX::autoEstCont(sc = sc, doPlot = T, forceAccept = T, verbose = F)
    grDevices::dev.off()
    adj.counts <- SoupX::adjustCounts(sc = sc, roundToInt = T)
    base::saveRDS(object = adj.counts, file = adj.count.file)
  }
  
  ### Run Seurat QC to remove low-quality genes and cells ------
  sobj <- Seurat::CreateSeuratObject(counts = adj.counts, project = curr.mouse, assay = "RNA", min.cells = min.cell, min.features = 0)
  sobj$mitoRatio <- Seurat::PercentageFeatureSet(object = sobj, pattern = mt.pattern)
  sobj$mitoRatio <- sobj$mitoRatio / 100
  sobj$nGenePerUMI <- sobj$nFeature_RNA / sobj$nCount_RNA
  sobj$cellID <- base::paste(sobj$orig.ident, base::rownames(x = sobj@meta.data), sep = "_")
  w <- 20
  h <- 5
  vp <- Seurat::VlnPlot(object = sobj, features = c("nCount_RNA", "nFeature_RNA", "mitoRatio", "nGenePerUMI"), ncol = 4, pt.size = 0)
  grDevices::pdf(file = base::file.path(fig.qc.dir, base::sprintf("QC_metrics_before_Seurat_QC_%s.pdf", curr.mouse)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = vp)
  grDevices::dev.off()
  vp <- Seurat::VlnPlot(object = sobj, features = c("nCount_RNA", "nFeature_RNA", "mitoRatio", "nGenePerUMI"), ncol = 4, pt.size = 0.4)
  grDevices::jpeg(file = base::file.path(fig.qc.dir, base::sprintf("QC_metrics_before_Seurat_QC_%s.jpg", curr.mouse)), 
                  width = w, height = h, units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = vp)
  grDevices::dev.off()
  sobj <- subset(x = sobj, subset = ((mitoRatio < max.mito) & (nCount_RNA > min.umi) & (nCount_RNA < max.umi) & 
                                       (nFeature_RNA > min.gene) & (nFeature_RNA < max.gene) & (nGenePerUMI > min.complex)))
  vp <- Seurat::VlnPlot(object = sobj, features = c("nCount_RNA", "nFeature_RNA", "mitoRatio", "nGenePerUMI"), ncol = 4, pt.size = 0)
  grDevices::pdf(file = base::file.path(fig.qc.dir, base::sprintf("QC_metrics_after_Seurat_QC_%s.pdf", curr.mouse)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = vp)
  grDevices::dev.off()
  vp <- Seurat::VlnPlot(object = sobj, features = c("nCount_RNA", "nFeature_RNA", "mitoRatio", "nGenePerUMI"), ncol = 4, pt.size = 0.4)
  grDevices::jpeg(file = base::file.path(fig.qc.dir, base::sprintf("QC_metrics_after_Seurat_QC_%s.jpg", curr.mouse)), 
                  width = w, height = h, units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = vp)
  grDevices::dev.off()
  
  sobj.before.df <- sobj
  
  ### Run DoubletFinder to remove doublets ------
  ## Normalize the data
  sobj <- Seurat::NormalizeData(object = sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  ## Identify highly variable genes (feature selection)
  sobj <- Seurat::FindVariableFeatures(object = sobj, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = F)
  ## Scale the data
  sobj <- Seurat::ScaleData(object = sobj, assay = "RNA", verbose = F)
  ## Perform linear dimensional reduction
  sobj <- Seurat::RunPCA(object = sobj, assay = "RNA", npcs = 50, verbose = F)
  ## Cluster the cells
  ndims <- 10
  sobj <- Seurat::FindNeighbors(object = sobj, reduction = "pca", dims = 1:ndims, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.8, verbose = F)
  ## Run non-linear dimensional reduction
  sobj <- Seurat::RunUMAP(object = sobj, dims = 1:ndims, reduction = "pca", umap.method = "uwot", metric = "cosine", verbose = F)
  ## pK Identification (no ground-truth)
  sweep.res.list <- DoubletFinder::paramSweep(seu = sobj, PCs = 1:ndims, sct = F)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.list = sweep.res.list, GT = F)
  grDevices::jpeg(file = base::file.path(fig.qc.dir, base::sprintf("DoubletFinder_BCmetric_pK_%s.jpg", curr.mouse)), 
                  width = 7, height = 4, units = "in", pointsize = 12, quality = 100, res = 300)
  bcmvn <- DoubletFinder::find.pK(sweep.stats = sweep.stats)
  grDevices::dev.off()
  base::sink(file = base::file.path(res.qc.dir, base::sprintf(fmt = "DoubletFinder_BCmvn_%s.txt", curr.mouse)), append = F, split = F)
  base::print(x = bcmvn)
  base::sink()
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations = sobj$seurat_clusters)
  nCell <- base::ncol(x = adj.counts)
  base::cat(base::sprintf(fmt = "\nNumber of Cells=%s\n", nCell))
  doublet.rate <- 0.004 * (nCell / 500)
  base::cat(base::sprintf(fmt = "Doublet Rate=%s\n", doublet.rate))
  nExp_poi <- base::round(doublet.rate * base::ncol(x = sobj))
  nExp_poi.adj <- base::round(nExp_poi * (1 - homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies
  select.pK <- base::as.numeric(x = base::as.character(x = bcmvn$pK[base::which.max(x = bcmvn$BCmetric)]))
  base::cat(base::sprintf(fmt = "\npK=%s\n", select.pK))
  sobj <- DoubletFinder::doubletFinder(seu = sobj, PCs = 1:ndims, pN = 0.25, pK = select.pK, nExp = nExp_poi, reuse.pANN = F, sct = F)
  pANN.to.reuse <- base::grep(pattern = "^pANN", x = base::colnames(x = sobj@meta.data), value = T)
  init.class <- base::grep(pattern = "^DF.class", x = base::colnames(x = sobj@meta.data), value = T)
  base::cat("\n")
  base::print(x = base::colnames(x = sobj@meta.data))
  base::cat(base::sprintf(fmt = "\npANN.to.reuse=%s\n", pANN.to.reuse))
  base::cat(base::sprintf(fmt = "init.class=%s\n", init.class))
  sobj <- DoubletFinder::doubletFinder(seu = sobj, PCs = 1:ndims, pN = 0.25, pK = select.pK, nExp = nExp_poi.adj, 
                                       reuse.pANN = pANN.to.reuse, sct = F)
  ## Identify doublets
  final.class <- base::grep(pattern = "^DF.class", x = base::colnames(x = sobj@meta.data), value = T)
  final.class <- final.class[final.class != init.class]
  base::cat("\n")
  base::print(x = base::colnames(x = sobj@meta.data))
  base::cat(base::sprintf(fmt = "\nfinal.class=%s\n", final.class))
  base::print(x = base::table(sobj[[final.class]]))
  doublets <- sobj[[final.class]]
  sobj.before.df <- Seurat::AddMetaData(object = sobj.before.df, metadata = doublets, col.name = "doublet")
  base::cat("\n")
  base::print(x = base::table(sobj.before.df$doublet))
  
  split.sobj[[curr.mouse]] <- sobj.before.df
}
sobj <- merge(x = split.sobj[[1]], y = split.sobj[2:base::length(x = split.sobj)], 
              add.cell.ids = curr.mice, merge.data = T, project = proj.name)
split.sobj <- NULL
base::gc()
if (base::ncol(x = sobj) != base::sum(base::rownames(x = sobj@meta.data) == sobj$cellID)) {
  base::stop("cell IDs are not the same.")
}
base::cat("\n")
base::print(x = sobj)
base::print(x = utils::head(x = sobj[[]], n = 3L))
base::print(x = base::table(sobj$doublet))
base::print(x = base::table(sobj$orig.ident))
sobj <- subset(x = sobj, subset = (doublet == "Singlet"))
base::cat("\n")
base::print(x = sobj)
base::print(x = utils::head(x = sobj[[]], n = 3L))
base::print(x = base::table(sobj$doublet))
base::print(x = base::table(sobj$orig.ident))
sobj$doublet <- NULL
base::cat("\n")
base::print(x = utils::head(x = sobj[[]], n = 3L))
models <- sobj$orig.ident
base::cat("\n")
base::print(x = base::table(models))
base::print(x = utils::head(x = models))
cell.names <- base::names(x = models)
base::cat("\n")
base::print(x = utils::head(x = cell.names))
models <- base::as.character(x = models)
base::cat("\n")
base::print(x = base::table(models))
base::print(x = utils::head(x = models))
models[models %in% base::c("C1")] <- "0h"
models[models %in% base::c("ST8")] <- "36h"
base::cat("\n")
base::print(x = base::table(models))
base::print(x = utils::head(x = models))
base::names(x = models) <- cell.names
base::cat("\n")
base::print(x = utils::head(x = models))
sobj <- Seurat::AddMetaData(object = sobj, metadata = models, col.name = "model")
base::cat("\n")
base::print(x = utils::head(x = sobj[[]], n = 3L))
base::print(x = base::table(sobj$model))
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_QC.RDS"))

base::cat("\n\"1_QC.R\" completed successfully!\n")
