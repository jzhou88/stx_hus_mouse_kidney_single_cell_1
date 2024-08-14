datasets.dir <- "/home/zlab/projects/datasets" # change here at different servers
cwd <- "/home/zlab/projects/kidney" # change here at different servers
src.dir <- "/home/zlab/projects/kidney/src" # change here at different servers
base::source(file = base::file.path(src.dir, "common.R"))

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_QC.RDS"))
sobj
head(sobj[[]], n=3L)
table(sobj$orig.ident)

sample.key <- "model"
curr.samples <- base::c("0h", "36h")
Seurat::Idents(object = sobj) <- sample.key
n.samples <- base::length(x = base::levels(x = sobj))
sample.colors <- scales::hue_pal()(n.samples)
base::names(x = sample.colors) <- curr.samples
features.to.plot <- base::c("nCount_RNA", "nFeature_RNA", "mitoRatio", "nGenePerUMI")
plots <- base::list()
for (feature.to.plot in features.to.plot) {
  plots[[feature.to.plot]] <- Seurat::VlnPlot(object = sobj, features = feature.to.plot, pt.size = 0, split.by = sample.key, ncol = 1) + 
    ggplot2::scale_fill_manual(values = sample.colors) + 
    ggplot2::theme(
      legend.position = "none", 
      axis.title.x = ggplot2::element_blank(), 
      plot.title = ggplot2::element_blank()
    )
}
metadata <- sobj@meta.data
metadata[[sample.key]] <- base::factor(x = metadata[[sample.key]], levels = curr.samples)
plots$nCell <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = sample.key, fill = sample.key)) + 
  ggplot2::geom_bar(color = "black", show.legend = F) + 
  ggplot2::theme_classic() + 
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.title = ggplot2::element_blank(), 
    plot.title = ggplot2::element_blank(), 
    text = ggplot2::element_text(size = 15, family = "sans", color = "black")
  )
w <- 25
h <- 5
base::options(repr.plot.width = w, repr.plot.height = h)
base::print(x = cowplot::plot_grid(plotlist = plots, nrow = 1, ncol = 5))
grDevices::pdf(file = base::file.path(fig.pub.dir, "QC_metrics.pdf"), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = cowplot::plot_grid(plotlist = plots, nrow = 1, ncol = 5))
grDevices::dev.off()

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_QC.RDS"))

Seurat2AnnData <- function(seurat.obj, save.dir, file.base, assay.to.use = "RNA") {
  trans.file = base::file.path(save.dir, base::sprintf(fmt = "%s.h5Seurat", file.base))
  if (base::file.exists(trans.file)) {
    base::file.remove(trans.file)
  }
  SeuratDisk::SaveH5Seurat(object = seurat.obj, filename = trans.file, overwrite = T, verbose = T)
  SeuratDisk::Convert(source = trans.file, dest = "h5ad", assay = assay.to.use, overwrite = T, verbose = T)
  if (base::file.exists(base::file.path(save.dir, base::sprintf(fmt = "%s.h5ad", file.base)))) {
    base::file.remove(trans.file)
  }
}

Seurat2AnnData(seurat.obj = sobj, save.dir = rds.dir, file.base = "after_QC", assay.to.use = "RNA")

sele.resol <- 1.1
orig.rds <- "/home/zlab/projects/kidney/2024-05-14/MouseStxKidneyCellRnaSeurat3.1/RDS/after_QC.RDS"
act.assay <- "RNA"
cluster.key <- base::sprintf(fmt = "leiden_res.%s", sele.resol)
celltype.key <- base::sprintf(fmt = "celltypes_leiden_res.%s", sele.resol)

adata <- anndata::read_h5ad(filename = base::file.path(rds.dir, base::sprintf(fmt = "after_DETesting_celltypes_res%s.h5ad", sele.resol)))
sobj <- base::readRDS(file = orig.rds)
sobj <- AnnData2Seurat(adata = adata, seurat.obj = sobj, act.assay = act.assay, cluster.key = cluster.key, celltype.key = celltype.key)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, base::sprintf(fmt = "after_annotation_res%s.RDS", sele.resol)))

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_annotation_res1.1.RDS"))
sobj
head(sobj[[]], n=3L)
sobj <- Seurat::NormalizeData(object = sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)

all.celltypes <- base::c(
  "Endo",
  "Podo",
  "PEC",
  "JG",
  "PT",
  "LOH",
  "DCT",
  "CD PC",
  "CD IC",
  "TEC",
  "Uro",
  "SMC",
  "Fib",
  "Mono",
  "Macro",
  "Neutro",
  "B lymph",
  "T lymph",
  "NK"
)
celltype.colors <- base::c(
  'B lymph'='#1f77b4',
  'CD IC'='#ff7f0e',
  'CD PC'='#279e68',
  'DCT'='#d62728',
  'Endo'='#aa40fc',
  'Fib'='#8c564b',
  'JG'='#e377c2',
  'LOH'='#b5bd61',
  'Macro'='#17becf',
  'Mono'='#aec7e8',
  'NK'='#ffbb78',
  'Neutro'='#98df8a',
  'PEC'='#ff9896',
  'PT'='#c5b0d5',
  'Podo'='#c49c94',
  'SMC'='#f7b6d2',
  'T lymph'='#dbdb8d',
  'TEC'='#9edae5',
  'Uro'='#ad494a'
)

sample.key <- "model"
split.sobj <- Seurat::SplitObject(object = sobj, split.by = sample.key)

features.to.plot <- base::c("Gsdmd", "Ripk3", "Mlkl")
fps <- base::list()
for (feature.to.plot in features.to.plot) {
  for (curr.sample in base::names(x = split.sobj)) {
    curr.sobj <- split.sobj[[curr.sample]]
    curr.name <- base::sprintf(fmt = "%s in %s", feature.to.plot, curr.sample)
    fps[[curr.name]] <- Seurat::FeaturePlot(object = curr.sobj, features = feature.to.plot, order = T, reduction = "umap", slot = "data", label = F) + 
      Seurat::NoAxes() + 
      viridis::scale_color_viridis(direction = -1, option = "D") + 
      ggplot2::labs(title = NULL) + 
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 32),
        legend.text = ggplot2::element_text(size = 18),
      )
  }
}
w <- 11
h <- 15
base::options(repr.plot.width = w, repr.plot.height = h)
base::print(x = cowplot::plot_grid(plotlist = fps, nrow = 3, ncol = 2))
grDevices::pdf(file = base::file.path(fig.pub.dir, base::sprintf(fmt = "featureplot_notitle_split%s_%s.pdf", sample.key, base::paste(features.to.plot, collapse = "_"))), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = cowplot::plot_grid(plotlist = fps, nrow = 3, ncol = 2))
grDevices::dev.off()

celltype.key <- "celltypes_leiden_res.1.1"
Seurat::Idents(object = sobj) <- celltype.key
Seurat::Idents(object = sobj) <- base::factor(x = Seurat::Idents(object = sobj), levels = all.celltypes)
levels(sobj)

markers.to.plot <- base::c(
  "Flt1", #"Emcn", # Endo
  "Nphs2", #"Nphs1", # Podo
  "Ncam1", # PEC
  "Ren1", # JG
  "Slc34a1", #"Slc13a1", "Slc4a4", PCT
  #"Slc7a12", # PST
  #"Slc27a2", "Lrp2", # PT
  #"Slc5a2", "Snhg11", "Slc5a12", # PT S1
  #"Slc22a6", "Slc13a3", # PT S2
  #"Slc7a13", "Atp11a", "Slc22a30", # PT S3
  "Umod", #"Slc12a1", # LOH
  "Slc12a3",# "Pvalb", # DCT
  "Hsd11b2", # "Aqp2", "H4c8", #"Krt18", "Cldn4", "Gdf15", "S100g", # "Scnn1g", "Aqp3", # CD PC
  "Atp6v1g3", #"Atp6v0d2", # CD IC
  "Slc14a2", #"Cryab", #"Cp", "Ncam1", #"Aqp1", # TEC
  "Upk1b", #"Upk3a", # Uro
  "Myh11", #"Acta2", "Tagln", "Rgs5", #"Rergl", "Map3k7cl", # SMC
  "Fbln5", #"Dcn", "Pdgfra", "Bgn", "Col1a2", "Col1a1", "Col3a1", # Fib
  #"Mki67",# "Stmn1", "Pclaf", "Top2a", "Birc5", "Cdca8", # Prolif
  #"Flt3",# "Ccl22", "Il4i1", "Cacnb3", "Ccr7", "Mreg", "Siglecg", "Fscn1", # "Cd209a", "Ccl17", "Clec10a", "Xcr1", # DC
  #"Plac8", "S100a4", # Mono
  "F13a1", #"Chil3", # Mono Ly6Chi
  #"Ace", "Treml4", # Mono Ly6Clo
  "C1qb", #"C1qa", "C1qc", # Macro
  "S100a9", #"S100a8", #"Retnlg", "Ly6g", "Cxcr2", "Ngp", "Mmp8", # Neutro
  "Cd79a", #"Ighm", "Bank1", "Cd79b", # "Igkc", "Ebf1", "Bach2", # B
  "Cd3g", #"Themis", "Itk", "Skap1", "Cd3d", "Cd3e", # T lymph
  #"Lef1",# "Tcf7", # naive T
  #"Cd8a",# "Cd8b1", "Eomes", # CD8 T
  #"Cd4", # CD4 Th
  #"Il12rb2",# "Txk", # Th1
  #"Il1rl1",# "Il17rb", "Areg", "Il9r", "Ccr8", "Gata3", # Th2
  #"Il23r",# "Il17a", # Th17
  #"Trgc2",# "Trgc4", "Trgv2", # GD T
  "Klrk1"# "Ncr1", "Nkg7", "Klrb1c", "Klre1"# "Klrb1a", "Klra8", "Klra4", "Klra7", "Klra9", "Klra1", "Klra3", "Klri2", "Klrb1f", "Gzma" # NK
)
c <- base::c("lightgrey", "blue")
#c <- base::c("white", "blue")
d <- 12
f <- 18
dp <- Seurat::DotPlot(object = sobj, assay = "RNA", features = markers.to.plot, cols = c, #col.min = 0,
                      dot.scale = d, cluster.idents = F, scale = T) +
  Seurat::RotatedAxis() +
  Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) +
  ggplot2::theme(
    text = ggplot2::element_text(size = f),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold"))
w <- 12
h <- 8.5
base::options(repr.plot.width = w, repr.plot.height = h)
base::print(x = dp)
grDevices::pdf(file = base::file.path(fig.pub.dir, base::sprintf(fmt = "dotplot_%s.pdf", celltype.key)), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

Seurat::Idents(object = sobj) <- "model"
dp <- Seurat::DimPlot(object = sobj, reduction = "umap", label = F) + 
  Seurat::NoAxes()
w <- 6
h <- 5
base::options(repr.plot.width = w, repr.plot.height = h)
base::print(x = dp)
grDevices::pdf(file = base::file.path(fig.pub.dir, "umap_model.pdf"), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

msdb.t2g <- msigdbr::msigdbr(species = "Mus musculus", category = "C5")

msdb.gobp.t2g <- msdb.t2g %>% dplyr::filter(gs_subcat == "GO:BP") %>% dplyr::select(gs_name, gene_symbol)

group.degs <- Xlsx2List(xlsx.file = base::file.path(res.deg.dir, "celltype_modelVS0h_wilcoxon_DEGs_leiden_res.1.1.xlsx"))

curr.group <- "36h in Endo"
curr.degs <- base::as.data.frame(x = group.degs[[curr.group]])
gene.logfc <- curr.degs$logfoldchanges
base::names(x = gene.logfc) <- curr.degs$names
gene.logfc <- base::sort(x = gene.logfc, decreasing = T)
length(gene.logfc)
head(gene.logfc)

gs <- clusterProfiler::GSEA(
  geneList = gene.logfc,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  TERM2GENE = msdb.gobp.t2g,
  nPermSimple = 100000
)

w <- 7
h <- 5
base::options(repr.plot.width = w, repr.plot.height = h)
gp <- enrichplot::gseaplot2(gs, geneSetID = which(as.data.frame(gs)$ID == "GOBP_PYROPTOSIS"), title = gs$Description[which(as.data.frame(gs)$ID == "GOBP_PYROPTOSIS")])
base::print(x = gp)
grDevices::pdf(file = base::file.path(fig.pub.dir, "gseaplot_pyroptosis.pdf"), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = gp)
grDevices::dev.off()

w <- 7
h <- 5
base::options(repr.plot.width = w, repr.plot.height = h)
gp <- enrichplot::gseaplot2(gs, geneSetID = which(as.data.frame(gs)$ID == "GOBP_NECROPTOTIC_SIGNALING_PATHWAY"), title = gs$Description[which(as.data.frame(gs)$ID == "GOBP_NECROPTOTIC_SIGNALING_PATHWAY")])
base::print(x = gp)
grDevices::pdf(file = base::file.path(fig.pub.dir, "gseaplot_necroptosis.pdf"), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = gp)
grDevices::dev.off()
