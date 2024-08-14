# common.R
base::cat("\nRunning \"common.R\" ...\n")

### Global settings ------
base::options(max.print = 99999, warn = 1)
base::suppressPackageStartupMessages({
  base::library(package = "anndata")
  base::library(package = "AnnotationHub")
  base::library(package = "Biobase")
  base::library(package = "BiocGenerics")
  base::library(package = "CellChat")
  base::library(package = "clusterProfiler")
  base::library(package = "ComplexHeatmap")
  base::library(package = "COSG")
  base::library(package = "cowplot")
  base::library(package = "data.table")
  base::library(package = "destiny")
  base::library(package = "DoubletFinder")
  base::library(package = "dplyr")
  base::library(package = "ensembldb")
  base::library(package = "future")
  base::library(package = "ggplot2")
  base::library(package = "grDevices")
  base::library(package = "hdWGCNA")
  base::library(package = "MAST")
  base::library(package = "methods")
  base::library(package = "monocle")
  base::library(package = "msigdbr")
  base::library(package = "openxlsx")
  base::library(package = "org.Mm.eg.db")
  base::library(package = "patchwork")
  base::library(package = "presto")
  base::library(package = "reshape2")
  base::library(package = "reticulate")
  base::library(package = "S4Vectors")
  base::library(package = "scales")
  base::library(package = "scITD")
  base::library(package = "Seurat")
  base::library(package = "SeuratDisk")
  base::library(package = "SoupX")
  base::library(package = "stringr")
  base::library(package = "tibble")
  base::library(package = "utils")
  base::library(package = "velocyto.R")
  base::library(package = "VGAM")
  base::library(package = "viridis")
  base::library(package = "WGCNA")
})
base::set.seed(seed = 1)

### Markers ------
known.mouse.markers <- base::c(
  "Nrp1", "Kdr", # GEC, Endo
  "Ehd3", # GEC
  "Igfbp3", # Endo
  "Nphs1", "Nphs2", # Podo
  "Slc27a2", "Lrp2", # PT, "Aqp1" is also a good marker
  "Slc34a1", # PCT
  "Slc7a12", # PST
  "Slc5a2", "Snhg11", "Slc5a12", # PT S1
  "Slc22a6", "Slc13a3", # PT S2
  "Slc7a13", "Atp11a", "Slc22a30", # PT S3
  "Fxyd5", "Il1b", "Cxcl2", "Ccl3", "Tyrobp", # Maladaptive PT (DOI: 10.1038/s41467-022-31772-9, Fig. 3b)
  "C3", "Vcam1", "Sparc", # Injured PT (DOI: 10.1038/s41467-022-31772-9, Fig. 3b)
  "Havcr1", "Gsta1", "Nqo1", "Krt20", # Injured PT (DOI: 10.1038/s41467-022-31772-9, Fig. 3b)
  "Nupr1", "Akap12", "Ankrd1", "Ifit3", # Injured PT (DOI: 10.1038/s41467-022-31772-9, Fig. 3b)
  "Ptgds", # Injured PT / DLOH (DOI: 10.1038/s41467-022-31772-9, Fig. 3b)
  "Slc4a11", "Cp", "Cryab", # DLOH
  "Fst", # DTL
  "Atp10b", # tAL
  "Ppp1r1b", "Slc12a1", # ALOH, TAL
  "Nos1", # Macular densa (MD)
  "Pvalb", "Slc12a3", # DCT, CNT
  "Trpv5", "Slc8a1", # CNT
  "Aqp2", "Hsd11b2", # PC
  "Insrr", "Rhbg", # Trans
  "Atp6v1g3", "Atp6v0d2", # IC
  "Slc4a1", "Aqp6", "Kit", # A IC
  "Slc26a4", "Hmx2", # B IC
  "Col1a1", "Col3a1", "Vim", "Fn1", # Fib
  "Plac8", "S100a4", # Mono
  "F13a1", "Chil3", # Mono Ly6Chi
  "Ace", "Treml4", # Mono Ly6Clo
  "Clec10a", # DC
  "Cd209a", # DC 11b+
  "Xcr1", # DC 11b-
  "Siglech", "Ccr9", # pDC
  "Fcer1a", "Mcpt8", "Csrp3", "Cd200r3", "Cyp11a1", "Ms4a2", # Baso
  "C1qa", "C1qb", "C1qc", # Macro
  "S100a8", "S100a9", "Hp", # Granul
  "Cd79a", "Cd79b", "Igkc", # B lymph
  "Ebf1", "Ighm", # B1
  "Igha", "Jchain", "Derl3", # B2
  "Ltb", "Cxcr6", # T lymph
  "Lef1", # T naive
  "Cd28", # T mem
  "Icos", "Rora", "Actn2", "Ly6g5b", # T gd
  "Cd3g", "Cd3d", "Cd8a", # CD8 effector
  "Ncr1", "Ccl5", "Nkg7", # NK
  "Gzma", # NK2
  "Mki67", "Cdca3", "Stmn1", "Lockd", # Novel
  "Top2a", # Novel2
  "Pdpn", "Wt1", "Mafb", "Synpo", "Cdkn1c", "Ptpro", # podocyte (DOI: 10.1681/asn.2020020220)
  "Cldn1", "Pax8", # + podocyte markers, parietal epithelial cell (PEC)
  "Ncam1", # PEC
  "Flt1", "Tie1", "Pecam1", "Emcn", # endothelial (DOI: 10.1681/asn.2020020220)
  "Gata3", "Des", "Itga8", # mesangial (DOI: 10.1681/asn.2020020220)
  "Plvap", "Prkca", "Art3", "Nt5e", # mesangial (DOI: 10.1681/asn.2020020220, Fig. 1)
  "Pdgfrb", # Fib, mesangial, juxtaglomerular apparatus (JGA)
  "Ren1", # JG, SMC
  "Acta2", "Myh11", "Tagln", # smooth muscle cell (SMC) (DOI: 10.1681/asn.2020020220)
  "Akr1b7", "Rgs5", "Rergl", "Map3k7cl", # SMC (DOI: 10.1681/asn.2020020220, Fig. 1)
  "Ptprc", "Lyz1", "Csf1r", "Itgam", "Ms4a1", # immune (DOI: 10.1681/asn.2020020220)
  "Fxyd2", "Slc14a2", "Aqp1", "Umod" # tubular epithelial cell (TEC) (DOI: 10.1681/asn.2020020220)
)

### Local settings ------
picked.name <- proj.name # change here to analyze different data
if (!(picked.name %in% base::names(x = all.mice))) {
  base::stop(base::sprintf("could not identify the picked name \"%s\".", picked.name))
}
curr.mice <- all.mice[[picked.name]]
out.dir <- base::file.path(cwd, time.stamp, picked.name)
base::dir.create(path = out.dir, showWarnings = T, recursive = T, mode = "0755")

fig.dir <- base::file.path(out.dir, "figures")
fig.ccc.dir <- base::file.path(fig.dir, "CCC")
fig.clu.dir <- base::file.path(fig.dir, "clusters")
fig.deg.dir <- base::file.path(fig.dir, "DEG")
fig.des.dir <- base::file.path(fig.dir, "trajectory", "destiny")
fig.mo2.dir <- base::file.path(fig.dir, "trajectory", "monocle2")
fig.pc.dir <- base::file.path(fig.dir, "PC")
fig.pub.dir <- base::file.path(fig.dir, "publication")
fig.qc.dir <- base::file.path(fig.dir, "QC")
base::dir.create(path = fig.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.ccc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.clu.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.deg.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.des.dir, showWarnings = T, recursive = T, mode = "0755")
base::dir.create(path = fig.mo2.dir, showWarnings = T, recursive = T, mode = "0755")
base::dir.create(path = fig.pc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.pub.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.qc.dir, showWarnings = T, recursive = F, mode = "0755")

rds.dir <- base::file.path(out.dir, "RDS")
base::dir.create(path = rds.dir, showWarnings = T, recursive = F, mode = "0755")

res.dir <- base::file.path(out.dir, "results")
res.clu.dir <- base::file.path(res.dir, "clusters")
res.deg.dir <- base::file.path(res.dir, "DEG")
res.mo2.dir <- base::file.path(res.dir, "trajectory", "monocle2")
res.pat.dir <- base::file.path(res.dir, "pathway")
res.pub.dir <- base::file.path(res.dir, "publication")
res.qc.dir <- base::file.path(res.dir, "QC")
res.wgc.dir <- base::file.path(res.dir, "WGCNA")
base::dir.create(path = res.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = res.clu.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = res.deg.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = res.mo2.dir, showWarnings = T, recursive = T, mode = "0755")
base::dir.create(path = res.pat.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = res.pub.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = res.qc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = res.wgc.dir, showWarnings = T, recursive = F, mode = "0755")
