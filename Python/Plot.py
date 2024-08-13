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
import scvelo
import seaborn
import string
import struct
import sys
import time

import AuxFunc
import Common

dirOut = os.path.join(cwd, stampTime, nameProj)
Common.mkDir(dirOut)

dirFig = os.path.join(dirOut, 'figures')
dirFigClu = os.path.join(dirFig, 'clusters')
dirFigDEG = os.path.join(dirFig, 'DEG')
dirFigPC = os.path.join(dirFig, 'PC')
dirFigPub = os.path.join(dirFig, 'publication')
Common.mkDir(dirFig)
Common.mkDir(dirFigClu)
Common.mkDir(dirFigDEG)
Common.mkDir(dirFigPC)
Common.mkDir(dirFigPub)

dirRDS = os.path.join(dirOut, 'RDS')
Common.mkDir(dirRDS)

dirRes = os.path.join(dirOut, 'results')
dirResCCC = os.path.join(dirRes, 'CCC')
dirResClu = os.path.join(dirRes, 'clusters')
dirResDEG = os.path.join(dirRes, 'DEG')
dirResPub = os.path.join(dirRes, 'publication')
Common.mkDir(dirRes)
Common.mkDir(dirResCCC)
Common.mkDir(dirResClu)
Common.mkDir(dirResDEG)
Common.mkDir(dirResPub)

scanpy.settings.verbosity = 3
scanpy.logging.print_header()
scanpy.settings.set_figure_params(dpi=150, dpi_save=400, vector_friendly=False, fontsize=14, facecolor='white')

fileInput = os.path.join(dirRDS, 'after_DETesting_celltypes_res1.1.h5ad')
adata = scanpy.read(fileInput, cache=True)
print(adata)
print(adata.obs.head(3))

keyCell = 'celltypes_leiden_res.1.1'
scanpy.pl.umap(adata, color=keyCell, use_raw=True, frameon=False, legend_loc='on data', legend_fontsize=7.5, legend_fontoutline=1, add_outline=True, title='')

adataEndo = adata[adata.obs[keyCell] == 'Endo'].copy()
with matplotlib.pyplot.rc_context({"figure.figsize": (3, 4)}):
    scanpy.pl.violin(adataEndo, ["Gsdmd"], groupby='model', use_raw=True, size=3, rotation=0)
    scanpy.pl.violin(adataEndo, ["Ripk3"], groupby='model', use_raw=True, size=3, rotation=0)
    scanpy.pl.violin(adataEndo, ["Mlkl"], groupby='model', use_raw=True, size=3, rotation=0)

