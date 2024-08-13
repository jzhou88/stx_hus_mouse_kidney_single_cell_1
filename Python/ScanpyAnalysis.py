#!/usr/bin/env python

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

import AuxFunc
import Common


def scanpyAnalysis(args):
	try:
		nameProj = args.p
		stampTime = args.t
		stepCurr = args.s
		resolCurr = args.r

		print('\nproject name: %s' % (nameProj))
		print('time stamp: %s' % (stampTime))
		print('current step: %s' % (stepCurr))
		print('current resolution: %s\n' % (resolCurr))

		cwd = '/home/zlab/projects/kidney'
		dirSrc = '/home/zlab/projects/kidney/src'

		dirOut = os.path.join(cwd, stampTime, nameProj)
		Common.mkDir(dirOut)

		dirFig = os.path.join(dirOut, 'figures')
		dirFigClu = os.path.join(dirFig, 'clusters')
		dirFigDEG = os.path.join(dirFig, 'DEG')
		dirFigPC = os.path.join(dirFig, 'PC')
		Common.mkDir(dirFig)
		Common.mkDir(dirFigClu)
		Common.mkDir(dirFigDEG)
		Common.mkDir(dirFigPC)

		dirRDS = os.path.join(dirOut, 'RDS')
		Common.mkDir(dirRDS)

		dirRes = os.path.join(dirOut, 'results')
		dirResClu = os.path.join(dirRes, 'clusters')
		dirResDEG = os.path.join(dirRes, 'DEG')
		Common.mkDir(dirRes)
		Common.mkDir(dirResClu)
		Common.mkDir(dirResDEG)

		scanpy.settings.verbosity = 3
		scanpy.logging.print_header()
		scanpy.settings.set_figure_params(dpi_save=400, vector_friendly=False, fontsize=14, facecolor='white')

		if (1 == stepCurr):
			pass

		elif (2 == stepCurr):
			with open(os.path.join(dirSrc, 'Integrate.py'), 'r') as f:
				exec(f.read(), globals(), locals())

		elif (3 == stepCurr):
			with open(os.path.join(dirSrc, 'Neighbor.py'), 'r') as f:
				exec(f.read(), globals(), locals())

		elif (4 == stepCurr):
			with open(os.path.join(dirSrc, 'Cluster.py'), 'r') as f:
				exec(f.read(), globals(), locals())
			
		elif (5 == stepCurr):
			with open(os.path.join(dirSrc, 'Annotate.py'), 'r') as f:
				exec(f.read(), globals(), locals())

		elif (6 == stepCurr):
			with open(os.path.join(dirSrc, 'DEG.py'), 'r') as f:
				exec(f.read(), globals(), locals())

		else:
			Common.error('could not identify Step \'%s\'' % (stepCurr))

		return 0
	except Exception as e:
		Common.exception(e, 'could not run Scanpy analysis')


def main():
	parser = argparse.ArgumentParser(description = 'run Scanpy analysis', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-p', required = True, type = str, help = 'specify project name')
	parser.add_argument('-t', required = True, type = str, help = 'specify time stamp')
	parser.add_argument('-s', required = True, type = int, help = 'select which step to run')
	parser.add_argument('-r', type = float, default = 0.0, help = 'optional: specify what resolution to use for clustering')
	args = parser.parse_args()
	scanpyAnalysis(args)
	return 0


if __name__ == '__main__':
	sys.exit(main())

