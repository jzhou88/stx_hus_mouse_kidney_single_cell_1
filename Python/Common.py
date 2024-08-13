## Common.py

import os
import sys


markersKnownMouse = {'GEC/Endo': ['Nrp1', 'Kdr', 'Ehd3', 'Igfbp3'], \
					 'Podo': ['Nphs1', 'Nphs2'], \
					 'PT(S1/S2/S3)': ['Slc27a2', 'Lrp2', 'Slc5a2', 'Snhg11', 'Slc5a12', 'Slc7a13', 'Atp11a', 'Slc22a30'], \
					 'DLOH': ['Slc4a11', 'Ptgds', 'Cp'], \
					 'ALOH': ['Slc12a1', 'Ppp1r1b'], \
					 'DCT/CNT': ['Slc12a3', 'Pvalb', 'Trpv5'], \
					 'CD PC': ['Aqp2', 'Hsd11b2'], \
					 'CD Trans': ['Insrr', 'Rhbg'], \
					 'CD IC(A-IC/B-IC)': ['Atp6v1g3', 'Atp6v0d2', 'Slc4a1', 'Aqp6', 'Slc26a4', 'Hmx2'], \
					 'Fib': ['Col1a1', 'Col3a1', 'Vim', 'Fn1'], \
					 'Mono(Ly6Chi/Ly6Clo)': ['Plac8', 'S100a4', 'F13a1', 'Chil3', 'Ace', 'Treml4'], \
					 'DC(11b+/11b-)/pDC': ['Clec10a', 'Cd209a', 'Xcr1', 'Siglech', 'Ccr9'], \
					 'Baso': ['Fcer1a', 'Mcpt8', 'Csrp3', 'Cd200r3', 'Cyp11a1', 'Ms4a2'], \
					 'Macro': ['C1qa', 'C1qb', 'C1qc'], \
					 'Granul': ['S100a8', 'S100a9', 'Hp'], \
					 'B lymph(B1/B2)': ['Cd79a', 'Cd79b', 'Igkc', 'Ebf1', 'Ighm', 'Igha', 'Jchain', 'Derl3'], \
					 'T lymph': ['Ltb', 'Cxcr6', 'Lef1', 'Cd28', 'Icos', 'Rora', 'Actn2', 'Ly6g5b', 'Cd3g', 'Cd3d', 'Cd8a'], \
					 'NK': ['Ncr1', 'Ccl5', 'Nkg7', 'Gzma'], \
					 'Proliferating': ['Mki67', 'Cdca3', 'Stmn1', 'Lockd', 'Top2a'], \
					 'Podocyte': ['Pdpn', 'Wt1', 'Mafb', 'Synpo', 'Cdkn1c', 'Ptpro'], \
					 'PEC': ['Cldn1', 'Pax8'], \
					 'Mesangial': ['Pdgfrb', 'Gata3', 'Des', 'Itga8', 'Plvap', 'Prkca', 'Art3', 'Nt5e'], \
					 'Endothelial': ['Flt1', 'Tie1', 'Pecam1', 'Emcn'], \
					 'SMC': ['Acta2', 'Myh11', 'Tagln', 'Ren1', 'Akr1b7', 'Rgs5', 'Rergl', 'Map3k7cl'], \
					 'Immune': ['Ptprc', 'Lyz1', 'Csf1r', 'Itgam', 'Ms4a1'], \
					 'TEC': ['Fxyd2', 'Slc14a2', 'Aqp1', 'Umod']}


def dirExist(dirname):
	if not(os.path.isdir(dirname)):
		error('dir \'%s\' does not exist' % (dirname))


def error(msg):
	print('Error: %s.' % (msg))
	sys.exit(-1)


def exception(e, msg):
	print('Exception:', e)
	error(msg)


def fileExist(filename):
	if not(os.path.isfile(filename)):
		error('file \'%s\' does not exist' % (filename))


def mkDir(pathDir):
	try:
		if not(os.path.isdir(pathDir)):
			os.mkdir(pathDir)
		return 0
	except Exception as e:
		exception(e, 'could not make directory \'%s\'' % (pathDir))


def readLinesFromFile(filename, chkLine = True):
	try:
		fileExist(filename)
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()
		if chkLine and (len(lines) <= 0):
			error('file \'%s\' contains nothing' % (filename))
		return lines
	except Exception as e:
		exception(e, 'could not read file \'%s\'' % (filename))

