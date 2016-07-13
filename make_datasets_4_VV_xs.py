from ROOT import *

import os
from optparse import OptionParser

ROOT.gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load("PDFs/Util_cxx.so")
ROOT.gSystem.Load("PDFs/hyperg_2F1_c.so")
ROOT.gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

parser	= OptionParser()
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ, defines signal region')
parser.add_option('--ch', dest='ch', default='el', help='channel, el or mu')
(options,args) = parser.parse_args()


ch = options.ch
cat = options.cat



fileInWS	= TFile.Open('cards_%s_HP%s/mj_closurepost.root'%(ch,cat[1]))
w		= fileInWS.Get('w4fitpost')
fileInWS.Close()
fileInWS	= TFile.Open('cards_%s_HP%s/mj_workspace.root'%(ch,cat[1]))
w_fits		= fileInWS.Get('workspace4fit_')
fileInWS.Close()
rrv_x		= w.var('rrv_mass_j')
VV_pdf		= w.pdf('model_VV_xww_%s_mj'%ch)
STop_pdf	= w.pdf('model_STop_xww_%s_mj'%ch)
TTbar_pdf	= w.pdf('model_TTbar_xww_%s_mj'%ch)
WJets_pdf	= w.pdf('model_WJets0_xww_%s_mj'%ch)
VV_norm		= w.var('rrv_number_VV_xww_%s_mj'%ch)
STop_norm	= w.var('rrv_number_STop_xww_%s_mj'%ch)
TTbar_norm	= w.var('rrv_number_TTbar_xww_%s_mj'%ch)
WJets_norm	= w.var('rrv_number_WJets0_xww_%s_mj'%ch)

VV_norm_val	= VV_norm.getVal()
STop_norm_val	= STop_norm.getVal()
TTbar_norm_val	= TTbar_norm.getVal()
WJets_norm_val	= WJets_norm.getVal()

TTbar_norm_err	= TTbar_norm.getError()

random		= TRandom2()

fileOut		= TFile.Open("VV_xs/datasets_%s_HP%s.root"%(ch,cat[1]),"recreate")

for i in range(500):
	VV_norm_tmp	= random.Poisson(VV_norm_val)
	STop_norm_tmp	= random.Poisson(STop_norm_val)
	TTbar_norm_tmp	= random.Poisson(TTbar_norm_val)
	WJets_norm_tmp	= random.Poisson(WJets_norm_val)
	#VV_norm_tmp	= int(VV_norm_val)
	#STop_norm_tmp	= int(STop_norm_val)
	#TTbar_norm_tmp = int(TTbar_norm_val)
	#WJets_norm_tmp	= int(WJets_norm_val)

	VV_dataset	= VV_pdf.generate(RooArgSet(rrv_x),VV_norm_tmp)
	STop_dataset	= STop_pdf.generate(RooArgSet(rrv_x),STop_norm_tmp)
	TTbar_dataset	= TTbar_pdf.generate(RooArgSet(rrv_x),TTbar_norm_tmp)
	WJets_dataset	= WJets_pdf.generate(RooArgSet(rrv_x),WJets_norm_tmp)

	allbkg_dataset = RooDataSet('allbkg_dataset_%s'%i,'allbkg_dataset_%s'%i,RooArgSet(rrv_x))
	allbkg_dataset.append(VV_dataset)
	allbkg_dataset.append(STop_dataset)
	allbkg_dataset.append(TTbar_dataset)
	allbkg_dataset.append(WJets_dataset)
	WJetsonly	= RooDataSet('Wjets_data_%s'%i,'Wjets_data_%s'%i,RooArgSet(rrv_x))
	VVonly		= RooDataSet('VV_data_%s'%i,'VV_data_%s'%i,RooArgSet(rrv_x))
	WJetsonly.append(WJets_dataset)
	VVonly.append(VV_dataset)
	if i%10==0:
		allbkg_dataset.Print()
	allbkg_dataset.Write()
	WJetsonly.Write()
	VVonly.Write()
fileOut.Close()
