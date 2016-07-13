from ROOT import *

import os
from optparse import OptionParser

ROOT.gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load("PDFs/Util_cxx.so")
ROOT.gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

parser	= OptionParser()
parser.add_option('--cat', dest='cat', default='WW', help='category, WW or WZ, defines signal region')
parser.add_option('--ch', dest='ch', default='el', help='channel, el or mu')
parser.add_option('--nuis', dest='nuis', action='store_true', default=False)
parser.add_option('--wjets', dest='onlywjets', action='store_true', default=False)
parser.add_option('--ttbar', dest='onlyttbar', action='store_true', default=False)
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

if options.nuis:
	fileOut		= TFile.Open("mj_closure/datasetsnuis_%s_HP%s.root"%(ch,cat[1]),"recreate")
elif options.onlywjets:
	fileOut		= TFile.Open("mj_closure_wjetsonly/datasetwjet_%s_HP%s.root"%(ch,cat[1]),"recreate")
elif options.onlyttbar:
	fileOut		= TFile.Open("mj_closure_wjetsonly/datasetttbar_%s_HP%s.root"%(ch,cat[1]),"recreate")
else:
	fileOut		= TFile.Open("mj_closure/datasets_%s_HP%s.root"%(ch,cat[1]),"recreate")
if options.onlywjets:
	VV_norm_tmp	= int(VV_norm_val)
	STop_norm_tmp	= int(STop_norm_val)
	TTbar_norm_tmp	= int(TTbar_norm_val)
	for i in range(100):
		WJets_norm_tmp  = random.Poisson(WJets_norm_val * 2)

		VV_dataset	= VV_pdf.generate(RooArgSet(rrv_x),VV_norm_tmp)
		STop_dataset	= STop_pdf.generate(RooArgSet(rrv_x),STop_norm_tmp)
		TTbar_dataset	= TTbar_pdf.generate(RooArgSet(rrv_x),TTbar_norm_tmp)
		WJets_dataset	= WJets_pdf.generate(RooArgSet(rrv_x),WJets_norm_tmp)

		allbkg_dataset = RooDataSet('allbkg_dataset_2wjets%s'%i,'allbkg_dataset_2wjets%s'%i,RooArgSet(rrv_x))
		allbkg_dataset.append(VV_dataset)
		allbkg_dataset.append(STop_dataset)
		allbkg_dataset.append(TTbar_dataset)
		allbkg_dataset.append(WJets_dataset)
		allbkg_dataset.Print()
		allbkg_dataset.Write()

		WJets_norm_tmp  = random.Poisson(WJets_norm.getVal() * 0.5)

		VV_dataset	= VV_pdf.generate(RooArgSet(rrv_x),VV_norm_tmp)
		STop_dataset	= STop_pdf.generate(RooArgSet(rrv_x),STop_norm_tmp)
		TTbar_dataset	= TTbar_pdf.generate(RooArgSet(rrv_x),TTbar_norm_tmp)
		WJets_dataset	= WJets_pdf.generate(RooArgSet(rrv_x),WJets_norm_tmp)

		allbkg_dataset = RooDataSet('allbkg_dataset_05wjets%s'%i,'allbkg_dataset_05wjets%s'%i,RooArgSet(rrv_x))
		allbkg_dataset.append(VV_dataset)
		allbkg_dataset.append(STop_dataset)
		allbkg_dataset.append(TTbar_dataset)
		allbkg_dataset.append(WJets_dataset)
		allbkg_dataset.Print()
		allbkg_dataset.Write()
elif options.onlyttbar:
	VV_norm_tmp	= int(VV_norm_val)
	STop_norm_tmp	= int(STop_norm_val)
	WJets_norm_tmp	= int(WJets_norm_val)
	for i in range(100):
		TTbar_norm_tmp	= random.Poisson(TTbar_norm_val + TTbar_norm_err)

		VV_dataset	= VV_pdf.generate(RooArgSet(rrv_x),VV_norm_tmp)
		STop_dataset	= STop_pdf.generate(RooArgSet(rrv_x),STop_norm_tmp)
		TTbar_dataset	= TTbar_pdf.generate(RooArgSet(rrv_x),TTbar_norm_tmp)
		WJets_dataset	= WJets_pdf.generate(RooArgSet(rrv_x),WJets_norm_tmp)

		allbkg_dataset = RooDataSet('allbkg_dataset_plusttbar%s'%i,'allbkg_dataset_plusttbar%s'%i,RooArgSet(rrv_x))
		allbkg_dataset.append(VV_dataset)
		allbkg_dataset.append(STop_dataset)
		allbkg_dataset.append(TTbar_dataset)
		allbkg_dataset.append(WJets_dataset)
		allbkg_dataset.Print()
		allbkg_dataset.Write()

		TTbar_norm_tmp	= random.Poisson(TTbar_norm_val - TTbar_norm_err)

		VV_dataset	= VV_pdf.generate(RooArgSet(rrv_x),VV_norm_tmp)
		STop_dataset	= STop_pdf.generate(RooArgSet(rrv_x),STop_norm_tmp)
		TTbar_dataset	= TTbar_pdf.generate(RooArgSet(rrv_x),TTbar_norm_tmp)
		WJets_dataset	= WJets_pdf.generate(RooArgSet(rrv_x),WJets_norm_tmp)

		allbkg_dataset = RooDataSet('allbkg_dataset_minusttbar%s'%i,'allbkg_dataset_minusttbar%s'%i,RooArgSet(rrv_x))
		allbkg_dataset.append(VV_dataset)
		allbkg_dataset.append(STop_dataset)
		allbkg_dataset.append(TTbar_dataset)
		allbkg_dataset.append(WJets_dataset)
		allbkg_dataset.Print()
		allbkg_dataset.Write()
else:
	for i in range(100):

		VV_norm_tmp	= random.Poisson(VV_norm_val)
		STop_norm_tmp	= random.Poisson(STop_norm_val)
		TTbar_norm_tmp	= random.Poisson(TTbar_norm_val)
		WJets_norm_tmp	= random.Poisson(WJets_norm_val)

		if options.nuis:
			VV_pars 	= w_fits.genobj('fitresult_model_VV_xww_%s_mj_rdataset4fit_VV_xww_%s_mj'%(ch,ch)).randomizePars()
			for j in range(VV_pars.getSize()):
				if 'number' not in VV_pars[j].GetName():
					w.var(VV_pars[j].GetName()).setVal(VV_pars[j].getVal())
			STop_pars	= w_fits.genobj('fitresult_model_STop_xww_%s_mj_rdataset4fit_STop_xww_%s_mj'%(ch,ch)).randomizePars()
			for j in range(STop_pars.getSize()):
				if 'number' not in STop_pars[j].GetName():
					w.var(STop_pars[j].GetName()).setVal(STop_pars[j].getVal())
			TTbar_pars	= w_fits.genobj('fitresult_model_TTbar_xww_%s_mj_rdataset4fit_TTbar_xww_%s_mj'%(ch,ch)).randomizePars()
			for j in range(TTbar_pars.getSize()):
				if 'number' not in TTbar_pars[j].GetName():
					w.var(TTbar_pars[j].GetName()).setVal(TTbar_pars[j].getVal())
			WJets_pars	= w_fits.genobj('fitresult_model_WJets0_xww_%s_mj_rdataset4fit_WJets0_xww_%s_mj'%(ch,ch)).randomizePars()
			for j in range(WJets_pars.getSize()):
				if 'number' not in WJets_pars[j].GetName():
					w.var(WJets_pars[j].GetName()).setVal(WJets_pars[j].getVal())
		VV_dataset	= VV_pdf.generate(RooArgSet(rrv_x),VV_norm_tmp)
		STop_dataset	= STop_pdf.generate(RooArgSet(rrv_x),STop_norm_tmp)
		TTbar_dataset	= TTbar_pdf.generate(RooArgSet(rrv_x),TTbar_norm_tmp)
		WJets_dataset	= WJets_pdf.generate(RooArgSet(rrv_x),WJets_norm_tmp)

		allbkg_dataset = RooDataSet('allbkg_dataset_%s'%i,'allbkg_dataset_%s'%i,RooArgSet(rrv_x))
		allbkg_dataset.append(VV_dataset)
		allbkg_dataset.append(STop_dataset)
		allbkg_dataset.append(TTbar_dataset)
		allbkg_dataset.append(WJets_dataset)
		allbkg_dataset.Print()
		allbkg_dataset.Write()
fileOut.Close()
