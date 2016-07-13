from ROOT import *
from optparse import OptionParser

gSystem.Load('PDFs/PdfDiagonalizer_cc.so')
gSystem.Load('PDFs/Util_cxx.so')
gSystem.Load('PDFs/hyperg_2F1_c.so')
gSystem.Load('PDFs/HWWLVJRooPdfs_cxx.so')

from ROOT import Calc_error_extendPdf, RooErfExpDecoPdf, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

parser = OptionParser()
parser.add_option('--cat', dest='cat', default='WW')
parser.add_option('--ch', dest='ch', default='el')
parser.add_option('-b', action='store_true', default=False)
(options,args) = parser.parse_args()



#ch 	= options.ch
#cat	= options.cat
#channel	= cat+'_'+ch
gROOT.SetBatch(options.b)

for ch in ['el','mu']:
	for cat in ['WW']:
		channel = cat+'_'+ch
		fileIn	= TFile.Open('cards_%s_HP%s/mj_closure.root'%(ch,cat[1]))
		w	= fileIn.Get('w4fit')
		fileIn.Close()

		rrv_x		= w.var('rrv_mass_j')
		VV_pdf		= w.pdf('model_VV_xww_%s_mj'%ch)
		STop_pdf	= w.pdf('model_STop_xww_%s_mj'%ch)
		TTbar_pdf	= w.pdf('model_TTbar_xww_%s_mj'%ch)
		WJets_pdf	= w.pdf('model_WJets0_xww_%s_mj'%ch)
		VV_norm		= w.var('rrv_number_VV_xww_%s_mj'%ch)
		STop_norm	= w.var('rrv_number_STop_xww_%s_mj'%ch)
		TTbar_norm	= w.var('rrv_number_TTbar_xww_%s_mj'%ch)
		WJets_norm	= w.var('rrv_number_WJets0_xww_%s_mj'%ch)
		c		= w.var('rrv_c_ErfExp_WJets0_xww_%s'%ch)
		width		= w.var('rrv_width_ErfExp_WJets0_xww_%s'%ch)
		offset		= w.var('rrv_offset_ErfExp_WJets0_xww_%s'%ch)


		VV_norm_val	= VV_norm.getVal()
		STop_norm_val	= STop_norm.getVal()
		TTbar_norm_val	= TTbar_norm.getVal()
		WJets_norm_val	= WJets_norm.getVal()
		c_val		= c.getVal()

		WJets_norm.setConstant(0)
		VV_norm.setConstant(0)
		TTbar_norm.setConstant(1)
		c.setConstant(0)
		width.setConstant(1)
		offset.setConstant(1)


		model_pdf4fit		= w.pdf('model_data_xww_%s_mj_4fit'%ch)

		norms		= []
		generated	= []
		WJetsnorms	= []
		ttbar		= []
		errors		= []
		errorsVV	= []
		fresults	= []
		#norm_nom 	= {'elWW' : 124.46, 'elWZ' : 102.76, 'muWW' : 181.28, 'muWZ' : 153.87}
		norm_nom 	= {'el' : 532.14, 'mu' : 772}
		ttbar_nom	= {'elWW' :68.87, 'elWZ' : 58.18, 'muWW' : 90.66, 'muWZ' : 76.63} 

		rrv_x.setRange('sb_lo',40,65)
		rrv_x.setRange('sb_hi',105,150)
		rrv_x.setRange('WW',65,85)
		rrv_x.setRange('WZ',85,105)
		rrv_x.setRange('all',40,150)


		fileInData	= TFile.Open('VV_xs/datasets_%s_HP%s.root'%(ch,cat[1]))
		for i in range(500):
			data 		= fileInData.Get('allbkg_dataset_%s'%i)
			wjetsdata	= fileInData.Get('Wjets_data_%s'%i)
			VVdata		= fileInData.Get('VV_data_%s'%i)
			WJets_norm.setVal(WJets_norm_val)
			TTbar_norm.setVal(TTbar_norm_val)
			c.setVal(c_val)
			VV_norm.setVal(VV_norm_val)

			#fresult		= model_pdf4fit.fitTo(data, RooFit.Range('all'), RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.PrintLevel(-1), RooFit.Warnings(kFALSE), RooFit.Verbose(kFALSE))
			fresult		= model_pdf4fit.fitTo(data, RooFit.Range('all'), RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Minimizer('Minuit2'), RooFit.PrintLevel(-1), RooFit.Warnings(kFALSE), RooFit.Verbose(kFALSE))
			fresults.append(fresult)

			WJetsnorms.append(WJets_norm.getVal())
			norms.append(VV_norm.getVal())

			WJets_norm.setError(Calc_error_extendPdf(data, WJets_pdf, fresult, 'all'))
			VV_norm.setError(Calc_error_extendPdf(data, VV_pdf, fresult, 'all'))
			errors.append(WJets_norm.getError())
			errorsVV.append(VV_norm.getError())
			generated.append(VVdata.sumEntries())





		norm_hist	= TH1F('norm_hist','norm_hist',30,0,150)
		for i in range(len(norms)):
			norm_hist.Fill(norms[i])
		#c1	= TCanvas('c1','c1',1)
		#c1.cd()
		#norm_hist.Draw()
		#c1.Draw()
		
		nom	= norm_nom[ch]
		pull1	= TH1F('pull1','',20,-3,3)
		pull2	= TH1F('pull2','',20,-3,3)
		hist2	= TH1F('hist1','',20,0,5)
		normshist	= TH1F('normshist','normshist',100,0,100)
		for i in range(len(WJetsnorms)):
			pull2.Fill((norms[i]-VV_norm_val)/errorsVV[i])
			pull1.Fill((norms[i]-generated[i])/errors[i])
			hist2.Fill((norms[i])/errorsVV[i])
			normshist.SetBinContent(i,norms[i])


		

		gStyle.SetTitle('')

		c2	= TCanvas('c2','MC-generated',1)
		c2.cd()
		pull1.GetXaxis().SetTitle('(N_{fit}^{WV}-N_{gen}^{WV})/#sigma_{fit}')
		pull1.GetYaxis().SetTitle('counts')
		pull1.Draw()
		c2.Draw()
		c2.Update()
		c3	= TCanvas('c3','MX-xs',1)
		c3.cd()
		pull2.GetXaxis().SetTitle('(N_{fit}^{WV}-N_{0}^{WV})/#sigma_{fit}')
		pull2.GetYaxis().SetTitle('counts')
		pull2.Draw()
		c3.Draw()
		c3.Update()
		c4	= TCanvas('c4','WV_norm/WV_error',1)
		c4.cd()
		hist2.GetXaxis().SetTitle('N_{fit}^{WV}/#sigma_{fit}')
		hist2.GetYaxis().SetTitle('counts')
		hist2.Draw()
		c4.Draw()
		c4.Update()
		c5	= TCanvas('norms','norms',1)
		c5.cd()
		normshist.Draw()
		c5.Draw()
		c5.Update()
		
		
		c4.SaveAs('VV_xs/MJ-VVoversigma_%s.pdf'%ch)
		c3.SaveAs('VV_xs/pullplot_MJ-nom_%s.pdf'%ch)
		c4.SaveAs('VV_xs/MJ-VVoversigma_%s.png'%ch)
		c3.SaveAs('VV_xs/pullplot_MJ-nom_%s.png'%ch)

		if not options.b:
			print VV_norm_val
			raw_input(channel)
			#c1.Close()
		fileInData.Close()






