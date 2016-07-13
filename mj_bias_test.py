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
parser.add_option('--nuis', action='store_true', default=False)
parser.add_option('--scale', dest='scale', default='0')
parser.add_option('-b', action='store_true', default=False)
(options,args) = parser.parse_args()



#ch 	= options.ch
#cat	= options.cat
#channel	= cat+'_'+ch
gROOT.SetBatch(options.b)

for ch in ['el','mu']:
	for cat in ['WW','WZ']:
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
		TTbar_norm.setConstant(0)
		c.setConstant(0)
		width.setConstant(1)
		offset.setConstant(1)



		model_pdf4fit		= w.pdf('model_data_xww_%s_mj_4fit'%ch)
		model_pdf4fit_alt	= w2.pdf('model_data_xww_%s_mj_4fit'%ch)

		norms		= []
		ttbar		= []
		errors		= []
		fresults	= []
		norm_nom 	= {'elWW' : 124.46, 'elWZ' : 102.76, 'muWW' : 181.28, 'muWZ' : 153.87}
		ttbar_nom	= {'elWW' :68.87, 'elWZ' : 58.18, 'muWW' : 90.66, 'muWZ' : 76.63} 

		rrv_x.setRange('sb_lo',40,65)
		rrv_x.setRange('sb_hi',105,150)
		rrv_x.setRange('WW',65,85)
		rrv_x.setRange('WZ',85,105)
		rrv_x_alt.setRange('sb_lo',40,65)
		rrv_x_alt.setRange('sb_hi',105,150)
		rrv_x_alt.setRange('WW',65,85)
		rrv_x_alt.setRange('WZ',85,105)

		if options.nuis:
			fileInData	= TFile.Open('mj_closure/datasetsnuis_%s_HP%s.root'%(ch,cat[1]))
		elif options.scale=='2' or options.scale=='05':
			fileInData	= TFile.Open('mj_closure_wjetsonly/datasetwjet_%s_HP%s.root'%(ch,cat[1]))
		elif options.scale=='plus' or options.scale=='minus':
			fileInData	= TFile.Open('mj_closure_wjetsonly/datasetttbar_%s_HP%s.root'%(ch,cat[1]))
		else:
			fileInData	= TFile.Open('mj_closure/datasets_%s_HP%s.root'%(ch,cat[1]))
		for i in range(100):
			if options.scale!='0':
				if options.scale == '2':
					data	= fileInData.Get('allbkg_dataset_2wjets%s'%i)
				elif options.scale == '05':
					data	= fileInData.Get('allbkg_dataset_05wjets%s'%i)
				elif options.scale == 'plus':
					data	= fileInData.Get('allbkg_dataset_plusttbar%s'%i)
				elif options.scale == 'minus':
					data	= fileInData.Get('allbkg_dataset_minusttbar%s'%i)
			else:
				data 	= fileInData.Get('allbkg_dataset_%s'%i)
			#WJets_norm.setVal(WJets_norm_val)
			#TTbar_norm.setVal(TTbar_norm_val)
			#c.setVal(c_val)

			fresult		= model_pdf4fit.fitTo(data, RooFit.Range('sb_lo,sb_hi'), RooFit.Save(kTRUE), RooFit.Extended(kTRUE), RooFit.Minimizer('Minuit2'), RooFit.PrintLevel(-1), RooFit.Warnings(kFALSE), RooFit.Verbose(kFALSE))

			fullIntWJets			= WJets_pdf.createIntegral(RooArgSet(rrv_x), RooArgSet(rrv_x))
			signalIntWJets			= WJets_pdf.createIntegral(RooArgSet(rrv_x), RooArgSet(rrv_x), (cat))
			signalInt_valWJets		= signalIntWJets.getVal()/fullIntWJets.getVal()
			norms.append(WJets_norm.getVal()*signalInt_valWJets)

			WJets_norm.setError(Calc_error_extendPdf(data, WJets_pdf, fresult, cat))

			errors.append(WJets_norm.getError())
			fresults.append(fresult)
			fullIntTTbar			= TTbar_pdf.createIntegral(RooArgSet(rrv_x), RooArgSet(rrv_x))
			signalIntTTbar			= TTbar_pdf.createIntegral(RooArgSet(rrv_x), RooArgSet(rrv_x), (cat))
			signalInt_valTTbar		= signalIntTTbar.getVal()/fullIntTTbar.getVal()
			ttbar.append(TTbar_norm.getVal()*signalInt_valTTbar)

		if options.scale=='2':
			nominal = norm_nom[ch+cat]*2
		elif options.scale=='05':
			nominal = norm_nom[ch+cat]*0.5
		else :
			nominal = norm_nom[ch+cat]
		ttbar_nominal	= ttbar_nom[ch+cat]


		hist		= TH1F('hist','hist',20,nominal-40,nominal+40)
		hist_pull	= TH1F('hist_pull','hist_pull',20,-3.5,3.5)
		hist_ttbar	= TH1F('ttbar','ttbar',20,ttbar_nominal-20,ttbar_nominal+20)
		for i in range(len(errors)):
			hist_pull.Fill((norms[i]-nominal)/errors[i])
			hist.Fill(norms[i])
		for i in range(len(ttbar)):
			hist_ttbar.Fill(ttbar[i])

		hist.GetXaxis().SetTitle('N_{WJets}')
		hist.GetYaxis().SetTitle('count')
		hist.SetTitle('')
		line = TLine(nominal,0,nominal,hist.GetMaximum())
		line.SetLineStyle(kDashed)
		hist_pull.GetXaxis().SetTitle("(N_{fit}^{WJets}-N_{0}^{WJets}) / #sigma_{fit}")
		hist_pull.SetTitle('')

		c1	= TCanvas('c1','c1',1)
		hist_pull.Draw()
		c1.Draw()
		c2	= TCanvas('c2','c2',1)
		c2.cd()
		hist.Draw()
		line.Draw("SAME")
		c2.Draw()




		if options.nuis:
			c1.SaveAs('FINALPLOTS/mj_closure_%s_%s_nuis_pull.pdf'%(cat,ch))
			c2.SaveAs('FINALPLOTS/mj_closure_%s_%s_nuis.pdf'%(cat,ch))
		elif options.scale!='0':
			c2.SaveAs('FINALPLOTS/mj_closure_%s_%s_%swjets.pdf'%(cat,ch,options.scale))
		else:
			c2.SaveAs("FINALPLOTS/mj_closure_%s_%s.pdf"%(cat,ch))

		print hist_ttbar.GetMean(),' +- ',hist_ttbar.GetRMS()
		N_ttbar	= RooRealVar('N_ttbar','N_ttbar',hist_ttbar.GetMean())
		N_ttbar.setError(hist_ttbar.GetRMS())
		if options.nuis:
			name='nuis'
		elif options.scale != '0':
			name=options.scale
		else:
			name='random'
		fileOut		= TFile.Open('fitbyhandroots/'+name+channel+'.root','recreate')
		N_ttbar.Write()
		fileOut.Close()
		if not options.b:
			raw_input(channel)






