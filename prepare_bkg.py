#python g1_exo_doFit_class.py -b -c mu > test.log
#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser
import CMS_lumi, tdrstyle
from array import array

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooNLLVar, RooAddition, RooProduct, RooConstraintSum, RooCustomizer, RooMinuit, RooArgSet, RooAbsData, RooAbsPdf, RooAbsReal, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite


############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('-s','--simple', action='store', dest='simple', default=False, help='pre-limit in simple mode')
parser.add_option('-b', action='store_true', dest='noX', default=True, help='no X11 windows')
parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--category', '--cat', action="store",type="string",dest="category",default="")
parser.add_option('--jetalgo', action="store",type="string",dest="jetalgo",default="Mjpruned")
parser.add_option('--hi', action='store', dest='mlvj_hi', type='float', default=3500, help='dont change atm!')
parser.add_option('--lo', action='store', dest='mlvj_lo', type='float', default=900, help='dont change atm!')
parser.add_option('--prefit', action='store_true', dest='prefitplot', default=False)
parser.add_option('--readtrees', action='store_true', dest='read_trees', default=False)
parser.add_option('--simfit', action='store_true', dest='simfit', default=False, help='fit simultaneously in mpruned and mww-sb')
parser.add_option('--noplots', action='store_true', dest='noplots', default=False, help='dont make any plots')


(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/hyperg_2F1_c.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf, RooErfExpDecoPdf

fitresultsmj	= []
fitresultsmlvj	= []
fitresultsfinal	= []

class doFit_wj_and_wlvj:

    def __init__(self, in_channel, jetalgo, in_mj_min=30, in_mj_max=140, in_mlvj_min=400., in_mlvj_max=1400., fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", input_workspace=None):

        tdrstyle.setTDRStyle()
        TGaxis.SetMaxDigits(3)
	
        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

	###save all fitresults

        ### set the channel type --> electron or muon
        self.channel=in_channel;
        self.leg = TLegend(); 
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;
	self.jetalgo = jetalgo
                
        print "########################################################################################"
        print "######## define class: binning, variables, cuts, files and nuissance parameters ########"
        print "########################################################################################"

        ### Set the mj binning for plots
        self.BinWidth_mj=5.;

        ### Set the binning for mlvj plots as a function of the model
        self.BinWidth_mlvj=100.;
            
        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.narrow_factor=1.

        ## correct the binning of mj 
        self.BinWidth_mj=self.BinWidth_mj/self.narrow_factor;
        nbins_mj=int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max=in_mj_min+nbins_mj*self.BinWidth_mj;
                   
        ## correct the binning of mlvj 
        self.BinWidth_mlvj=self.BinWidth_mlvj/self.narrow_factor;
        nbins_mlvj=int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        in_mlvj_max=in_mlvj_min+nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
	varname = "Softdrop jet mass"
	if self.jetalgo == "Mjpruned": varname = "M_{pruned}"
        rrv_mass_j = RooRealVar("rrv_mass_j",varname,(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV");
        rrv_mass_j.setBins(nbins_mj);

        ## define invariant mass WW variable
        rrv_mass_lvj= RooRealVar("rrv_mass_lvj","M_{WV}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV");
        rrv_mass_lvj.setBins(nbins_mlvj);

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_");
        else:
            self.workspace4fit_ = input_workspace;
        getattr(self.workspace4fit_,"import")(rrv_mass_j);
        getattr(self.workspace4fit_,"import")(rrv_mass_lvj);

        #prepare workspace for unbin-Limit -> just fo the stuff on which running the limit 
        self.workspace4limit_ = RooWorkspace("workspace4limit_","workspace4limit_");

        #define sidebands
        self.mj_sideband_lo_min = in_mj_min;
        self.mj_sideband_lo_max = 65
        self.mj_sideband_hi_min = 105
        self.mj_sideband_hi_max = in_mj_max;
        self.mj_signal_min = 65
        self.mj_signal_max = 105
        if (options.category).find('W') != -1:
         self.mj_signal_min = 65
         self.mj_signal_max = 85
	elif (options.category).find('Z') != -1:
         self.mj_signal_min = 85
         self.mj_signal_max = 105



        ## zone definition in the jet mass 
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
	rrv_mass_j.setRange("total_signal_region",self.mj_sideband_lo_max,self.mj_sideband_hi_min)
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

	## one zone for MWW
	rrv_mass_lvj.setRange("total_region",in_mlvj_min,in_mlvj_max)
	rrv_mass_lvj.setRange('over3500',3500,5000)
	rrv_mass_lvj.setRange('signal_region',900,3500)

        #prepare the data and mc files --> set the working directory and the files name
        self.file_Directory="AnaSigTree_25ns/"+self.channel+"/";
                 
        #prepare background data and signal samples            

        self.file_data		= ("treeEDBR_data_xww.root");#keep blind!!!!
        self.file_WJets0_mc 	= ("treeEDBR_WJets_xww.root");
        self.file_VV_mc     	= ("treeEDBR_VV_xww.root");# WW+WZ
        self.file_TTbar_mc   	= ("treeEDBR_TTBARpowheg_xww.root");
        self.file_STop_mc    	= ("treeEDBR_SingleTop_xww.root");
        #self.file_data 		= ("tree_data_%s.root"%self.channel);#keep blind!!!!
        #self.file_WJets0_mc  	= ("tree_WJets_%s.root"%self.channel);
        #self.file_VV_mc      	= ("tree_VV_%s.root"%self.channel);# WW+WZ
        #self.file_TTbar_mc   	= ("tree_TTBARpowheg_%s.root"%self.channel);
        #self.file_STop_mc    	= ("tree_SingleTop_%s.root"%self.channel);

        
        #self.mean_shift = -0.8
        #self.sigma_scale=1.086
        self.mean_shift = -1.294
        self.sigma_scale=0.958

	self.wtagger_label	= options.category
        
	if options.mlvj_hi!=3500:
		raw_input('cuts different from M_WV<3500 not supported atm!')

	self.plotsDir = 'plots_%s_%s_%s_%s' %(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi)
	if not os.path.isdir("cards_%s_%s_%s_%s"%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi)):
        	os.system("mkdir cards_%s_%s_%s_%s"%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi));
	self.rlt_DIR="cards_%s_%s_%s_%s/"%(self.channel,self.wtagger_label,options.mlvj_lo,options.mlvj_hi)

        #result files: The event number, parameters and error write into a txt file. The dataset and pdfs write into a root file


        ## extra text file
        self.file_rlt_txt = self.rlt_DIR+"other_wwlvj_%s_%s.txt"%(self.channel,self.wtagger_label)
        ## workspace for limit
        self.file_rlt_root = self.rlt_DIR+"wwlvj_%s_%s_workspace.root"%(self.channel,self.wtagger_label)
        ## datacard for the ubninned limit
        self.file_datacard_unbin = self.rlt_DIR+"wwlvj_%s_%s_unbin.txt"%(self.channel,self.wtagger_label)
        ## workspace for the binned limit
        self.file_datacard_counting = self.rlt_DIR+"wwlvj_%s_%s_counting.txt"%(self.channel,self.wtagger_label)
        
        self.file_out=open(self.file_rlt_txt,"w");
        self.file_out.write("Welcome:\n");
        self.file_out.close()
        self.file_out=open(self.file_rlt_txt,"a+");

        ## color palette 
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : kGreen+1,
            'VV' : kRed,
            'STop' : kBlue,
            'TTbar' : kOrange,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        # parameters of data-driven method to get the WJets background event number.
        self.number_WJets_insideband=-1;
        self.datadriven_alpha_WJets_unbin=-1;
        self.datadriven_alpha_WJets_counting=-1;

        #### Set systematic on the Wjets shape   and TTbar due to PS, fitting function etc..
        self.shape_para_error_WJets0 = 1.4;
        self.shape_para_error_alpha  = 1.4;
        self.shape_para_error_TTbar = 2.0;
        self.shape_para_error_VV    = 1.;
        self.shape_para_error_STop  = 1.;
                                                                  
        # shape parameter uncertainty
        self.FloatingParams=RooArgList("floatpara_list");

      
    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1., isalpha=False, ismj=False):
        print "############### draw the legend ########################"
        if left==-1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetLineColor(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.37+x_offset_low, 0.50+y_offset_low, 0.72+x_offset_high, 0.82+y_offset_high, "", "NDC");            
            theLeg.SetName("theLegend");
            if ismj: theLeg = TLegend(0.3715365+x_offset_low,0.505+y_offset_low,0.8526448+x_offset_high,0.845+y_offset_high, "", "NDC"); 
            if TwoCoulum :
                theLeg.SetNColumns(2);
            if isalpha: theLeg = TLegend(0.3944724+x_offset_low,0.4370629+y_offset_low,0.7650754+x_offset_high,0.8374126+y_offset_high, "", "NDC");  
            
        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.05);
        theLeg.SetTextFont(42);

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        if   self.channel == 'el': legHeader="e#nu";
        elif self.channel == 'mu': legHeader="#mu#nu";
        
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
	  #if objName.find("TPave") != -1: continue
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
  	    elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);  objName_before=objName;
            elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "Data, W#rightarrow"+legHeader,"PE");  objName_before=objName;                 
            else: objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        
                   
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="WJets" : objName_before=objName;                 
            else:  objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
	    if objName.find("TPave") != -1: continue
            if objName == "errorband" : objName = "Uncertainty";
            print objName;
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data()=="STop" : theLeg.AddEntry(theObj, "Single Top","F");
                    #elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ","F");
                    elif TString(objName).Data()=="data" :  objName_before=objName; entryCnt = entryCnt+1; continue ;
                    elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F"); entryCnt = entryCnt+1; continue;
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");

                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
           theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;

    def get_canvas(self,cname,isalpha=False):

       #tdrstyle.setTDRStyle()
       CMS_lumi.lumi_13TeV = "2.3 fb^{-1}"
       CMS_lumi.writeExtraText = 1
       if isalpha:
	       	CMS_lumi.extraText = "simulation\n preliminary"
       else:
       		CMS_lumi.extraText = "preliminary"

       iPos = 11
       if( iPos==0 ): CMS_lumi.relPosX = 0.15

       H_ref = 600; 
       W_ref = 800; 
       W = W_ref
       H  = H_ref

       T = 0.08*H_ref
       B = 0.12*H_ref 
       L = 0.12*W_ref
       R = 0.06*W_ref

       canvas = ROOT.TCanvas(cname,cname,W,H)
       canvas.SetFillColor(0)
       canvas.SetBorderMode(0)
       canvas.SetFrameFillStyle(0)
       canvas.SetFrameBorderMode(0)
       canvas.SetLeftMargin( L/W )
       canvas.SetRightMargin( R/W )
       canvas.SetTopMargin( T/H )
       canvas.SetBottomMargin( B/H+0.03 )
       canvas.SetTickx()
       canvas.SetTicky()
       if isalpha:
        canvas.SetTicky(0)
       
       return canvas

    #### just drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0, frompull=0, isalpha=0, fix_axis=0, force_plots=0):
      
	if options.noplots and not force_plots:
	  return 0

        print "############### draw the canvas without pull ########################"
        cMassFit = self.get_canvas("cMassFit",isalpha)#TCanvas("cMassFit","cMassFit", 600,600);

	if fix_axis == 0:
		if frompull and logy :
		    in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/200)
		elif not frompull and logy :
		    in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())

        if is_range:
            h2=TH2D("h2","",100,400,1400,4,0.00001,4);
            h2.Draw();
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045);
        in_obj.GetXaxis().SetTitleOffset(1.15);
        in_obj.GetXaxis().SetLabelSize(0.04);

        in_obj.GetYaxis().SetTitleSize(0.04);
        in_obj.GetYaxis().SetTitleOffset(1.40);
        in_obj.GetYaxis().SetLabelSize(0.04);

        self.leg.SetTextSize(0.031); 
        if isalpha: self.leg.SetTextSize(0.038)

	cMassFit.Update()
	cMassFit.cd()
	if isalpha:
		CMS_lumi.CMS_lumi(cMassFit, 0, 11, 0.075)
	else:
		CMS_lumi.CMS_lumi(cMassFit, 4, 11)
	cMassFit.cd()
	cMassFit.Update()
	cMassFit.RedrawAxis()
	frame = cMassFit.GetFrame()
	frame.Draw()   
	cMassFit.cd()
	cMassFit.Update()
	        
        Directory=TString(in_directory); #+"_%02d_%02d/"%(options.cprime,options.BRnew));
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
        else:
            rlt_file.ReplaceAll(".root","");
            rlt_file = rlt_file.Append(".png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".pdf",".root");
        #cMassFit.SaveAs(rlt_file.Data());

        if logy:
	    if not isalpha and fix_axis == 0:
                in_obj.GetYaxis().SetRangeUser(1e-3,in_obj.GetMaximum()*200);
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            #cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());
	
	cMassFit.Delete()

    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, rrv_x, datahist, mplot, mplot_pull,ndof,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0,ismj=0, fix_axis=0, force_plots=0):# mplot + pull

	if options.noplots and not force_plots:
	  return 0

        print "############### draw the canvas with pull ########################" 
	#@#mplot.GetXaxis().SetTitle(rrv_x.GetTitle() + " (GeV)")
	mplot.GetXaxis().SetTitle("")
        mplot.GetYaxis().SetTitleSize(0.07)
        mplot.GetYaxis().SetTitleOffset(0.9)
        mplot.GetYaxis().SetLabelSize(0.06)
	mplot.GetXaxis().SetLabelSize(0);
	
        cMassFit = self.get_canvas("cMassFit")
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
         pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
         pad3=TPad("pad3","pad3",0.8,0.,1,1);
         pad1.Draw();
         pad2.Draw();
         pad3.Draw();
        else:
         pad1=TPad("pad1","pad1",0.,0. ,1,0.30); #pad1 - pull
         pad2=TPad("pad2","pad2",0.,0.3,1.,1. ); #pad0

   	 pad2.SetRightMargin(0.1);
   	 pad2.SetTopMargin(0.1);
   	 pad2.SetBottomMargin(0.0001);
   	 pad1.SetRightMargin(0.1)
   	 pad1.SetTopMargin(0)
   	 pad1.SetBottomMargin(0.4)   
         pad1.Draw();
         pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
	
        mplot.Draw();

        pad1.cd();
        mplot_pull.Draw("AP");

        medianLine = TLine(mplot.GetXaxis().GetXmin(),0.,mplot.GetXaxis().GetXmax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
	medianLine.Draw()
	mplot_pull.Draw("Psame");
	
        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

	cMassFit.Update()
	pad2.cd()
	CMS_lumi.CMS_lumi(pad2, 4, 11)	
	pad2.cd()
	pad2.Update()
	pad2.RedrawAxis()
	frame = pad2.GetFrame()
	frame.Draw()   
	cMassFit.cd()
	cMassFit.Update()
			
        ## create the directory where store the plots
        Directory = TString(in_directory);
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
        else:
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","");
            rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());
        
        rlt_file.ReplaceAll(".pdf",".root");
        #cMassFit.SaveAs(rlt_file.Data());

        string_file_name = TString(in_file_name);
        if string_file_name.EndsWith(".root"):
            string_file_name.ReplaceAll(".root","_"+in_model_name);
        else:
            string_file_name.ReplaceAll(".root","");
            string_file_name.Append("_"+in_model_name);

        if logy:
	    if fix_axis == 0:
		    mplot.GetYaxis().SetRangeUser(0.002,mplot.GetMaximum()*200);
	    pad2.SetLogy() ;
	    pad2.Update();
	    cMassFit.Update();
	    rlt_file.ReplaceAll(".root","_log.root");
	    #cMassFit.SaveAs(rlt_file.Data());
	    rlt_file.ReplaceAll(".root",".pdf");
	    cMassFit.SaveAs(rlt_file.Data());
	    rlt_file.ReplaceAll(".pdf",".png");
	    cMassFit.SaveAs(rlt_file.Data());

        #self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1,0,fix_axis);

    def calculate_chi2(self,hist,rrv_x,mplot_orig,ndof,ismj):

        pulls = array('d',[])
        print "############### calculate chi2 (new) ########################"
        hpull = mplot_orig.pullHist();
	bins = 0
	bins_ = 0
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
	  hist.get(bins_)
	  hist.weightError(RooAbsData.SumW2)
	  print x,y,bins_,hist.get(bins_).getRealValue(rrv_x.GetName()),hist.weight(),hist.weightError(RooAbsData.SumW2)
	  if hist.weight() != 0: pulls.append(y)
	  bins_+=1
	  
	chi2 = 0
	for p in pulls:
	 chi2+=(p*p)
	 
	print "Chi2/ndof = %f/%f = %f" %(chi2,ndof,chi2/ndof)
	return chi2,ndof
	       
    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist();
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
	  #print x,y
          if(y == 0):
           hpull.SetPoint(ipoint,x,10)
       
   	gt = ROOT.TH1F("gt","gt",int(rrv_x.getBins()/self.narrow_factor),rrv_x.getMin(),rrv_x.getMax());
   	gt.SetMinimum(-3.999);
   	gt.SetMaximum(3.999);
   	gt.SetDirectory(0);
   	gt.SetStats(0);
   	gt.SetLineStyle(0);
   	gt.SetMarkerStyle(20);
   	gt.GetXaxis().SetTitle(rrv_x.GetTitle() + " (GeV)");
   	gt.GetXaxis().SetLabelFont(42);
   	gt.GetXaxis().SetLabelOffset(0.02);
   	gt.GetXaxis().SetLabelSize(0.15);
   	gt.GetXaxis().SetTitleSize(0.15);
   	gt.GetXaxis().SetTitleOffset(1.2);
   	gt.GetXaxis().SetTitleFont(42);
   	gt.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
   	gt.GetYaxis().CenterTitle(True);
   	gt.GetYaxis().SetNdivisions(205);
   	gt.GetYaxis().SetLabelFont(42);
   	gt.GetYaxis().SetLabelOffset(0.007);
   	gt.GetYaxis().SetLabelSize(0.15);
   	gt.GetYaxis().SetTitleSize(0.15);
   	gt.GetYaxis().SetTitleOffset(0.35);
   	gt.GetYaxis().SetTitleFont(42);
   	#gt.GetXaxis().SetNdivisions(505)
	hpull.SetHistogram(gt)
	return hpull
  

    def getData_PoissonInterval(self,data_obs,mplot):
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        datahist   = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone");
        data_histo = datahist.createHistogram("histo_data",rrv_x) ;
        data_histo.SetName("data");
        data_plot  = RooHist(data_histo);
        data_plot.SetMarkerStyle(20);
        data_plot.SetMarkerSize(1);
        
        alpha = 1 - 0.6827;
        for iPoint  in range(data_plot.GetN()):
          	N = data_plot.GetY()[iPoint];
          	if N==0 : 
			L = 0;
         	else : 
			L = (ROOT.Math.gamma_quantile(alpha/2,N,1.));
          	U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1);
          	data_plot.SetPointEYlow(iPoint, N-L);
          	data_plot.SetPointEYhigh(iPoint,U-N);
         	data_plot.SetPointEXlow(iPoint,0);        
          	data_plot.SetPointEXhigh(iPoint,0);        
        
        mplot.addPlotable(data_plot,"PE");

    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self,label, mlvj_region):
        return self.workspace4fit_.pdf("model"+label+mlvj_region+"_"+self.channel+"_mlvj");

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region="_signal_region"):
        print "########### Fixing a general mlvj model  ############"
        rdataset_General_mlvj = self.workspace4fit_.data("rdataset%s%s_%s_mlvj"%(label, mlvj_region,self.channel))
        model_General = self.get_mlvj_Model(label,mlvj_region);
        rdataset_General_mlvj.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_mlvj);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return self.get_mlvj_Model(label,mlvj_region);

    ###### get TTbar model mlvj in a region 
    def get_TTbar_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing TTbar mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_TTbar_xww",mlvj_region);

    ###### get Single Top model mlvj in a region 
    def get_STop_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing Stop mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_STop_xww",mlvj_region);

    ###### get VV mlvj in a region 
    def get_VV_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing VV mlvj for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_VV_xww",mlvj_region);
	
    ### get an mj model from the workspace givin the label
    def get_mj_Model(self,label):
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_mj_Model(self, label ):
        print "########### Fixing a general mj model  ############"
        rdataset_General_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_General = self.get_mj_Model(label);
        rdataset_General_mj.Print();
        model_General.Print();
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_mj);
        par=parameters_General.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                param.Print();
            param.setConstant(kTRUE);
            param=par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    ### fix only the ttbar component using the default label --> for extended pdf
    def get_TTbar_mj_Model(self,label="_TTbar_xww"):
        print "########### Fixing only the TTbar mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the stop component using the default label --> for extended pdf
    def get_STop_mj_Model(self,label="_STop_xww"):
        print "########### Fixing only the Stop mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the VV component using the default label --> for extended pdf
    def get_VV_mj_Model(self,label="_VV_xww"):
        print "########### Fixing only the VV mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the WJets model --> for extended pdf (just fix shape parameters of width, offset of ErfExp and p1 of User1 function
    def get_WJets_mj_Model(self,label):
        print "########### Fixing only the WJets mj Shape --> just the printed parameters  ############"
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_WJets = self.get_mj_Model(label);
        rdataset_WJets_mj.Print();
        model_WJets.Print();
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mj);
        par=parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            if ( paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or paraName.Contains("rrv_p1_User1_WJets")):
             param.setConstant(kTRUE);
             param.Print();
            else:
             param.setConstant(0);
            param=par.Next()
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="total_region",mass_spectrum="_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj  ############"        
        rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
        model = self.get_mlvj_Model(label,mlvj_region);
        model.Print();
        rdataset.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param=par.Next()

    ### fix a pdf in a different way --> for RooAbsPdf 
    def fix_Pdf(self,model_pdf,argset_notparameter):
        print "########### Fixing a RooAbsPdf for mlvj or mj  ############"        
        parameters = model_pdf.getParameters(argset_notparameter);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()


    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters          
    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):
        if TString(mass_spectrum).Contains("_mj"): rrv_x = self.workspace4fit_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4fit_.var("rrv_mass_lvj");

        ## ExpN pdf for W+jets bkg fit
        if in_model_name == "ExpN":
            
            print "########### ExpN funtion for W+jets mlvj ############"
            rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-2e-3,-1e-1,-1e-5);
	    rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 4e3, 0, 1e4);
	    if rrv_x.getMin() == 700:
	       rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 0, -10000, 10000);	       
            #if(ismc==1):
            #  rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 1e3, -1e2, 1e4);	
            #else :	    
            #  if self.channel == "el" :
            #     rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 1e3, -1e2, 1e4);
            #  elif self.wtagger_label == "LP" :
            #     rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 1e3, -1e2, 1e4);
            #  else:
            #     rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, 5e2, 0, 1e3);		 

            model_pdf = ROOT.RooExpNPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ExpN, rrv_n_ExpN);
                                                                                                             
        ## levelled exp for W+jets bkg fit
        if in_model_name == "ExpTail":
            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if self.wtagger_label.find("LP") != -1:
             rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
             rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1e-1,-1.e-2,1e6);
            else:
                if self.channel == "el" :
                 if ismc == 1 and label_tstring.Contains("sb_lo"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 139,0.,355);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2e-2,-1.e-2,5.5e-2);                     
                 elif ismc == 1 and label_tstring.Contains("signal_region"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 162,18,395);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1.6e-2,-1.e-2,5.5e-2);
                 elif ismc == 0 :  
                     rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,70,240);
                     rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);
                           
                if self.channel == "mu" or self.channel == "em":
                 if ismc == 1 and label_tstring.Contains("sb_lo"):
                   #rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 99,10,255);
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-1e-2,7.5e-2);                   
                   #rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-2e-2,7.5e-2);                   
                 elif ismc == 1 and label_tstring.Contains("signal_region"):
                  # rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 110,20,242);
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 110,20,500);
                   #rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1e-2,7.5e-2);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1,7.5e-2);
                 elif ismc == 0 :  
                     rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,40,280);
                     rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);    
      
      	    print "HERE I AM"
	    rrv_s_ExpTail.Print()     
	    rrv_a_ExpTail.Print()     
            model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

        ## sum of two exponential 
        if in_model_name == "Exp" or in_model_name == "Exp_sr":
            print "########### Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,-0.05,-0.1,0.);
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp);

        ## Erf times Exp for mj spectrum
        if in_model_name == "ErfExp" :
            print "########### Erf*Exp for mj fit  ############"
            #rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.05,-0.1,-1e-4);
	    rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,15.,0.,30);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,60.,30.,150);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,50.,30, 100.);
            #model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);
            model_pdf         = ROOT.RooErfExpDecoPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v1" :
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.006,-0.1,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,450.,400.,550.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,70.,10,100.);
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v2" : 
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.005,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,450.,400.,500.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, 50.,10,100.);
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.channel,"rrv_residue_ErfExp"+label+"_"+self.channel,0.,0.,1.);
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)*(1.+TMath::Erf((%s-%s)/%s))/2. "%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(), rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v3" : #different init-value and range
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.005,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,450.,400,500.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, 50.,10,100.);
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.channel,"rrv_residue_ErfExp"+label+"_"+self.channel,0.,0.,1.);
            rrv_high_ErfExp    = RooRealVar("rrv_high_ErfExp"+label+"_"+self.channel,"rrv_high_ErfExp"+label+"_"+self.channel,1.,0.,400);
            rrv_high_ErfExp.setConstant(kTRUE);
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)* TMath::Power( ((1+TMath::Erf((%s-%s)/%s))/2.), %s )"%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(),rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName(), rrv_high_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_high_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

        ## User1 function 
        if in_model_name == "User1":
            print "########### User 1 Pdf  for mlvj fit ############"
            rrv_p0     = RooRealVar("rrv_p0_User1"+label+"_"+self.channel,"rrv_p0_User1"+label+"_"+self.channel, 12, 10, 30);
            if self.wtagger_label=="HP": #change this!
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel,"rrv_p1_User1"+label+"_"+self.channel, -4, -9, -2);
            else:
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel,"rrv_p1_User1"+label+"_"+self.channel, -2.5, -4, 0.);
            model_pdf=RooUser1Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_p0,rrv_p1);

        ## Exp+Gaus or mj spectrum
        if in_model_name == "ExpGaus":
            print "########### Exp + Gaus for mj  fit  ############"
            rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,0.05,-0.2,0.2);
            exp             = ROOT.RooExponential("exp"+label+"_"+self.channel,"exp"+label+"_"+self.channel,rrv_x,rrv_c_Exp);

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,78,90);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,10);
            rrv_high        = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.56,0.,1.);
            gaus            = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            model_pdf       = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))


        ## Erf*Exp + 2Gaus for mj spectrum
	###TTBARMJFIT
        if in_model_name == "2Gaus_ErfExp":

            print "########### 2Gaus + Erf*Exp for mj fit  ############"
            mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            #frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;
	    frac_tmp       = 0; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp, deltamean_tmp-deltamean_tmp_err*10, deltamean_tmp+deltamean_tmp_err*10);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*10, scalesigma_tmp+scalesigma_tmp_err*10);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel,"rrv_frac_2gaus"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            c0_tmp     = -2.9893e-02 ; c0_tmp_err     = 6.83e-03;
            offset_tmp = 7.9350e+01  ; offset_tmp_err = 9.35e+00;
            width_tmp  = 3.3083e+01  ; width_tmp_err  = 2.97e+00;

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-10, width_tmp+10);
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum,"erfexp"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel, 0.5,0.,1.);
	    #only one gauss
            #model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)
	    model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1),RooArgList(rrv_frac),1)
        ## 2Gaus+2Gaus for VV mj spectrum -> WZ and WW
        if in_model_name == "2_2Gaus":

            print "########### 2Gaus +2Gaus for mj fit  ############"
            mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
            #mean1_tmp      = 7.3141e+01; mean1_tmp_err      = 1.63e-01; @@@ JEN            
	    deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_shift = RooRealVar("rrv_shift"+label+"_"+self.channel,"rrv_shift"+label+"_"+self.channel,10.8026) # Z mass: 91.1876; shift=91.1876-80.385=10.8026

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-6, mean1_tmp+6);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,0.,-8,10);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac1 = RooRealVar("rrv_frac1"+label+"_"+self.channel,"rrv_frac1"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*10, frac_tmp+frac_tmp_err*10);
            gausguas_1 =RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum,"gausguas_1"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)

            rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_shift));
            rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean2_gaus, rrv_shift));
            gaus3 = RooGaussian("gaus3"+label+"_"+self.channel,"gaus3"+label+"_"+self.channel, rrv_x,rrv_mean3_gaus,rrv_sigma1_gaus);
            gaus4 = RooGaussian("gaus4"+label+"_"+self.channel,"gaus4"+label+"_"+self.channel, rrv_x,rrv_mean4_gaus,rrv_sigma2_gaus);
            gausguas_2 = RooAddPdf("gausguas_2"+label+"_"+self.channel+mass_spectrum,"gausguas_2"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus3,gaus4),RooArgList(rrv_frac1),1)

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,0.74)#,0.5,1.0)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gausguas_1,gausguas_2),RooArgList(rrv_frac),1)
            
        ## return the pdf
        getattr(self.workspace4fit_,"import")(model_pdf)
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)


    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500):

      ##### define an extended pdf from a standard Roofit One
      print " "
      print "###############################################"
      print "## Make model : ",label," ",in_model_name,"##";
      print "###############################################"
      print " "

      ## call the make RooAbsPdf method
      rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum,"rrv_number"+label+"_"+self.channel+mass_spectrum,area_init_value,0.,1e7);
      model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList,ismc_wjet)
      print "######## Model Pdf ########"        
      model_pdf.Print();
      rrv_number.Print();
      
      ## create the extended pdf
      model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
      print "######## Model Extended Pdf ########"        

      #### put all the parameters ant the shape in the workspace
      getattr(self.workspace4fit_,"import")(rrv_number)
      getattr(self.workspace4fit_,"import")(model) 
      self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();
      ## return the total extended pdf

      return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum);

    ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0):

        print "############### Fit mlvj in mj sideband: ",label," ",mlvj_region,"  ",mlvj_model," ##################"
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset_data_mlvj = self.workspace4fit_.data("rdataset_data_xww%s_%s_mlvj"%(mlvj_region,self.channel))

        ## get the minor component shapes in the sb low
        model_VV_backgrounds    = self.get_VV_mlvj_Model("_sb_lo");
        number_VV_sb_lo_mlvj    = self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel))
        model_TTbar_backgrounds = self.get_TTbar_mlvj_Model("_sb_lo");
        number_TTbar_sb_lo_mlvj = self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel))
        model_STop_backgrounds  = self.get_STop_mlvj_Model("_sb_lo");
        number_STop_sb_lo_mlvj  = self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel))

	
        self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).Print();

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_sb_lo_from_fitting"%(label), mlvj_model,"_mlvj");
        model_pdf_WJets.Print();
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_sb_lo = self.workspace4fit_.var("rrv_number%s_sb_lo_%s_mlvj"%(label,self.channel)).clone("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel));
        model_WJets =RooExtendPdf("model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),"model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),model_pdf_WJets,number_WJets_sb_lo);
        model_pdf_WJets.Print();
        number_WJets_sb_lo.Print()

        ## Add the other bkg component fixed to the total model
        model_data = RooAddPdf("model_data%s%s_mlvj"%(label,mlvj_region),"model_data%s%s_mlvj"%(label,mlvj_region),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds));
        
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
	fitresultsmlvj.append(rfresult)
        getattr(self.workspace4fit_,"import")(model_data)
        self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).Print();

        model_WJets.Print();
        model_WJets.getParameters(rdataset_data_mlvj).Print("v");
        self.workspace4fit_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label,self.channel)).getParameters(rdataset_data_mlvj).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),"rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),
                                                 self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getVal() );

        rrv_number_data_sb_lo_mlvj.setError( TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getError()));

        getattr(self.workspace4fit_,"import")(rrv_number_data_sb_lo_mlvj)

        ### plot for WJets default + default shape
        if TString(label).Contains("_WJets0"):

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0) );

            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent)) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent)) ;

            model_data.plotOn(mplot, RooFit.Components("model_VV_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent));

            model_data.plotOn(mplot, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent));


            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(), RooAbsReal.NumEvent)) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(), RooAbsReal.NumEvent)) ;

            model_data.plotOn(mplot, RooFit.Components("model_VV_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(), RooAbsReal.NumEvent));

            model_data.plotOn(mplot, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(), RooAbsReal.NumEvent));
 
            ### draw the error band 
            draw_error_band(rdataset_data_mlvj, model_data,self.workspace4fit_.var("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible());
            self.getData_PoissonInterval(rdataset_data_mlvj,mplot);

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);  

            ### Add the legend to the plot 
            #self.leg=self.legend4Plot(mplot,0,1,0., 0.16, 0.16, 0.1);
	    self.leg=self.legend4Plot(mplot,0,1,0.25,0.,0.1,0.,0);

            mplot.addObject(self.leg)

            ### calculate the chi2
            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            #print mplot.chiSquare();
            #print "#################### nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            ### write the result in the output
            #self.file_out.write("\n fit_mlvj_in_Mj_sideband: nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof ) );

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_mass_lvj,mplot);
            parameters_list = model_data.getParameters(rdataset_data_mlvj);
            datahist = rdataset_data_mlvj.binnedClone( rdataset_data_mlvj.GetName()+"_binnedClone",rdataset_data_mlvj.GetName()+"_binnedClone" )
            
	    mplot.GetXaxis().SetTitle(rrv_mass_lvj.GetTitle() + " (GeV)")
	    mplot.GetYaxis().SetRangeUser(2e-3,2e5)
            #@#self.draw_canvas_with_pull( rrv_mass_lvj,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), "m_lvj_sb_lo%s"%(label),"",1,1)
	    self.draw_canvas_with_pull( rrv_mass_lvj,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), "m_lvj_sb_lo%s"%(label),"",1,1,0,1)
	    

        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sb_lo_from_fitting_mlvj"%(label));
	purity = self.wtagger_label[0]+self.wtagger_label[1]
        Deco      = PdfDiagonalizer("Deco%s_sb_lo_from_fitting_%s_%s_mlvj_13TeV"%(label,self.channel,purity),wsfit_tmp,rfresult);
        print"#################### diagonalize data sideband fit "
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets);
        print "##################### workspace for decorrelation ";
        wsfit_tmp.Print("v");
        print "##################### original  parameters ";
        model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
        print "##################### original  decorrelated parameters ";
        model_pdf_WJets_deco.getParameters(rdataset_data_mlvj).Print("v");
        print "##################### original  pdf ";
        model_pdf_WJets.Print();
        print "##################### decorrelated pdf ";
        model_pdf_WJets_deco.Print();


        getattr(self.workspace4fit_,"import")(model_pdf_WJets_deco);


        #### Call the alpha evaluation in automatic
        self.get_WJets_mlvj_correction_sb_lo_to_signal_region(label,mlvj_model);

        ### Fix the pdf of signal, TTbar, STop and VV in the signal region 
        self.fix_Model("_TTbar_xww","_signal_region","_mlvj")
        self.fix_Model("_STop_xww","_signal_region","_mlvj")
        self.fix_Model("_VV_xww","_signal_region","_mlvj")

        ### Call the evaluation of the normalization in the signal region for TTbar, VV, STop, and WJets after the extrapolation via alpha
        self.get_mlvj_normalization_insignalregion("_TTbar_xww");
        self.get_mlvj_normalization_insignalregion("_STop_xww");
        self.get_mlvj_normalization_insignalregion("_VV_xww");
        self.get_mlvj_normalization_insignalregion(label,"model_pdf%s_signal_region_%s_after_correct_mlvj"%(label,self.channel));    

    ##### Function that calculate the normalization inside the mlvj signal region (mass window around the resonance in order to fill datacards)
    def get_mlvj_normalization_insignalregion(self, label, model_name=""):
        
        print "############### get mlvj normalization inside SR ",label," ",model_name," ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_signal_region"+"_"+self.channel+"_mlvj");
        else:
            model = self.workspace4fit_.pdf(model_name);

        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj");
        fullInt   = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj) );
        signalInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("signal_region"));
        highMassInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("high_mass"));
	over3500Int = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("over3500"))
        
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        highMassInt_val = highMassInt.getVal()/fullInt_val 
	over3500Int_val = over3500Int.getVal()/fullInt_val

        ## integral in the signal region
        print "######### integral in SR: ",label+"signalInt=%s"%(signalInt_val)

        print "####### Events Number in MC Dataset:"
        self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").Print();
	if 'WJets01' not in label:
        	print "########## Events Number get from fit:"
        	rrv_tmp		= self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_postfit_signal_region");
		rrv_tmp_pre	= self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_prefit_signal_region");

        	print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)
		rrv_tmp.setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_postfit_signal_region").getVal()*(signalInt_val-over3500Int_val))
		rrv_tmp_pre.setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_prefit_signal_region").getVal()*(signalInt_val-over3500Int_val))
	else:
		rrv_tmp		= self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj")

        #### store the info in the output file
        self.file_out.write( "\n%s++++++++++++++++++get_mlvj_normalization_insignalregion++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in All Region from dataset : %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset: %s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getVal()) )
        self.file_out.write( "\nRatio signal_region/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal() ) )
        self.file_out.write( "\nEvents Number in All Region from fitting : %s\n"%(rrv_tmp.getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting: %s\n"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nEvents Number in High Mass Region from fitting: %s\n"%(rrv_tmp.getVal()*highMassInt_val) )
        self.file_out.write( "\nRatio signal_region/all_range from fitting :%s"%(signalInt_val ) )
#LUCA
        if not self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj"):
            rrv_number_fitting_signal_region_mlvj = RooRealVar("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj","rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj", rrv_tmp.getVal()*signalInt_val );
            getattr(self.workspace4fit_,"import")(rrv_number_fitting_signal_region_mlvj);
        else :
            self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj").setVal(rrv_tmp.getVal()*signalInt_val);


        self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj").Print();

    ### method to get the alpha function to extrapolate the wjets in the signal region
    def get_WJets_mlvj_correction_sb_lo_to_signal_region(self,label, mlvj_model):

        print" ############# get the extrapolation function alpha from MC : ",label,"   ",mlvj_model," ###############";          

        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here 
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        rdataset_WJets_sb_lo_mlvj = self.workspace4fit_.data("rdataset4fit%s_sb_lo_%s_mlvj"%(label,self.channel))
        rdataset_WJets_signal_region_mlvj = self.workspace4fit_.data("rdataset4fit%s_signal_region_%s_mlvj"%(label,self.channel))

        ### create a frame for the next plots 
        mplot = rrv_x.frame(RooFit.Title("correlation_pdf"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor))) ;
        mplot.GetYaxis().SetTitle("F_{W+jets}^{SR,MC},F_{W+jets}^{SB,MC} (Arbitrary units)");

        if mlvj_model=="Exp":
            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Exp%s_sb_lo_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Exp%s_%s"%(label,self.channel),"rrv_delta_c_Exp%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal(),
                                      self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                                      self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )

            correct_factor_pdf = RooExponential("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);
            
        if mlvj_model=="Pow":

            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Pow%s_sb_lo_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Pow%s_%s"%(label,self.channel),"rrv_delta_c_Pow%s_%s"%(label,self.channel),0., -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            correct_factor_pdf = RooPowPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);
		
        if mlvj_model=="ExpN":
            rrv_c_sb  = self.workspace4fit_.var("rrv_c_ExpN%s_sb_lo_%s"%(label,self.channel));
            rrv_n_sb  = self.workspace4fit_.var("rrv_n_ExpN%s_sb_lo_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_ExpN%s_%s"%(label,self.channel),"rrv_delta_c_ExpN%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal(),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )
            rrv_delta_n = RooRealVar("rrv_delta_n_ExpN%s_%s"%(label,self.channel),"rrv_delta_n_ExpN%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal(),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal()-4*rrv_n_sb.getError(),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal()+4*rrv_n_sb.getError() )

            correct_factor_pdf = RooExpNPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c, rrv_delta_n);
 
        if mlvj_model=="ExpTail":
            rrv_s_sb =self.workspace4fit_.var("rrv_s_ExpTail%s_sb_lo_%s"%(label,self.channel));
            rrv_a_sb =self.workspace4fit_.var("rrv_a_ExpTail%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_s = RooRealVar("rrv_delta_s_ExpTail%s_%s"%(label,self.channel),"rrv_delta_s_ExpTail%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal(),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal()-4*rrv_s_sb.getError(),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal()+4*rrv_s_sb.getError() )
            rrv_delta_a = RooRealVar("rrv_delta_a_ExpTail%s_%s"%(label,self.channel),"rrv_delta_a_ExpTail%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal(),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal()-4*rrv_a_sb.getError(),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal()+4*rrv_a_sb.getError() )
                     
            rrv_a_sr = RooFormulaVar("rrv_a_ExpTail_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_a_sb, rrv_delta_a ) );
            rrv_s_sr = RooFormulaVar("rrv_s_ExpTail_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_s_sb, rrv_delta_s ) );

            correct_factor_pdf = RooAlpha4ExpTailPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_s_sr, rrv_a_sr, rrv_s_sb, rrv_a_sb);

        if mlvj_model=="ErfExp_v1":

            rrv_c_sb       = self.workspace4fit_.var("rrv_c_ErfExp%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb  = self.workspace4fit_.var("rrv_offset_ErfExp%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb   = self.workspace4fit_.var("rrv_width_ErfExp%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfExp%s_%s"%(label,self.channel),"rrv_delta_c_ErfExp%s_%s"%(label,self.channel),0.,
                                           -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfExp%s_%s"%(label,self.channel),"rrv_delta_offset_ErfExp%s_%s"%(label,self.channel),0.,
                                           -100*rrv_offset_sb.getError(),100*rrv_offset_sb.getError());
            rrv_delta_width = RooRealVar("rrv_delta_width_ErfExp%s_%s"%(label,self.channel),"rrv_delta_width_ErfExp%s_%s"%(label,self.channel),0.,
                                          -100*rrv_width_sb.getError(),100*rrv_width_sb.getError());

            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c_sb, rrv_delta_c ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr,rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb, rrv_x.getMin(), rrv_x.getMax());
            
        if mlvj_model=="ErfPow":

            rrv_c_sb      = self.workspace4fit_.var("rrv_c_ErfPow%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfPow%s_%s"%(label,self.channel),"rrv_delta_c_ErfPow%s_%s"%(label,self.channel),0.,
                                          -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPow%s_%s"%(label,self.channel),"rrv_delta_offset_ErfPow%s_%s"%(label,self.channel),0.,
                                          -100*rrv_offset_sb.getError(),100*rrv_offset_sb.getError());
            rrv_delta_width  = RooRealVar("rrv_delta_width_ErfPow%s_%s"%(label,self.channel),"rrv_delta_width_ErfPow%s_%s"%(label,self.channel),0.,
                                          -100*rrv_width_sb.getError(),100*rrv_width_sb.getError());
            
            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c_sb, rrv_delta_c ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPowPdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr,rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="ErfPow2":

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow2%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c0      = RooRealVar("rrv_delta_c0_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_c0_ErfPow2%s_%s"%(label,self.channel),-8, -20 ,0);
            rrv_delta_c1      = RooRealVar("rrv_delta_c1_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_c1_ErfPow2%s_%s"%(label,self.channel),0., -5, 5);
            rrv_delta_offset  = RooRealVar("rrv_delta_offset_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_offset_ErfPow2%s_%s"%(label,self.channel),30., 1.,80 );
            rrv_delta_width   = RooRealVar("rrv_delta_width_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_width_ErfPow2%s_%s"%(label,self.channel),15,1.,100*rrv_width_sb.getError());
            
            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPow2Pdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr,rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="ErfPowExp": ## take initial value from what was already fitted in the SR

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPowExp%s_sb_lo_%s"%(label,self.channel));
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPowExp%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPowExp%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPowExp%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c0  = RooRealVar("rrv_delta_c0_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_c0_ErfPowExp%s_%s"%(label,self.channel),
                                        self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal(),
                                        self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal()-4*rrv_c0_sb.getError(),
                                        self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal()+4*rrv_c0_sb.getError() )

            rrv_delta_c1 = RooRealVar("rrv_delta_c1_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_c1_ErfPowExp%s_%s"%(label,self.channel),
                                       self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal(),
                                       self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal()-4*rrv_c1_sb.getError(),
                                       self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal()+4*rrv_c1_sb.getError() )

            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_offset_ErfPowExp%s_%s"%(label,self.channel),
                                       self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal(),
                                       self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal()-4*rrv_offset_sb.getError(),
                                       self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal()+4*rrv_offset_sb.getError())

            rrv_delta_width = RooRealVar("rrv_delta_width_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_width_ErfPowExp%s_%s"%(label,self.channel),
                                         self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal(),
                                         self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal()-4*rrv_width_sb.getError(),
                                         self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal()+4*rrv_width_sb.getError() )

            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPowExpPdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr,rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb);
 
        ### define the category and do the simultaneous fit taking the combined dataset of events in mlvj sb and sr
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");
        combData4fit = self.workspace4fit_.data("combData4fit%s_%s"%(label,self.channel));


        model_pdf_sb_lo_WJets         = self.workspace4fit_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label,self.channel));
        model_pdf_signal_region_WJets = RooProdPdf("model_pdf%s_signal_region_%s_mlvj"%(label,self.channel),"model_pdf%s_signal_region_%s_mlvj"%(label,self.channel) ,model_pdf_sb_lo_WJets,correct_factor_pdf);

        simPdf = RooSimultaneous("simPdf","simPdf",data_category);
        simPdf.addPdf(model_pdf_sb_lo_WJets,"sideband");
        simPdf.addPdf(model_pdf_signal_region_WJets,"signal_region");
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE));
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
	fitresultsfinal.append(rfresult)
	#raw_input('...')

        ### Decorrelate the parameters in the alpha shape
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sim_mlvj"%(label));
        print "############### diagonalizer alpha ";
        Deco      = PdfDiagonalizer("Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label),wsfit_tmp,rfresult);
        correct_factor_pdf_deco = Deco.diagonalize(correct_factor_pdf);
        print "##################### workspace for decorrelation ";
        wsfit_tmp.Print("v");
        print "##################### original  parameters ";
        correct_factor_pdf.getParameters(rdataset_WJets_signal_region_mlvj).Print("v");
        print "##################### original  decorrelated parameters ";
        correct_factor_pdf_deco.getParameters(rdataset_WJets_signal_region_mlvj).Print("v");
        print "##################### original  pdf ";
        correct_factor_pdf.Print();
        print "##################### decorrelated pdf ";
        correct_factor_pdf_deco.Print();

        getattr(self.workspace4fit_,"import")(correct_factor_pdf_deco);
                     
        ## in case of default Wjets with default shape
        if TString(label).Contains("_WJets0"):

            ### only mc plots in the SB region
            mplot_sb_lo = rrv_x.frame(RooFit.Title("WJets sb low"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
            
            rdataset_WJets_sb_lo_mlvj.plotOn(mplot_sb_lo, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_sb_lo_WJets.plotOn(mplot_sb_lo);
            mplot_pull_sideband = self.get_pull(rrv_x,mplot_sb_lo);
            parameters_list     = model_pdf_sb_lo_WJets.getParameters(rdataset_WJets_sb_lo_mlvj);
            mplot_sb_lo.GetYaxis().SetRangeUser(1e-2,mplot_sb_lo.GetMaximum()*1.2);

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot_sb_lo.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            #print mplot_sb_lo.chiSquare();
            #print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot_sb_lo.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            datahist = rdataset_WJets_sb_lo_mlvj.binnedClone( rdataset_WJets_sb_lo_mlvj.GetName()+"_binnedClone",rdataset_WJets_sb_lo_mlvj.GetName()+"_binnedClone" )

            self.draw_canvas_with_pull( rrv_x,datahist,mplot_sb_lo, mplot_pull_sideband,ndof,parameters_list,"%s/other/"%(self.plotsDir), "m_lvj%s_sb_lo_sim"%(label),"",1,1)
            ### only mc plots in the SR region
            mplot_signal_region = rrv_x.frame(RooFit.Title("WJets sr"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));

            rdataset_WJets_signal_region_mlvj.plotOn(mplot_signal_region, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_signal_region_WJets.plotOn(mplot_signal_region);
            mplot_pull_signal_region = self.get_pull(rrv_x, mplot_signal_region);
            parameters_list = model_pdf_signal_region_WJets.getParameters(rdataset_WJets_signal_region_mlvj);
            mplot_signal_region.GetYaxis().SetRangeUser(1e-2,mplot_signal_region.GetMaximum()*1.2);

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot_signal_region.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            #print mplot_signal_region.chiSquare();
            #print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot_signal_region.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );

            self.draw_canvas_with_pull( rrv_x,datahist,mplot_signal_region, mplot_pull_signal_region,ndof,parameters_list,"%s/other/"%(self.plotsDir), "m_lvj%s_signal_region_sim"%(label),"",1,1);

        ### Total plot shape in sb_lo, sr and alpha
	alpha_norm_factor = 1
        model_pdf_sb_lo_WJets.plotOn(mplot,RooFit.Name("Sideband"));
        model_pdf_signal_region_WJets.plotOn(mplot, RooFit.LineColor(kRed), RooFit.Name("Signal Region"));
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("Transfer function #alpha"), RooFit.Normalization(alpha_norm_factor,RooAbsReal.Relative));

        ### plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
        if TString(label).Contains("_WJets0"):
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha: Alternate PS") );

            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha: Alternate Function") , RooFit.Normalization(alpha_norm_factor,RooAbsReal.Relative));

        paras=RooArgList();

        if mlvj_model=="ErfExp_v1" or mlvj_model=="ErfPow" or mlvj_model=="2Exp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig4"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig5"%(label,self.channel, self.wtagger_label) ));

        if mlvj_model=="ErfPow2" or mlvj_model=="ErfPowExp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig4"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig5"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig6"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig7"%(label,self.channel, self.wtagger_label) ));

        if mlvj_model=="Exp" or mlvj_model=="Pow":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));

        if mlvj_model=="ExpN" or mlvj_model=="ExpTail" or mlvj_model=="Pow2":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
        
        if TString(label).Contains("_WJets0") or TString(label).Contains("_WJets1"): ### draw error band at 1 and 2 sigma using the decorrelated shape
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha #pm",20,400,alpha_norm_factor);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,2 ,mplot,kGreen+2,"F",3002,"#alpha #pm",20,400,alpha_norm_factor);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha_invisible #pm",20,400,alpha_norm_factor);
            
        ### plot on the same canvas
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha_invisible"), RooFit.Normalization(alpha_norm_factor,RooAbsReal.Relative))

        if TString(label).Contains("_WJets0") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") , RooFit.Normalization(alpha_norm_factor,RooAbsReal.Relative));

        elif TString(label).Contains("_WJets01") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        ### Add the legend
        self.leg=self.legend4Plot(mplot,1,0, 0, 0., 0., -0.1, 0., True);
        mplot.addObject(self.leg);
        
        ## set the Y axis in arbitrary unit 
	tmp_y_max=0.125
	tmp_y_min=1e-5
        mplot.GetYaxis().SetRangeUser(tmp_y_min,tmp_y_max);

        #### Draw another axis with the real value of alpha
        model_pdf_sb_lo_WJets.getVal(RooArgSet(rrv_x))
        model_pdf_signal_region_WJets.getVal(RooArgSet(rrv_x))
        correct_factor_pdf_deco.getVal(RooArgSet(rrv_x))
        tmp_alpha_ratio = ( model_pdf_signal_region_WJets.getVal(RooArgSet(rrv_x))/model_pdf_sb_lo_WJets.getVal(RooArgSet(rrv_x)) );
        tmp_alpha_pdf   = correct_factor_pdf_deco.getVal(RooArgSet(rrv_x)) * mplot.getFitRangeBinW(); ## value of the pdf in each point
        tmp_alpha_scale = tmp_alpha_ratio/tmp_alpha_pdf;

        #add alpha scale axis
        axis_alpha=TGaxis( rrv_x.getMax(), 0, rrv_x.getMax(), tmp_y_max, tmp_y_min * tmp_alpha_scale, tmp_y_max*tmp_alpha_scale, 510, "+L" ); #-,-+,+,L
        axis_alpha.SetTitle("Transfer function #alpha");
        axis_alpha.SetTitleOffset(0.5);
        axis_alpha.SetTitleSize(0.05);
        axis_alpha.SetLabelSize(0.045);
        axis_alpha.SetTitleFont(42);
        axis_alpha.SetLabelFont(42);
        #axis_alpha.RotateTitle(1);
        mplot.addObject(axis_alpha);

        self.draw_canvas(mplot,"%s/other/"%(self.plotsDir),"correction_pdf%s_%s_M_lvj_signal_region_to_sideband"%(label,mlvj_model),0,1,0,1);
	#@#also with log scale
	self.draw_canvas(mplot,"%s/other/"%(self.plotsDir),"correction_pdf%s_%s_M_lvj_signal_region_to_sideband"%(label,mlvj_model),0,0,0,1)
	#@#

	#@#make nice plots with log-scale
	#tmp='''
	
	#'''
	#@#

        correct_factor_pdf_deco.getParameters(rdataset_WJets_sb_lo_mlvj).Print("v");
        model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco = self.workspace4fit_.pdf("model_pdf%s_sb_lo_from_fitting_%s_mlvj_Deco%s_sb_lo_from_fitting_%s_%s_mlvj_13TeV"%(label,self.channel,label, self.channel,self.wtagger_label[0]+self.wtagger_label[1]));
        model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco.Print("v");

        ### Wjets shape in the SR correctedfunction * sb 
        model_pdf_WJets_signal_region_after_correct_mlvj = RooProdPdf("model_pdf%s_signal_region_%s_after_correct_mlvj"%(label,self.channel),"model_pdf%s_signal_region_%s_after_correct_mlvj"%(label,self.channel),model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco,self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label)) );
        model_pdf_WJets_signal_region_after_correct_mlvj.Print("v")
        ### fix the parmaters and import in the workspace
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_signal_region_after_correct_mlvj)

        ##### calculate the normalization and alpha for limit datacard
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).Print();
        self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)).Print();
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).setVal(self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)).getVal());
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).setError(self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)).getError());

        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).setConstant(kTRUE);

    ##### Counting of the events of each component in the signal region taking the label for the model
    def get_mj_normalization_insignalregion(self, label, forlimit):
        print "################## get mj normalization ",label," ################## ";
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        model      = self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj");

        fullInt   	= model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        sb_loInt  	= model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("sb_lo"));
        signalInt 	= model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signal_region"));
        sb_hiInt  	= model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("sb_hi"));

        
        fullInt_val   	= fullInt.getVal()
        sb_loInt_val  	= sb_loInt.getVal()/fullInt_val
        sb_hiInt_val  	= sb_hiInt.getVal()/fullInt_val
        signalInt_val 	= signalInt.getVal()/fullInt_val


        print "########### Events Number in MC Dataset: #############"
        self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").Print();
        self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").Print();
        self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").Print();

        print "########### Events Number get from fit: ##############"
	if forlimit==1 and 'data' not in label and 'WJets01' not in label:
	        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_postfit");

		rrv_tmp_signal	= rrv_tmp.clone("rrv_number"+label+"_"+self.channel+"_mj_postfit_signal_region")
		rrv_tmp_signal.setVal(rrv_tmp_signal.getVal()*signalInt_val)
		rrv_tmp_signal.setError(rrv_tmp.getError()/rrv_tmp.getVal()*rrv_tmp_signal.getVal())
		getattr(self.workspace4fit_,'import')(rrv_tmp_signal)
		rrv_tmp_pre = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj_prefit");
		rrv_tmp_signal_pre	= rrv_tmp_pre.clone("rrv_number"+label+"_"+self.channel+"_mj_prefit_signal_region")
		rrv_tmp_signal_pre.setVal(rrv_tmp_signal_pre.getVal()*signalInt_val)
		rrv_tmp_signal_pre.setError(rrv_tmp_pre.getError()/rrv_tmp_pre.getVal()*rrv_tmp_signal_pre.getVal())
		getattr(self.workspace4fit_,'import')(rrv_tmp_signal_pre)
	else:
	        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj");
		rrv_tmp_signal	= rrv_tmp.clone("rrv_number"+label+"_"+self.channel+"_mj_signal_region")
		rrv_tmp_signal.setVal(rrv_tmp_signal.getVal()*signalInt_val)
        rrv_tmp.Print();



        print "Events Number in sideband_low :%s"%(rrv_tmp.getVal()*sb_loInt_val)
        print "Events Number in Signal Region:%s"%(rrv_tmp.getVal()*signalInt_val)
        print "Events Number in sideband_high:%s"%(rrv_tmp.getVal()*sb_hiInt_val)
        print "Total Number in sidebands :%s"%(rrv_tmp.getVal()*(sb_loInt_val+sb_hiInt_val) )
        print "Ratio signal_region/sidebands :%s"%(signalInt_val/(sb_loInt_val+sb_hiInt_val) )

	if forlimit==1:
        ##### Save numbers in the output text file
		self.file_out.write( "\n%s++++++++++++++++get_mj_normalization_insignalregion++++++++++++++++++++"%(label) )
		self.file_out.write( "\nEvents Number in sideband_low from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal() ) )
		self.file_out.write( "\nEvents Number in Signal Region from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal() ) )
		self.file_out.write( "\nEvents Number in sideband_high from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal() ) )
		self.file_out.write( "\nTotal Number in sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal() ) )
		self.file_out.write( "\nTotal number of events:%s"%rrv_tmp.getVal())
		if (self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()) != 0:
		   self.file_out.write( "\nRatio signal_region/sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()/(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()) ) )

		self.file_out.write( "\nEvents Number in sideband_low from fitting:%s"%(rrv_tmp.getVal()*sb_loInt_val) )
		self.file_out.write( "\nEvents Number in Signal Region from fitting:%s"%(rrv_tmp.getVal()*signalInt_val) )
		self.file_out.write( "\nEvents Number in sideband_high from fitting:%s"%(rrv_tmp.getVal()*sb_hiInt_val) )
		self.file_out.write( "\nTotal Number in sidebands from fitting:%s"%(rrv_tmp.getVal()*(sb_loInt_val+sb_hiInt_val) ) )
		self.file_out.write( "\nRatio signal_region/sidebands from fitting:%s"%(signalInt_val/(sb_loInt_val+sb_hiInt_val) ) )

    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0): # to get the normalization of WJets in signal_region

	if options.prefitplot:
		self.mj_prefit_plot("_WJets0_xww")
		exit(0)

        print "############### Fit mj Normalization ##################"
        ## fit the two version of pdf for Wjets shape if available
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww");
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets01_xww");

        ## in case fit also the scaled jet mass distributions in order to have the jet mass scale sys included
        if scaleJetMass :
         self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww_massup","massup");
         self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww_massdn","massdn");
         self.fit_WJetsNormalization_in_Mj_signal_region("_WJets1_xww");

        ## take the normalization numbers
        rrv_WJets0  = self.workspace4fit_.var("rrv_number_WJets0_xww_in_mj_signal_region_from_fitting_%s"%(self.channel));
        rrv_WJets01 = self.workspace4fit_.var("rrv_number_WJets01_xww_in_mj_signal_region_from_fitting_%s"%(self.channel));
        rrv_WJets0.Print();
        rrv_WJets01.Print();
        if scaleJetMass :
         rrv_WJets1 = self.workspace4fit_.var("rrv_number_WJets1_xww_in_mj_signal_region_from_fitting_%s"%(self.channel));
         rrv_WJets1.Print();
         rrv_WJets0massup.Print();
         rrv_WJets0massdn.Print();

        ### total uncertainty combining the result with two different shapes
	#@#increase W+Jets uncertainty from fit by 20% since not all shape parameters were floating
        total_uncertainty = TMath.Sqrt( TMath.Power(rrv_WJets0.getError()*1.2,2) + TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(),2) );
	additional_uncert	= RooRealVar("alt_func_uncert_%s_%s"%(self.channel,self.wtagger_label),"alt_func_uncert_%s_%s"%(self.channel,self.wtagger_label),rrv_WJets01.getVal()-rrv_WJets0.getVal())
	stat_uncert		= RooRealVar("stat_uncert_%s_%s"%(self.channel,self.wtagger_label),"stat_uncert_%s_%s"%(self.channel,self.wtagger_label),1.2*rrv_WJets0.getError())
	getattr(self.workspace4limit_,'import')(additional_uncert)
	getattr(self.workspace4limit_,'import')(stat_uncert)
        rrv_WJets0.setError(total_uncertainty);


	wjets_tmp_uncert = self.workspace4fit_.var("rrv_number_WJets0_xww_%s_mj_postfit_signal_region"%self.channel).clone("WJets_norm_noly_fit_uncertainty")
	getattr(self.workspace4limit_,'import')(wjets_tmp_uncert)
	self.workspace4fit_.var("rrv_number_WJets0_xww_%s_mj_postfit_signal_region"%self.channel).setError(total_uncertainty)
        rrv_WJets0.Print();

        ##jet mass uncertainty on WJets normalization and the other bkg component
        if self.workspace4fit_.var("rrv_number_WJets0_xww_massup_in_mj_signal_region_from_fitting_%s"%(self.channel)) and self.workspace4fit_.var("rrv_number_WJets0_xww_massdn_in_mj_signal_region_from_fitting_%s"%(self.channel)):            
          rrv_WJets0massup = self.workspace4fit_.var("rrv_number_WJets0_xww_massup_in_mj_signal_region_from_fitting_%s"%(self.channel));
          rrv_WJets0massdn = self.workspace4fit_.var("rrv_number_WJets0_xww_massdn_in_mj_signal_region_from_fitting_%s"%(self.channel));
          self.WJets_normlization_uncertainty_from_jet_mass= ( TMath.Abs(rrv_WJets0massup.getVal()-rrv_WJets0.getVal())+TMath.Abs(rrv_WJets0massdn.getVal()-rrv_WJets0.getVal() ) )/2./rrv_WJets0.getVal();

        rrv_STop  = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww__%s_mj"%(self.channel));

        if self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massup_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massdn_%s_mj"%(self.channel)) :
         rrv_STopmassup = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massup_%s_mj"%(self.channel));
         rrv_STopmassdn = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massdn_%s_mj"%(self.channel));
         self.STop_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_STopmassup.getVal()-rrv_STop.getVal())+TMath.Abs(rrv_STopmassdn.getVal()-rrv_STop.getVal() ) )/2./rrv_STop.getVal();

        rrv_TTbar = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww__%s_mj"%(self.channel));
        if self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massup_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massdn_%s_mj"%(self.channel)):
         rrv_TTbarmassup = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massup_%s_mj"%(self.channel));
         rrv_TTbarmassdn = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massdn_%s_mj"%(self.channel));
         self.TTbar_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_TTbarmassup.getVal()-rrv_TTbar.getVal())+TMath.Abs(rrv_TTbarmassdn.getVal()-rrv_TTbar.getVal() ) )/2./rrv_TTbar.getVal();


    def simultaneous_fit_mj_and_mlvj_sb(self,label,mlvj_region,mlvj_model):
	###
	###
	#mj model
	###
	###
	massscale=''
        print "############### Fit mj Normalization: ",label," ",massscale," ##################"
        rrv_mass_j 	= self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj	= self.workspace4fit_.var("rrv_mass_lvj")

        ## get real data in mj distribution --> mass up and down have only an effect on Wjets shape -> effect on the normalization -> evaluated in the MC and fit data
        rdataset_data_mj	= self.workspace4fit_.data("rdataset_data_xww_%s_mj"%(self.channel))
	rdataset_data_mlvj	= self.workspace4fit_.data('rdataset_data_xww_sb_lo_%s_mlvj'%self.channel)

        ### Fix TTbar, VV and STop
        model_TTbar_mj = self.get_TTbar_mj_Model("_TTbar_xww"+massscale);
        model_STop_mj  = self.get_STop_mj_Model("_STop_xww"+massscale);
        model_VV_mj    = self.get_VV_mj_Model("_VV_xww"+massscale);

        ## only two parameters are fix, offset and width, while the exp is floating , otherwise if shape different User1 or ErfExp everything is flaoting
        model_WJets_mj = self.get_WJets_mj_Model(label);
	model_WJets_mj.Print()

        ## Total Pdf and fit 
        model_mj	 = RooAddPdf("model_data_xww%s_%s_mj"%(massscale,self.channel),"model_data_xww%s_%s_mj"%(massscale,self.channel),RooArgList(model_WJets_mj,model_VV_mj,model_TTbar_mj,model_STop_mj))
	#@#save value of rrv_number_XX from MC-fit so fit with alternative function has same mean+constraint 
	if '_WJets0_' in label:
		meanval_ttbar	= RooRealVar('meanval_ttbar','meanval_ttbar',self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%(self.channel)).getVal())
		meanval_VV	= RooRealVar('meanval_VV','meanval_VV',self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%(self.channel)).getVal())
		meanval_STop	= RooRealVar('meanval_STop','meanval_STop',self.workspace4fit_.var("rrv_number_STop_xww_%s_mj"%(self.channel)).getVal())
		ttbar_error	= RooRealVar('ttbar_error','ttbar_error',self.workspace4fit_.var('rrv_number_TTbar_xww_%s_mj'%(self.channel)).getError())
		VV_error	= RooRealVar('VV_error','VV_error',self.workspace4fit_.var('rrv_number_VV_xww_%s_mj'%(self.channel)).getError())
		STop_error	= RooRealVar('STop_error','STop_error',self.workspace4fit_.var('rrv_number_STop_xww_%s_mj'%(self.channel)).getError())
		meanval_ttbar.setError(ttbar_error.getVal())
		meanval_VV.setError(VV_error.getVal())
		meanval_STop.setError(STop_error.getVal())
		getattr(self.workspace4fit_,'import')(meanval_ttbar)
		getattr(self.workspace4fit_,'import')(ttbar_error)
		getattr(self.workspace4fit_,'import')(meanval_VV)
		getattr(self.workspace4fit_,'import')(VV_error)
		getattr(self.workspace4fit_,'import')(meanval_STop)
		getattr(self.workspace4fit_,'import')(STop_error)
		#@# save prefit values
		N_STop_pre	= self.workspace4fit_.var("rrv_number_STop_xww_%s_mj"%self.channel).clone("rrv_number_STop_xww_%s_mj_prefit"%self.channel)
		N_VV_pre	= self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%self.channel).clone("rrv_number_VV_xww_%s_mj_prefit"%self.channel)
		N_TTbar_pre	= self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%self.channel).clone("rrv_number_TTbar_xww_%s_mj_prefit"%self.channel)
		N_WJets_pre	= self.workspace4fit_.var("rrv_number_WJets0_xww_%s_mj"%self.channel).clone("rrv_number_WJets0_xww_%s_mj_prefit"%self.channel)
		getattr(self.workspace4fit_,'import')(N_STop_pre)
		getattr(self.workspace4fit_,'import')(N_VV_pre)
		getattr(self.workspace4fit_,'import')(N_TTbar_pre)
		getattr(self.workspace4fit_,'import')(N_WJets_pre)
		getattr(self.workspace4limit_,'import')(N_VV_pre)
	#@#set value and error of rrv_number_ttbar/VV to prefit for WJets01-fit
	elif '_WJets01_' in label:
		self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var('meanval_ttbar').getVal())
		self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setError(self.workspace4fit_.var('ttbar_error').getVal())
		self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var('meanval_VV').getVal())
		self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setError(self.workspace4fit_.var('VV_error').getVal())

	mean_ttbar	= RooRealVar('ttbar_norm_mean_%s'%label,'ttbar_norm_mean_%s'%label,self.workspace4fit_.var('meanval_ttbar').getVal())
	mean_VV		= RooRealVar('VV_norm_mean_%s'%label,'VV_norm_mean_%s'%label,self.workspace4fit_.var('meanval_VV').getVal())
	ttbar_uncert	= {'el' : 0.2, 'mu' : 0.2}
	VV_uncert	= {'el' : 1., 'mu' : 1.}
	ttbar_xs_uncertainty	= ttbar_uncert[self.channel]
	VV_xs_uncertainty	= VV_uncert[self.channel]
	sigma_ttbar		= RooRealVar('ttbar_norm_sigma_%s'%label,'ttbar_norm_sigma_%s'%label,ttbar_xs_uncertainty * mean_ttbar.getVal())
	sigma_VV		= RooRealVar('VV_norm_sigma_%s'%label,'VV_norm_sigma_%s'%label,VV_xs_uncertainty * mean_VV.getVal())
	ttbar_gauss		= RooGaussian("ttbar_gauss_uncert_%s"%label,"ttbar_gauss_uncert_%s"%label,self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)),mean_ttbar,sigma_ttbar)
	VV_gauss		= RooGaussian("VV_gauss_uncert_%s"%label,"VV_gauss_uncert_%s"%label,self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)),mean_VV,sigma_VV)
	
	###
	###
	#mlvj model
	###
	###
        ## get the minor component shapes in the sb low
        model_VV_backgrounds   	= self.workspace4fit_.pdf("model_pdf_VV_xww_sb_lo_%s_mlvj"%self.channel)
        self.fix_Pdf(model_VV_backgrounds,RooArgSet(rrv_mass_lvj))
        model_TTbar_backgrounds	= self.workspace4fit_.pdf("model_pdf_TTbar_xww_sb_lo_%s_mlvj"%self.channel);
        self.fix_Pdf(model_TTbar_backgrounds,RooArgSet(rrv_mass_lvj))
        model_STop_backgrounds 	= self.get_STop_mlvj_Model("_sb_lo");
        self.fix_Pdf(model_STop_backgrounds,RooArgSet(rrv_mass_lvj))
        number_STop_sb_lo_mlvj 	= self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel))
         ### Make the Pdf for the WJets
        model_pdf_WJets 	= self.make_Pdf("%s_sb_lo_from_fitting"%(label), mlvj_model,"_mlvj");


	###make norm dependent on mjpruned fit result
	signalIntWJets		= model_WJets_mj.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), 'total_signal_region')
	number_WJets_sb_lo 	= RooFormulaVar("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),"rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),'(1-@0)*@1',RooArgList(signalIntWJets,self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel))))
        model_WJets 		= RooExtendPdf("model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),"model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),model_pdf_WJets,number_WJets_sb_lo);
        getattr(self.workspace4fit_,'import')(number_WJets_sb_lo)

        signalIntTTbar		= model_TTbar_mj.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), 'total_signal_region')
        number_TTbar_sb_lo	= RooFormulaVar("rrv_number_TTbar_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),"rrv_number_TTbar_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),'(1-@0)*@1',RooArgList(signalIntTTbar,self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%(self.channel))))
        model_TTbar 		= RooExtendPdf("model_TTbar_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),"model_TTbar_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),model_TTbar_backgrounds,number_TTbar_sb_lo);
        getattr(self.workspace4fit_,'import')(number_TTbar_sb_lo)
        
        signalIntVV		= model_VV_mj.createIntegral(RooArgSet(rrv_mass_j), RooArgSet(rrv_mass_j), 'total_signal_region')
        number_VV_sb_lo		= RooFormulaVar("rrv_number_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),"rrv_number_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),'(1-@0)*@1',RooArgList(signalIntVV,self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%(self.channel))))
        model_VV 		= RooExtendPdf("model_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),"model_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel),model_VV_backgrounds,number_VV_sb_lo);
        getattr(self.workspace4fit_,'import')(number_VV_sb_lo)     
        
        ## Add the other bkg component fixed to the total model
        model_mlvj = RooAddPdf("model_data%s%s_mlvj"%(label,mlvj_region),"model_data%s%s_mlvj"%(label,mlvj_region),RooArgList(model_WJets,model_VV, model_TTbar, model_STop_backgrounds));

	#set which normalizations float in addition to W+Jets
	self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kFALSE)
	self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kFALSE)
	
	ttbar_constraint= RooConstraintSum('ttbar_constraint','ttbar_constraint',RooArgSet(ttbar_gauss),RooArgSet(self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel))))
	VV_constraint	= RooConstraintSum('VV_constraint','VV_constraint',RooArgSet(VV_gauss),RooArgSet(self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel))))	
	NLL_mj		= RooNLLVar('NLL_mj','NLL_mj', model_mj, rdataset_data_mj, RooFit.Extended(kTRUE), RooFit.Range("sb_lo_to_sb_hi"), RooFit.Minimizer('Minuit2'))
	NLL_mlvj	= RooNLLVar('NLL_mlvj','NLL_mlvj', model_mlvj, rdataset_data_mlvj, RooFit.Extended(kTRUE), RooFit.Minimizer('Minuit2'))
	NLL		= RooAddition('NLL','NLL',RooArgSet(NLL_mj,NLL_mlvj,ttbar_constraint,VV_constraint))

	m	= RooMinuit(NLL)
	m.migrad()
	m.setStrategy(2)
	m.migrad()
	m.hesse()
	rfresult2 = m.save()
	rfresult2.Print()


	###
	###
	#plot mj
	###
	###
        ## Total numver of event 
        rrv_number_data_mj = RooRealVar("rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),"rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),
                                         self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal());


        rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()));
        getattr(self.workspace4fit_,"import")(rrv_number_data_mj);
        
	## make the final plot
	mplot_mj = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(int(rrv_mass_j.getBins()/self.narrow_factor)));
	rdataset_data_mj.plotOn(mplot_mj, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

	## plot solid style 

	model_mj.plotOn(mplot_mj,RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"), RooFit.VLines());  
	model_mj.plotOn(mplot_mj,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"), RooFit.VLines());
	model_mj.plotOn(mplot_mj,RooFit.Name("STop"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"), RooFit.VLines());
	model_mj.plotOn(mplot_mj,RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"), RooFit.VLines());
	
	### solid line
	model_mj.plotOn( mplot_mj,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
	model_mj.plotOn( mplot_mj,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
	model_mj.plotOn( mplot_mj,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) , RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
	model_mj.plotOn( mplot_mj,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Normalization(rrv_number_data_mj.getVal(),RooAbsReal.NumEvent), RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
	
	rdataset_data_mj.plotOn(mplot_mj, RooFit.Name("data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

	### draw the error band using the sum of all the entries component MC + fit         
	draw_error_band(rdataset_data_mj, model_mj, rrv_number_data_mj,rfresult2,mplot_mj,self.color_palet["Uncertainty"],"F");
	rdataset_data_mj.plotOn(mplot_mj, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

	### Get the pull and plot it 
	mplot_pull_mj=self.get_pull(rrv_mass_j,mplot_mj);
	### legend of the plot
	self.leg = self.legend4Plot(mplot_mj,0,1,0.3,0.,0.,0.,0,0,1);
	mplot_mj.addObject(self.leg);
	mplot_mj.GetYaxis().SetRangeUser(1e-2,mplot_mj.GetMaximum()*1.1);

	datahist = rdataset_data_mj.binnedClone( rdataset_data_mj.GetName()+"_binnedClone",rdataset_data_mj.GetName()+"_binnedClone" )
	parameters_list = model_mj.getParameters(rdataset_data_mj);
	self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot_mj, mplot_pull_mj,1,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), "m_j_sideband%s_simfit"%(label),"",1,0,1,0,1)

	###
	###
	#plot mlvj
	###
	###

	rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),"rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),
                                                 self.workspace4fit_.function("rrv_number_TTbar_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.function("rrv_number_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.function("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getVal() );
	
        rrv_number_data_sb_lo_mlvj.setError( TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()*(1-signalIntWJets.getVal())*
                                                        self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()*(1-signalIntWJets.getVal())+
                                                        self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%(self.channel)).getError()*(1-signalIntTTbar.getVal())*
                                                        self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%(self.channel)).getError()*(1-signalIntTTbar.getVal())+
                                                        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%(self.channel)).getError()*(1-signalIntVV.getVal())*
                                                        self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%(self.channel)).getError()*(1-signalIntVV.getVal())));


        getattr(self.workspace4fit_,"import")(rrv_number_data_sb_lo_mlvj)
	mplot_mlvj = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));

	#rdataset_data_mlvj.plotOn( mplot_mlvj , RooFit.Invisible(), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0) );

	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_from_fitting_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_from_fitting_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent)) ;
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model_TTbar_xww_sb_lo_from_fitting_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent)) ;
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model_VV_xww_sb_lo_from_fitting_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent));
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent));

	model_STop_backgrounds.plotOn(mplot_mlvj, RooFit.LineColor(kCyan), RooFit.Normalization(number_STop_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent))



	#solid line
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_from_fitting_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_from_fitting_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent)) ;
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model_TTbar_xww_sb_lo_from_fitting_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_from_fitting_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent)) ;
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model_VV_xww_sb_lo_from_fitting_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent));
	model_mlvj.plotOn(mplot_mlvj, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent) );

	### draw the error band 
	draw_error_band(rdataset_data_mlvj, model_mlvj, self.workspace4fit_.var("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel)) ,rfresult2,mplot_mlvj,self.color_palet["Uncertainty"],"F");
	model_mlvj.plotOn( mplot_mlvj , RooFit.Invisible(), RooFit.Normalization(rrv_number_data_sb_lo_mlvj.getVal(),RooAbsReal.NumEvent));
	self.getData_PoissonInterval(rdataset_data_mlvj,mplot_mlvj);

	mplot_mlvj.GetYaxis().SetRangeUser(1e-2,mplot_mlvj.GetMaximum()*1.2);  

	### Add the legend to the plot 
	self.leg=self.legend4Plot(mplot_mlvj,0,1,0.25,0.,0.1,0.,0);
	mplot_mlvj.addObject(self.leg)

	### calculate the chi2
	self.nPar_float_in_fitTo = rfresult2.floatParsFinal().getSize();
	nBinX = mplot_mlvj.GetNbinsX();
	ndof  = nBinX-self.nPar_float_in_fitTo;

	### get the pull plot and store the canvas
	mplot_mlvj_pull = self.get_pull(rrv_mass_lvj,mplot_mlvj);
	parameters_list = model_mlvj.getParameters(rdataset_data_mlvj);
	datahist = rdataset_data_mlvj.binnedClone( rdataset_data_mlvj.GetName()+"_binnedClone",rdataset_data_mlvj.GetName()+"_binnedClone" )
	
	mplot_mlvj.GetXaxis().SetTitle(rrv_mass_lvj.GetTitle() + " (GeV)")
	mplot_mlvj.GetYaxis().SetRangeUser(2e-3,2e5)
	#@#self.draw_canvas_with_pull( rrv_mass_lvj,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), "m_lvj_sb_lo%s"%(label),"",1,1)
	self.draw_canvas_with_pull( rrv_mass_lvj,datahist,mplot_mlvj, mplot_mlvj_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), "m_lvj_sb_lo_simfit%s"%(label),"",1,1,0,1,1)
	    
	rfresult2.Print()
	exit(0)

    #### make the mj sideband fit on data to get the Wjets normaliztion 
    def fit_WJetsNormalization_in_Mj_signal_region(self,label,massscale=""): 

	if options.simfit:
	        self.simultaneous_fit_mj_and_mlvj_sb(label,'sb_lo',self.MODEL_4_mlvj)

        print "############### Fit mj Normalization: ",label," ",massscale," ##################"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        ## get real data in mj distribution --> mass up and down have only an effect on Wjets shape -> effect on the normalization -> evaluated in the MC and fit data
        rdataset_data_mj=self.workspace4fit_.data("rdataset_data_xww_%s_mj"%(self.channel))

        ### Fix TTbar, VV and STop
        model_TTbar = self.get_TTbar_mj_Model("_TTbar_xww"+massscale);
        model_STop  = self.get_STop_mj_Model("_STop_xww"+massscale);
        model_VV    = self.get_VV_mj_Model("_VV_xww"+massscale);

	model_TTbar.Print()
	
	#@#errorTTbar 
        ## only two parameters are fix, offset and width, while the exp is floating , otherwise if shape different User1 or ErfExp everything is flaoting
        model_WJets = self.get_WJets_mj_Model(label);
	model_WJets.Print()


        ## Total Pdf and fit 
        model_data = RooAddPdf("model_data_xww%s_%s_mj"%(massscale,self.channel),"model_data_xww%s_%s_mj"%(massscale,self.channel),RooArgList(model_WJets,model_VV,model_TTbar,model_STop))
	#@#let ttbar float with gaussian constraint
	#@#save value of rrv_number_XX from MC-fit so fit with alternative function has same mean+constraint 
	if '_WJets0_' in label:
		meanval_ttbar	= RooRealVar('meanval_ttbar','meanval_ttbar',self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getVal())
		meanval_VV	= RooRealVar('meanval_VV','meanval_VV',self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getVal())
		meanval_STop	= RooRealVar('meanval_STop','meanval_STop',self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getVal())
		ttbar_error	= RooRealVar('ttbar_error','ttbar_error',self.workspace4fit_.var('rrv_number_TTbar_xww%s_%s_mj'%(massscale,self.channel)).getError())
		VV_error	= RooRealVar('VV_error','VV_error',self.workspace4fit_.var('rrv_number_VV_xww%s_%s_mj'%(massscale,self.channel)).getError())
		STop_error	= RooRealVar('STop_error','STop_error',self.workspace4fit_.var('rrv_number_STop_xww%s_%s_mj'%(massscale,self.channel)).getError())
		meanval_ttbar.setError(ttbar_error.getVal())
		meanval_VV.setError(VV_error.getVal())
		meanval_STop.setError(STop_error.getVal())
		getattr(self.workspace4fit_,'import')(meanval_ttbar)
		getattr(self.workspace4fit_,'import')(ttbar_error)
		getattr(self.workspace4fit_,'import')(meanval_VV)
		getattr(self.workspace4fit_,'import')(VV_error)
		getattr(self.workspace4fit_,'import')(meanval_STop)
		getattr(self.workspace4fit_,'import')(STop_error)
	#@#set value and error of rrv_number_ttbar to prefit for WJets01-fit
	elif '_WJets01_' in label:
		self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var('meanval_ttbar').getVal())
		self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setError(self.workspace4fit_.var('ttbar_error').getVal())
		self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var('meanval_VV').getVal())
		self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setError(self.workspace4fit_.var('VV_error').getVal())

	mean_ttbar	= RooRealVar('ttbar_norm_mean_%s'%label,'ttbar_norm_mean_%s'%label,self.workspace4fit_.var('meanval_ttbar').getVal())
	mean_VV		= RooRealVar('VV_norm_mean_%s'%label,'VV_norm_mean_%s'%label,self.workspace4fit_.var('meanval_VV').getVal())
	#mean_STop	= RooRealVar('STop_norm_mean_%s'%label,'STop_norm_mean_%s'%label,self.workspace4fit_.var('meanval_STop').getVal())
	ttbar_uncert	= {'el' : 0.20, 'mu' : 0.20}
	VV_uncert	= {'el' : 1., 'mu' : 1.}
	#STop_uncert	= {'el' : 0.03, 'mu' : 0.05}
	ttbar_xs_uncertainty	= ttbar_uncert[self.channel]
	VV_xs_uncertainty	= VV_uncert[self.channel]
	#STop_xs_uncertainty	= STop_uncert[self.channel]
	sigma_ttbar		= RooRealVar('ttbar_norm_sigma_%s'%label,'ttbar_norm_sigma_%s'%label,ttbar_xs_uncertainty * mean_ttbar.getVal())
	sigma_VV		= RooRealVar('VV_norm_sigma_%s'%label,'VV_norm_sigma_%s'%label,VV_xs_uncertainty * mean_VV.getVal())
	#sigma_STop		= RooRealVar('STop_norm_sigma_%s'%label,'STop_norm_sigma_%s'%label,STop_xs_uncertainty * mean_STop.getVal())
	ttbar_gauss		= RooGaussian("ttbar_gauss_uncert_%s"%label,"ttbar_gauss_uncert_%s"%label,self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)),mean_ttbar,sigma_ttbar)
	VV_gauss		= RooGaussian("VV_gauss_uncert_%s"%label,"VV_gauss_uncert_%s"%label,self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)),mean_VV,sigma_VV)
	#STop_gauss		= RooGaussian("STop_gauss_uncert_%s"%label,"STop_gauss_uncert_%s"%label,self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)),mean_STop,sigma_STop)
	
	#ttbar and VV constrained
	model_data4fit	= RooProdPdf("model_data_xww%s_%s_mj_4fit"%(massscale,self.channel),"model_data_xww%s_%s_mj_4fit"%(massscale,self.channel),RooArgList(model_data,ttbar_gauss,VV_gauss))
	model_data4fit.Print()
	#raw_input('...')

	#set which normalizations float in addition to W+Jets
	self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kFALSE)
	self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kFALSE)
	#self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kFALSE)

	#@# save prefit values
	if 'WJets0_' in label:
		N_STop_pre	= self.workspace4fit_.var("rrv_number_STop_xww_%s_mj"%self.channel).clone("rrv_number_STop_xww_%s_mj_prefit"%self.channel)
		N_VV_pre	= self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%self.channel).clone("rrv_number_VV_xww_%s_mj_prefit"%self.channel)
		N_TTbar_pre	= self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%self.channel).clone("rrv_number_TTbar_xww_%s_mj_prefit"%self.channel)
		N_WJets_pre	= self.workspace4fit_.var("rrv_number_WJets0_xww_%s_mj"%self.channel).clone("rrv_number_WJets0_xww_%s_mj_prefit"%self.channel)
		getattr(self.workspace4fit_,'import')(N_STop_pre)
		getattr(self.workspace4fit_,'import')(N_VV_pre)
		getattr(self.workspace4fit_,'import')(N_TTbar_pre)
		getattr(self.workspace4fit_,'import')(N_WJets_pre)
		getattr(self.workspace4limit_,'import')(N_VV_pre)


	#@#test
	#self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getVal()*2)
	#@#
	###
	###use mj signal region as well for fit
        rfresult = model_data4fit.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo_to_sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2) );
       	rfresult = model_data4fit.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo_to_sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2), RooFit.Minimizer("Minuit2") );
	##4DEBUG
	#rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2) );
        #rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2), RooFit.Minimizer("Minuit2") );
	##
	fitresultsfinal.append(rfresult)

	self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kTRUE)
	self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kTRUE)

	#model_data.SetName("model_data_xww%s_%s_mj"%(massscale,self.channel))
        getattr(self.workspace4fit_,"import")(model_data)
	getattr(self.workspace4fit_,'import')(rfresult)
	getattr(self.workspace4limit_,'import')(rfresult)
	if 'WJets0_' in label:
		tmpfile=TFile.Open('cards_%s_%s_900_3500/workspace_mj_%s.root'%(self.channel,self.wtagger_label,self.channel),'recreate')
		self.workspace4fit_.Write()
		tmpfile.Close()
	

        rfresult.Print();
	print 'realtive error of ttbar-normalization: ' + str(self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError() / self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getVal())


        rfresult.covarianceMatrix().Print();

        ## Total numver of event 
        rrv_number_data_mj = RooRealVar("rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),"rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),
                                         self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal());

        rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()));
        getattr(self.workspace4fit_,"import")(rrv_number_data_mj);
        
        ## if fit on Wjets default with the default shape
        if TString(label).Contains("_WJets0"):

            ## make the final plot
            mplot = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(int(rrv_mass_j.getBins()/self.narrow_factor)));
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ## plot solid style 
            model_data.plotOn(mplot,RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
              
	    model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));

            model_data.plotOn(mplot,RooFit.Name("STop"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));

            model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
 	    
    
	    ### solid line
	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
	
	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));

	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
		   
	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));

 
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### draw the error band using the sum of all the entries component MC + fit         
            draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );
	    print rdataset_data_mj.sumEntries()
	    #raw_input(label)

            ### Get the pull and plot it 
            mplot_pull=self.get_pull(rrv_mass_j,mplot);

            ### signal window zone with vertical lines
	    signal_mid	= 85
            lowerLine 	= TLine(65,0.,65,mplot.GetMaximum()*1); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(7);
            middleLine 	= TLine(signal_mid,0.,signal_mid,mplot.GetMaximum()*1); middleLine.SetLineWidth(2); middleLine.SetLineColor(kBlack); middleLine.SetLineStyle(7);
            upperLine 	= TLine(105,0.,105,mplot.GetMaximum()*1); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(7);
            mplot.addObject(lowerLine);
	    mplot.addObject(middleLine);
            mplot.addObject(upperLine);

	    pt = ROOT.TPaveText(0.3592965,0.02847153,0.5728643,0.1008991,"NDC")
	    pt.SetTextFont(42)
	    pt.SetTextSize(0.04995005)
	    #pt.SetTextAlign(12)
	    pt.SetFillColor(self.color_palet["WJets"])
	    pt.SetBorderSize(0)
	    text = pt.AddText("#leftarrow signal region #rightarrow")
	    text.SetTextFont(62)

	    pt1 = ROOT.TPaveText(0.3592965,0.1183816,0.4535176,0.1908092,"NDC")
	    pt1.SetTextFont(42)
	    pt1.SetTextSize(0.037)
	    #pt.SetTextAlign(12)
	    pt1.SetFillColor(self.color_palet["WJets"])
	    #pt1.SetFillStyle(0)
	    pt1.SetBorderSize(0)
	    text = pt1.AddText("WW")
	    text.SetTextFont(62)
	    text = pt1.AddText("category")
	    text.SetTextFont(62)

	    pt2 = ROOT.TPaveText(0.4723618,0.1183816,0.5728643,0.1908092,"NDC")
	    pt2.SetTextFont(42)
	    pt2.SetTextSize(0.037)
	    #pt.SetTextAlign(12)
	    pt2.SetFillColor(self.color_palet["WJets"])
	    pt2.SetBorderSize(0)
	    #pt2.SetFillStyle(0)
	    text = pt2.AddText("WZ")
	    text.SetTextFont(62)
	    text = pt2.AddText("category")
	    text.SetTextFont(62)

	    mplot.addObject(pt)
	    mplot.addObject(pt1)
	    mplot.addObject(pt2)
	    	    
            ### legend of the plot
            self.leg = self.legend4Plot(mplot,0,1,0.3,0.,0.,0.,0,0,1);
            #self.leg = self.legend4Plot(mplot,0,1,-0.10,-0.01,0.10,0.01);

            mplot.addObject(self.leg);
            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.1);

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            #print mplot.chiSquare();
            #print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            datahist = rdataset_data_mj.binnedClone( rdataset_data_mj.GetName()+"_binnedClone",rdataset_data_mj.GetName()+"_binnedClone" )

            parameters_list = model_data.getParameters(rdataset_data_mj);
            self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_j_fitting/"%(self.plotsDir), "m_j_sideband%s"%(label),"",1,0,1)

            ### call the function for getting the normalizatio in signal region for data, TTbar, STop, VV and W+jets = label -> store in a output txt file
	    if 'WJets0_' in label:
		N_STop_post	= self.workspace4fit_.var("rrv_number_STop_xww_%s_mj"%self.channel).clone("rrv_number_STop_xww_%s_mj_postfit"%self.channel)
		N_VV_post	= self.workspace4fit_.var("rrv_number_VV_xww_%s_mj"%self.channel).clone("rrv_number_VV_xww_%s_mj_postfit"%self.channel)
		N_TTbar_post	= self.workspace4fit_.var("rrv_number_TTbar_xww_%s_mj"%self.channel).clone("rrv_number_TTbar_xww_%s_mj_postfit"%self.channel)
		N_WJets_post	= self.workspace4fit_.var("rrv_number_WJets0_xww_%s_mj"%self.channel).clone("rrv_number_WJets0_xww_%s_mj_postfit"%self.channel)
		getattr(self.workspace4fit_,'import')(N_STop_post)
		getattr(self.workspace4fit_,'import')(N_VV_post)
		getattr(self.workspace4fit_,'import')(N_TTbar_post)
		getattr(self.workspace4fit_,'import')(N_WJets_post)
           	self.get_mj_normalization_insignalregion("_data_xww",1);
	        self.get_mj_normalization_insignalregion("_TTbar_xww",1);
        	self.get_mj_normalization_insignalregion("_STop_xww",1);
            	self.get_mj_normalization_insignalregion("_VV_xww",1);
	        self.get_mj_normalization_insignalregion(label,1);
	    elif 'WJets01' in label:
            	self.get_mj_normalization_insignalregion("_data_xww",-1);
	        self.get_mj_normalization_insignalregion("_TTbar_xww",-1);
        	self.get_mj_normalization_insignalregion("_STop_xww",-1);
            	self.get_mj_normalization_insignalregion("_VV_xww",-1);
	        self.get_mj_normalization_insignalregion(label,-1);

	#### to calculate the WJets's normalization and error in M_J signal_region. The error must contain the shape error: model_WJets have new parameters fitting data

        fullInt   	= model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        signalInt 	= model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signal_region"));
 	signalIntTot 	= model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("total_signal_region"));
        fullInt_val	= fullInt.getVal()
        signalInt_val 	= signalInt.getVal()/fullInt_val
	signalIntTot_val= signalIntTot.getVal()/fullInt_val

	
        ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
        rrv_number_WJets_in_mj_signal_region_from_fitting = RooRealVar("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),"rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal()*signalInt_val);
	#@#WJETSNORMSB
	rrv_number_WJets_in_mj_sb_region_from_fitting = RooRealVar("rrv_number%s_in_mj_sb_region_from_fitting_%s"%(label,self.channel),"rrv_number%s_in_mj_sb_region_from_fitting_%s"%(label,self.channel),self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal()*(1-signalIntTot_val)) 
	#@#

        #### Error on the normalization --> from a dedicated function taking into account shape uncertainty
	print rrv_number_WJets_in_mj_signal_region_from_fitting.getError()
        rrv_number_WJets_in_mj_signal_region_from_fitting.setError( Calc_error_extendPdf(rdataset_data_mj, model_WJets, rfresult,"signal_region") );

	if 'WJets0_' in label:
		self.workspace4fit_.var('rrv_number_WJets0_xww_%s_mj_postfit_signal_region'%self.channel).setError(Calc_error_extendPdf(rdataset_data_mj, model_WJets, rfresult,"signal_region") )
	print rrv_number_WJets_in_mj_signal_region_from_fitting.getError()
	print rrv_number_WJets_in_mj_signal_region_from_fitting.getError()/rrv_number_WJets_in_mj_signal_region_from_fitting.getVal()
        print "########## error on the normalization due to shape + norm = %s"%(rrv_number_WJets_in_mj_signal_region_from_fitting.getError());
        getattr(self.workspace4fit_,"import")(rrv_number_WJets_in_mj_signal_region_from_fitting);
	getattr(self.workspace4fit_,"import")(rrv_number_WJets_in_mj_sb_region_from_fitting);
        rrv_number_WJets_in_mj_signal_region_from_fitting.Print();
	
	#@# save postfit values
	if 'WJets01_' in label:
		self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var("meanval_ttbar").getVal())
		self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setError(self.workspace4fit_.var("meanval_ttbar").getError())
		self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var("meanval_VV").getVal())
		self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).setVal(self.workspace4fit_.var("meanval_VV").getError())
	
	if 'WJets0_' in label:
	  #@#prepare functions for simultaneous fit in combine
	  model_WJets_pdf	= self.workspace4fit_.pdf("model_pdf_WJets0_xww_%s_mj"%self.channel)
	  model_TTbar_pdf	= self.workspace4fit_.pdf("model_pdf_TTbar_xww_%s_mj"%self.channel)
	  model_VV_pdf		= self.workspace4fit_.pdf("model_pdf_VV_xww_%s_mj"%self.channel)
	  model_STop_pdf	= self.workspace4fit_.pdf('model_pdf_STop_xww_%s_mj'%self.channel)
	  bkgs			= ['WJets','TTbar','VV','STop']
	  mj_bkg_pdfs	= {'WJets':model_WJets_pdf, 'TTbar':model_TTbar_pdf, 'VV':model_VV_pdf, 'STop':model_STop_pdf}
	  for bkg in bkgs:
	    getattr(self.workspace4limit_,'import')(mj_bkg_pdfs[bkg].clone('m_j_%s_%s'%(bkg,self.channel)))
	  for region in ['W'+self.wtagger_label[2]+'_sig','sb_lo','sb_hi']:
	    for bkg in bkgs:
	      custom_mj		= RooCustomizer(mj_bkg_pdfs[bkg],'%s_mj_%s'%(bkg,region))
	      custom_mj.replaceArg(self.workspace4fit_.var('rrv_mass_j'),self.workspace4fit_.var('mj_%s'%region))
	      m_pruned_pdf	= custom_mj.build()
	      m_pruned_pdf.Print()
	      getattr(self.workspace4limit_,'import')(m_pruned_pdf.clone('m_j_%s_%s_pdf_%s'%(bkg,region,self.channel)),RooFit.RecycleConflictNodes())
	self.workspace4limit_.Print()




    def mj_prefit_plot(self,label,massscale=""): 

        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        ## get real data in mj distribution --> mass up and down have only an effect on Wjets shape -> effect on the normalization -> evaluated in the MC and fit data
        rdataset_data_mj=self.workspace4fit_.data("rdataset_data_xww_%s_mj"%(self.channel))

        ### Fix TTbar, VV and STop
        model_TTbar = self.get_TTbar_mj_Model("_TTbar_xww"+massscale);
        model_STop  = self.get_STop_mj_Model("_STop_xww"+massscale);
        model_VV    = self.get_VV_mj_Model("_VV_xww"+massscale);

	model_TTbar.Print()
	

        ## only two parameters are fix, offset and width, while the exp is floating , otherwise if shape different User1 or ErfExp everything is flaoting
        model_WJets = self.get_WJets_mj_Model(label);

        ## Total Pdf and fit only in sideband 
        model_data = RooAddPdf("model_data_xww%s_%s_mj"%(massscale,self.channel),"model_data_xww%s_%s_mj"%(massscale,self.channel),RooArgList(model_WJets,model_VV,model_TTbar,model_STop))

	self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kFALSE)



	self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).setConstant(kTRUE)
        getattr(self.workspace4fit_,"import")(model_data)


        ## Total numver of event 
        rrv_number_data_mj = RooRealVar("rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),"rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),
                                         self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal());

        rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()));
        getattr(self.workspace4fit_,"import")(rrv_number_data_mj);
        
        ## if fit on Wjets default with the default shape
        if TString(label).Contains("_WJets0"):

            ## make the final plot
            mplot = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(int(rrv_mass_j.getBins()/self.narrow_factor)));
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ## plot solid style 
            model_data.plotOn(mplot,RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
              
	    model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));

            model_data.plotOn(mplot,RooFit.Name("STop"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));

            model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo_to_sb_hi"));
 	    
	    	    
	    ### solid line
	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));

	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));

	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));
	   
	    model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo_to_sb_hi"), RooFit.Range("sb_lo_to_sb_hi"));

 
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### draw the error band using the sum of all the entries component MC + fit         
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### Get the pull and plot it 

            ### signal window zone with vertical lines
    	    signal_mid	= 85
    	    lowerLine 	= TLine(65,0.,65,mplot.GetMaximum()); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(7);
    	    middleLine 	= TLine(signal_mid,0.,signal_mid,mplot.GetMaximum()); middleLine.SetLineWidth(2); middleLine.SetLineColor(kBlack); middleLine.SetLineStyle(7);
    	    upperLine 	= TLine(105,0.,105,mplot.GetMaximum()); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(7);
    	    mplot.addObject(lowerLine);
    	    mplot.addObject(middleLine);
    	    mplot.addObject(upperLine);

	    	    
            ### legend of the plot
            self.leg = self.legend4Plot(mplot,0,1,0.3,0.,0.,0.,0,0,1);
            #self.leg = self.legend4Plot(mplot,0,1,-0.10,-0.01,0.10,0.01);

            mplot.addObject(self.leg);
            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.1);

            datahist = rdataset_data_mj.binnedClone( rdataset_data_mj.GetName()+"_binnedClone",rdataset_data_mj.GetName()+"_binnedClone" )

            parameters_list = model_data.getParameters(rdataset_data_mj);
            mplot_pull=self.get_pull(rrv_mass_j,mplot);


            self.draw_canvas_with_pull( rrv_mass_j,rdataset_data_mj,mplot,mplot_pull,0,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), "m_j_prefit_%s"%(label),'',1,0,1)


    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_mlvj_model_single_MC(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",in_file_name," ",label,"  ",mlvj_model,"  ",in_range," ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.channel+"_mlvj");
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,mlvj_model,"_mlvj",constrainslist,ismc);


        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
	fitresultsmlvj.append(rfresult)

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,6,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        nPar = rfresult.floatParsFinal().getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-nPar;
        #print mplot.chiSquare();
        #print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(nPar ,mplot.chiSquare(nPar)*ndof, ndof );
        datahist = rdataset.binnedClone( rdataset.GetName()+"_binnedClone",rdataset.GetName()+"_binnedClone" )
	rdataset.Print()
	
        ## get the pull 
        mplot_pull      = self.get_pull(rrv_mass_lvj,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
        
        self.draw_canvas_with_pull( rrv_mass_lvj, datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, logy);
	
	#@# make missing plot
	if 'TTBAR' in in_file_name and 'sb' in in_range:
		self.draw_canvas_with_pull( rrv_mass_lvj, datahist,mplot, mplot_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, 1)

         
        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_mlvj"%(label,in_range,self.channel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.channel+"_mlvj");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.channel+"_"+self.wtagger_label+"_mlvj_13TeV",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
            print "##################### diagonalize ";
            model_pdf_deco = Deco.diagonalize(model_pdf); ## diagonalize            
            print "##################### workspace for decorrelation ";
            wsfit_tmp.Print("v");
            print "##################### original  parameters ";
            model_pdf.getParameters(rdataset).Print("v");
            print "##################### original  decorrelated parameters ";
            model_pdf_deco.getParameters(rdataset).Print("v");
            print "##################### original  pdf ";
            model_pdf.Print();
            print "##################### decorrelated pdf ";
            model_pdf_deco.Print();


            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4fit_,"import")(model_pdf_deco);

            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame( RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));
            
            if label=="_TTbar_xww" and in_range=="_signal_region":
                
                rdataset.plotOn(mplot_deco, RooFit.Name("Powheg Sample"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name("TTbar_Powheg"),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset = RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## draw the error band with the area
                self.workspace4fit_.var("rrv_number_TTbar_xww_signal_region_%s_mlvj"%(self.channel)).Print();
            else:
                rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.leg = self.legend4Plot(mplot_deco,0); ## add the legend                
            mplot_deco.addObject(self.leg);

            self.draw_canvas( mplot_deco, "%s/other/"%(self.plotsDir), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print();



    ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

        print "############### Fit mj single MC sample",in_file_name," ",label,"  ",in_model_name," ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.channel+"_mj");
        rdataset_mj.Print();

        ## make the extended model
	model = self.make_Model(label,in_model_name);
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.Extended(kTRUE) );
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
	getattr(self.workspace4fit_,'import')(rfresult)
	fitresultsmj.append(rfresult)
	#if 'WJets0_' in label:
	#	raw_input(label)
	
        ## Plot the result
        mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins()/self.narrow_factor)) );
        rdataset_mj.plotOn( mplot, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,6,"L");
        ## re-draw the dataset
        rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## draw the function
        model.plotOn( mplot );# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = self.get_pull(rrv_mass_j, mplot); 
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        nPar = rfresult.floatParsFinal().getSize();
        #nBinX = mplot.GetNbinsX();
	nBinX = 22
        ndof  = nBinX-nPar;
        #print "#################### nPar=%s, nBinX=%s , chiSquare=%s/%s"%(nPar,nBinX,mplot.chiSquare('model%s_%s_mj_Norm[rrv_mass_j]'%(label,self.channel),'h_rdataset4fit%s_%s_mj'%(label,self.channel),nPar)*ndof,ndof);
	print self.calculate_chi2(rdataset_mj,rrv_mass_j, mplot, ndof,1)
	mplot.Print()

        datahist = rdataset_mj.binnedClone( rdataset_mj.GetName()+"_binnedClone",rdataset_mj.GetName()+"_binnedClone" )
        parameters_list = model.getParameters(rdataset_mj);
	#@here@#
	#@# add name to plot
	pt1 = ROOT.TPaveText(0.5680905,0.7,0.75,0.75,"NDC")
	pt1.SetTextFont(42)
	pt1.SetTextSize(0.1)
	#pt.SetTextAlign(12)
	pt1.SetFillColor(0)
	pt1.SetFillStyle(0)
	pt1.SetBorderSize(0)
	if '_STop' in label:
		text = pt1.AddText('Single Top')
	if '_TTbar' in label:
		#text = pt1.AddText('t#bar{t}')
		legx = TLegend(0.7,0.6,0.9,0.85)
		legx.SetBorderSize(0)
		legx.SetFillColor(0)
		legx.AddEntry(mplot.getObject(3),'MC data','p')
		legx.AddEntry(mplot.getObject(4),'fit','l')
		legx.AddEntry(mplot.getObject(1),'error band','l')
		legx.SetHeader('t#bar{t}')
		#mplot.addObject(legx)
	if '_VV' in label:
		text = pt1.AddText('Diboson')
	if '_WJets' in label:
		text = pt1.AddText('W+Jets')
	#@#mplot.addObject(pt1)
        self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_j_fitting/"%(self.plotsDir), label+in_file_name, in_model_name)

	#@#more plots
	if 'TTbar' in label or 'VV' in label:
	    lowerLine = TLine(65,0.,65,mplot.GetMaximum()*0.8); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
            middleLine = TLine(85,0.,85,mplot.GetMaximum()*0.9); middleLine.SetLineWidth(2); middleLine.SetLineColor(kBlack); middleLine.SetLineStyle(9);
            upperLine = TLine(105,0.,105,mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
            mplot.addObject(lowerLine);
	    mplot.addObject(middleLine);
            mplot.addObject(upperLine);

	    pt = ROOT.TPaveText(0.3592965,0.02847153,0.5728643,0.1008991,"NDC")
	    pt.SetTextFont(42)
	    pt.SetTextSize(0.04995005)
	    #pt.SetTextAlign(12)
	    pt.SetFillColor(0)
	    pt.SetBorderSize(0)
	    text = pt.AddText("#leftarrow signal region #rightarrow")
	    text.SetTextFont(62)
	    mplot.addObject(pt);

	    pt1 = ROOT.TPaveText(0.3555276,0.1183816,0.4535176,0.1908092,"NDC")
	    pt1.SetTextFont(42)
	    pt1.SetTextSize(0.037)
	    #pt.SetTextAlign(12)
	    pt1.SetFillColor(0)
	    pt1.SetFillStyle(0)
	    pt1.SetBorderSize(0)
	    text = pt1.AddText("WW")
	    text.SetTextFont(62)
	    text = pt1.AddText("category")
	    text.SetTextFont(62)
	    mplot.addObject(pt1);

	    pt2 = ROOT.TPaveText(0.4723618,0.1183816,0.5678392,0.1908092,"NDC")
	    pt2.SetTextFont(42)
	    pt2.SetTextSize(0.037)
	    #pt.SetTextAlign(12)
	    pt2.SetFillColor(0)
	    pt2.SetBorderSize(0)
	    pt2.SetFillStyle(0)
	    text = pt2.AddText("WZ")
	    text.SetTextFont(62)
	    text = pt2.AddText("category")
	    text.SetTextFont(62)
	    mplot.addObject(pt2);
	    self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), label+in_file_name, in_model_name)
	#@#more plots
	if 'WJets' in label:
	    lowerLine = TLine(65,0.,65,mplot.GetMaximum()*0.8); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
            middleLine = TLine(85,0.,85,mplot.GetMaximum()*0.9); middleLine.SetLineWidth(2); middleLine.SetLineColor(kBlack); middleLine.SetLineStyle(9);
            upperLine = TLine(105,0.,105,mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
            mplot.addObject(lowerLine);
	    mplot.addObject(middleLine);
            mplot.addObject(upperLine);

	    pt = ROOT.TPaveText(0.3592965,0.02847153,0.5728643,0.1008991,"NDC")
	    pt.SetTextFont(42)
	    pt.SetTextSize(0.04995005)
	    #pt.SetTextAlign(12)
	    pt.SetFillColor(0)
	    pt.SetBorderSize(0)
	    text = pt.AddText("#leftarrow signal region #rightarrow")
	    text.SetTextFont(62)
	    mplot.addObject(pt);
	    self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir), label+in_file_name, in_model_name)

	self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();



	self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print();
        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();


	fullInt   	= model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        signalInt 	= model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signal_region"));
        fullInt_val	= fullInt.getVal()
        signalInt_val 	= signalInt.getVal()/fullInt_val

	
        ##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV 
	
        par=parameters_list.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                #param.Print();
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
            param=par.Next()

    def IsGoodEvent(self,tree,label):
       	#@#always True for new trees -> only HP
	keepEvent = True
	return keepEvent
       	#@#
              
    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass ):# to get the shape of m_lvj,jet_mass="jet_mass_pr"

	if not options.read_trees:
		self.get_mj_and_mlvj_dataset_from_file(in_file_name, label, jet_mass )
	else:

		print "################### get_mj_and_mlvj_dataset : ",in_file_name,"  ",label,"  ##################";

		fileIn_name = TString(options.inPath+"/"+self.file_Directory+in_file_name);
		fileIn = TFile(fileIn_name.Data());
		treeIn = fileIn.Get("tree");
		
		rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
		rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
		rrv_weight   = RooRealVar("rrv_weight","rrv_weight",-1000. ,10000000.)
		rrv_mass_j_sb_lo	= rrv_mass_j.clone("mj_sb_lo")
		rrv_mass_j_sb_hi	= rrv_mass_j.clone("mj_sb_hi")
		rrv_mass_j_sig		= rrv_mass_j.clone("mj_W%s_sig"%self.wtagger_label[2])
		rrv_mass_j_sb_lo.setRange(self.mj_sideband_lo_min,self.mj_sideband_lo_max)
		rrv_mass_j_sb_hi.setRange(self.mj_sideband_hi_min,self.mj_sideband_hi_max)
		rrv_mass_j_sig.setRange(self.mj_signal_min,self.mj_signal_max)

		self.BinWidth_mlvj=self.BinWidth_mlvj/self.narrow_factor;
		nbins_mlvj=int((rrv_mass_lvj.getMax()-rrv_mass_lvj.getMin())/self.BinWidth_mlvj);
		rrv_mass_lvj.setBins(nbins_mlvj);
		   
		##### dataset of m_j -> scaled and not scaled to lumi 
		rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
		rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
		##### dataset of m_lvj -> scaled and not scaled to lumi in different region
		rdataset_sb_lo_mlvj = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

		rdataset_signal_region_mlvj = RooDataSet("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

		rdataset_sb_hi_mlvj = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );
		

		rdataset4fit_sb_lo_mlvj = RooDataSet("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

		rdataset4fit_signal_region_mlvj = RooDataSet("rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

		rdataset4fit_sb_hi_mlvj = RooDataSet("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );
		
		if 'data' in label:
		  ### datasets for simultaneous fit with combine
		  dataset_mj_sb_lo	= RooDataSet('dataset_mj_sb_lo','dataset_mj_sb_lo',RooArgSet(rrv_mass_j_sb_lo, rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))
		  dataset_mj_sb_hi	= RooDataSet('dataset_mj_sb_hi','dataset_mj_sb_hi',RooArgSet(rrv_mass_j_sb_hi, rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))
		  dataset_mj_sig	= RooDataSet('dataset_mj_W%s_sig'%self.wtagger_label[2],'dataset_mj_W%s_sig'%self.wtagger_label[2],RooArgSet(rrv_mass_j_sig, rrv_mass_lvj, rrv_weight), RooFit.WeightVar(rrv_weight))


		### categorize the event in sideband and signal region --> combined dataset 

		data_category = RooCategory("data_category","data_category");
		data_category.defineType("sideband");
		data_category.defineType("signal_region");
		combData = RooDataSet("combData"+label+"_"+self.channel,"combData"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
		combData4fit = RooDataSet("combData4fit"+label+"_"+self.channel,"combData4fit"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
		
		print "###### N entries: ", treeIn.GetEntries()

		hnum_4region=TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
		hnum_2region=TH1D("hnum_2region"+label+"_"+self.channel,"hnum_2region"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total
		

		for i in range(treeIn.GetEntries()):

		    if i % 10000 == 0: print "iEntry: ",i
		    treeIn.GetEntry(i);
		    	
		    if i==0:
			if TString(label).Contains('data'):
				tmp_scale_to_lumi = 1
			else:
	      			tmp_scale_to_lumi = treeIn.totEventWeight / (treeIn.puweight*(treeIn.genweight/abs(treeIn.genweight)))


		    tmp_jet_mass=getattr(treeIn, jet_mass);

		    self.isGoodEvent = 0;   
		    if self.IsGoodEvent(treeIn,label) and treeIn.MWW> rrv_mass_lvj.getMin() and treeIn.MWW<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax():
		        self.isGoodEvent = 1;  
	 
	    	    if self.isGoodEvent == 1:
		        ### weigh MC events                                                        
			#@#change weight for new trees, lumi is already taken care of in totEventWeight
		        if TString(label).Contains('data'):
		            tmp_event_weight=1.;
		            tmp_event_weight4fit=1.; 
			else:
		       	    tmp_event_weight 		= treeIn.totEventWeight
			    tmp_event_weight4fit 	= tmp_event_weight/tmp_scale_to_lumi
			#@#	

                   
		        
		        rrv_mass_lvj.setVal(treeIn.MWW);
			#whole sideband region (still named sb_lo), added upper sideband
		        if (tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max) or (tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max):
		            rdataset_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
		            rdataset4fit_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

		            data_category.setLabel("sideband");
		            combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
		            combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
			    	
			#signal region, blind data if needed
			if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
			    rdataset_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
			    rdataset4fit_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
			    
			    data_category.setLabel("signal_region");
			    combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
			    combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
			    hnum_2region.Fill(1,tmp_event_weight);

			    hnum_2region.Fill(0,tmp_event_weight);
			#sideband hi only, for sim-fit
		        if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
		            rdataset_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
		            rdataset4fit_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
		    
		            
		        rrv_mass_j.setVal( tmp_jet_mass );
			#mj spectrum, cut out data in signal region if needed
			if 'data' in label:
			  rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight )
			  rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit )
			  #datasets for sim. combine fit
			  #sideband lo only, for sim-fit
			  if tmp_jet_mass>= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
			    rrv_mass_j_sb_lo.setVal(tmp_jet_mass)
			    dataset_mj_sb_lo.add(RooArgSet(rrv_mass_j_sb_lo,rrv_mass_lvj), 1)
			  #sideband hi
			  if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
			    rrv_mass_j_sb_hi.setVal(tmp_jet_mass)
			    dataset_mj_sb_hi.add(RooArgSet(rrv_mass_j_sb_hi,rrv_mass_lvj), 1)	
			  #signal
			  if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
			    rrv_mass_j_sig.setVal(tmp_jet_mass)
			    dataset_mj_sig.add(RooArgSet(rrv_mass_j_sig,rrv_mass_lvj), 1)	
			else:
			  rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight )
			  rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit )

		        if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max:
		            hnum_4region.Fill(-1,tmp_event_weight );
		        if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
		            hnum_4region.Fill(0,tmp_event_weight);
		        if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max:
		            hnum_4region.Fill(1,tmp_event_weight);

		        hnum_4region.Fill(2,tmp_event_weight);


		### scale to lumi for MC in 4fit datasets
		rrv_scale_to_lumi=RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi)
		rrv_scale_to_lumi.Print()
		getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)

		### prepare m_lvj dataset to be compared with the fit results
		rrv_number_dataset_signal_region_mlvj=RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(1));
		rrv_number_dataset_AllRange_mlvj=RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(2));
		
		getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj)
		getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj)

		### import the dataset       
		getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj);
		getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj);
		getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj);
		getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj);
		getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj);
		getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj);
		getattr(self.workspace4fit_,"import")(combData);
		getattr(self.workspace4fit_,"import")(combData4fit);
		if 'data' in label:
		  getattr(self.workspace4fit_,'import')(dataset_mj_sb_lo)
		  getattr(self.workspace4fit_,'import')(dataset_mj_sb_hi)
		  getattr(self.workspace4fit_,'import')(dataset_mj_sig)

		### write in the output 
		self.file_out.write("\n%s events number in m_lvj from dataset: %s"%(label,rdataset_signal_region_mlvj.sumEntries()))

		### prepare m_j dataset
		rrv_number_dataset_sb_lo_mj=RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
		rrv_number_dataset_signal_region_mj=RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));
		rrv_number_dataset_sb_hi_mj=RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));
		getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj)
		getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj)
		getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj)
		        
		getattr(self.workspace4fit_,"import")(rdataset_mj)
		getattr(self.workspace4fit_,"import")(rdataset4fit_mj)

		#@#write datasets to file
		w_tmp	= RooWorkspace('w_tmp_%s'%label)
		getattr(w_tmp,'import')(rrv_scale_to_lumi)
		getattr(w_tmp,'import')(rrv_number_dataset_signal_region_mlvj)
		getattr(w_tmp,'import')(rrv_number_dataset_AllRange_mlvj)
		getattr(w_tmp,'import')(rdataset_sb_lo_mlvj)
		getattr(w_tmp,'import')(rdataset_signal_region_mlvj)
		getattr(w_tmp,'import')(rdataset_sb_hi_mlvj)
		getattr(w_tmp,'import')(rdataset4fit_sb_lo_mlvj)
		getattr(w_tmp,'import')(rdataset4fit_signal_region_mlvj)
		getattr(w_tmp,'import')(rdataset4fit_sb_hi_mlvj)
		getattr(w_tmp,'import')(combData)
		getattr(w_tmp,'import')(combData4fit)
		getattr(w_tmp,'import')(rrv_number_dataset_sb_lo_mj)
		getattr(w_tmp,'import')(rrv_number_dataset_signal_region_mj)
		getattr(w_tmp,'import')(rrv_number_dataset_sb_hi_mj)
		getattr(w_tmp,'import')(rdataset_mj)
		getattr(w_tmp,'import')(rdataset4fit_mj)
		if 'data' in label:
		  getattr(w_tmp,'import')(dataset_mj_sb_lo)
		  getattr(w_tmp,'import')(dataset_mj_sb_hi)
		  getattr(w_tmp,'import')(dataset_mj_sig)
		if 'WJets0_' in label:
			fileOut	= TFile.Open('cards_%s_%s_900_3500/datasets_%s_%s.root'%(self.channel,self.wtagger_label,self.channel,self.wtagger_label),'recreate')
		else:
			fileOut	= TFile.Open('cards_%s_%s_900_3500/datasets_%s_%s.root'%(self.channel,self.wtagger_label,self.channel,self.wtagger_label),'update')
		w_tmp.Write()
		fileOut.Close()
		#@#

		#### print everything
		
		rdataset_sb_lo_mlvj.Print();
		rdataset_signal_region_mlvj.Print();
		rdataset_sb_hi_mlvj.Print();
		rdataset_mj.Print();
		rdataset4fit_sb_lo_mlvj.Print();
		rdataset4fit_signal_region_mlvj.Print();
		rdataset4fit_sb_hi_mlvj.Print();
		rdataset4fit_mj.Print();
		rrv_number_dataset_signal_region_mlvj.Print()
		rrv_number_dataset_AllRange_mlvj.Print()
		rrv_number_dataset_sb_lo_mj.Print()
		rrv_number_dataset_signal_region_mj.Print()
		rrv_number_dataset_sb_hi_mj.Print()
		combData.Print();
		combData4fit.Print();

    def get_mj_and_mlvj_dataset_from_file(self,in_file_name, label, jet_mass ):
	fileIn	= TFile.Open('cards_%s_%s_900_3500/datasets_%s_%s.root'%(self.channel,self.wtagger_label,self.channel,self.wtagger_label))
	w_tmp	= fileIn.Get('w_tmp_%s'%label)
	getattr(self.workspace4fit_,'import')(w_tmp.var('rrv_scale_to_lumi'+label+'_'+self.channel))
	getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("combData"+label+"_"+self.channel))
	getattr(self.workspace4fit_,'import')(w_tmp.data("combData4fit"+label+"_"+self.channel))
	getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj"))
	getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj"))
	getattr(self.workspace4fit_,'import')(w_tmp.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset"+label+"_"+self.channel+"_mj"))
	getattr(self.workspace4fit_,'import')(w_tmp.data("rdataset4fit"+label+"_"+self.channel+"_mj"))
	if 'data' in label:
	  getattr(self.workspace4fit_,'import')(w_tmp.data("dataset_mj_sb_lo"))
	  getattr(self.workspace4fit_,'import')(w_tmp.data("dataset_mj_sb_hi"))
	  getattr(self.workspace4fit_,'import')(w_tmp.data("dataset_mj_W%s_sig"%self.wtagger_label[2]))
	fileIn.Close()

	
        
    #### function to run the selection on data to build the datasets 
    def get_data(self):
        print "############### get_data ########################"
        self.get_mj_and_mlvj_dataset(self.file_data,"_data_xww", self.jetalgo)
        #self.get_mj_and_mlvj_dataset(self.file_data,"_data_xww", "Mjsoftdrop")
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_dataset_signal_region_data_xww_%s_mlvj"%(self.channel)).clone("observation_for_counting_xww"))

    #### Define the steps to fit STop MC in the mj and mlvj spectra
    def fit_STop(self):
        print "############################## fit_STop  #################################"
        self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop_xww", self.jetalgo)
        self.fit_mj_single_MC(self.file_STop_mc,"_STop_xww","ExpGaus");
        self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop_xww","_sb_lo","Exp", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop_xww","_signal_region","Exp", 1, 0, 1);
        print "________________________________________________________________________"

    ##### Define the steps to fit VV MC in the mj and mlvj spectra
    def fit_VV(self):
        print "############################# fit_VV ################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV_xww", self.jetalgo)
        ### fitting shape as a function of the mlvj region -> signal mass
        if self.wtagger_label.find("LP") != -1:
           self.fit_mj_single_MC(self.file_VV_mc,"_VV_xww","ExpGaus");
        else:
           self.fit_mj_single_MC(self.file_VV_mc,"_VV_xww","2_2Gaus");       
        self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV_xww","_sb_lo","Exp", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV_xww","_signal_region",self.MODEL_4_mlvj, 1, 0, 1);
        print "________________________________________________________________________"

    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_TTbar(self):
        print "################################ fit_TTbar #########################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar_xww", self.jetalgo)# to get the shape of m_lvj
        if self.wtagger_label.find("LP") != -1: self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_xww","ExpGaus");
        else:                          self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_xww","2Gaus_ErfExp");
	self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar_xww","_sb_lo","ExpN");
        self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar_xww","_signal_region","ExpN",1, 0, 1);
        print "________________________________________________________________________"

    ##### Define the steps to fit WJets MC in the mj and mlvj spectra
    def fit_WJets(self):
        print "######################### fit_WJets ########################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0_xww", self.jetalgo)# to get the shape of m_lvj
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets01_xww", self.jetalgo)# to get the shape of m_lvj
	
	### Fit in mj depends on the mlvj lower limit -> fitting the turn on at low mass or not
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_xww","ErfExp");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","ErfExp")        
	self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","User1");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_xww","User1");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","ErfExp");

            
        #### Fit the mlvj in sb_lo, signal region using two different model as done in the mj
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0_xww","_sb_lo",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0_xww","_signal_region",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets01_xww","_sb_lo",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets01_xww","_signal_region",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);       
        print "________________________________________________________________________"


    ##### Fit of all the MC in both mj and mlvj : Signal, TTbar, STop, VV and Wjets
    def fit_AllSamples_Mj_and_Mlvj(self):
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
#	sys.exit()
        self.fit_WJets()
        self.fit_TTbar()
        self.fit_VV()
        self.fit_STop()
        print "________________________________________________________________________"

    ##### Analysis with sideband alpha correction 
    def analysis_sideband_correction_method1(self):
        print "##################### Start sideband correction full analysis ##############";
        ### Fit all MC components in both mj and mlvj
        self.fit_AllSamples_Mj_and_Mlvj();
        ### take the real data
        self.get_data()
        ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
        self.fit_WJetsNorm();
        ### fit data in the mlvj low sideband with two different models
        self.fit_mlvj_in_Mj_sideband("_WJets01_xww","_sb_lo",self.MODEL_4_mlvj_alter,1)
        self.fit_mlvj_in_Mj_sideband("_WJets0_xww","_sb_lo",self.MODEL_4_mlvj,1)

        ### Prepare the workspace and datacards     
        self.prepare_limit("sideband_correction_method1",1,0,0)
        ### finale plot and check of the workspace
        self.read_workspace(1)
        
    ##### Analysis with no shape uncertainty on alpha
    def analysis_sideband_correction_method1_without_shape_and_psmodel_systematic(self):
        #### fit all the MC samples 
        self.fit_AllSamples_Mj_and_Mlvj()
        #### take the real data
        self.get_data()
        #### fit WJets just with one shape parametrization
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww");
        #### fit sb lo with just one parametrization
        self.fit_mlvj_in_Mj_sideband("_WJets0_xww","_sb_lo", self.MODEL_4_mlvj,1)
        #### prepare limit 
        self.prepare_limit("sideband_correction_method1",1,0,0)
        #### read the workspace
        self.read_workspace(1)

    ##### Prepare the workspace for the limit and to store info to be printed in the datacard
    def prepare_limit(self,mode, isTTbarFloating=0, isVVFloating=0, isSTopFloating=0):
	
	#get correct normalizations(mj signal region, 900 < mlvj < 3500), save pre- and post-fit values
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_STop_xww_%s_mj_postfit_signal_region'%self.channel).clone("rate_STop_xww_for_unbin"))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_VV_xww_%s_mj_prefit_signal_region'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_VV_xww_%s_mj_postfit_signal_region'%self.channel).clone("rate_VV_xww_for_unbin"))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_TTbar_xww_%s_mj_postfit_signal_region'%self.channel).clone("rate_TTbar_xww_for_unbin"))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_WJets0_xww_%s_mj_postfit_signal_region'%self.channel).clone("rate_WJets_xww_for_unbin"))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_STop_xww_%s_mj_postfit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_VV_xww_%s_mj_postfit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_TTbar_xww_%s_mj_postfit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_WJets0_xww_%s_mj_postfit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_STop_xww_%s_mj_prefit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_VV_xww_%s_mj_prefit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_TTbar_xww_%s_mj_prefit'%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.var('rrv_number_WJets0_xww_%s_mj_prefit'%self.channel))
	self.workspace4fit_.var('rrv_number_STop_xww_%s_mj_postfit_signal_region'%self.channel).Print()
	self.workspace4fit_.var('rrv_number_VV_xww_%s_mj_prefit_signal_region'%self.channel).Print()
	self.workspace4fit_.var('rrv_number_VV_xww_%s_mj_postfit_signal_region'%self.channel).Print()
	self.workspace4fit_.var('rrv_number_TTbar_xww_%s_mj_postfit_signal_region'%self.channel).Print()
	self.workspace4fit_.var('rrv_number_WJets0_xww_%s_mj_postfit_signal_region'%self.channel).Print()
	
	###prepare simultaneaous fit
	#save sb pdfs
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf("model_data_WJets0_xww_sb_lo_mlvj").clone("model_data_sb_%s_mlvj"%(self.channel)))
	#get m_pruned and mlvj_sb datasets
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.data("rdataset_data_xww_sb_lo_%s_mlvj"%self.channel))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('dataset_mj_sb_lo'))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('dataset_mj_sb_hi'))
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.data('dataset_mj_W%s_sig'%self.wtagger_label[2]))	
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.data("rdataset_data_xww_%s_mj"%self.channel))
	#get mj-pdfs
	getattr(self.workspace4limit_,'import')(self.workspace4fit_.pdf('model_data_xww_%s_mj'%self.channel))
	

        print "####################### prepare_limit for %s method ####################"%(mode);

        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_mass_lvj"));


        ### Get the dataset for data into the signal region
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.data("rdataset_data_xww_signal_region_%s_mlvj"%(self.channel)).Clone("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)))
        ### Take the corrected pdf from the alpha method for the WJets
        if mode=="sideband_correction_method1":
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WJets0_xww_signal_region_%s_after_correct_mlvj"%(self.channel)).clone("WJets_xww_%s_%s"%(self.channel, self.wtagger_label)));

        if isTTbarFloating:
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_xww_signal_region_%s_mlvj_Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV"%(self.channel, self.channel, self.wtagger_label)).clone("TTbar_xww_%s_%s"%(self.channel,self.wtagger_label)))
        else :
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_xww_signal_region_%s_mlvj"%(self.channel)).clone("TTbar_xww_%s_%s"%(self.channel,self.wtagger_label)))

        if isSTopFloating :     
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_STop_xww_signal_region_%s_mlvj_Deco_xww_STop_signal_region_%s_%s_mlvj_13TeV"%(self.channel, self.channel, self.wtagger_label)).clone("STop_xww_%s_%s"%(self.channel,self.wtagger_label)))
        else :
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_STop_xww_signal_region_%s_mlvj"%(self.channel)).clone("STop_xww_%s_%s"%(self.channel,self.wtagger_label)))

        if isVVFloating :    
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_xww_signal_region_%s_mlvj_Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV"%(self.channel, self.channel, self.wtagger_label)).clone("VV_%s_%s"%(self.channel,self.wtagger_label)))
        else:
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_xww_signal_region_%s_mlvj"%(self.channel)).clone("VV_xww_%s_%s"%(self.channel,self.wtagger_label)))

	### Fix all the Pdf parameters 
        rrv_x = self.workspace4limit_.var("rrv_mass_lvj");

        self.fix_Pdf(self.workspace4limit_.pdf("TTbar_xww_%s_%s"%(self.channel,self.wtagger_label)), RooArgSet(rrv_x) );            
        self.fix_Pdf(self.workspace4limit_.pdf("STop_xww_%s_%s"%(self.channel,self.wtagger_label)), RooArgSet(rrv_x));
        self.fix_Pdf(self.workspace4limit_.pdf("VV_xww_%s_%s"%(self.channel,self.wtagger_label)), RooArgSet(rrv_x));
        self.fix_Pdf(self.workspace4limit_.pdf("WJets_xww_%s_%s"%(self.channel,self.wtagger_label)), RooArgSet(rrv_x));
        
        print " ############## Workspace for limit ";
        parameters_workspace = self.workspace4limit_.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next()
        
        params_list = [];
        ### main modality for the alpha function method
        if mode=="sideband_correction_method1":

            if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));

                
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)))

                if isTTbarFloating !=0 :
                 self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                 params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));

            if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));

                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)))


                ### TTbar use expN
                if isTTbarFloating !=0:
                    print "##################### TTbar will float in the limit procedure + final plot ######################";
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
		    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
		    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)));
		    

                ### VV use ExpTail:
                if isVVFloating !=0:
                  print "##################### VV will float in the limit procedure + final plot ######################";
                  self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_VV);
                  self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_VV);
                  params_list.append(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
                  params_list.append(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)));
 
                ### STop use Exp:
                if isSTopFloating !=0:
                  print "##################### STop will float in the limit procedure + final plot ######################";
                  self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_STop);
                  params_list.append(self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
                                       
        
        ### calculate the shape uncertainty for cut-and-counting
        self.rrv_counting_uncertainty_from_shape_uncertainty = RooRealVar("rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.channel),"rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.channel),0);
        self.rrv_counting_uncertainty_from_shape_uncertainty.setError( Calc_error("WJets_xww_%s_%s"%(self.channel,self.wtagger_label), "rrv_mass_lvj" ,self.FloatingParams,self.workspace4limit_,"signal_region") );
        self.rrv_counting_uncertainty_from_shape_uncertainty.Print();

        print " param list ",params_list ;
            
        ### Print the datacard for unbin and couting analysis
        self.print_limit_datacard("unbin",params_list);
        self.print_limit_datacard("counting");

        if mode=="sideband_correction_method1":

          if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)) );

                if isTTbarFloating!=0:
                 self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel,self.wtagger_label)));


          if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );


                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)) );

                if isTTbarFloating!=0:
                    self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel,self.wtagger_label)));

                if isVVFloating!=0:     
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)));

                if isSTopFloating!=0:
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel,self.wtagger_label)));
                    

        ### Add the floating list to the combiner --> the pdf which are not fixed are floating by default
        getattr(self.workspace4limit_,"import")(self.FloatingParams);

        ### Save the workspace
        self.save_workspace_to_file();

    #### Method used to print the general format of the datacard for both counting and unbinned analysis
    def print_limit_datacard(self, mode, params_list=[]):
        print "############## print_limit_datacard for %s ################"%(mode)
        if not (mode == "unbin" or mode == "counting"):
            print "print_limit_datacard use wrong mode: %s"%(mode);raw_input("ENTER");

        ### open the datacard    
        datacard_out = open(getattr(self,"file_datacard_%s"%(mode)),"w");

        ### start to print inside 
        datacard_out.write( "imax 1" )
        datacard_out.write( "\njmax 4" )
        datacard_out.write( "\nkmax *" )
        datacard_out.write( "\n--------------- ")

        if mode == "unbin":
            fnOnly = ntpath.basename(self.file_rlt_root) ## workspace for limit --> output file for the workspace
                
            datacard_out.write("\nshapes WJets_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes TTbar_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes STop_xww   CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes VV_xww     CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes data_obs   CMS_xww_%s1J%s  %s %s:$PROCESS_xww_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write( "\n--------------- ")
            
        datacard_out.write( "\nbin CMS_xww_%s1J%s "%(self.channel,self.wtagger_label));    
        if mode == "unbin":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)).sumEntries()) )
        if mode == "counting":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.var("observation_for_counting_xww").getVal()) )
            
        datacard_out.write( "\n------------------------------" );

        datacard_out.write( "\nbin CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label));

        datacard_out.write( "\nprocess -1 1 2 3 4" );

        datacard_out.write( "\n-------------------------------- " )

 
        ### WJets Normalization from data fit -> data driven
        if self.number_WJets_insideband >0:
            datacard_out.write( "\nCMS_xww_WJ_norm_13TeV gmN %0.3f %0.3f - - -"%(self.number_WJets_insideband, getattr(self, "datadriven_alpha_WJets_xww_%s"%(mode)) ) )
        else:
            datacard_out.write( "\nCMS_xww_WJ_norm_%s_%s_13TeV lnN - %0.3f - - -"%(self.channel, self.wtagger_label, 1+ self.workspace4limit_.var('rate_WJets_xww_for_unbin').getError()/self.workspace4limit_.var('rate_WJets_xww_for_unbin').getVal() ) );

        if self.channel == "mu":
            self.channel_short = "m"
        elif self.channel =="el":
            self.channel_short = "e"
        elif self.channel =="em":
            self.channel_short = "em"

        ### print shapes parameter to be taken int account
        if mode == "unbin":
            for ipar in params_list:
              print "Name %s",ipar.GetName();
              if TString(ipar.GetName()).Contains("Deco_TTbar_xww_signal_region"):
               datacard_out.write( "\n%s param %0.1f %0.1f "%( ipar.GetName(), ipar.getVal(), ipar.getError() ) )
              else:
               datacard_out.write( "\n%s param %0.1f %0.1f "%( ipar.GetName(), ipar.getVal(), ipar.getError() ) )

	#@# write WJets normalisation 
	datacard_out.write("\n #WJets number of events sideband region: %s"%self.workspace4fit_.var("rrv_number_WJets0_xww_in_mj_sb_region_from_fitting_%s"%self.channel).getVal())

    #### Method used in order to save the workspace in a output root file
    def save_workspace_to_file(self):
        self.workspace4limit_.writeToFile(self.file_rlt_root);
        self.file_out.close()

    #### Read the final workspace and produce the latest plots 
    def read_workspace(self, logy=0):

        ### Taket the workspace for limits  
        file = TFile(self.file_rlt_root) ;
        workspace = file.Get("workspace4limit_") ;
        workspace.Print()

        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------";
        parameters_workspace = workspace.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next()
        print "---------------------------------------------";

        workspace.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)).Print()

        print "----------- Pdf in the Workspace -------------";
        pdfs_workspace = workspace.allPdfs();
        par = pdfs_workspace.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            param.Print();
            param = par.Next()
        print "----------------------------------------------";

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label));
            
        model_pdf_WJets  = workspace.pdf("WJets_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_VV     = workspace.pdf("VV_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_TTbar  = workspace.pdf("TTbar_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_STop   = workspace.pdf("STop_xww_%s_%s"%(self.channel,self.wtagger_label));
	
        model_pdf_WJets.Print();
        model_pdf_VV.Print();
        model_pdf_TTbar.Print();
        model_pdf_STop.Print();
	
        rrv_number_WJets  = workspace.var("rate_WJets_xww_for_unbin");
        rrv_number_VV     = workspace.var("rate_VV_xww_for_unbin");
        rrv_number_TTbar  = workspace.var("rate_TTbar_xww_for_unbin");
        rrv_number_STop   = workspace.var("rate_STop_xww_for_unbin");

        rrv_number_WJets.Print();
        rrv_number_VV.Print();
        rrv_number_TTbar.Print();
        rrv_number_STop.Print();

        #### Prepare the final plot starting from total background 
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww","rrv_number_Total_background_MC_xww",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal()+
                rrv_number_TTbar.getVal()+
                rrv_number_STop.getVal());

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets.getError()* rrv_number_WJets.getError()+
                rrv_number_VV.getError()* rrv_number_VV.getError()+
                rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
                rrv_number_STop.getError() *rrv_number_STop.getError()
                ));

        #### Total pdf 
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww","model_Total_background_MC_xww",RooArgList(model_pdf_WJets,model_pdf_VV,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_TTbar,rrv_number_STop));

        if data_obs.sumEntries() != 0:
         #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
         scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
        else:
	 scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()   
        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());
        
	model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV"), RooFit.Components("VV_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines());


	model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop"), RooFit.Components("STop_xww_%s_%s"%(self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());

        #solid line
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());

	model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar_line_invisible"), RooFit.Components("VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_%s_%s,STop_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());
	
	model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop_line_invisible"), RooFit.Components("STop_xww_%s_%s"%(self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(1), RooFit.VLines());


        #### plot the observed data using poissonian error bar
	self.getData_PoissonInterval(data_obs,mplot);
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Invisible());

        mplot_pull=self.get_pull(rrv_x,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC,self.FloatingParams,workspace ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
	self.leg=self.legend4Plot(mplot,0,1,0.25,0.,0.1,0.,0);

        #self.leg.SetTextSize(0.036);
        mplot.addObject(self.leg);
	pt1 = ROOT.TPaveText(0.6180905,0.515644,0.8291457,0.537992,"NDC")
	pt1.SetTextFont(42)
	pt1.SetTextSize(0.05)
	#pt.SetTextAlign(12)
	pt1.SetFillColor(0)
	pt1.SetFillStyle(0)
	pt1.SetBorderSize(0)
	text = pt1.AddText("")
	if options.category.find('Z') != -1: text = pt1.AddText("WZ category")
	elif options.category.find('W') != -1: text = pt1.AddText("WW category")
	text.SetTextFont(62)
	#mplot.addObject(pt1)
	            
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
	mplot.GetXaxis().SetTitle('M_{WV}')
            
        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal());
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-self.nPar_float_in_fitTo;
        #print "nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof );
        datahist = data_obs.binnedClone( data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone" )
	
        parameters_list = RooArgList();
        self.draw_canvas_with_pull( rrv_x,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir),"check_workspace_for_limit","",0,1);
	#self.draw_canvas_with_pull( rrv_x,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/ExtraPlots/"%(self.plotsDir),"check_workspace_for_limit","",0,1);
	if self.channel == 'el':
		mplot.GetYaxis().SetRangeUser(1e-3,2e5)
	mplot.GetXaxis().SetTitle('M_{WV} (GeV)')
	mplot.GetYaxis().SetRangeUser(2e-3,1e5)
	self.draw_canvas(mplot,"%s/ExtraPlots/"%(self.plotsDir),"check_workspace_for_limit",0,logy,0,0,1);
	mplot.Print()

		
        

### funtion to run the complete alpha analysis
def pre_limit_sb_correction(method, channel, jetalgo="Mjsoftdrop", in_mj_min=30, in_mj_max=140, in_mlvj_min=400, in_mlvj_max=1400, fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1"): 

    print "#################### pre_limit_sb_correction: channel %s, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    
    boostedW_fitter=doFit_wj_and_wlvj(channel, jetalgo,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    getattr(boostedW_fitter,"analysis_sideband_correction_%s"%(method) )();

                                                
### funtion to run the analysis without systematics
def pre_limit_sb_correction_without_systematic( channel, in_mj_min=30, in_mj_max=140, in_mlvj_min=400, in_mlvj_max=1400, fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1"):

    print "#################### pre_limit_sb_correction_without_systermatic: channel %s, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    boostedW_fitter=doFit_wj_and_wlvj(channel,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    boostedW_fitter.analysis_sideband_correction_method1_without_shape_and_psmodel_systematic()


def pre_limit_simple(channel):
    print "######################### pre_limit_simple for %s sampel"%(channel)
    pre_limit_sb_correction_without_systematic(channel,40,130,700,5000,"ExpN","ExpTail")

#### Main Code
if __name__ == '__main__':

    channel=options.channel;
            
    if options.simple:
        print '################# simple mode for %s sample'%(channel)
        pre_limit_simple(channel);
    else:
        pre_limit_sb_correction("method1",channel,options.jetalgo,40,150,options.mlvj_lo,5000,"ExpN","ExpTail")

    for i in fitresultsmj:
	i.Print()
    for i in fitresultsmlvj:
	i.Print()
    for i in fitresultsfinal:
	i.Print()




