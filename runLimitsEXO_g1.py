#! /usr/bin/env python 
#python runLimitsEXO_g1.py --computeLimits --plotLimits --channel mu --category HP --signal Wprime_WZ --dir cards_mu
import os
import glob
import math
import array
import sys
import time
import random

from array import array

import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TH1D, TString, TPaveText, TGaxis
import subprocess
from subprocess import Popen
from optparse import OptionParser
import CMS_lumi, tdrstyle

############################################
#             Job steering                 #
############################################

parser = OptionParser()

### to compute limits
parser.add_option('--computeLimits', action='store_true', dest='computeLimits', default=False, help='compute limits')

### to compute p-value
parser.add_option('--computePvalue', action='store_true', dest='computePvalue', default=False, help='compute p-value')

### to plot limits
parser.add_option('--plotLimits', action='store_true', dest='plotLimits', default=False, help='plot limits')

### to plot p-value
parser.add_option('--plotPvalue', action='store_true', dest='plotPvalue', default=False, help='plot p-value')

### other options 
parser.add_option('--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--signal', action="store",type="string",dest="signal",default="BulkG_WW")
parser.add_option('--dir', action="store",type="string",dest="dir",default="results")
parser.add_option('--batch', '-b',action="store_true",dest="batch",default=False)

(options, args) = parser.parse_args()

if options.batch: gROOT.SetBatch(ROOT.kTRUE)

### Mass points ###

#mass      = [800,1000,1200,1400,1600,1800,2000,2500,3000,4000]
mass      = [800,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500]
#mass = [2000]

### cross-sections for BulkG ###

#xsBulkG = {800:0.0030423317346004058,1000:0.00081998770967753234,1200:0.00027358240573576689,1400:0.00010454961631658292,1600:0.0000932846,1800:0.0000420266,
#           2000:9.5807525950148653e-06,2500:1.7942831152359325e-06,3000:3.9299443484538374e-07,3500:0.4787424623e-07,4000:0.0652212882e-07,
#	   4500:0.0088853962e-07}

xsBulkG = {800:0.076058293365,1000:0.0204996927419,1200:0.00683956014339,1400:0.00261374040791,1600:0.0011483744,1800:0.0004882930,
           2000:0.000239518814875,2500:4.48570778809e-05,3000:9.82486087113e-06,3500:41.9758278604e-07,4000:24.3824379692e-07,
     	   4500:20.8967256514e-07}

### cross-sections for RSG ###

xsRSG = {800:1.16691,1000:0.377865,1200:0.144482,1400:0.0616708,1600:0.0288651,1800:0.0141334,
         2000:0.00751431,2500:0.00167726,3000:0.000443483,3500:0.000133915,4000:0.0000424117,4500:0.0000130705}

### cross-sections for Wprime ###

xsWprime = {800:1.587885*0.428731, 1000:0.986533*0.460359, 1200:0.535394*0.467605, 1400:0.2955239*0.470248, 1600:0.1681478*0.471413,
            1800:0.0984325*0.471991, 2000:0.058998*0.472301, 2500:0.01771031*0.472619, 3000:0.00567529*0.472715, 3500:0.001878491*0.472748,
	    4000:6.0983956656e-08*0.472759,4500:168.1479623273e-07}

TGaxis.SetMaxDigits(3)
	 
### methods ###
def get_canvas(text):

   tdrstyle.setTDRStyle()
   CMS_lumi.lumi_13TeV = "2.1 fb^{-1}," + text
   CMS_lumi.writeExtraText = 1
   CMS_lumi.extraText = ""

   iPos = 0
   if( iPos==0 ): CMS_lumi.relPosX = 0.12

   H_ref = 600 
   W_ref = 800 
   W = W_ref
   H  = H_ref

   T = 0.08*H_ref
   B = 0.12*H_ref 
   L = 0.12*W_ref
   R = 0.04*W_ref

   canvas = ROOT.TCanvas("c2","c2",50,50,W,H)
   canvas.SetFillColor(0)
   canvas.SetBorderMode(0)
   canvas.SetFrameFillStyle(0)
   canvas.SetFrameBorderMode(0)
   canvas.SetLeftMargin( L/W+0.01 )
   canvas.SetRightMargin( R/W+0.03 )
   canvas.SetTopMargin( T/H )
   canvas.SetBottomMargin( B/H+0.03 )
   canvas.SetGrid()
   canvas.SetLogy()
   
   return canvas
   
def getAsymLimits(file):
    
    f = ROOT.TFile(file)
    t = f.Get("limit")
    entries = t.GetEntries()
    
    lims = [0,0,0,0,0,0]
    
    for i in range(entries):
        
        t.GetEntry(i)
        t_quantileExpected = t.quantileExpected
        t_limit = t.limit
        
        #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected
        
        if t_quantileExpected == -1.: lims[0] = t_limit
        elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit
        elif t_quantileExpected >= 0.15  and t_quantileExpected <= 0.17:  lims[2] = t_limit
        elif t_quantileExpected == 0.5: lims[3] = t_limit
        elif t_quantileExpected >= 0.83  and t_quantileExpected <= 0.85:  lims[4] = t_limit
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit
        else: print "Unknown quantile!"
    
    return lims

####################################################
### Get PValue from combine -M ProfileLikelihood ###
####################################################

def getPValueFromCard( file ):

    f = ROOT.TFile(file)
    t = f.Get("limit")
    entries = t.GetEntries()
    
    lims = 1
    
    for i in range(1):
        
        t.GetEntry(i)
        lims = t.limit
    
    return lims


##########################################
### Make Limit Plot --> Brazilian Plot ###
##########################################

def doLimitPlot( suffix ):
    
    print "**********************************"
    print suffix
    xbins     = array('d', [])
    xbins_env = array('d', [])
    ybins_exp = array('d', [])
    ybins_obs = array('d', [])
    ybins_1s  = array('d', [])
    ybins_2s  = array('d', [])
    ybins_xs  = array('d', [])

    br_lvjj = 2*0.322*0.6760
    if options.signal.find('Wprime') != -1: br_lvjj = 0.322*0.6991
    
    print br_lvjj
    if TString(suffix).Contains("_el_") :
        text = "W#rightarrow e#nu"
    #    br_lvjj = 2*0.1075*0.6760
    elif TString(suffix).Contains("_mu_") :
        text = "W#rightarrow #mu#nu"
    #    br_lvjj = 2*0.1057*0.6760
    elif TString(suffix).Contains("_combo_") :    
        text = "W#rightarrow l#nu"
    #    br_lvjj = 2*0.108*0.6760
	
    if TString(suffix).Contains("HP"):
        text+=" HP"
    elif TString(suffix).Contains("LP") and not(TString(suffix).Contains("ALLP")):
        text+=" LP"		
    elif TString(suffix).Contains("ALLP"):
        text+=" ALLP"	
	        
    signal = "G_{Bulk}"
    xsDict = xsBulkG
    if options.signal.find("RS1") != -1: 
       signal = "G_{RS1}"
       xsDict = xsRSG
    if options.signal.find("Wprime") != -1: 
       signal = "W'"
       xsDict = xsWprime
    
    for i in range(len(mass)):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3]*0.01/br_lvjj );
        ybins_obs.append( curAsymLimits[0]*0.01/br_lvjj );
        ybins_2s.append( curAsymLimits[1]*0.01/br_lvjj );
        ybins_1s.append( curAsymLimits[2]*0.01/br_lvjj );
        ybins_xs.append(xsDict[mass[i]]);
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5]*0.01/br_lvjj );
        ybins_1s.append( curAsymLimits[4]*0.01/br_lvjj );
    
    
    canv = get_canvas(text)
    canv.cd()
        
    curGraph_exp    = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp)
    curGraph_obs    = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs)
    curGraph_xs     = ROOT.TGraph(nPoints,xbins,ybins_xs)
    curGraph_1s     = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s)
    curGraph_2s     = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s)
    
    curGraph_obs.SetMarkerStyle(20)
    curGraph_obs.SetLineWidth(3)
    curGraph_obs.SetLineStyle(1)
    curGraph_obs.SetMarkerSize(1.6)
    curGraph_exp.SetMarkerSize(1.3)
    curGraph_exp.SetMarkerColor(ROOT.kBlack)

    curGraph_exp.SetLineStyle(2)
    curGraph_exp.SetLineWidth(3)
    curGraph_exp.SetMarkerSize(2)
    curGraph_exp.SetMarkerStyle(24)
    curGraph_exp.SetMarkerColor(ROOT.kBlack)

    curGraph_xs.SetLineStyle(ROOT.kSolid)
    curGraph_xs.SetFillStyle(3344)
    curGraph_xs.SetLineWidth(2)
    curGraph_xs.SetMarkerSize(2)
    curGraph_xs.SetLineColor(ROOT.kRed)
    
    curGraph_1s.SetFillColor(ROOT.kGreen)
    curGraph_1s.SetFillStyle(1001)
    curGraph_1s.SetLineStyle(ROOT.kDashed)
    curGraph_1s.SetLineWidth(3)

    curGraph_2s.SetFillColor(ROOT.kYellow)
    curGraph_2s.SetFillStyle(1001)
    curGraph_2s.SetLineStyle(ROOT.kDashed)
    curGraph_2s.SetLineWidth(3)

    hrl_SM = canv.DrawFrame(750,10e-8, 4050, 10)
    ytitle = ""
    print "suffix is %s" %(suffix)
    if suffix.find("_el_") != -1:
       ytitle = "#sigma_{95%} #times BR_{"+signal+"#rightarrow WW} (pb)"
    if suffix.find("_mu_") != -1:
       ytitle = "#sigma_{95%} #times BR_{"+signal+"#rightarrow WW} (pb)"
    if suffix.find("_combo_") != -1:
       ytitle = "#sigma_{95%} #times BR_{"+signal+"#rightarrow WW} (pb)"
    if options.signal.find('Wprime') != -1: ytitle = ytitle.replace('WW','WZ')
    hrl_SM.GetYaxis().SetTitle(ytitle)
    hrl_SM.GetYaxis().CenterTitle()
    hrl_SM.GetYaxis().SetTitleSize(0.06)
    hrl_SM.GetXaxis().SetTitleSize(0.06)
    hrl_SM.GetXaxis().SetLabelSize(0.045)
    hrl_SM.GetYaxis().SetLabelSize(0.045)
    hrl_SM.GetYaxis().SetTitleOffset(1)
    hrl_SM.GetXaxis().SetTitleOffset(1.1)
    hrl_SM.GetXaxis().SetTitle("M_{"+signal+"} (GeV)")
    hrl_SM.GetXaxis().CenterTitle()
    hrl_SM.SetMinimum(0.0001)
    if options.signal.find('Wprime') != -1: hrl_SM.SetMinimum(0.001)
    hrl_SM.SetMaximum(100)
    hrl_SM.GetYaxis().SetNdivisions(505)
        
    curGraph_2s.Draw("F")
    curGraph_1s.Draw("Fsame")
    curGraph_exp.Draw("Lsame")
    curGraph_obs.Draw("LPsame")
    curGraph_xs.Draw("Csame")

    pt = ROOT.TPaveText(0.6595477,0.222028,0.8944724,0.3653846,"NDC")
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.SetTextFont(62)
    pt.SetTextSize(0.038)
    pt.SetTextAlign(31)
    pt.SetFillColor(0)
    
    if TString(suffix).Contains("W"):
     text = pt.AddText("WW category")
    if TString(suffix).Contains("Z"):
     text = pt.AddText("WZ category")
    if not(TString(suffix).Contains("Z")) and not(TString(suffix).Contains("W")): 
     text = pt.AddText("WW+WZ categories combined")
    #text.SetTextFont(62)
    if TString(suffix).Contains("HP"):
     text = pt.AddText("#tau_{21} < 0.6")
    if TString(suffix).Contains("LP") and not(TString(suffix).Contains("ALLP")):
     text = pt.AddText("0.6 < #tau_{21} < 0.75")
    #text.SetTextFont(62)
       
    leg2 = ROOT.TLegend(0.4334171,0.6695804,0.9183417,0.8706294)
    leg2.SetLineWidth(2)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.038)
    leg2.SetTextAlign(12)

    leg2.AddEntry(curGraph_1s,"Asympt. CL_{S}  Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s,"Asympt. CL_{S}  Expected #pm 2#sigma","LF")
    theoleg = "#sigma_{TH} #times BR_{"+signal+"#rightarrow WW}"
    if options.signal.find('Wprime') != -1: theoleg = theoleg.replace('WW','WZ')
    if options.signal.find('Bulk') != -1: theoleg+=" #tilde{k}=0.5"
    if options.signal.find('Wprime') != -1: theoleg+=" , HVT_{B}"
    leg2.AddEntry(curGraph_xs,theoleg,"L")     
    leg2.Draw()
    pt.Draw()
             
    canv.Update()   
    canv.cd()
    CMS_lumi.CMS_lumi(canv, 4, 0)	   
    canv.cd()
    canv.Update()
    canv.RedrawAxis()
    canv.RedrawAxis("g")
    frame = canv.GetFrame()
    frame.Draw()   
    canv.cd()
    canv.Update()    
    
    os.system('mkdir LimitResult')
    
    canv.SaveAs("./LimitResult/Lim%s%s.png"%(suffix,options.signal))
    canv.SaveAs("./LimitResult/Lim%s%s.pdf"%(suffix,options.signal))
    canv.SaveAs("./LimitResult/Lim%s%s.root"%(suffix,options.signal))
    canv.SaveAs("./LimitResult/Lim%s%s.C"%(suffix,options.signal))
    
def doLimitPlot2( suffix ):
        
    xbins     = array('d', [])
    xbins_env = array('d', [])
    ybins_exp = array('d', [])
    ybins_obs = array('d', [])
    ybins_1s  = array('d', [])
    ybins_2s  = array('d', [])
    ybins_xs  = array('d', [])

    br_lvjj = 2*0.322*0.6760
    factor = 1.
    if options.signal.find('Wprime') != -1:
       br_lvjj = 0.322*0.6991
       factor = 1.
        
    if suffix.find("_el_") != -1:
        text = "W#rightarrow e#nu"
    #    br_lvjj = 0.1075*0.6760
    if suffix.find("_mu_") != -1:
        text = "W#rightarrow #mu#nu"
    #    br_lvjj = 0.1057*0.6760
    if suffix.find("_combo_") != -1:    
        text = "W#rightarrow l#nu"
    #    br_lvjj = 0.108*0.6760

    if TString(suffix).Contains("HP"):
        text+=" HP"
    elif TString(suffix).Contains("LP") and not(TString(suffix).Contains("ALLP")):
        text+=" LP"	
    elif TString(suffix).Contains("ALLP"):
        text+=" ALLP"	
		
    signal = "G_{Bulk}"
    xsDict = xsBulkG
    if options.signal.find("RS1") != -1: 
       signal = "G_{RS1}"
       xsDict = xsRSG
    if options.signal.find("Wprime") != -1: 
       signal = "W'"
       xsDict = xsWprime
                 
    for i in range(len(mass)):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3]*0.01 );
        ybins_obs.append( curAsymLimits[0]*0.01 );
        ybins_2s.append( curAsymLimits[1]*0.01  );
        ybins_1s.append( curAsymLimits[2]*0.01  );
        ybins_xs.append( xsDict[mass[i]]*br_lvjj );
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "higgsCombine_lim_%03d%s.Asymptotic.mH%03d.root"%(mass[i],suffix,mass[i]);
        curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5]*0.01 );
        ybins_1s.append( curAsymLimits[4]*0.01 );
    
    
    canv = get_canvas(text)
    canv.cd()
        
    curGraph_exp    = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp)
    curGraph_obs    = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs)
    curGraph_xs     = ROOT.TGraph(nPoints,xbins,ybins_xs)
    curGraph_1s     = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s)
    curGraph_2s     = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s)
    
    curGraph_obs.SetMarkerStyle(20)
    curGraph_obs.SetLineWidth(3)
    curGraph_obs.SetLineStyle(1)
    curGraph_obs.SetMarkerSize(1.6)
    curGraph_exp.SetMarkerSize(1.3)
    curGraph_exp.SetMarkerColor(ROOT.kBlack)

    curGraph_exp.SetLineStyle(2)
    curGraph_exp.SetLineWidth(3)
    curGraph_exp.SetMarkerSize(2)
    curGraph_exp.SetMarkerStyle(24)
    curGraph_exp.SetMarkerColor(ROOT.kBlack)

    curGraph_xs.SetLineStyle(ROOT.kSolid)
    curGraph_xs.SetFillStyle(3344)
    curGraph_xs.SetLineWidth(2)
    curGraph_xs.SetMarkerSize(2)
    curGraph_xs.SetLineColor(ROOT.kRed)
    
    curGraph_1s.SetFillColor(ROOT.kGreen)
    curGraph_1s.SetFillStyle(1001)
    curGraph_1s.SetLineStyle(ROOT.kDashed)
    curGraph_1s.SetLineWidth(3)

    curGraph_2s.SetFillColor(ROOT.kYellow)
    curGraph_2s.SetFillStyle(1001)
    curGraph_2s.SetLineStyle(ROOT.kDashed)
    curGraph_2s.SetLineWidth(3)

    hrl_SM = canv.DrawFrame(750,10e-8, 4050, 10)
    ytitle = ""
    print "suffix is %s" %(suffix)
    if suffix.find("_el_") != -1:
       ytitle = "#sigma_{95%} #times BR_{"+signal+"#rightarrow WW} #times BR_{WW #rightarrow e#nuqq} (pb)"
    if suffix.find("_mu_") != -1:
       ytitle = "#sigma_{95%} #times BR_{"+signal+"#rightarrow WW} #times BR_{WW #rightarrow #mu#nuqq} (pb)"
    if suffix.find("_combo_") != -1:
       ytitle = "#sigma_{95%} #times BR_{"+signal+"#rightarrow WW} #times BR_{WW #rightarrow l#nuqq} (pb)"
    if options.signal.find('Wprime') != -1: ytitle.replace('WW','WZ')
    hrl_SM.GetYaxis().SetTitle(ytitle)
    hrl_SM.GetYaxis().CenterTitle()
    hrl_SM.GetYaxis().SetTitleSize(0.06)
    hrl_SM.GetXaxis().SetTitleSize(0.06)
    hrl_SM.GetXaxis().SetLabelSize(0.045)
    hrl_SM.GetYaxis().SetLabelSize(0.045)
    hrl_SM.GetYaxis().SetTitleOffset(1)
    hrl_SM.GetXaxis().SetTitleOffset(1.1)
    hrl_SM.GetXaxis().SetTitle("M_{"+signal+"} (GeV)")
    hrl_SM.GetXaxis().CenterTitle()
    hrl_SM.SetMinimum(0.0001)
    if options.signal.find('Wprime') != -1: hrl_SM.SetMinimum(0.001)
    hrl_SM.SetMaximum(100)
    hrl_SM.GetYaxis().SetNdivisions(505)
        
    curGraph_2s.Draw("F")
    curGraph_1s.Draw("Fsame")
    curGraph_exp.Draw("Lsame")
    #curGraph_obs.Draw("LPsame")
    curGraph_xs.Draw("Csame")

    pt = ROOT.TPaveText(0.6595477,0.222028,0.8944724,0.3653846,"NDC")
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.SetTextFont(62)
    pt.SetTextSize(0.038)
    pt.SetTextAlign(31)
    pt.SetFillColor(0)
    
    if TString(suffix).Contains("W"):
     text = pt.AddText("WW category")
    if TString(suffix).Contains("Z"):
     text = pt.AddText("WZ category")
    if not(TString(suffix).Contains("Z")) and not(TString(suffix).Contains("W")): 
     text = pt.AddText("WW+WZ categories combined")
    #text.SetTextFont(62)
    if TString(suffix).Contains("HP"):
     text = pt.AddText("#tau_{21} < 0.6")
    if TString(suffix).Contains("LP") and not(TString(suffix).Contains("ALLP")):
     text = pt.AddText("0.6 < #tau_{21} < 0.75")
    #text.SetTextFont(62)
           
    leg2 = ROOT.TLegend(0.3555276,0.6643357,0.8404523,0.8653846)
    leg2.SetLineWidth(2)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.038)
    leg2.SetTextAlign(12)

    leg2.AddEntry(curGraph_1s,"Asympt. CL_{S}  Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s,"Asympt. CL_{S}  Expected #pm 2#sigma","LF")
    theoleg = "#sigma_{TH} #times BR_{"+signal+"#rightarrow WW} #times BR_{WW #rightarrow e#nuqq}"
    if suffix.find("_mu_") != -1: theoleg = theoleg.replace("e#nuqq","#mu#nuqq")
    if suffix.find("_combo_") != -1: theoleg = theoleg.replace("e#nuqq","l#nuqq")
    if options.signal.find('Wprime') != -1: theoleg = theoleg.replace('WW','WZ')
    if options.signal.find('Bulk') != -1: theoleg+=" #tilde{k}=0.5"
    if options.signal.find('Wprime') != -1: theoleg+=" , HVT_{B}"
    leg2.AddEntry(curGraph_xs,theoleg,"L")  
         
    leg2.Draw()
    pt.Draw()
             
    canv.Update()   
    canv.cd()
    CMS_lumi.CMS_lumi(canv, 4, 0)	   
    canv.cd()
    canv.Update()
    canv.RedrawAxis()
    canv.RedrawAxis("g")
    frame = canv.GetFrame()
    frame.Draw()   
    canv.cd()
    canv.Update()    
    
    os.system('mkdir LimitResult')
    canv.SaveAs("./LimitResult/Lim%s%s_lvjj.png"%(suffix,options.signal))
    canv.SaveAs("./LimitResult/Lim%s%s_lvjj.pdf"%(suffix,options.signal))
    canv.SaveAs("./LimitResult/Lim%s%s_lvjj.root"%(suffix,options.signal))
    canv.SaveAs("./LimitResult/Lim%s%s_lvjj.C"%(suffix,options.signal))

#################
### Main Code ###    
#################
    
if __name__ == '__main__':

    
    CHAN = options.channel
    DIR = options.dir
    CAT = options.category
    SIG = options.signal
    
    nPoints = len(mass)
    mLo = 0
    mHi = nPoints
    
    if CAT != 'ALLP' and (CHAN == 'mu' or CHAN == 'el') and (CAT.find('W') != -1 or CAT.find('Z') != -1): DIR = DIR+"/cards_%s_%s/" %(CHAN,CAT)
    
    os.chdir(DIR)

    ### Compute Limits
    if options.computeLimits:

     for i in range(mLo,mHi):
	
      print "##################"+str(mass[i])+"#####################"
      time.sleep(0.3)

      if options.category.find('W') != -1 or options.category.find('Z') != -1:
	  
       runCmmd = "combine -M Asymptotic -m %03d -n _lim_%03d_%s_%s_ -d wwlvj_%s_lvjj_M%03d_%s_%s_unbin.txt --run both -H ProfileLikelihood"%(mass[i],mass[i],CHAN,CAT,SIG,mass[i],CHAN,CAT);
       print runCmmd
       os.system(runCmmd)

       time.sleep(0.1)                	 

      if options.category.find('W') == -1 and options.category.find('Z') == -1 and CHAN != 'em':
	
       if options.category != 'ALLP':
       
        cmd = 'combineCards.py '
	cmd += 'cards_%s_%sW/wwlvj_%s_lvjj_M%03d_%s_%sW_unbin.txt ' %(CHAN,CAT,SIG,mass[i],CHAN,CAT)
	cmd += 'cards_%s_%sZ/wwlvj_%s_lvjj_M%03d_%s_%sZ_unbin.txt ' %(CHAN,CAT,SIG,mass[i],CHAN,CAT)
	cmd += '> wwlvj_%s_lvjj_M%03d_%s_%s_unbin.txt'%(SIG,mass[i],CHAN,CAT)
	
	print cmd
	os.system(cmd)

        runCmmd = "combine -M Asymptotic -m %03d -n _lim_%03d_%s_%s_ -d wwlvj_%s_lvjj_M%03d_%s_%s_unbin.txt --run both -H ProfileLikelihood"%(mass[i],mass[i],CHAN,CAT,SIG,mass[i],CHAN,CAT);
        print runCmmd
        os.system(runCmmd)

        time.sleep(0.1)
	
       elif options.category == 'ALLP':
       	 
        cmd = 'combineCards.py '
	cmd += 'cards_%s_HPW/wwlvj_%s_lvjj_M%03d_%s_HPW_unbin.txt ' %(CHAN,SIG,mass[i],CHAN)
	cmd += 'cards_%s_HPZ/wwlvj_%s_lvjj_M%03d_%s_HPZ_unbin.txt ' %(CHAN,SIG,mass[i],CHAN)
	cmd += 'cards_%s_LPW/wwlvj_%s_lvjj_M%03d_%s_LPW_unbin.txt ' %(CHAN,SIG,mass[i],CHAN)
	cmd += 'cards_%s_LPZ/wwlvj_%s_lvjj_M%03d_%s_LPZ_unbin.txt ' %(CHAN,SIG,mass[i],CHAN)
	cmd += '> wwlvj_%s_lvjj_M%03d_%s_%s_unbin.txt'%(SIG,mass[i],CHAN,CAT)
	
	print cmd
	os.system(cmd)

        runCmmd = "combine -M Asymptotic -m %03d -n _lim_%03d_%s_%s_ -d wwlvj_%s_lvjj_M%03d_%s_%s_unbin.txt --run both -H ProfileLikelihood"%(mass[i],mass[i],CHAN,CAT,SIG,mass[i],CHAN,CAT);
        print runCmmd
        os.system(runCmmd)

        time.sleep(0.1)

      elif CHAN == 'em' and options.category == 'ALLP':
	       	 	 
        cmd = 'combineCards.py '
	cmd += 'cards_mu_HPW/wwlvj_%s_lvjj_M%03d_mu_HPW_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_mu_HPZ/wwlvj_%s_lvjj_M%03d_mu_HPZ_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_mu_LPW/wwlvj_%s_lvjj_M%03d_mu_LPW_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_mu_LPZ/wwlvj_%s_lvjj_M%03d_mu_LPZ_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_el_HPW/wwlvj_%s_lvjj_M%03d_el_HPW_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_el_HPZ/wwlvj_%s_lvjj_M%03d_el_HPZ_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_el_LPW/wwlvj_%s_lvjj_M%03d_el_LPW_unbin.txt ' %(SIG,mass[i])
	cmd += 'cards_el_LPZ/wwlvj_%s_lvjj_M%03d_el_LPZ_unbin.txt ' %(SIG,mass[i])
	cmd += '> wwlvj_%s_lvjj_M%03d_combo_%s_unbin.txt'%(SIG,mass[i],CAT)
	
	print cmd
	os.system(cmd)

        runCmmd = "combine -M Asymptotic -m %03d -n _lim_%03d_combo_%s_ -d wwlvj_%s_lvjj_M%03d_combo_%s_unbin.txt --run both -H ProfileLikelihood"%(mass[i],mass[i],CAT,SIG,mass[i],CAT);
        print runCmmd
        os.system(runCmmd)

        time.sleep(0.1)
			
    ### Make the limit plots    
    if options.plotLimits:

       if CHAN != "em":    
         doLimitPlot("_%s_%s_"%(CHAN,CAT))
         doLimitPlot2("_%s_%s_"%(CHAN,CAT))
       elif CHAN == "em":	 
	 doLimitPlot("_combo_%s_"%(CAT))
	 doLimitPlot2("_combo_%s_"%(CAT))

    ### Make p-value plots
    if options.plotPvalue:

       if CHAN != "em":    
         doPvaluePlot("_%s_"%(CHAN))
       elif CHAN == "em":	 
	 doPvaluePlot("_combo_")

    #time.sleep(10)
