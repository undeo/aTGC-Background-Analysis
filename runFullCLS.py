import os,commands
import sys
import random
from array import array
from optparse import OptionParser
from optparse import OptionGroup
import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TH1D, TString, TPaveText, TGaxis
import subprocess
from subprocess import Popen
from optparse import OptionParser
import CMS_lumi, tdrstyle
import time

#$1: outfile
#$2: jobid
#$3: cmsswdir
#$4: datacard
#$5: singlepoint
#$6: seed
#$7: mass
#$8: suffix
#$9: seoutdir
#$10: toys
#$11: iterations
#$12: subid

### methods ###
def get_canvas(text):

   tdrstyle.setTDRStyle()
   CMS_lumi.lumi_13TeV = "2.2 fb^{-1}," + text
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
        
        print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected
        
        if t_quantileExpected == -1.: lims[0] = t_limit
        elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit
        elif t_quantileExpected >= 0.15  and t_quantileExpected <= 0.17:  lims[2] = t_limit
        elif t_quantileExpected == 0.5: lims[3] = t_limit
        elif t_quantileExpected >= 0.83  and t_quantileExpected <= 0.85:  lims[4] = t_limit
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit
        else: print "Unknown quantile!"
    
    return lims

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

    nPoints = len(mass)
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
        curFile = "m%i-fullCLS-%s/higgsCombine_%s_.HybridNew.mH%i.all.root"%(mass[i],signal_,signal_,mass[i])
        curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
        xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3]*0.01/br_lvjj );
        ybins_obs.append( curAsymLimits[0]*0.01/br_lvjj );
        ybins_2s.append( curAsymLimits[1]*0.01/br_lvjj );
        ybins_1s.append( curAsymLimits[2]*0.01/br_lvjj );
        ybins_xs.append(xsDict[mass[i]]);
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "m%i-fullCLS-%s/higgsCombine_%s_.HybridNew.mH%i.all.root"%(mass[i],signal_,signal_,mass[i])
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
       
parser = OptionParser()
parser.add_option("--outdir", dest="outdir", action="store",
                  help="SE output subdir (default = None)", metavar="OUTDIR", type="string", 
                  default="")
parser.add_option("--cmsswdir", dest="cmsswdir", action="store",
                  help="cmssw dir", metavar="CMSSWDIR", type="string", 
                  default="")
parser.add_option("--datadir", dest="datadir", action="store",
                  help="datacard dir", metavar="DATADIR", type="string", 
                  default="")
parser.add_option("--channel", dest="channel", action="store",
                  help="channel (el,mu,combo)", metavar="CHANNEL", type="string", 
                  default="")
parser.add_option('--signal', action="store",type="string",dest="signal",default="BulkG_WW")
parser.add_option("--toys", dest="toys", action="store",
                  help="number of toys", metavar="TOYS", type="int", 
                  default=1)
parser.add_option("--i", dest="iterations", action="store",
                  help="number of iterations", metavar="TOYS", type="int", 
                  default=3)
parser.add_option("--merge", dest="mergefiles", action="store_true",
                  help="merge job outputs (default = False)", metavar="MERGEFILES", 
                  default=False)
parser.add_option("--runtoys", dest="runtoys", action="store_true",
                  help="run toys (default = False)", metavar="MERGEFILES", 
                  default=False)	
parser.add_option("--runfullcls", dest="runfullcls", action="store_true",
                  help="run full cls (default = False)", metavar="MERGEFILES", 
                  default=False)
parser.add_option("--plotlimits", dest="plotlimits", action="store_true",
                  help="plot limits (default = False)", metavar="PLOTLIMITS", 
                  default=False)
		  		  		  	  		  		  		  
(options, args) = parser.parse_args()

cmsswdir_ = options.cmsswdir
outdir_ = options.outdir
datadir_ = options.datadir
channel_ = options.channel
toys_ = options.toys
it_ = options.iterations
signal_ = options.signal

points = []		  
for p in range(1,10):
   points+=[float(p/10.)]
   points+=[float(p/10.+0.05)]
   points+=[float(p/1.)]
   points+=[float(p/1.+0.5)]
   points+=[float(p*10.)]
   points+=[float(p*10.+5.)]

#mass = [800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500]
#mass = [1000,2000,3000,4000]
mass = [1000,2000,3000]

xsBulkG = {800:0.076058293365,1000:0.0204996927419,1200:0.00683956014339,1400:0.00261374040791,1600:0.0011483744,1800:0.0004882930,
           2000:0.000239518814875,2500:4.48570778809e-05,3000:9.82486087113e-06,3500:41.9758278604e-07,4000:24.3824379692e-07,
     	   4500:20.8967256514e-07}

xsWprime = {800:1.587885*0.428731, 1000:0.986533*0.460359, 1200:0.535394*0.467605, 1400:0.2955239*0.470248, 1600:0.1681478*0.471413,
            1800:0.0984325*0.471991, 2000:0.058998*0.472301, 2500:0.01771031*0.472619, 3000:0.00567529*0.472715, 3500:0.001878491*0.472748,
	    4000:6.0983956656e-08*0.472759,4500:168.1479623273e-07}
	  
if options.runtoys:

   for m in mass:

      for p in range(len(points)):

         point = points[p]
	 for i in range(1): #range(10)
            seed = int(random.random()*pow(10,7))
            outfile = "higgsCombine_%s_.HybridNew.mH%i.%i.root" %(signal_,m,seed)
            datacard = "%s/wwlvj_%s_lvjj_M%i_combo_ALLP_unbin.txt" %(datadir_,signal_,m)
            cmd = "qsub -q all.q submitJobsOnT3batch.sh %s %s %s %s %f %i %i %s %s %i %i %i" %(outfile,p,cmsswdir_,datacard,point,seed,m,signal_,outdir_,toys_,it_,i)
            print cmd
            os.system(cmd)

if options.mergefiles:

   user = os.popen('whoami').read()
   user = user.split('\n')[0]
   
   for m in mass: 
        
      status,ls_la = commands.getstatusoutput( 'ls -l m%i-fullCLS-%s'%(m,signal_) )  													      
      if status:																				      
         os.system('mkdir m%i-fullCLS-%s'%(m,signal_) )
      #else:
         #os.system('rm -rf m%i-fullCLS'%(m) )
	 #os.system('mkdir m%i-fullCLS'%(m) )	 
      
      for p in range(len(points)):#xrange(45,len(points)):
      
         for i in range(1): #range(10)
            outdir = "jobtmp"
            cmd = "uberftp t3se01.psi.ch 'ls /pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i-%s'" %(user,outdir_,m,p,i,signal_)
            status,ls_la = commands.getstatusoutput( cmd )
            filename = ''
            list_ = ls_la.split(os.linesep)

            for a in list_:
              b = a.split(" ")
              if b[-1:][0].find("root") != -1:
            	 filename = b[len(b)-1].split('\r')[0]
       
            cmd = "lcg-cp -bD srmv2 srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i-%s/%s m%i-fullCLS-%s/%s" %(user,outdir_,m,p,i,signal_,filename,m,signal_,filename) 
            print cmd
	    os.system(cmd)
	    cmd = "gfal-rm srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i-%s/%s" %(user,outdir_,m,p,i,signal_,filename)
	    print cmd
	    os.system(cmd)
	    cmd = "lcg-del -d srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i-%s/" %(user,outdir_,m,p,i,signal_)
            print cmd
	    os.system(cmd)
	    
         hadd = "hadd -f m%i-fullCLS-%s/grid_mX%i_p%i_%s.root m%i-fullCLS-%s/higgsCombine*" %(m,signal_,m,p,signal_,m,signal_)
         print hadd
         os.system(hadd)
	 rm = "rm m%i-fullCLS-%s/higgsCombine*" %(m,signal_)
         print rm
	 os.system(rm)
	    
      hadd = "hadd -f m%i-fullCLS-%s/grid_mX%i_%s.root m%i-fullCLS-%s/grid_mX*" %(m,signal_,m,signal_,m,signal_)
      print hadd
      os.system(hadd)

if options.runfullcls:

   for m in mass:
   
      datacard = "%s/wwlvj_%s_lvjj_M%i_combo_ALLP_unbin.txt" %(datadir_,signal_,m)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS-%s/grid_mX%i_%s.root -m %i -n _%s_" %(datacard,m,signal_,m,signal_,m,signal_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS-%s/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.5" %(datacard,m,signal_,m,signal_,m,signal_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS-%s/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.16" %(datacard,m,signal_,m,signal_,m,signal_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS-%s/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.84" %(datacard,m,signal_,m,signal_,m,signal_)
      print cmd
      os.system(cmd)                  
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS-%s/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.025" %(datacard,m,signal_,m,signal_,m,signal_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS-%s/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.975" %(datacard,m,signal_,m,signal_,m,signal_)
      print cmd
      os.system(cmd)
      
      cmd = "mv higgsCombine* m%i-fullCLS-%s/."%(m,signal_)
      print cmd
      os.system(cmd)

      cmd = "hadd m%i-fullCLS-%s/higgsCombine_%s_.HybridNew.mH%i.all.root m%i-fullCLS-%s/higgsCombine*"%(m,signal_,signal_,m,m,signal_)
      print cmd
      os.system(cmd)
            
if options.plotlimits:

   doLimitPlot("_combo_ALLP_")      
      
#python runFullCLS.py --channel combo --signal BulkG_WW --cmsswdir /shome/jngadiub/EXOVVAnalysisRunII/CMSSW_7_1_5 --datadir /shome/jngadiub/EXOVVAnalysisRunII/CMSSW_7_1_5/src/HiggsAnalysis/CombinedLimit/13TeV_datacards_Spring15/cards_BulkG_unblind/ --outdir jobtmp --runtoys --toys 100 --i 30
#python runFullCLS.py --channel combo --outdir jobtmp --merge
#python runFullCLS.py --channel combo --datadir /shome/jngadiub/EXOVVAnalysisRunII/CMSSW_7_1_5/src/HiggsAnalysis/CombinedLimit/13TeV_datacards_Spring15/cards_BulkG_unblind/ --runfullcls
#python runFullCLS.py --channel combo --plotlimits
