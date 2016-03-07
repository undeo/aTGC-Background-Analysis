import ROOT
from ROOT import *
import time
import sys
from array import array
import tdrstyle

#BulkG
#xsBulkG = {800:0.0030423317346004058,1000:0.00081998770967753234,1200:0.00027358240573576689,1400:0.00010454961631658292,
#           2000:9.5807525950148653e-06,2500:1.7942831152359325e-06,3000:3.9299443484538374e-07}
#Wprime
xs = {800:1.587885*0.428731, 1000:0.986533*0.460359, 1200:0.535394*0.467605, 1400:0.2955239*0.470248, 1600:0.1681478*0.471413,
            1800:0.0984325*0.471991, 2000:0.058998*0.472301, 2500:0.01771031*0.472619, 3000:0.00567529*0.472715, 3500:0.001878491*0.472748}

x = array('d',[])
y = array('d',[])	   
for m,xs in xs.iteritems():
   x.append(m)
   y.append(TMath.Log10(xs*0.5*0.5/(0.1*0.1)))
   print m,xs*0.5*0.5/(0.1*0.1)	   
   
gr = ROOT.TGraph(len(x),x,y)
gr.SetMarkerColor(kRed)
gr.SetMarkerStyle(20)

func =  TF1("func","pol2",800,5000)
#func =  TF1("func","expo(0)",800,5000)
#func =  TF1("func","[0]/(x+[1])",800,5000)
func.SetLineColor(kBlue)

c = TCanvas()
c.cd()
gr.Draw('AP') 
gr.Fit(func,'R')

print "m = 1600 xsec = %.10f" %(TMath.Power(10,func.Eval(1600)))
print "m = 1800 xsec = %.10f" %(TMath.Power(10,func.Eval(1800)))  
print "m = 3500 xsec = %.10f" %(TMath.Power(10,func.Eval(3500))*10000000.)  
print "m = 4000 xsec = %.10f" %(TMath.Power(10,func.Eval(4000))*10000000.) 
print "m = 4500 xsec = %.10f" %(TMath.Power(10,func.Eval(4500))*10000000.) 

print func.Eval(1600)
print func.Eval(1800)
print func.Eval(3500)
print func.Eval(4000)
print func.Eval(4500)
x1 = array('d',[1600,1800,3500,4000,4500])
y1 = array('d',[func.Eval(1600),func.Eval(1800),func.Eval(3500),func.Eval(4000),func.Eval(4500)])

#c1 = TCanvas()
#c1.cd()
gr1 = TGraph(len(x1),x1,y1)
gr1.SetMarkerStyle(20)
gr1.SetMarkerColor(kBlack)
gr1.Draw('Psame')
gr1.GetXaxis().SetRangeUser(800,5000)

time.sleep(1000)
