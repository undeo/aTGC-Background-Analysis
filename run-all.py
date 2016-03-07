import os,commands
import sys
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('-s', action="store",type="string",dest="sample",default="BulkG_WW")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--jetalgo', action="store",type="string",dest="jetalgo",default="Mjpruned")
parser.add_option('--interpolate', action="store_true",dest="interpolate",default=False)
parser.add_option('--batchMode', action="store_true",dest="batchMode",default=False)
(options, args) = parser.parse_args()

currentDir = os.getcwd();

#masses = [800,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500]
masses = [800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
          3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500] #use this one only for the interpolation mode!
#masses = [3000,3500,4000,4500]
for m in masses:
   if (options.interpolate==True and options.batchMode==True):

      fn = "Job/Job_%s_%s_%d.sh"%(options.channel,options.category,m)
      outScript = open(fn+".sh","w");
 
      outScript.write('#!/bin/bash');
      outScript.write("\n"+'cd '+currentDir);
      outScript.write("\n"+'eval `scram runtime -sh`');
      outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
      outScript.write("\n"+'echo ${PATH}');
      outScript.write("\n"+'ls');
#      cmd = "python g1_exo_doFit_class.py -b -c %s --mass %i --category %s --sample %s_lvjj --jetalgo %s --interpolate True > log/%s_M%i_%s_%s.log" %(options.channel,m,options.category,options.sample,options.jetalgo,options.sample,m,options.channel,options.category)
      cmd = "python g1_exo_doFit_class.py -b -c %s --mass %i --category %s --sample %s_lvjj --jetalgo %s --interpolate True" %(options.channel,m,options.category,options.sample,options.jetalgo)
      outScript.write("\n"+cmd);
#      outScript.write("\n"+'rm *.out');
      outScript.close();

      os.system("chmod 777 "+currentDir+"/"+fn+".sh");
      os.system("bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh");

   elif (options.interpolate==True and not options.batchMode==True):
      cmd = "python g1_exo_doFit_class.py -b -c %s --mass %i --category %s --sample %s_lvjj --jetalgo %s --interpolate True > log/%s_M%i_%s_%s.log" %(options.channel,m,options.category,options.sample,options.jetalgo,options.sample,m,options.channel,options.category)
      print cmd
      os.system(cmd)

   else:   
      cmd = "python g1_exo_doFit_class.py -b -c %s --mass %i --category %s --sample %s_lvjj --jetalgo %s > log/%s_M%i_%s_%s.log" %(options.channel,m,options.category,options.sample,options.jetalgo,options.sample,m,options.channel,options.category)
      print cmd
      os.system(cmd)

#python run-all.py --channel mu -s Wprime_WZ --jetalgo Mjsoftdrop --category HP
#python run-all.py -c mu -s BulkG_WW --category HPW
