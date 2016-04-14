from ROOT import *
from array import array
import math
from struct import *



def merge_trees(list_):    

  mergedtree = TChain("BasicTree")
 
  for file_ in list_:    
    mergedtree.Add(file_)

  return mergedtree
  


def prepare_trees_with_cuts(filenameOut, list_of_trees, channel, is_data = False):
  
  if is_data:
    oldtree = TChain("treeDumper/BasicTree")
    for file_ in list_of_trees:    
      oldtree.Add(file_)
  else:
    oldtree = merge_trees(list_of_trees)
  
  
  nEntries = oldtree.GetEntries()
  newtree = oldtree.CloneTree(0)  
  
  ## create new branches
  if is_data:
    var_weight = array('f',[1])
    branch_weight = newtree.Branch("totEventWeight2", var_weight, "weight")
    var_weight[0] = 1
  
  ##fill new branches  
  used_events = 0
  print list_of_trees
  print "Number of Entries: " + str(nEntries) + " , " "Number of used events: "
  for i in range(nEntries):    
  	use_event = False
  	oldtree.GetEntry(i)
    	##apply cuts
  	if oldtree.jet_pt>200 and oldtree.jet_tau2tau1<0.6 and oldtree.Mjpruned<150\
	and oldtree.Mjpruned>40 and oldtree.W_pt>200 and abs(oldtree.deltaR_LeptonWJet)>math.pi/2\
	and abs(oldtree.deltaPhi_WJetMet)>2 and abs(oldtree.deltaPhi_WJetWlep)>2 and oldtree.nbtag==0:
      		use_event = True

	
	#if ch == 'ele' and not unpack('?',oldtree.bit_HLT_Ele_105)[0] :
		#use_event = False

    	if use_event==True:
      		used_events += 1
      		if used_events%1000 == 0:
			print str(used_events) + " / " + str(i) + " / " + str(nEntries)    
		newtree.Fill()

  print str(used_events) + " / " + str(nEntries)
  print "____"
  
  ##set name to 'tree'
  newtree.SetName("tree")
  basepath = '../AnaSigTree_25ns/%s/'%channel[:2]
  fileOut = TFile(basepath + filenameOut, "recreate")
  newtree.Write()
  fileOut.Close()



chan = ["ele", "mu"]
#chan = ['ele']
for ch in chan:
  ###merge and rename trees
  data_list_of_trees = ["data-RunD_%s.root"%(ch)]
  prepare_trees_with_cuts("treeEDBR_data_xww.root", data_list_of_trees,ch, True)

  W_list_of_trees = ["WJets_Ht100To200_%s.root"%(ch), "WJets_Ht200To400_%s.root"%(ch), "WJets_Ht400To600_%s.root"%(ch),\
			 "WJets_Ht600To800_%s.root"%(ch), "WJets_Ht800To1200_%s.root"%(ch), "WJets_Ht1200To2500_%s.root"%(ch),\
			 "WJets_Ht2500ToInf_%s.root"%(ch)]
  prepare_trees_with_cuts("treeEDBR_WJets_xww.root", W_list_of_trees, ch)

  T_list_of_trees = ["s-ch_%s.root"%(ch), "tW-ch-antitop_%s.root"%(ch), "tW-ch-top_%s.root"%(ch), "t-ch_%s.root"%(ch)]
  prepare_trees_with_cuts("treeEDBR_SingleTop_xww.root", T_list_of_trees, ch)

  VV_list_of_trees = ["WW_%s.root"%(ch), "WZ_%s.root"%(ch)]
  prepare_trees_with_cuts("treeEDBR_VV_xww.root", VV_list_of_trees, ch)

  TTBAR_list_of_trees = ["ttbar_%s.root"%(ch)]
  prepare_trees_with_cuts("treeEDBR_TTBARpowheg_xww.root", TTBAR_list_of_trees, ch)

