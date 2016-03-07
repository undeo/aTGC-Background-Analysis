//c++ -o interpolate interpolate.cpp `root-config --cflags --glibs`
//to run:
//./interpolate BulkG el_HPW WW /afs/cern.ch/user/l/lbrianza/work/PHD/WW_ANALYSIS_RUN2/WWSEMILEP/ANALISI_NOVEMBER/TEST/CMSSW_7_1_5/src/13TeV_datacards_Spring15/ log/
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include <sstream>
#include "TMath.h"
#include <cmath>

using namespace std;

double DoubleCB(Double_t *x,Double_t *par)
{

  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2

  double xx = x[0];
  double mean   = par[2] ; // mean
  double sigma  = par[6] ; // sigma of the positive side of the gaussian
  // double sigmaN = par[3] ; // sigma of the negative side of the gaussian
  double alpha  = par[0] ; // junction point on the positive side of the gaussian
  double n      = par[3] ; // power of the power law on the positive side of the gaussian
  double alpha2 = par[1] ; // junction point on the negative side of the gaussian
  double n2     = par[4] ; // power of the power law on the negative side of the gaussian
  double N      = par[5] ;

  if ((xx-mean)/sigma > fabs(alpha))
    {
      double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
      double B = n/fabs(alpha) - fabs(alpha);
    
      return N * A * pow(B + (xx-mean)/sigma, -1.*n);
    }
  
  else if ((xx-mean)/sigma < -1.*fabs(alpha2))
    {
      double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
      double B = n2/fabs(alpha2) - fabs(alpha2);
    
      return N * A * pow(B - (xx-mean)/sigma, -1.*n2);
    }
  
  else if ((xx-mean) > 0)
    {
      return N * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
    }
  
  else
    {
      return N * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
    }
   
}


int main(int argc, char** argv)
{
  std::string signal =argv[1];
  std::string category =argv[2];
  std::string type =argv[3];
  std::string folderDatacard = argv[4];
  std::string folderLogFile = argv[5];

  string line;
  string line2;
  ifstream logFile[11];
  ifstream datacard[11];
  int mass[11] = {800,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
  bool goOut=false;
  int nPar = 8, nMass=11;
  string namePar[8];
  float val[8][11];
  float err[8][11];
  string parName[8]={"alpha1","alpha2","mean","n1","n2","number","sigma","Ndatacard"}; 
  //  string folderDatacard = "/afs/cern.ch/user/l/lbrianza/work/PHD/WW_ANALYSIS_RUN2/WWSEMILEP/ANALISI_NOVEMBER/TEST/CMSSW_7_1_5/src/13TeV_datacards_Spring15/";
  //  string folderLogFile = "log/";

  //  Char %20louer%20furnished%20studio%20in%20St%20Genis&FolderCTID=0x0120020014A4FDC150F34C48869500E594CCB052&SiteMapTitle=Housing&SiteMapUrl=https%3A%2F%2Fsocial%2Ecern%2Ech%2Fcommunity%2Fcern-market%2FSitePages%2FCategory%2Easpx%3FCategoryID%3D3%26SiteMapTitle%3DHousingsignal[200]="BulkG_WW";
  //  char category[200]="el_HPW";
  //  char signal[200]={signalString.c_str()};
  //char category[200]={categoryString.c_str()};

  for(int iMass=0; iMass<nMass; iMass++)
    {
      TString NameFile = Form("%s%s_%s_M%d_%s.log",folderLogFile.c_str(),signal.c_str(),type.c_str(),mass[iMass],category.c_str());
      logFile[iMass].open(NameFile.Data());
      if (logFile[iMass]) 
	{
	  goOut = false;
	  while (getline(logFile[iMass],line) && !goOut) 
	    {
	      if (TString(line).Contains("  --------------------  --------------------------")) 
		{
		  std::cout<<"file: "<<NameFile.Data()<<" content: "<<line<<std::endl;
		  for (int par=0; par<nPar-1; par++) 
		    {
		      string pm;
		      //		      getline(logFile[i],line);
		      //std::cout<<"file: "<<NameFile.Data()<<" content: "<<line<<std::endl;
		      //		      std::stringstream linestream(line);
		      //		      std::cout<<linestream<<std::endl;
		      logFile[iMass] >> namePar[par] >> val[par][iMass] >> pm >> err[par][iMass];		      
		      std::cout<<"par: "<<namePar[par]<<" val: "<<val[par][iMass]<<" err: "<<err[par][iMass]<<std::endl;
		    }		  
		  goOut = true;
		}	    
	    }
	}
      logFile[iMass].close();

      TString NameCard = Form("%scards_%s_unblind/cards_%s/wwlvj_%s_%s_lvjj_M%d_%s_unbin.txt",folderDatacard.c_str(),signal.c_str(),category.c_str(),signal.c_str(),type.c_str(),mass[iMass],category.c_str());
      std::cout<<"open datacard: "<<NameCard.Data()<<std::endl;

      datacard[iMass].open(NameCard.Data());
      if (datacard[iMass]) 
	{
	  goOut = false;
	  while (getline(datacard[iMass],line2) && !goOut) 
	    {
	      //	      std::cout<<line2<<std::endl;
	      if (TString(line2).Contains("process")) 
		{
		  string rate;
		  getline(datacard[iMass],line2);
		  std::cout<<"file: "<<NameFile.Data()<<" content: "<<line2<<std::endl;
		  datacard[iMass] >> rate >> val[nPar-1][iMass];
		  err[nPar-1][iMass]=0;		      
		  std::cout<<"par: "<<namePar[nPar-1]<<" val: "<<val[nPar-1][iMass]<<" err: "<<err[nPar-1][iMass]<<std::endl;
		  goOut = true;
		}	    
	    }
	}
      
    }

  TString NameOutput = Form("interpolationFiles/%s_%s_lvjj_%s.root",signal.c_str(),type.c_str(),category.c_str());
  TFile *outputFile = new TFile(NameOutput.Data(),"RECREATE");
  outputFile->cd();

  TGraphErrors *gPar[8];
  TF1 *fitPar[8];

  ofstream param;
  param.open((string("param_")+NameOutput.Data()+string(".txt")).c_str());

  for (int iPar=0; iPar<nPar; iPar++)
    {
      TString NameFile = Form("%s",parName[iPar].c_str());
      TString NameFunc = Form("fit_%s",parName[iPar].c_str());
      //gPar[iPar] = new TGraphErrors(NameFile,NameFile);
      gPar[iPar] = new TGraphErrors();
      //      fitPar[iPar] = new TF1("pol1","pol1",2,1000,3600);

      for (int iMass=0; iMass<nMass; iMass++)
	{
	  gPar[iPar]->SetPoint(iMass,mass[iMass],val[iPar][iMass]);
	  gPar[iPar]->SetPointError(iMass,0,err[iPar][iMass]);
	}
      /*
      param << parName[iPar]<<" : ";
      if (TString(parName[iPar]).Contains("number") || TString(parName[iPar]).Contains("n1") || TString(parName[iPar]).Contains("n2")) {
	gPar[iPar]->Fit("pol3");
	TF1 *func = gPar[iPar]->GetFunction("pol3");
	//	func->SetRange(0,3600);
	gPar[iPar]->Fit(func);
	if (func->GetParameter(1)>0 && func->GetParameter(2)>0 && func->GetParameter(3)>0 )
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x+"<<func->GetParameter(2)<<"*x*x+"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)>0 && func->GetParameter(2)>0 && func->GetParameter(3)<0 )
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x+"<<func->GetParameter(2)<<"*x*x"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)>0 && func->GetParameter(2)<0 && func->GetParameter(3)>0 )
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x"<<func->GetParameter(2)<<"*x*x+"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)>0 && func->GetParameter(2)<0 && func->GetParameter(3)<0 )
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x"<<func->GetParameter(2)<<"*x*x"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)<0 && func->GetParameter(2)>0 && func->GetParameter(3)>0 )
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x+"<<func->GetParameter(2)<<"*x*x+"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)<0 && func->GetParameter(2)>0 && func->GetParameter(3)<0 )
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x+"<<func->GetParameter(2)<<"*x*x"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)<0 && func->GetParameter(2)<0 && func->GetParameter(3)>0 )
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x"<<func->GetParameter(2)<<"*x*x+"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
	else if (func->GetParameter(1)<0 && func->GetParameter(2)<0 && func->GetParameter(3)<0 )
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x"<<func->GetParameter(2)<<"*x*x"<<func->GetParameter(3)<<"*x*x*x"<<std::endl;
      }
      else if (TString(parName[iPar]).Contains("n1") || TString(parName[iPar]).Contains("n2")) {
	gPar[iPar]->Fit("pol2");
	TF1 *func = gPar[iPar]->GetFunction("pol2");
	//	func->SetRange(0,3600);
	gPar[iPar]->Fit(func);
	if (func->GetParameter(1)>0 && func->GetParameter(2)>0 )
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x+"<<func->GetParameter(2)<<"*x*x"<<std::endl;
	else if (func->GetParameter(1)>0 && func->GetParameter(2)<0 )
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x"<<func->GetParameter(2)<<"*x*x"<<std::endl;
	else if (func->GetParameter(1)<0 && func->GetParameter(2)>0 )
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x+"<<func->GetParameter(2)<<"*x*x"<<std::endl;
	else if (func->GetParameter(1)<0 && func->GetParameter(2)<0 )
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x"<<func->GetParameter(2)<<"*x*x"<<std::endl;
      }
      else {
	gPar[iPar]->Fit("pol1");
	TF1 *func = gPar[iPar]->GetFunction("pol1");
	//func->SetRange(0,3600);
	gPar[iPar]->Fit(func);
	if (func->GetParameter(1)>0)
	  param<<func->GetParameter(0)<<"+"<<func->GetParameter(1)<<"*x"<<std::endl;
	else
	  param<<func->GetParameter(0)<<func->GetParameter(1)<<"*x"<<std::endl;
      }
      */
      gPar[iPar]->GetXaxis()->SetTitle("mass (GeV)");
      gPar[iPar]->GetYaxis()->SetTitle(parName[iPar].c_str());
      gPar[iPar]->Write(NameFile.Data());
    }

  outputFile->Close();
  
  
  return 0;
}
