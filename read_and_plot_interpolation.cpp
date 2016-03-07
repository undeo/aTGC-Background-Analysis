//c++ -o read_and_plot read_and_plot_interpolation.cpp `root-config --cflags --glibs`
// ./read_and_plot_interpolation BulkG el_HPW WW interpolationFiles/
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
  double N      = par[7] ;

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
  std::string inputFolder = argv[4];

  string line;
  ifstream logFile[11];
  bool goOut=false;
  int nPar = 8;
  string namePar[8];
  Double_t val[8][11];
  Double_t err[8][11];
  string parName[8]={"alpha1","alpha2","mean","n1","n2","number","sigma","Ndatacard"}; 

  TString *readFile = new TString (Form("%s%s_%s_lvjj_%s.root",inputFolder.c_str(),signal.c_str(),type.c_str(),category.c_str()));
  TFile* inputFile = new TFile(readFile->Data());

  std::cout<<"open: "<<readFile->Data()<<std::endl;
  //  inputFile->cd();

  TGraph *gPar[8];
  Double_t paramVec[8];
  TF1* fTest[100];

  for (int iPar=0; iPar<nPar; iPar++)
    {
      TString NameFile = Form("%s",parName[iPar].c_str());
      gPar[iPar] = (TGraph*)inputFile->Get(NameFile.Data());
    }

  int i=0;
  for (int mass=800; mass<5001; mass+=200,i++) {

    TString NameFunc = Form("f_%d",mass);
    fTest[i] = new TF1 (NameFunc.Data(), DoubleCB, 0, 10000, 8) ;

    for (int iPar=0; iPar<nPar; iPar++) {
      paramVec[iPar] = gPar[iPar]->Eval(mass);
      std::cout<<"mass: "<<mass<<" par: "<<paramVec[iPar]<<std::endl;
    }
    //    if (mass==4200) std::cout<<paramVec[8]<<std::endl;
    fTest[i]->SetParameters (paramVec) ;
    fTest[i]->SetLineColor (kBlue-9+(i%5));
    //fTest[i]->SetLineColor (kGreen);
  }

  inputFile->Close();

  TCanvas *c = new TCanvas (Form("%s_%s_%s",signal.c_str(),category.c_str(),type.c_str()) ,"",500,525);
  c->cd();
  c->DrawFrame(0,3,3,5000);

  i=0;
  for (int mass=800; mass<5001; mass+=200,i++) {
    //    if (mass!=4200) continue;
    //else fTest[i]->Draw();
    if (i==0)
      fTest[i]->Draw();
    else
      fTest[i]->Draw("same");
    
  }

  TString NameOutput = Form("extrapolation_%s_%s_%s.root",signal.c_str(),category.c_str(),type.c_str());
  TFile *outputFile = new TFile(NameOutput.Data(),"RECREATE");
  outputFile->cd();

  c->Write();
    
  
  return 0;
}
