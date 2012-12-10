//
// This scripts looks up the correlation histograms returned by combine
// and output the large correlation sources > 20%
// 

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "TRint.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCut.h"
#include "THStack.h"


void drawcorr(TString fileName = "mlfit.root", TString fitName="fit_s", double threshold = 0.2){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();
  
  
  TFile* file = new TFile(fileName, "READ");
  TH2F *corr; 
  RooFitResult *fit = (RooFitResult*)file->Get(fitName);
  corr = (TH2F*) fit->correlationHist();
  int nBinsX = corr->GetXaxis()->GetNbins();
  int nBinsY = corr->GetYaxis()->GetNbins();


  for ( int binx = 1; binx <= nBinsX; binx++) 
    for ( int biny = 1; biny <= nBinsY+1-binx; biny++) {
      
      double bincontent = corr->GetBinContent(binx, biny); 
      
      if ( ( binx + biny) == (nBinsX+1) ) continue;

      if ( TMath::Abs(bincontent) > threshold ) {
	std::cout << Form("Correlation between %s and %s [%i][%i] is %.4f\n", 
			  corr->GetXaxis()->GetBinLabel(binx), 
			  corr->GetYaxis()->GetBinLabel(biny), binx, biny, bincontent); 
      }
    }
  
}
