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


void drawcorr(){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();
  
  
  TFile* file = new TFile("mlfit_18fb.root", "READ");
  TH2F *corr; 
  RooFitResult *fit = (RooFitResult*)file->Get("fit_s");
  corr = (TH2F*) fit->correlationHist();
  int nBinsX = corr->GetXaxis()->GetNbins();
  int nBinsY = corr->GetYaxis()->GetNbins();


  for ( int binx = 1; binx <= nBinsX; binx++) 
    for ( int biny = 1; biny <= nBinsY+1-binx; biny++) {
      
      double bincontent = corr->GetBinContent(binx, biny); 
      
      if ( ( binx + biny) == (nBinsX+1) ) continue;

      if ( TMath::Abs(bincontent) > 0.2 ) {
	std::cout << Form("Correlation between %s and %s [%i][%i] is %.4f\n", 
			  corr->GetXaxis()->GetBinLabel(binx), 
			  corr->GetYaxis()->GetBinLabel(biny), binx, biny, bincontent); 
      }
    }
  
}
