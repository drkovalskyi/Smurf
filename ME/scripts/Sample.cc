#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"

#include <vector>
#include <math.h>
#include <iostream>

// Matrix Element related
#include "../TVar.hh"
#include "Sample.h"

Sample::Sample(TVar::Process LRProcess, TString datasample, TChain *chain, Color_t color, bool stack, float lumi) {
  // initialize data members
  LRProcess_ = LRProcess;
  datasample_ = datasample;
  chain_ = chain;
  color_ = color;
  stack_ = stack;
  lumi_ = lumi;
}

Sample::~Sample() {
}

TH1F * Sample::GetHistogram(TString var, Int_t nBins, Float_t binMin, Float_t binMax) {
  
  TH1F *tmp = new TH1F(Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), nBins, binMin, binMax);
  // VV contains only the non-peaking 
  // Use the lepton mother id to distinguish lep1McId and lep2McId 
  // This variable is 11 for leptons from W and 21 for leptons from Z (and 0 for fakes)
  if (datasample_ == "VV") {
    chain_->Project(Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), var,  "scale1fb*(lep1McId!=lep2McId)");
  }
  
  else if (datasample_ == "Zjets") {
    chain_->Project(Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), var, "scale1fb*(lep1McId==lep2McId&&lep1McId==21)");
  }
  
  else
    chain_->Project(Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), var, "scale1fb");
  
  tmp->Scale(lumi_/1000.0);
  
  if (stack_) {
    tmp->SetFillColor(color_);
    tmp->SetLineColor(color_);
  }
  else {
    tmp->SetLineColor(color_);
    tmp->SetMarkerColor(color_);
    tmp->SetLineWidth(3);
  }
  return tmp;
}

  

