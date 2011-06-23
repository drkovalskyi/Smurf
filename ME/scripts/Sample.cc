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
#include "TCut.h"

#include <vector>
#include <math.h>
#include <iostream>

// Matrix Element related
#include "../TVar.hh"
#include "Sample.h"
#include "../../Core/SmurfTree.h"


Sample::Sample(TVar::Process LRProcess, TString datasample, TChain *chain, Color_t color, bool stack, float lumi) {
  // initialize data members
  LRProcess_ = LRProcess;
  datasample_ = datasample;
  chain_ = chain;
  color_ = color;
  stack_ = stack;
  lumi_ = lumi;
  sf_ = 1.0;
}

Sample::~Sample() {
}

TH1F * Sample::GetHistogram(TString var, Int_t nBins, Float_t binMin, Float_t binMax) {
  
  TH1F *tmp = new TH1F(Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), nBins, binMin, binMax);
  tmp->Sumw2();

  TCut cut; 
  // use only the wz/zz events where the two leptons are from different mothers....
  
  if (datasample_ == "Data") 
    cut = TCut("mycut", "scaleME");
  // Very bad hard coded stuff...fixme
  // 49: wz 50 :zz
  // requiring the two leptons from the same mother if the Zjets samples are wz or zz
  else if (datasample_ == "VV")
    cut = TCut("mycut", Form("%f*scale1fb*scaleME*(lep1MotherMcId!=lep2MotherMcId)", sf_*lumi_/1000.0));
  else if (datasample_ == "Zjets") 
    cut = TCut("mycut", Form("%f*scale1fb*scaleME*((dstype!=%i&& dstype!=%i)||(lep1MotherMcId==lep2MotherMcId&&lep1MotherMcId==23))",sf_*lumi_/1000.0, 49, 50));
  else
    cut = TCut("mycut", Form("%f*scale1fb*scaleME", sf_*lumi_/1000.0));
  chain_->Project(Form("hist_%s_%s", TVar::SmurfProcessName(LRProcess_).Data(), datasample_.Data()), var, cut);
  
  /*
    // for future developement
  float scale1fb_ = 0;
  double scaleME_ = 0;
  double dPhi_ = 0;
  int  lep1MotherMcId_;                                                                                        
  int  lep2MotherMcId_;   
  SmurfTree::DataType dstype_; 
  double var_;
  
  chain_->SetBranchAddress( "scale1fb", &scale1fb_); 
  chain_->SetBranchAddress( "scaleME", &scaleME_); 
  chain_->SetBranchAddress( "dPhi", &dPhi_); 
  chain_->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId_); 
  chain_->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId_); 
  chain_->SetBranchAddress( "dstype", &dstype_); 
  chain_->SetBranchAddress( var.Data(), &var_);
  
  for(int ievt=0; ievt < chain_->GetEntries();ievt++){
    chain_->GetEntry(ievt);
    // skipping the events where both leptons come from the same parent
    if (datasample_ == "VV" && (lep1MotherMcId_ == lep2MotherMcId_)  ) continue;
    if (datasample_ == "Zjets") {
      if (dstype_ == SmurfTree::wz || dstype_ == SmurfTree::zz) {
	if (lep1MotherMcId_ != lep2MotherMcId_ || lep1MotherMcId_ != 23) 
	  continue;
      }
    }
    tmp->Fill(var_, scale1fb_*scaleME_*lumi_/1000.0);
  }

  */
    
  if (stack_) {
    tmp->SetFillColor(color_);
    tmp->SetLineColor(color_);
    tmp->SetMarkerColor(color_);
  }
  else {
    tmp->SetLineColor(color_);
    tmp->SetMarkerColor(color_);
    tmp->SetLineWidth(3);
  }
  return tmp;
}

  

