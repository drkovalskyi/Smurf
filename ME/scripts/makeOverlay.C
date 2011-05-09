//
// Plotting script to overlay various Higgs signal with the SM background
// D. Evans and Y. Gao
// ==================
// Instructions
// 1. Create a soft link called data where the full statistics smurfNtuples live
// 2. By default the smurfNtuples containing the ME output is the current directory
// 3. Specify the higgs samples to consider in the main function makeOverlay()
// 4. Do root -l makeOverlay.C+
// 

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

const unsigned int kNDilep = 4;


class Sample {

public:
  
  Sample(TVar::Process process, Color_t color, bool stack, TString smurfFDir, TString outputDir, TString meFSuffix, Float_t massCut) { 
    
    // initialize data members
    process_ = process;
    stack_ = stack;
    color_ = color;
    file_me_ = TFile::Open(outputDir+TVar::SmurfProcessName(process)+meFSuffix);
    if (file_me_ == 0x0) {
      std::cout << "Sample (" << outputDir+TVar::SmurfProcessName(process_) << ") ME file is BAD" << std::endl;
    }
    else {
      tree_me_ = (TTree*)file_me_->Get("tree");
    }
    
    // get the yields for the normalisations
    InitializeYields(smurfFDir, massCut);
    
  }
  ~Sample() {
    file_me_->Close();
  };
  
 
  TVar::Process GetProcess() {
    return process_;
  }
  
  void SetYield(float yield, unsigned int j) { 
    yield_[j] = yield; 
  }
  const float &GetYield(unsigned int j) { 
    return yield_[j]; 
  }
  float GetTotalYield() {
    float total = 0.0;
    for (unsigned int j = 0; j < kNDilep; ++j) total += yield_[j];
    return total;
  }
  
  TH1F *GetHistogram(TVar::Process processLR, Int_t nBins, Float_t binMin, Float_t binMax) {
    TH1F *tmp = new TH1F(Form("h1_%sLR_%s",  TVar::SmurfProcessName(processLR).Data(), TVar::SmurfProcessName(GetProcess()).Data()), "LR", nBins, binMin, binMax);
    tree_me_->Project(Form("h1_%sLR_%s", TVar::SmurfProcessName(processLR).Data(), TVar::SmurfProcessName(GetProcess()).Data()), Form("LR[%i]", processLR));
    tmp->Scale(GetTotalYield()/tmp->Integral(0, nBins+1));
    if (stack_) {
      tmp->SetFillColor(color_);
      tmp->SetLineColor(color_);
      // tmp->SetFillStyle(1);
    }
    else {
      tmp->SetLineColor(color_);
      tmp->SetMarkerColor(color_);
      tmp->SetLineWidth(3);
    }
    return tmp;
  }
  
  bool Stack() { return stack_; }
  
private:
  
  //
  // member functions
  //
  
  void InitializeYields(TString inputDir, Float_t massCut) 
  {
    
    TFile *file_norm = TFile::Open(inputDir+TVar::SmurfProcessName(GetProcess())+".root");
    TTree *tree_norm = (TTree*)file_norm->Get("tree");
    
    if (file_norm == 0x0) {
      std::cout << "Sample (" << TVar::SmurfProcessName(GetProcess()) << ") NORM file is bad" << std::endl;
      return;
    }
    
    gROOT->cd();
    for (unsigned int j = 0; j < kNDilep; j++)
      {
	TH1F *tmp = new TH1F("tmp", "tmp", 20, 0, 4);

	
	if(j==0 || j==2)
	  tree_norm->Project("tmp", "dPhi", Form("scale1fb*(type==%i&&dilep.mass()<%f)", j, massCut));
	if(j==1 || j==3)
	  tree_norm->Project("tmp", "dPhi", Form("scale1fb*((type==%i&&lep2.pt()>15)&&dilep.mass()<%f)", j, massCut));
	if (tmp != 0x0) yield_[j] = tmp->Integral(0, 9999);
	else yield_[j] = 0.0;
	delete tmp;
      }
    
    file_norm->Close();
    
  }
  
  //
  // member data
  //
  
  TFile *file_me_;
  TTree *tree_me_;
  TVar::Process process_;
  bool stack_;
  Color_t color_;
  float yield_[kNDilep];
};

// 
// End of Definition of HiggsProcess Class
// 

// 
// draw the overlay of several LR distributions for signal and bkgs
// 

void drawLROverlay(Sample* & higgsSample, std::vector<Sample*> & bkgsamples, 
		   Int_t nBins, Float_t binMin, Float_t binMax, TString outputDir)
{   
  TH1F *hist_sig = (TH1F*)higgsSample->GetHistogram(higgsSample->GetProcess(), nBins, binMin, binMax);
  
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(size_t(bkgsamples.size()));
  
  TLegend *stacklg = new TLegend(0.62, 0.6, 0.88, 0.9);
  stacklg->SetBorderSize(0);
  stacklg->SetFillStyle(0);
  stacklg->SetShadowColor(0);

  TLegend *overlaylg = new TLegend(0.62, 0.6, 0.88, 0.9);
  overlaylg->SetBorderSize(0);
  overlaylg->SetFillStyle(0);
  overlaylg->SetShadowColor(0);
  
  stacklg->AddEntry(hist_sig, TVar::SmurfProcessName(higgsSample->GetProcess()), "l");
  overlaylg->AddEntry(hist_sig,TVar::SmurfProcessName(higgsSample->GetProcess()), "l");

  
  THStack *h1_stack = new THStack(); 
  for (unsigned int s = 0; s < bkgsamples.size(); ++s) {
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(higgsSample->GetProcess(), nBins, binMin, binMax));
    h1_stack->Add(hist_bkg[s]);
    stacklg->AddEntry(hist_bkg[s], TVar::SmurfProcessName(bkgsamples[s]->GetProcess()), "f");
    overlaylg->AddEntry(hist_bkg[s], TVar::SmurfProcessName(bkgsamples[s]->GetProcess()), "l");
  }
  
  // draw stacked plots
  TCanvas *c1 = new TCanvas();
  c1->cd();
  c1->SetLogy(0);
  h1_stack->Draw("HIST");
  h1_stack->GetXaxis()->SetTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsSample->GetProcess()).Data()));
  hist_sig->Draw("SAMEHIST");
  stacklg->Draw();
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR.eps");
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR.png");
  
  // draw stacked plots - log
  c1->SetLogy();
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR_log.eps");
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR_log.png");

  // draw overlay plots
  c1->Clear();
  c1->SetLogy(0);
  hist_sig->Draw("HIST");
  hist_sig->GetXaxis()->SetTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsSample->GetProcess()).Data()));
  float yMax = hist_sig->GetMaximum();
  for (unsigned int i = 0; i < hist_bkg.size(); i++) {
    hist_bkg[i]->SetLineWidth(3);
    hist_bkg[i]->Draw("SAMEHIST");
    hist_bkg[i]->SetFillStyle(0);
    if (hist_bkg[i]->GetMaximum() > yMax) yMax = hist_bkg[i]->GetMaximum();
  }
  hist_sig->GetYaxis()->SetRangeUser(0.001, yMax  + 2*sqrt(yMax));
  overlaylg->Draw();
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR_overlay.eps");
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR_overlay.png");
  
  // draw overlay plots - log
  c1->SetLogy(1);
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR_overlay_log.eps");
  c1->SaveAs( outputDir + TVar::SmurfProcessName(higgsSample->GetProcess()) + "_LR_overlay_log.png");  

  // tidy up
  delete hist_sig; 
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
  delete h1_stack;
  delete c1;
  delete stacklg;
  delete overlaylg;

}


double calcsig(double nsig, double nbkg, Int_t fom) {
  
  if(nsig == 0.0 || nbkg == 0.0) return 0.0;
  if(fom == 0) return nsig/sqrt(nsig+nbkg+pow(0.35*nbkg, 2));
  if(fom == 1) return nsig/sqrt(nsig+nbkg);
  if(fom == 2) return nsig/nbkg;
  return 0.;
}


// 
// Calculate the significance for a signal histogram and a set of background histograms
//
void drawSignificance(Sample* & higgsSample, std::vector<Sample*> & bkgsamples, Int_t nBins, Float_t binMin, Float_t binMax,
		      Float_t nsig_CB, Float_t nbkg_CB, Float_t nsig_MVA, Float_t nbkg_MVA, TString outputDir, Int_t fom) 
{
  std::cout << "drawSignficance for " << TVar::SmurfProcessName(higgsSample->GetProcess()) << "....\n";
  // load the signal and bkg histograms 
  TH1F *hist_sig = (TH1F*)higgsSample->GetHistogram(higgsSample->GetProcess(), nBins, binMin, binMax);
  
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(bkgsamples.size());
  
  for (unsigned int s = 0; s < bkgsamples.size(); ++s) 
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(higgsSample->GetProcess(), nBins, binMin, binMax));
  
  // define the SOverB histogram
  TH1F *hist_SOverB = new TH1F("SOverB", "SOverB", nBins, binMin, binMax);
  hist_SOverB->GetXaxis()->SetTitle(Form("LR(%s) Mininum",  TVar::SmurfProcessName(higgsSample->GetProcess()).Data()));
  TString yTitle; 

  for (int i = 1;i < hist_SOverB->GetNbinsX()+1; i++) 
    {
      double nsig (0.0), nbkg(0.0);
      nsig = hist_sig->Integral(i, 100000);

      for(unsigned int j = 0; j < bkgsamples.size(); j++) 
	nbkg += hist_bkg[j]->Integral(i, 100000);
      
      if (nbkg <= 0) continue;
      Float_t SOverB = calcsig(nsig, nbkg, fom);
      hist_SOverB->SetBinContent(i, SOverB);
    }
  
  
  Float_t sign_CB = calcsig(nsig_CB, nbkg_CB, fom);
  Float_t sign_MVA = calcsig(nsig_MVA, nbkg_MVA, fom); 
  
  // get the ranges of the histogram right
  Float_t maxBin =  hist_SOverB->GetMaximum(); 
  maxBin = sign_MVA > maxBin ? sign_MVA : maxBin;
  maxBin = sign_CB > maxBin ? sign_CB : maxBin;

  Float_t minBin =  hist_SOverB->GetMinimum(); 
  minBin = sign_MVA < minBin ? sign_MVA : minBin;
  minBin = sign_CB < minBin ? sign_CB : minBin;
  
  hist_SOverB->GetYaxis()->SetRangeUser(minBin*0.9, maxBin*1.1);
  
  if(fom == 0)  hist_SOverB->GetYaxis()->SetTitle("N_{S}/#sqrt{(N_{S}+N_{B}}+#sigmaN_{B}^{2})");
  if(fom == 1)  hist_SOverB->GetYaxis()->SetTitle("N_{S}/#sqrt{(N_{S}+N_{B}})");
  if(fom == 2)  hist_SOverB->GetYaxis()->SetTitle("N_{S}/N_{B}");
  
  TArrow *arr_CB = new TArrow(binMin, sign_CB, binMax, sign_CB);
  arr_CB->SetLineStyle(7);
  arr_CB->SetLineWidth(3);
  arr_CB->SetLineColor(kBlue);
  
  TArrow *arr_MVA = new TArrow(binMin, sign_MVA, binMax, sign_MVA);
  arr_MVA->SetLineStyle(7);
  arr_MVA->SetLineWidth(3);
  arr_MVA->SetLineColor(kRed);
  
  TLegend *lg = new TLegend(0.4, 0.15, 0.7, 0.25);
  if (sign_CB != 0.0)
    lg->AddEntry(arr_CB, "Cut-Based Analysis", "l");
  if (sign_MVA != 0.0)
    lg->AddEntry(arr_MVA, "MVA Analysis", "l");
  
  lg->SetTextSize(0.03);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetShadowColor(0);
  
  TCanvas *c1 = new TCanvas();
  c1->cd();
  hist_SOverB->Draw();
  if(sign_CB != 0.0)
    arr_CB->Draw("same");
  if(sign_MVA != 0.0)
    arr_MVA->Draw("same");
  lg->Draw("same");
  c1->SaveAs(outputDir+ Form("SOverB_%s_fom%i.eps", TVar::SmurfProcessName(higgsSample->GetProcess()).Data(), fom));
  c1->SaveAs(outputDir + Form("SOverB_%s_fom%i.png", TVar::SmurfProcessName(higgsSample->GetProcess()).Data(), fom));
  
  delete hist_sig;
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
  delete hist_SOverB;
  delete c1;
  delete arr_CB;
  delete arr_MVA;

}




void makeOverlay(TVar::Process higgsProcess, TString smurfFDir, TString outputDir, TString meFSuffix, Float_t massCut, 
		 Float_t nsig_CB, Float_t nbkg_CB, Float_t nsig_MVA, Float_t nbkg_MVA)
{

  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");

  Int_t LRBins = 20;
  Float_t LRMin = 0.0;
  Float_t LRMax = 1.0;

  Sample *higgsSample = new Sample(higgsProcess, kRed, false, smurfFDir, outputDir, meFSuffix, massCut);
  
  // now define background samples
  std::vector<Sample*> bkgSamples;
  bkgSamples.push_back(new Sample(TVar::ttbar, kMagenta, true, smurfFDir, outputDir, meFSuffix, massCut));
  bkgSamples.push_back(new Sample(TVar::Wp_1jet, kCyan, true,  smurfFDir, outputDir, meFSuffix, massCut));
  bkgSamples.push_back(new Sample(TVar::WW, kYellow+2, true,  smurfFDir, outputDir, meFSuffix, massCut));
  
  // now make the overlay plots
  drawLROverlay(higgsSample, bkgSamples, LRBins, LRMin, LRMax, outputDir);
  
  // HWW130
  drawSignificance(higgsSample, bkgSamples, LRBins, LRMin, LRMax, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, outputDir, 0);
  drawSignificance(higgsSample, bkgSamples, LRBins, LRMin, LRMax, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, outputDir, 1);
  drawSignificance(higgsSample, bkgSamples, LRBins, LRMin, LRMax, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, outputDir, 2);
  
}

