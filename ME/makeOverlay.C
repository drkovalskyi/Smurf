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
#include "TVar.hh"

const unsigned int kNDilep = 4;


class Sample {

public:
  
  Sample(TVar::Process process, Color_t color, bool stack) { 
    
    // initialize data members
    process_ = process;
    stack_ = stack;
    color_ = color;
    file_me_ = TFile::Open(TVar::SmurfProcessName(process)+"_ME.root");
    if (file_me_ == 0x0) {
      std::cout << "Sample (" << TVar::SmurfProcessName(process_) << ") ME file is BAD" << std::endl;
    }
    else {
      tree_me_ = (TTree*)file_me_->Get("tree");
    }
    
    // get the yields for the normalisations
    InitializeYields();
    
  }
  ~Sample() {
    file_me_->Close();
  };
  
  TString GetProcessName() {
    return TVar::SmurfProcessName(process_);
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
    TH1F *tmp = new TH1F(Form("h1_%sLR_%s",  TVar::SmurfProcessName(processLR).Data(), GetProcessName().Data()), "LR", nBins, binMin, binMax);
    tree_me_->Project(Form("h1_%sLR_%s", TVar::SmurfProcessName(processLR).Data(), GetProcessName().Data()), Form("LR[%i]", processLR));
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
  
  void InitializeYields() 
  {
    
    TFile *file_norm = TFile::Open("data/"+GetProcessName()+".root");
    TTree *tree_norm = (TTree*)file_norm->Get("tree");
    
    if (file_norm == 0x0) {
      std::cout << "Sample (" << GetProcessName() << ") NORM file is bad" << std::endl;
      return;
    }
    
    gROOT->cd();
    for (unsigned int j = 0; j < kNDilep; j++)
      {
	TH1F *tmp = new TH1F("tmp", "tmp", 20, 0, 4);
	if(j==0 || j==2)
	  tree_norm->Project("tmp", "dPhi", Form("scale1fb*(type==%i)", j));
	if(j==1 || j==3)
	  tree_norm->Project("tmp", "dPhi", Form("scale1fb*(type==%i&&lep2.pt()>15)", j));
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
// draw the overlay of several LR distributions for signal and bkgs
// 

void drawLROverlay(TVar::Process higgsProcess, std::vector<Sample*> & bkgsamples, Int_t nBins, Float_t binMin, Float_t binMax)
{  
  Sample higgsSample(higgsProcess, kRed, false);
  TH1F *hist_sig = (TH1F*)higgsSample.GetHistogram(higgsProcess, nBins, binMin, binMax);
  
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(size_t(bkgsamples.size()));
  
  TLegend *lg = new TLegend(0.62, 0.6, 0.88, 0.9);
  lg->AddEntry(hist_sig, TVar::SmurfProcessName(higgsProcess), "lp");
  
  THStack *h1_stack = new THStack(); 
  for (unsigned int s = 0; s < bkgsamples.size(); ++s) {
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(higgsProcess, nBins, binMin, binMax));
    h1_stack->Add(hist_bkg[s]);
    lg->AddEntry(hist_bkg[s], bkgsamples[s]->GetProcessName(), "f");
  }

  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetShadowColor(0);
  
  // LR histograms
  TCanvas *c1 = new TCanvas();
  c1->Clear();
  c1->cd();
  h1_stack->Draw("HIST");
  h1_stack->GetXaxis()->SetTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsProcess).Data()));
  hist_sig->Draw("SAMEHIST");
  lg->Draw();
  c1->Print( TVar::SmurfProcessName(higgsProcess)+"_LR.eps");
  c1->Print( TVar::SmurfProcessName(higgsProcess)+"_LR.png");
  
  // draw log plots
  c1->Clear();
  c1->cd();
  c1->SetLogy();
  h1_stack->Draw("HIST");
  h1_stack->GetXaxis()->SetTitle(Form("LR(%s)", TVar::SmurfProcessName(higgsProcess).Data()));
  hist_sig->Draw("SAMEHIST");
  lg->Draw();
  c1->Print( TVar::SmurfProcessName(higgsProcess)+"_LR_log.eps");
  c1->Print( TVar::SmurfProcessName(higgsProcess)+"_LR_log.png");
  
  // tidy up
  delete hist_sig; 
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
  delete h1_stack;
  delete c1;
  delete lg;

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

void drawSignificance(TVar::Process higgsProcess, std::vector<Sample*> & bkgsamples, Int_t nBins, Float_t binMin, Float_t binMax, Int_t fom, Float_t nsig_CB, Float_t nbkg_CB, Float_t nsig_MVA, Float_t nbkg_MVA) 
{
  std::cout << "drawSignficance for " << TVar::SmurfProcessName(higgsProcess) << "....\n";
  // load the signal and bkg histograms 

  Sample higgsSample(higgsProcess, kRed, false);
  TH1F *hist_sig = (TH1F*)higgsSample.GetHistogram(higgsProcess, nBins, binMin, binMax);
  
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(bkgsamples.size());
  
  for (unsigned int s = 0; s < bkgsamples.size(); ++s) 
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(higgsProcess, nBins, binMin, binMax));
  
  // define the SOverB histogram
  TH1F *hist_SOverB = new TH1F("SOverB", "SOverB", nBins, binMin, binMax);
  hist_SOverB->GetXaxis()->SetTitle(Form("LR(%s) Mininum",  TVar::SmurfProcessName(higgsProcess).Data()));
  TString yTitle; 

  for (int i = 1;i < hist_SOverB->GetNbinsX(); i++) 
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
  
  if(sign_MVA > hist_SOverB->GetMaximum()) 
    hist_SOverB->GetYaxis()->SetRangeUser(hist_SOverB->GetMinimum(), sign_MVA*1.1);
  
  
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
  
  TCanvas *c1 = new TCanvas ();
  c1->cd();
  hist_SOverB->Draw();
  if(sign_CB != 0.0)
    arr_CB->Draw("same");
  if(sign_MVA != 0.0)
    arr_MVA->Draw("same");
  lg->Draw("same");
  c1->SaveAs(Form("SOverB_%s_fom%i.eps", TVar::SmurfProcessName(higgsProcess).Data(), fom));
  c1->SaveAs(Form("SOverB_%s_fom%i.png", TVar::SmurfProcessName(higgsProcess).Data(), fom));
  
  delete hist_sig;
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
  delete hist_SOverB;
  delete c1;
  delete arr_CB;
  delete arr_MVA;
}

// 
// End of Definition of HiggsProcess Class
// 

void makeOverlay()
{

  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");


  Int_t LRBins = 20;
  Float_t LRMin = 0.0;
  Float_t LRMax = 1.0;

  // now define background samples
  std::vector<Sample*> bkgsamples;
  bkgsamples.push_back(new Sample(TVar::ttbar, kMagenta, true)); 
  bkgsamples.push_back(new Sample(TVar::Wp_1jet, kCyan, true));
  bkgsamples.push_back(new Sample(TVar::WW, kYellow+2, true)); 
  
  // now make the overlay plots
  drawLROverlay(TVar::HWW130, bkgsamples, LRBins, LRMin, LRMax);
  drawLROverlay(TVar::HWW160, bkgsamples, LRBins, LRMin, LRMax);
  drawLROverlay(TVar::HWW200, bkgsamples, LRBins, LRMin, LRMax);



  // mow make the significance plots
  // The last four numbers follows (nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA)
  // Read off the numbers from the Higgs Note

  // HWW130
  drawSignificance(TVar::HWW130, bkgsamples, LRBins, LRMin, LRMax, 0, 15.0, 1.9+16.8+43.1, 6.14, 9.79+1.12);
  drawSignificance(TVar::HWW130, bkgsamples, LRBins, LRMin, LRMax, 1, 15.0, 1.9+16.8+43.1, 6.14, 9.79+1.12);
  drawSignificance(TVar::HWW130, bkgsamples, LRBins, LRMin, LRMax, 2, 15.0, 1.9+16.8+43.1, 6.14, 9.79+1.12);


  // HWW160
  drawSignificance(TVar::HWW160, bkgsamples, LRBins, LRMin, LRMax, 0, 30.9, 1.3+17.3, 11.05, 1.79+0.14);
  drawSignificance(TVar::HWW160, bkgsamples, LRBins, LRMin, LRMax, 1, 30.9, 1.3+17.3, 11.05, 1.79+0.14);
  drawSignificance(TVar::HWW160, bkgsamples, LRBins, LRMin, LRMax, 2, 30.9, 1.3+17.3, 11.05, 1.79+0.14);
  
  // HWW200
  drawSignificance(TVar::HWW200, bkgsamples, LRBins, LRMin, LRMax, 0, 0.0, 0.0, 8.05, 8.32+0.94);
  drawSignificance(TVar::HWW200, bkgsamples, LRBins, LRMin, LRMax, 1, 0.0, 0.0, 8.05, 8.32+0.94);
  drawSignificance(TVar::HWW200, bkgsamples, LRBins, LRMin, LRMax, 2, 0.0, 0.0, 8.05, 8.32+0.94);
  
      
}

