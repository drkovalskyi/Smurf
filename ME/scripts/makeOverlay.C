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
#include <fstream>

// Matrix Element related
#include "../TVar.hh"
#include "Sample.cc"

// This function is the create the matching histogram
// name as in the BDT analysis when creating the 
// LR histogram root files for the LandS

// 
// draw the overlay of several LR distributions for signal and bkgs
// 

void getProcess(int mH, TVar::Process & k);
void drawLROverlay( Sample* & dataSample, Sample* & higgsSample, std::vector<Sample*> & bkgsamples, 
		    Int_t nBins, Float_t binMin, Float_t binMax, TString outputDir)
{   
  std::cout << "Drawing the LR Overlay for process " << TVar::SmurfProcessName(higgsSample->GetLRProcess()) << "\n";
  TString var = Form("LR[%i]",higgsSample->GetLRProcess());

  // get the MVA output histograms
  TH1F *hist_data = (TH1F*)dataSample->GetHistogram(var, nBins, binMin, binMax);
  TH1F *hist_sig = (TH1F*)higgsSample->GetHistogram(var, nBins, binMin, binMax);
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(size_t(bkgsamples.size()));
  
  // define the stack and overlay legends..
  TLegend *stacklg = new TLegend(0.0, 0.4, 1.0, 0.95);
  stacklg->SetBorderSize(0);
  stacklg->SetFillStyle(0);
  stacklg->SetShadowColor(0);
  stacklg->SetTextSize(0.16);

  TLegend *overlaylg = new TLegend(0.0, 0.4, 1.0, 0.95);
  overlaylg->SetBorderSize(0);
  overlaylg->SetFillStyle(0);
  overlaylg->SetShadowColor(0);
  overlaylg->SetTextSize(0.16);
  
  // Fill the stack plots...
  THStack *stack_bkg = new THStack(); 
  for (unsigned int s = 0; s < bkgsamples.size(); ++s) {
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(var, nBins, binMin, binMax));
    stack_bkg->Add(hist_bkg[s]);
    stacklg->AddEntry(hist_bkg[s], bkgsamples[s]->name(), "f");
    overlaylg->AddEntry(hist_bkg[s], bkgsamples[s]->name(), "l");
  }
  
  stacklg->AddEntry(hist_sig, higgsSample->name(), "l");
  if (dataSample->lumi() == higgsSample->lumi())
    stacklg->AddEntry(hist_data, dataSample->name().Data(), "lp");
  overlaylg->AddEntry(hist_sig, higgsSample->name(), "l");
  
  float yMax = stack_bkg->GetMaximum();
  if (dataSample->lumi() == higgsSample->lumi())
    yMax = yMax > hist_data->GetMaximum() ? yMax : hist_data->GetMaximum();
  stack_bkg->SetMaximum(yMax + 1.2 * sqrt(yMax));

  // == Draw histograms
  // draw stacked plots
  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.0, 0.82, 1.0);
  pad1->SetBottomMargin(0.13);
  pad1->SetRightMargin(0.07);
  pad1->Draw();
  c1->cd();
  
  TPad *pad2 = new TPad("p_leg", "p_leg", 0.82, 0.0, 1.0, 1.0);
  pad2->SetTopMargin(0.01);
  pad2->SetRightMargin(0.01);
  pad2->SetBottomMargin(0.13);
  pad2->Draw();
  
  pad1->cd();
  pad1->SetLogy(0);
  stack_bkg->Draw("HIST");
  stack_bkg->GetXaxis()->SetTitle(Form("LR(%s)",  TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data()));
  stack_bkg->GetYaxis()->SetTitle(Form("Number of Events / %.2f", stack_bkg->GetXaxis()->GetBinWidth(1)));
  hist_sig->Draw("SAMEHIST");
  if (dataSample->lumi() == higgsSample->lumi())
    hist_data->Draw("SAMEE1");
  pad2->cd();
  stacklg->Draw();
  
  // draw stacked plots - linear
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR.eps");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR.png");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR.pdf");
  // draw stacked plots - log
  pad1->SetLogy();
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_log.eps");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_log.png");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_log.pdf");

  // draw overlay plots
  c1->cd();
  pad1->Clear();
  pad1->cd();
  pad1->SetLogy(0);
  hist_sig->GetXaxis()->SetTitle(Form("LR(%s)", higgsSample->name().Data()));
  hist_sig->GetYaxis()->SetTitle(Form("Number of Events / %.2f", stack_bkg->GetXaxis()->GetBinWidth(1)));
  hist_sig->Draw("HIST"); 

  yMax = hist_sig->GetMaximum();
  for (unsigned int i = 0; i < hist_bkg.size(); i++) {
    hist_bkg[i]->SetLineWidth(3);
    hist_bkg[i]->SetFillStyle(0);
    hist_bkg[i]->Draw("SAMEHIST");
    yMax = yMax > hist_bkg[i]->GetMaximum() ? yMax : hist_bkg[i]->GetMaximum();
  }
  hist_sig->SetMaximum(yMax  + 1.2 *sqrt(yMax));

  
  pad2->Clear();
  pad2->cd();
  overlaylg->Draw();

  // draw overlay plots - linear
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_overlay.eps");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_overlay.png");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_overlay.pdf");
  
  // draw overlay plots - log
  pad1->SetLogy(1);
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_overlay_log.eps");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_overlay_log.png");
  c1->SaveAs( outputDir + "plots/" + TVar::SmurfProcessName(higgsSample->GetLRProcess()) + "_LR_overlay_log.pdf");  
  
  delete pad1;
  delete pad2;
  // end of drawing histograms
  
  // tidy up
  delete hist_sig; 
  delete hist_data;
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
  delete stack_bkg;
  delete c1;
  delete stacklg;
  delete overlaylg;
}

float getTotErr(TH1F* & tmp ) {
  float err = 0;
  for (int i = 0; i < tmp->GetNbinsX()+1 ; i++) {
    err += pow(tmp->GetBinError(i),2);
  }
  return sqrt(err);
}

void writeMVAOutput(int mH, TString method, Sample* & dataSample, Sample* & higgsSample, std::vector<Sample*> & bkgsamples, 
		    Int_t nBins, Float_t binMin, Float_t binMax, TString outputDir)
{   
  cout << "Writing out the " << method << " Output for process " << TVar::SmurfProcessName(higgsSample->GetLRProcess()) << "\n";
  TString var;
  TString xtitle; 
  if (method == "ME") { 
    var = Form("LR[%i]",higgsSample->GetLRProcess());
    xtitle = "ME output";
  }
  if (method == "BDT") {
    if (mH == 150 || mH == 160 || mH ==190 )
      var = Form("bdtg_%s_ww", TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data());
    else 
      var = Form("nn_%s_ww", TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data());
    xtitle = "MVA Output";
  }

  TH1F *hist_data = (TH1F*)dataSample->GetHistogram(var, nBins, binMin, binMax);
  hist_data->SetName("histo_Data");
  hist_data->SetXTitle(xtitle);
  hist_data->SetYTitle("Number of Events");

  TH1F *hist_sig = (TH1F*)higgsSample->GetHistogram(var, nBins, binMin, binMax);
  hist_sig->SetName("histo_Higgs");
  hist_sig->SetXTitle(xtitle);
  hist_sig->SetYTitle("Number of Events");
  
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(size_t(bkgsamples.size()));

  for(unsigned int s = 0; s < bkgsamples.size(); ++s) {
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(var, nBins, binMin, binMax));
 }
  
  // Write out the input files to the limits setting
  TFile *limits_histo = new TFile(outputDir + "limits/" + "histo_limits_ntuples_" + method + "_" +TVar::SmurfProcessName(higgsSample->GetLRProcess())+".root", "RECREATE");
  limits_histo->cd();
  
  // determine the overal normalization for each processes
  float nHiggs (0.0), nqqWW(0.0), nggWW(0.0), nVV(0.0), nTop(0.0), nZjets(0.0), nWjets(0.0), nWgamma(0.0); 
  float sigmanHiggs(0.0), sigmanqqWW(0.0), sigmanggWW(0.0), sigmanVV(0.0), sigmanTop(0.0), sigmanZjets(0.0), sigmanWjets(0.0), sigmanWgamma(0.0); 
  nHiggs = hist_sig->Integral(0, 10000);
  sigmanHiggs = getTotErr(hist_sig);
  Int_t nData (0);
  if (hist_data != 0) { 
    hist_data->Scale(higgsSample->lumi()/dataSample->lumi());
    nData = hist_data->Integral(0,1000);
  }
  cout << "For higgs mass = " << mH << ": Observed " << nData << ", expected Higgs yield " << nHiggs << " \t +/- " << sigmanHiggs << "\n"; 
  
  // for the SM backgrounds
  for(unsigned int s = 0; s < bkgsamples.size(); ++s) {
    hist_bkg[s]->SetName("histo_"+bkgsamples[s]->name());
    hist_bkg[s]->SetXTitle(xtitle);
    hist_bkg[s]->SetYTitle("Number of Events");
    if (bkgsamples[s]->name() == "qqWW")   {
      nqqWW = hist_bkg[s]->Integral(0,100);
      sigmanqqWW = getTotErr(hist_bkg[s]);
    }
    if (bkgsamples[s]->name() == "ggWW") {
      nggWW = hist_bkg[s]->Integral(0,100);
      sigmanggWW = getTotErr(hist_bkg[s]);
    }
    if (bkgsamples[s]->name() == "VV") {
      nVV = hist_bkg[s]->Integral(0,100);
      sigmanVV = getTotErr(hist_bkg[s]);
    }
    if (bkgsamples[s]->name() == "Top")  {
      nTop = hist_bkg[s]->Integral(0,100);
      sigmanTop = getTotErr(hist_bkg[s]);
    }
    if (bkgsamples[s]->name() == "Wjets") {
      nWjets = hist_bkg[s]->Integral(0,100);
      sigmanWjets = getTotErr(hist_bkg[s]);
    }
    if (bkgsamples[s]->name() == "Zjets") {
      nZjets = hist_bkg[s]->Integral(0,100);
      sigmanZjets = getTotErr(hist_bkg[s]);
    }
    if (bkgsamples[s]->name() == "Wgamma") {
      nWgamma = hist_bkg[s]->Integral(0,100);
      sigmanWgamma = getTotErr(hist_bkg[s]);
    }
    hist_bkg[s]->Write();
  }
  
  hist_sig->Write();
  hist_data->Write();
  limits_histo->Close();
  
  // now write out the cards..
  // some really not complete documentation on these awesome numbers...
  // http://www.t2.ucsd.edu/tastwiki/bin/view/Smurf/MVAShape
  ofstream text; 
  text.open(outputDir + "limits/" + TVar::SmurfProcessName(higgsSample->GetLRProcess())+"_"+method.Data()+".card");
  text << "imax 1 number of channels\n";
  text << "jmax 6 number of background\n";
  text << "kmax 25 number of nuisance parameters\n";
  text << Form("Observation %i\n", nData);
  text << Form("shapes *          * histo_limits_ntuples_%s_%s.root  histo_$PROCESS\n", method.Data(), TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data());
  text << Form("shapes data_obs   1 histo_limits_ntuples_%s_%s.root  histo_Data\n", method.Data(), TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data());
  text << "bin 1 1 1 1 1 1 1\n";
  text << "process Higgs qqWW ggWW VV Top Zjets Wjets Wgamma\n";
  text << "process 0 1 2 3 4 5 6 7\n";
  text << Form("rate %.3f  %.3f  %.3f %.3f %.3f %.3f %.3f %.3f\n", 
	       nHiggs, nqqWW, nggWW, nVV, nTop, nZjets, nWjets, nWgamma);
  text << "1 lnN 1.040 1.040 1.040 1.040 1.000 1.000 1.000 1.040\n";
  text << "2 lnN 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n";
  text << "3 lnN 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n";
  text << "4 lnN 1.025 1.025 1.025 1.025 1.000 1.000 1.000 1.025\n";
  text << "5 lnN 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n";
  text << "6 lnN 1.015 1.015 1.015 1.015 1.000 1.000 1.000 1.015\n";
  text << "7 lnN 1.005 1.005 1.005 1.005 1.000 1.000 1.000 1.005\n";
  text << "8 lnN 1.070 1.000 1.054 1.054 1.000 1.000 1.000 1.054\n";
  text << "9 lnN 1.030 1.030 1.030 1.030 1.000 1.000 1.000 1.030\n";
  text << "10 lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040 1.100\n";
  // fakerate
  text << "11 lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.360 1.000\n";
  // Higgs cross-section, currently constant throughout the mass points
  text << "12 lnN 1.150 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n";
  // VV cross-section
  text << "13 lnN 1.000 1.000 1.000 1.030 1.000 1.000 1.000 1.000\n";
  // WW normalization
  text << "14 lnN 1.000 1.200 1.200 1.000 1.000 1.000 1.000 1.000\n";
  // extra Xsection normalization for ggWW
  text << "15 lnN 1.000 1.000 1.500 1.000 1.000 1.000 1.000 1.000\n";
  // top normalization
  text << "16 lnN 1.000 1.000 1.000 1.000 1.700 1.000 1.000 1.000\n";
  // dy normalization
  text << "17 lnN 1.000 1.000 1.000 1.000 1.000 1.450 1.000 1.000\n";
  // Begin the statistical uncertainties...
  // Higgs
  text << Form("18 lnN %.3f 1.000 1.000 1.000 1.000 1.000 1.000 1.000\n", 1.0 + sigmanHiggs/nHiggs);
  // qqWW
  text << Form("19 lnN 1.000 %.3f 1.000 1.000 1.000 1.000 1.000 1.000\n", 1.0 + sigmanqqWW/nqqWW);
  // ggWW 
  text << Form("20 lnN 1.000 1.000 %.3f 1.000 1.000 1.000 1.000 1.000\n", 1.0 + sigmanggWW/nggWW);
  // VV 
  text << Form("21 lnN 1.000 1.000 1.000 %.3f 1.000 1.000 1.000 1.000\n", 1.0 + sigmanVV/nVV);
  // Top
  text << Form("22 lnN 1.000 1.000 1.000 1.000 %.3f 1.000 1.000 1.000\n", 1.0 + sigmanTop/nTop);
  // Zjets
  text << Form("23 lnN 1.000 1.000 1.000 1.000 1.000 %.3f 1.000 1.000\n", 1.0 + sigmanZjets/nZjets);
  // Wjets
  text << Form("24 lnN 1.000 1.000 1.000 1.000 1.000 1.000 %.3f 1.000\n", 1.0 + sigmanWjets/nWjets);
  // Wgamma
  text << Form("25 lnN 1.000 1.000 1.000 1.000 1.000 1.000 1.000 %.3f\n", 1.0 + sigmanWgamma/nWgamma);
  text.close();
  
  // tidy up
  delete hist_sig; 
  delete hist_data;
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
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
void drawSignificance(Sample* & higgsSample, std::vector<Sample*> & bkgsamples, Int_t nBins, Float_t binMin, Float_t binMax,  Float_t nsig_CB, Float_t nbkg_CB, Float_t nsig_MVA, Float_t nbkg_MVA, TString outputDir, Int_t fom) 
{
  std::cout << "drawSignficance for " << TVar::SmurfProcessName(higgsSample->GetLRProcess()) << "....\n";
  
  TString var = Form("LR[%i]",higgsSample->GetLRProcess());
  // load the signal and bkg histograms 
  TH1F *hist_sig = (TH1F*)higgsSample->GetHistogram(var, nBins, binMin, binMax);
  
  std::vector<TH1F*> hist_bkg;
  hist_bkg.reserve(bkgsamples.size());
  
  for (unsigned int s = 0; s < bkgsamples.size(); ++s) 
    hist_bkg.push_back(bkgsamples[s]->GetHistogram(var, nBins, binMin, binMax));
  
  // define the SOverB histogram
  TH1F *hist_SOverB = new TH1F("SOverB", "SOverB", nBins, binMin, binMax);
  hist_SOverB->GetXaxis()->SetTitle(Form("LR(%s) Mininum",  TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data()));
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
  c1->SaveAs(outputDir + "plots/" + Form("SOverB_%s_fom%i.eps", TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data(), fom));
  c1->SaveAs(outputDir + "plots/" + Form("SOverB_%s_fom%i.png", TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data(), fom));
  c1->SaveAs(outputDir + "plots/" + Form("SOverB_%s_fom%i.pdf", TVar::SmurfProcessName(higgsSample->GetLRProcess()).Data(), fom));
  
  delete hist_sig;
  for (unsigned int s = 0; s < bkgsamples.size(); s++) delete hist_bkg[s];
  delete hist_SOverB;
  delete c1;
  delete arr_CB;
  delete arr_MVA;

}


void makeOverlay(int mH, TString outputDir, Float_t nsig_CB, Float_t nbkg_CB, Float_t nsig_MVA, Float_t nbkg_MVA, Float_t lumi, Float_t datalumi)
{
  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");

  TVar::Process higgsProcess;
  getProcess(mH, higgsProcess);
  
  // Data 
  TChain *DataChain = new TChain("tree");
  DataChain->Add(outputDir + "data_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *DataSample = new Sample(higgsProcess, "Data", DataChain, kBlack, false, datalumi);
  
  // higgs 
  TChain *HiggsChain = new TChain("tree");
  HiggsChain->Add(outputDir + TVar::SmurfProcessName(higgsProcess) + "_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *HiggsSample = new Sample(higgsProcess, "Higgs", HiggsChain, kRed, false, lumi);
  HiggsSample->setSF(1.13);
  
  // VV
  TChain *VVChain = new TChain("tree");
  VVChain->Add(outputDir + "zz_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  VVChain->Add(outputDir + "wz_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *VVSample = new Sample(higgsProcess, "VV", VVChain, kGreen, true, lumi);

  // Zjets
  TChain *ZjetsChain = new TChain("tree");
  ZjetsChain->Add(outputDir + "zz_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  ZjetsChain->Add(outputDir + "wz_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  ZjetsChain->Add(outputDir + "dyee_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  ZjetsChain->Add(outputDir + "dymm_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  ZjetsChain->Add(outputDir + "dytt_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *ZjetsSample = new Sample(higgsProcess, "Zjets", ZjetsChain, kBlue, true, lumi);
  
  // top
  TChain *TopChain = new TChain("tree");
  TopChain->Add(outputDir + "ttbar_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  TopChain->Add(outputDir + "tw_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *TopSample = new Sample(higgsProcess, "Top", TopChain, kMagenta, true, lumi);

  // wjets
  TChain *WjetsChain = new TChain("tree");
  //WjetsChain->Add(outputDir + "wjets_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  WjetsChain->Add(outputDir + "wjets_pythia_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *WjetsSample = new Sample(higgsProcess, "Wjets", WjetsChain, kCyan, true, lumi);

  // wgamma
  TChain *WgammaChain = new TChain("tree");
  WgammaChain->Add(outputDir + "wgamma_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *WgammaSample = new Sample(higgsProcess, "Wgamma", WgammaChain, kOrange-3, true, lumi);
  
  // ggWW
  TChain *ggWWChain = new TChain("tree");
  ggWWChain->Add(outputDir + "ggww_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *ggWWSample = new Sample(higgsProcess, "ggWW", ggWWChain, kYellow+3, true, lumi);
  
  // WW
  TChain *WWChain = new TChain("tree");
  WWChain->Add(outputDir + "qqww_LR_"+ TVar::SmurfProcessName(higgsProcess) + ".root");
  Sample *WWSample = new Sample(higgsProcess, "qqWW", WWChain, kYellow+2, true, lumi);
  
  std::vector<Sample*> bkgSamples;
  bkgSamples.push_back(ZjetsSample);
  bkgSamples.push_back(VVSample);
  bkgSamples.push_back(TopSample);
  bkgSamples.push_back(WjetsSample);
  bkgSamples.push_back(WgammaSample);
  bkgSamples.push_back(ggWWSample);
  bkgSamples.push_back(WWSample);
  
  // now make the LR overlay plots 
  Int_t LRBins = 10;
  Float_t LRMin = 0.0;
  Float_t LRMax = 1.0;
  bool drawData = false;
  if (datalumi == lumi) drawData = true;
  drawLROverlay(DataSample, HiggsSample, bkgSamples, 10, 0, 1, outputDir);
  drawSignificance(HiggsSample, bkgSamples, LRBins, LRMin, LRMax, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, outputDir, 0);
  drawSignificance(HiggsSample, bkgSamples, LRBins, LRMin, LRMax, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, outputDir, 1);
  drawSignificance(HiggsSample, bkgSamples, LRBins, LRMin, LRMax, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, outputDir, 2);
  
  // now write out the mva outputs for the limits setting
  Int_t BDTBins = 40;
  Float_t BDTMin = -1.0;
  Float_t BDTMax = 1.0;
  writeMVAOutput(mH, "BDT", DataSample, HiggsSample, bkgSamples, BDTBins, BDTMin, BDTMax, outputDir);
  
  // use the same number of bins in the shape analysis
  LRBins = BDTBins; // this makes sure that we are using the same 
  writeMVAOutput(mH, "ME", DataSample, HiggsSample, bkgSamples, LRBins, LRMin, LRMax, outputDir);
  
  delete DataChain;
  delete DataSample;

  delete HiggsChain;
  delete HiggsSample;
  
  delete ZjetsChain;
  delete ZjetsSample;

  delete VVChain;
  delete VVSample;

  delete TopChain;
  delete TopSample;
  
  delete WjetsChain;
  delete WjetsSample;

  delete WgammaChain;
  delete WgammaSample;
  
  delete ggWWChain;
  delete ggWWSample;

  delete WWChain;
  delete WWSample;
    
}

void getProcess(int mH, TVar::Process & k)
{
  switch (mH) {
  case (115):
    k = TVar::HWW115;
    break;
  case (120):
    k = TVar::HWW120;
    break;
  case (130):
    k = TVar::HWW130;
    break;
  case (140):
    k = TVar::HWW140;
    break;
  case (150):
    k = TVar::HWW150;
    break;
  case (160):
    k = TVar::HWW160;
    break;
  case (170):
    k = TVar::HWW170;
    break;
  case (180):
    k = TVar::HWW180;
    break;
  case (190):
    k = TVar::HWW190;
    break;
  case (200):
    k = TVar::HWW200;
    break;
  case (210):
    k = TVar::HWW210;
    break;
  case (220):
    k = TVar::HWW220;
    break;
  case (230):
    k = TVar::HWW230;
    break;
  case (250):
    k = TVar::HWW250;
    break;
  case (300):
    k = TVar::HWW300;
    break;
  default:
    break;
  }
  

}
