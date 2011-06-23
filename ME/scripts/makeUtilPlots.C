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

// Yanyan's drawing tool box
#include "histTool.C"

void makeUtilPlots( TString UtilFile = "../Util.root") 
{
  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  TGaxis *gaxis = new TGaxis();
  gaxis->SetMaxDigits(3);
  
  drawBoost(UtilFile, "output/plots/");
  drawEff(UtilFile, "output/plots/", "Eta", "Lepton #eta", -2.5, 2.5);
  drawEff(UtilFile, "output/plots/", "Pt", "Lepton pT (GeV)", 0, 100);
  drawFR(UtilFile, "output/plots/", "wjets_heleGenFR");
  drawFR(UtilFile, "output/plots/", "wjets_hmuGenFR");
}


void drawFR(TString UtilFile, TString outputdir, TString histName)
{
  TFile *file = TFile::Open(UtilFile);
  assert(file);
  gROOT->cd();

  TH2F *fr = (TH2F*) file->Get(histName);
  fr->GetYaxis()->SetRangeUser(10, 40);
  fr->SetXTitle("Lepton |#eta|");
  fr->SetYTitle("Lepton pT (GeV)");
  fr->GetYaxis()->SetTitleOffset(1.6);
  fr->GetXaxis()->SetTitleOffset(1.2);
  
  TCanvas *c1 = new TCanvas();
  c1->SetRightMargin(0.16);
  fr->Draw("colztexte1");
  
  c1->SaveAs(outputdir+histName+".png");
  c1->SaveAs(outputdir+histName+".pdf");
  
  // tidy up
  delete c1;
  delete fr;
  file->Close();
}


void drawEff(TString UtilFile, TString outputdir, TString var, TString x_title, float xMin, float xMax)
{
  TFile *file = TFile::Open(UtilFile);
  assert(file);
  gROOT->cd();
  
  TString y_title = "Lepton Efficiency"; 
  TLegend *overlaylg = new TLegend(0.3, 0.3, 0.88, 0.45);
  overlaylg->SetBorderSize(0);
  overlaylg->SetFillStyle(0);
  overlaylg->SetShadowColor(0);

  TH1F *els_eff = (TH1F*) file->Get("ww_heleEff"+var);
  setStyle(els_eff, 1, false, kBlue, 20, x_title, y_title, xMin, xMax);
  els_eff->GetYaxis()->SetRangeUser(0,1.0);
  overlaylg->AddEntry(els_eff, "Electron", "lp");
  
  TH1F *mus_eff = (TH1F*) file->Get("ww_hmuEff"+var);
  setStyle(mus_eff, 1, false, kRed, 22, x_title, y_title, xMin, xMax);
  mus_eff->GetYaxis()->SetRangeUser(0,1.0);
  overlaylg->AddEntry(mus_eff, "Muon", "lp");
  
  TCanvas *c1 = new TCanvas();
  els_eff->Draw("he1");
  mus_eff->Draw("samehe1");
  overlaylg->Draw("SAME");
  
  c1->SaveAs(outputdir+"lepton_eff_"+var+".png");
  c1->SaveAs(outputdir+"lepton_eff_"+var+".pdf");
  
  // tidy up
  delete c1;
  delete overlaylg;
  delete els_eff;
  delete mus_eff;
  file->Close();
}

void drawBoost(TString UtilFile, TString outputdir) {
  
  std::vector<TString> process;
  process.push_back("ww");
  process.push_back("hww130");
  process.push_back("hww160");
  process.push_back("hww200");
  
  std::vector<TH1F*> hist;
  hist.reserve(size_t(process.size()));
  
  TFile *file = TFile::Open(UtilFile);
  assert(file);
  gROOT->cd();
  
  float yMax (0.);
  for(unsigned int i = 0; i < process.size(); i++) {
    hist[i] = new TH1F(Form("%s_kt", process[i].Data()),  Form("%s_kt", process[i].Data()),  20, 0, 100);
    hist[i] = (TH1F*)  file->Get(process[i]+"_kt");
    hist[i]->Scale(1.0/hist[i]->Integral(0, 10000));
    yMax = hist[i]->GetMaximum() > yMax ? hist[i]->GetMaximum() : yMax;
  }
  

  TLegend *overlaylg = new TLegend(0.62, 0.6, 0.88, 0.9);
  overlaylg->SetBorderSize(0);
  overlaylg->SetFillStyle(0);
  overlaylg->SetShadowColor(0);

  
  TCanvas *c1 = new TCanvas();
  for (unsigned int i = 0; i < process.size(); i++) {
    Color_t color (kBlack);
    Style_t style (1);
    setColor(process[i], color, style);
    hist[i]->GetYaxis()->SetRangeUser(0, yMax*1.1);
    hist[i]->SetXTitle("System Boost (GeV)");
    hist[i]->SetYTitle("A.U.");
    hist[i]->SetLineColor(color);
    hist[i]->SetMarkerColor(color);
    hist[i]->SetLineStyle(style);
    hist[i]->SetLineWidth(3);
    hist[i]->SetTitleOffset(1.2, "X");
    hist[i]->SetTitleOffset(2., "Y");
    overlaylg->AddEntry(hist[i], process[i], "l");
    
    if (i == 0)
      hist[i]->Draw();
    else 
      hist[i]->Draw("samehist");
    overlaylg->Draw("SAME");
  }
    

  c1->SaveAs(outputdir+"boost.pdf");
  c1->SaveAs(outputdir+"boost.png");

  // tidy up
  delete c1;
  for (unsigned int i = 0; i < process.size() ; i++) 
    delete hist[i];
  delete overlaylg;
  file->Close();
}

void setColor(TString process, Color_t & color, Style_t & style) 
{
  
  if (process == "ww") {
    color = kBlack;
    style = 1;
  }
  if (process == "hww130") {
    color = kRed;
    style = 1;
  }

  if (process == "hww160") {
    color = kCyan;
    style = 1;
  }

  if (process == "hww200") {
    color = kBlue;
    style = 1;
  }
 
}


void draw2D(TString fName, TString histName, TString plotName, TString option, TString x_title, TString y_title) {
  
  TFile *file = TFile::Open(fName.Data());
  gROOT->cd();

  TH2F *hist = (TH2F*) file->Get(histName);
  if (hist==0x0) return;
  setStyle(hist, kBlue, 20, x_title, y_title, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());
  
  TCanvas *c1 = new TCanvas("c1", "c1");
  hist->Draw(option);	   
  c1->SetRightMargin(0.15);
  c1->SaveAs(TString(plotName+".png"));
  c1->SaveAs(TString(plotName+".eps"));

  delete c1;
  delete hist;
  delete file;
}
