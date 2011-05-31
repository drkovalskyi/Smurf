// Make plots from histograms
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "THStack.h"
#include <iostream>
#include <stdlib.h>
#include <vector.h>


float GetEntries(const TH1F *h, int lowBin, int highBin) {
  
  if(lowBin < 0 || highBin > h->GetNbinsX()+1) {
    cout << "Bins out of range. Setting the lowBin to the underflow and the highBin to the overflow" << endl;
    lowBin = 0;
    highBin = h->GetNbinsX();
  }
  
  float nentries = 0;
  for(int i = lowBin; i < highBin+1; i++) 
    nentries = nentries + h->GetBinContent(i);
  
  return nentries;
}

void getEff(TH1F*  hist, TH1F* hist_eff) {
  for (int i=1;i<hist->GetNbinsX()+1;i++) {
    double denom = hist->Integral(0, hist->GetNbinsX()+1);
    Double_t num = hist->Integral(0, i);
    Double_t eff = num/denom;
    hist_eff->SetBinContent(i,eff);
    hist_eff->SetBinError( i, sqrt(eff*(1-eff)/hist->GetEntries()));
  }
}

void getMetEff(TH1F*  hist, TH1F* hist_eff) {
  for (int i=1;i<hist->GetNbinsX()+1;i++) {
    double denom = hist->Integral(0, hist->GetNbinsX()+1);
    Double_t num = hist->Integral(i+1,hist->GetNbinsX()+1);
    Double_t eff = num/denom;
    hist_eff->SetBinContent(i,eff);
    hist_eff->SetBinError( i, sqrt(eff*(1-eff)/hist->GetEntries()));
  }
}

void fixRangeY(TH1* r,TH1* s){
  double ymin = (r->GetBinContent(r->GetMinimumBin()) < s->GetBinContent(s->GetMinimumBin())) ? 
    r->GetBinContent(r->GetMinimumBin()) : s->GetBinContent(s->GetMinimumBin());
  double ymax = (r->GetBinContent(r->GetMaximumBin()) > s->GetBinContent(s->GetMaximumBin())) ?
    r->GetBinContent(r->GetMaximumBin()) : s->GetBinContent(s->GetMaximumBin());
  r->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
  s->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
}


void fixRangeY(TH1 **h, const int nHist) {
  double ymin = h[0]->GetBinContent(h[0]->GetMinimumBin());
  double ymax =  h[0]->GetBinContent(h[0]->GetMaximumBin());
  for(int i=1;i<nHist;i++) {
    ymin = ymin < h[i]->GetBinContent(h[i]->GetMinimumBin()) ? ymin : h[i]->GetBinContent(h[i]->GetMinimumBin());
    ymax = ymax > h[i]->GetBinContent(h[i]->GetMaximumBin()) ? ymax : h[i]->GetBinContent(h[i]->GetMaximumBin());
  }
  for (int i=0;i<nHist;i++) 
    h[i]->GetYaxis()->SetRangeUser(ymin*0.9,ymax*1.1);
}

void setStyle(TH1F *& hist, int rebin, bool scale, Color_t color, Style_t marker, TString x_title, TString y_title, float xMin, float xMax)
{
  hist->Rebin(rebin);
  if(scale)
    hist->Scale(1.0/rebin);
  hist->SetLineColor(color);
  hist->SetLineWidth(3);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->GetXaxis()->SetRangeUser(xMin,xMax);

  hist->GetXaxis()->SetTitle(x_title);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);

  hist->GetYaxis()->SetTitle(y_title);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.04);

  hist->SetTitleOffset(1.4,"Y");
  hist->SetTitleOffset(1.1,"X");
  hist->SetTitleFont(42);

  hist->SetStats(0);
  hist->SetFillColor(0);
  //hist->SetBinContent(1, hist->GetBinContent(0)+hist->GetBinContent(1));
}


void setStyle(TH2F *& hist, Color_t color, Style_t marker, TString x_title, TString y_title, float xMin, float xMax, float yMin, float yMax)
{
  hist->SetLineColor(color);
  hist->SetLineWidth(3);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.2);
  hist->GetXaxis()->SetRangeUser(xMin,xMax);
  hist->GetYaxis()->SetRangeUser(yMin,yMax);

  hist->GetXaxis()->SetTitle(x_title);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);

  hist->GetYaxis()->SetTitle(y_title);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.04);

  hist->SetTitleOffset(1.4,"Y");
  hist->SetTitleOffset(1.1,"X");
  hist->SetTitleFont(42);

  hist->SetStats(0);
  hist->SetFillColor(0);
  //hist->SetBinContent(1, hist->GetBinContent(0)+hist->GetBinContent(1));
}


TArrow *getArrow(const int nHist, TH1F **hist, float cutVal)
{
  float cutVal; 
  float max (0);
  for(int i=0; i<nHist;i++) 
    max = max > hist[i]->GetMaximum() ? max : hist[i]->GetMaximum();
    
  TArrow *arr_cut = new TArrow(cutVal, max, cutVal, 0, 0.08, "|>");
  arr_cut->SetAngle(25.0);
  arr_cut->SetLineWidth(3.0);
  arr_cut->SetLineColor(kBlack);
  arr_cut->SetFillColor(kBlack);
  arr_cut->SetLineStyle(kDashed);
  arr_cut->SetFillStyle(3001);
  return arr_cut;
}

TLegend *getLegend(const int nHist, TH1F**hist, TString* & label, double x1, double y1, double x2, double y2) {
  TLegend *leg = new TLegend(x1, y1, x2, y2, "", "brNDC");
  for(int i = 0;i<nHist;i++) 
    leg->AddEntry(hist[i], label[i], "lp");
  leg->SetFillColor(0);
  return leg;
}

unsigned int getBin(TH1F* h, float x) 
{
  for(int i=0;i<h->GetNbinsX();i++) {
    if(x>h->GetBinLowEdge(i) && x <= h->GetBinLowEdge(i) + h->GetBinWidth(i)) 
      return i;
  }
}


void drawOverlay(TString fileName, TString hist1Name, TString hist2Name, TString plotName, bool normalize, float xMin, float xMax, TString option)
{
  TFile *_file0 = TFile::Open(fileName);
  TH1F *h1 = (TH1F*)_file0->Get(hist1Name);
  TH1F* h2 = (TH1F*)_file0->Get(hist2Name);

  if(normalize) {
    h1->Scale(1.0/h1->Integral(0,h1->GetNbinsX()+1));
    h2->Scale(1.0/h2->Integral(0,h2->GetNbinsX()+1));
    setStyle(h1, 1, false, kBlue,  20, h1->GetXaxis()->GetTitle(), "A.U", xMin, xMax);
    setStyle(h2, 1, false, kRed,  20, h2->GetXaxis()->GetTitle(), "A.U",xMin, xMax);
  }

  setStyle(h1, 1, false, kBlue,  20, h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle(), xMin, xMax);
  setStyle(h2, 1, false, kRed,  20, h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle(),xMin, xMax);

  fixRangeY(h1,h2);

  TLegend *leg =  new TLegend(0.35, 0.7, 0.8, 0.9, "", "brNDC");
  leg->AddEntry(h1,hist1Name, "lp");
  leg->AddEntry(h2,hist2Name, "lp");
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1->Draw(option);
  h2->Draw(TString("sameh"+option));
  leg->Draw("same");
  c1->SaveAs(TString(plotName+".png"));
  c1->SaveAs(TString(plotName+".eps"));
}





void drawOverlay(TString file1Name, TString file2Name, TString histName, TString plotName, bool normalize, float xMin, float xMax, TString option)
{
  TFile *_file1 = TFile::Open(file1Name);
  TFile *_file2 = TFile::Open(file2Name);
  TH1F *h1 = (TH1F*)_file1->Get(histName);
  TH1F* h2 = (TH1F*)_file2->Get(histName);

  double maxInte = h1->Integral(0, 1000) > h2->Integral(0, 1000) ? h1->Integral(0, 100) : h2->Integral(0, 100);
  if(normalize) {
    h1->Scale(1.0/h1->Integral(0, 1000));
    h2->Scale(1.0/h2->Integral(0, 1000));
    setStyle(h1, 1, false, kBlue,  20, h1->GetXaxis()->GetTitle(), "A.U.", xMin, xMax);
    setStyle(h2, 1, false, kRed,  20, h2->GetXaxis()->GetTitle(), "A.U.",xMin, xMax);
  }
  
  setStyle(h1, 1, false, kBlue,  20, h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle(), xMin, xMax);
  setStyle(h2, 1, false, kRed,  20, h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle(),xMin, xMax);
  
  fixRangeY(h1,h2);

  TLegend *leg =  new TLegend(0.35, 0.7, 0.8, 0.9, "", "brNDC");
  leg->AddEntry(h1,file1Name, "lp");
  leg->AddEntry(h2,file2Name, "lp");
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);


  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1->Draw(option);
  h2->Draw(TString("sameh"+option));
  leg->Draw("same");
  c1->SaveAs(TString(plotName+".eps"));
}


void drawPlot(TH1F** hist, const int nHist, TString dir, TString plotName, bool logy, std::vector<TString> legend, TString legend_name, float* legend_range, bool norm)
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();  
  //gPad->SetGridy();
  //gPad->SetGridx();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);

  if(norm) {
    for(int i=0;i<nHist;i++)
      hist[i]->Scale(1.0/hist[i]->Integral(0,1000000));
  }

  double maxY = 0; 
  double minY = 0; 
  for(int i=0;i<nHist;i++) {
    int maxbin = hist[i]->GetMaximumBin();
    int minbin = hist[i]->GetMinimumBin();
    maxY = maxY < hist[i]->GetBinContent(maxbin) ? hist[i]->GetBinContent(maxbin) : maxY;
    minY = minY < hist[i]->GetBinContent(maxbin) ? minY : hist[i]->GetBinContent(minbin);
  }
  
  TLegend *leg = new TLegend(legend_range[0], legend_range[1], legend_range[2], legend_range[3], "brNDC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(legend_name);
  
  for(int i=0;i<nHist;i++) {
    if(!logy) // why is that!!
      hist[i]->GetYaxis()->SetRangeUser(minY, maxY*1.1);
    leg->AddEntry(hist[i], legend[i], "lp");
    if (i==0) hist[i]->Draw("hist");
    else
      hist[i]->Draw("sameh");
  }
  leg->Draw("SAME"); 

  if (logy) plotName+="_log";
  c1->SaveAs(TString(dir+"epsfiles/"+plotName+".eps"));
  c1->SaveAs(TString(dir+"pngfiles/"+plotName+".png"));

  delete leg;  
  delete c1;
}

void drawPlot(TH2F** hist, const int nHist, TString dir, TString plotName, std::vector<TString> legend, TString legend_name, float* legend_range, TString drawoption)
{
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();  
  //gPad->SetGridy();
  //gPad->SetGridx();

  double maxY = 0; 
  double minY = 0; 
  for(int i=0;i<nHist;i++) {
    int maxbin = hist[i]->GetMaximumBin();
    int minbin = hist[i]->GetMinimumBin();
    maxY = maxY < hist[i]->GetBinContent(maxbin) ? hist[i]->GetBinContent(maxbin) : maxY;
    minY = minY < hist[i]->GetBinContent(maxbin) ? minY : hist[i]->GetBinContent(minbin);
  }
  
  TLegend *leg = new TLegend(legend_range[0], legend_range[1], legend_range[2], legend_range[3], "brNDC");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(legend_name);
  
  for(int i=0;i<nHist;i++) {
    hist[i]->GetYaxis()->SetRangeUser(minY, maxY*1.1);
    leg->AddEntry(hist[i], legend[i], "lp");
    if (i==0) hist[i]->Draw(drawoption);
    else
      //hist[i]->Draw(TString("same"+drawoption));
      hist[i]->Draw("samebox");
  }
  leg->Draw("SAME");

  c1->SaveAs(TString(dir+"epsfiles/"+plotName+".eps"));
  c1->SaveAs(TString(dir+"pngfiles/"+plotName+".png"));

  delete leg;  
  delete c1;
}


void normalizeByArea(TH1F** hist, const int nHist) 
{
  // Assumes the first histogram is from the data
  
  for(int i=1;i<nHist;i++) {
    hist[i]->Scale(GetEntries(hist[0], 0, hist[0]->GetNbinsX()+1) / GetEntries(hist[i], 0, hist[i]->GetNbinsX()+1));
  }
}

// ===========
// This assumes that hist_sig and hist_bkg all have the same binning and ranges!!
// ===========
void getSignficance(std::vector<TH1F*> & hist_sig, std::vector<TH1F*> & hist_bkg, TH1F* & hist_SOverB, TH1F* & hist_nsig, TH1F* & hist_nbkg)
{
  if(hist_sig.size()<1 || hist_bkg.size() < 1) return;

  for (int i=1;i<hist_SOverB->GetNbinsX();i++) 
    {
      double nsig = 0;
      double nbkg = 0;
      for(int j=0;j<hist_sig.size();j++) 
	nsig += hist_sig[j]->Integral(i, 100000);

      for(int j=0;j<hist_bkg.size();j++) 
	nbkg  += hist_bkg[j]->Integral(i, 100000);

      double SOverB;
      if (nbkg > 0 ) SOverB = nsig/sqrt(nsig+1.35*nbkg);
      hist_SOverB->SetBinContent(i, SOverB);
      hist_nsig->SetBinContent(i, nsig);
      hist_nbkg->SetBinContent(i, nbkg);
      //cout << "Bin "<<i<<": nsig = "<<nsig << "; nbkg = "<<nbkg<< "; SOverB = "<<SOverB<<endl;
    }
}
