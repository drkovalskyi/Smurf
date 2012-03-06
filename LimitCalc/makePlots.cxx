#include "TROOT.h"
#include "TFile.h"
#include "TList.h"
#include <map>
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <iostream>
#include "TCanvas.h"
#include "TH1D.h"
#include <algorithm>
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"

#ifndef __CINT__
// typedef std::map<TString,TGraphAsymmErrors*> SubMap;
typedef std::map<TString,TH1*> SubMap;
typedef std::map<TString,SubMap > Map;
TH1* makeHistogram(const TGraphAsymmErrors* graph){
  TH1D* hist = new TH1D(Form("h_%s",graph->GetName()),Form("%s",graph->GetName()),20,-1,1);
  hist->SetDirectory(0);
  for (int i=0; i<graph->GetN();++i){
    Int_t bin = hist->FindBin(graph->GetX()[i]);
    hist->SetBinContent(bin,graph->GetY()[i]);
  }
  hist->SetStats(0);
  return hist;
}
Map loadPlots(const char* file){
  Map map;
  TFile* f = TFile::Open(file);
  assert(f);
  TList* list = f->GetListOfKeys();
  for( int i=0; i<list->GetSize(); ++i ){
    // list_sb->At(i)->Print();
    TString s(list->At(i)->GetName());
    TObjArray *subStrL = TPRegexp("^([^_]+)_([^_]+)$").MatchS(s);
    const Int_t nrSubStr = subStrL->GetLast()+1;
    if (nrSubStr > 2) {
      TString card = ((TObjString *)subStrL->At(1))->GetString();
      TString process  = ((TObjString *)subStrL->At(2))->GetString();
      TGraphAsymmErrors* obj = dynamic_cast<TGraphAsymmErrors*>(f->Get(s));
      assert(obj);
      map[card][process] = makeHistogram(obj);
    }
  }
  return map;
}
Map loadPlots(std::string file){
  return loadPlots(file.c_str());
}

typedef std::pair<TString,double> StringPair;

bool comp(const StringPair& a, const StringPair& b){
  return a.second>b.second;
}

TH1* getACopy(const TH1* iHist, Color_t fillColor){
  TH1* outHist = (TH1*)iHist->Clone();
  outHist->SetDirectory(0);
  outHist->SetFillColor(fillColor);
  outHist->SetLineColor(kBlack);
  outHist->SetLineStyle(1);
  return outHist;
}

//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
const Float_t _tsize   = 0.03;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;
 
void DrawLegend(Float_t x1,
		Float_t y1,
		const TH1*   hist,
		TString label,
		TString option)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + _xoffset,
				y1 + _yoffset);
  
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (_tsize);
  legend->AddEntry(hist, label.Data(), option.Data());
  legend->Draw();
}

void DrawLegends(SubMap& map, const TH1* data, bool showHiggs = false){
  float xPos[9] = {0.22,0.22,0.22,0.39,0.39,0.39,0.56,0.56,0.56}; 
  float yOff[9] = {0,1,2,0,1,2,0,1,2};
  size_t j=0;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, data,       " data",  "lp"); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["qqWW"],  " WW",      "f" ); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["Zjets"], " Z+jets",  "f" ); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["Top"],   " top",     "f" ); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["VV"],    " WZ/ZZ",   "f" ); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["Wjets"], " W+jets",  "f" ); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["Wgamma"]," W#gamma", "f" ); j++;
  if (showHiggs) {DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, map["ggH"], " Higgs", "f" ); j++;}
}

void DrawLegends2(const TH1* h1, const TH1* h2, const TH1* h3)
{  
  float xPos[3] = {0.17,0.39,0.61}; 
  float yOff[3] = {0,0,0};
  size_t j=0;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, h1,       " Sig+Bkg fit",  "l"); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, h2,       " Bkg fit",  "l"); j++;
  DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, h3,       " No fit",  "l"); j++;
}


#endif

void makePlots(unsigned int mass=120)
{
  std::string file_sb  = Form("limits_nj_shape-%u-cls_asym_fittedShape_floatMu.root",mass);
  std::string file_b   = Form("limits_nj_shape-%u-cls_asym_fittedShape_mu0.root",mass);
  std::string file_nom = Form("limits_nj_shape-%u-cls_asym_nominalShape.root",mass);
  std::string file_n0of = Form("Full2011-v2/%u/hwwof_0j.input.root",mass);
  std::string file_n0sf = Form("Full2011-v2/%u/hwwsf_0j.input.root",mass);
  std::string file_n1of = Form("Full2011-v2/%u/hwwof_1j.input.root",mass);
  std::string file_n1sf = Form("Full2011-v2/%u/hwwsf_1j.input.root",mass);
  
  gROOT->SetStyle("Plain");
  // gStyle->SetOptTitle(0);
  Map map_sb(loadPlots(file_sb));
  std::cout << "size: " << map_sb.size() << std::endl;
  Map map_b(loadPlots(file_b));
  assert(map_sb.size()==map_b.size());
  Map map_nom(loadPlots(file_nom));
  assert(map_sb.size()==map_nom.size());

  TFile* f_n0of = TFile::Open(file_n0of.c_str());
  assert(f_n0of);
  TH1* histo_Data_n0of = (TH1*)f_n0of->Get("histo_Data");
  histo_Data_n0of->SetMarkerStyle(20);
  histo_Data_n0of->SetDirectory(0);
  
  TFile* f_n0sf = TFile::Open(file_n0sf.c_str());
  assert(f_n0sf);
  TH1* histo_Data_n0sf = (TH1*)f_n0sf->Get("histo_Data");
  histo_Data_n0sf->SetMarkerStyle(20);
  histo_Data_n0sf->SetDirectory(0);
  
  TFile* f_n1of = TFile::Open(file_n1of.c_str());
  assert(f_n1of);
  TH1* histo_Data_n1of = (TH1*)f_n1of->Get("histo_Data");
  histo_Data_n1of->SetMarkerStyle(20);
  histo_Data_n1of->SetDirectory(0);
  
  TFile* f_n1sf = TFile::Open(file_n1sf.c_str());
  assert(f_n1sf);
  TH1* histo_Data_n1sf = (TH1*)f_n1sf->Get("histo_Data");
  histo_Data_n1sf->SetMarkerStyle(20);
  histo_Data_n1sf->SetDirectory(0);

  std::map<TString,Color_t> color_map;
  color_map["Top"] = kYellow;
  color_map["VV"] = kAzure-2;
  color_map["WH"] = kRed+1;
  color_map["Wgamma"] = kBlue;
  color_map["Wjets"] = kGray+1;
  color_map["ZH"] = kRed+1;
  color_map["Zjets"] = kGreen+2;
  color_map["ggH"] = kRed+1;
  color_map["ggWW"] = kAzure-9;
  color_map["qqH"] = kRed+1;
  color_map["qqWW"] = kAzure-9;

//   hbkg = (TH1*)histo_Data->Clone();
//   hbkg->SetTitle("Bkg");
//   hbkg->SetName("hbkg");
//   hbkg->Add(histo_Top,-1);
//   hbkg->Add(histo_VV,-1);
//   hbkg->Add(histo_Wgamma,-1);
//   hbkg->Add(histo_Wjets,-1);
//   hbkg->Add(histo_Zjets,-1);
//   hbkg->Add(histo_ggWW,-1);
//   hbkg->Add(histo_qqWW,-1);

//   new TCanvas("cc","cc",600,600);
//   histo_Data->Draw("e");
//   hbkg->Draw("hist same");
// )


  unsigned int nProcesses = 0;
  for(Map::iterator itr = map_sb.begin(); itr!=map_sb.end(); ++itr){
    if (nProcesses==0) nProcesses = itr->second.size();
    assert(nProcesses==itr->second.size());
    TString channel = itr->first;
    TCanvas* c = new TCanvas(Form("c_%u_%s",mass,channel.Data()),"",900,900);
    unsigned int sidex = int(sqrt(nProcesses));
    unsigned int sidey = sidex;
    if (sidex*sidey < nProcesses) sidex++;
    if (sidex*sidey < nProcesses) sidey++;
    c->Divide(sidex,sidey);
    double maxValue(0);
    double maxInt(0);
    std::vector<StringPair> allProcesses;
    std::vector<StringPair> bkgProcesses;
    for(SubMap::const_iterator itr_p = itr->second.begin(); itr_p!=itr->second.end(); ++itr_p){
      double value = itr_p->second->GetMaximum();
      double integral = itr_p->second->Integral();
      if (value>maxValue) maxValue = value;
      if (integral>maxInt) maxInt = integral;
    }
    // std::cout << "maxInt: " << maxInt << std::endl;
    // const double cutoff = 0.10; // ignore contributions that are less than the threshold, which is a fraction with respect to the biggest contribution
    for(SubMap::const_iterator itr_p = itr->second.begin(); itr_p!=itr->second.end(); ++itr_p){
      double integral = itr_p->second->Integral();
      if (itr_p->first=="ggH" || itr_p->first=="qqH" || itr_p->first=="ZH" || itr_p->first=="WH" )
	{ allProcesses.push_back(StringPair(itr_p->first,1e9)); //put them on top
	  continue;
	}
      allProcesses.push_back(StringPair(itr_p->first,integral));
      // if (integral<maxInt*cutoff) continue;
      bkgProcesses.push_back(StringPair(itr_p->first,integral));
    }
    std::sort(bkgProcesses.begin(),bkgProcesses.end(),comp);
    std::sort(allProcesses.begin(),allProcesses.end(),comp);
    std::reverse(allProcesses.begin(),allProcesses.end());
    // std::cout << "bkgProcesses: " << bkgProcesses.size() << std::endl;
    maxValue *= 1.2;
    unsigned int i=1;
    for(SubMap::iterator itr_p = itr->second.begin(); itr_p!=itr->second.end(); ++itr_p){
      c->cd(i);
      TString process = itr_p->first;
      TH1* hist_sb = itr_p->second;
      // hist_sb->SetTitle(Form("%s_%s_%u",channel.Data(),process.Data(),mass));
      hist_sb->SetTitle(process.Data());
      hist_sb->SetMaximum(maxValue);
      hist_sb->SetLineColor(kRed);
      hist_sb->Draw();
      TH1* hist_b = map_b[channel][process];
      hist_b->SetTitle(Form("%s_%s_%u",channel.Data(),process.Data(),mass));
      hist_b->SetLineColor(kBlack);
      hist_b->Draw("same");
      TH1* hist_nom = map_nom[channel][process];
      hist_nom->SetTitle(Form("%s_%s_%u",channel.Data(),process.Data(),mass));
      hist_nom->SetLineColor(kBlue);
      hist_nom->SetLineStyle(3);
      hist_nom->Draw("same");
      i++;
    }
    c->Print(Form("c_%u_%s.pdf",mass,channel.Data()));
    TCanvas* c2 = new TCanvas(Form("c2_%u_%s",mass,channel.Data()),"",1600,400);
    c2->Divide(4,1);
    for ( unsigned int j=0; j<4; ++j ){
      c2->cd(j+1);
      TString process = bkgProcesses.at(j).first;
      map_sb[channel][process]->Draw();
      map_b[channel][process]->Draw("same");
      map_nom[channel][process]->Draw("same");
      DrawLegends2(map_sb[channel][process],map_b[channel][process],map_nom[channel][process]);
    }
    c2->Print(Form("bdt2_%u_%s.pdf", mass, channel.Data()));
    TCanvas* c3 = new TCanvas(Form("c3_%u_%s",mass, channel.Data()),"",1200,400);
    c3->Divide(3,1);
    THStack* h_sb  = new THStack(Form("hsb_%u_%s",  mass, channel.Data()), "Signal+background fit");
    THStack* h_b   = new THStack(Form("hb_%u_%s",   mass, channel.Data()), "Background only fit");
    THStack* h_nom = new THStack(Form("hnom_%u_%s", mass, channel.Data()), "Nominal");
    SubMap submap_nom, submap_b, submap_sb;
    for (std::vector<StringPair>::const_iterator iProcess = allProcesses.begin();
	 iProcess != allProcesses.end(); ++iProcess){
      Color_t color = color_map[iProcess->first];
      TH1* hist = getACopy(map_sb[channel][iProcess->first],color); 
      h_sb->Add(hist);
      submap_sb[iProcess->first] = hist;
      if ( iProcess->first=="ggH" || iProcess->first=="qqH" || iProcess->first=="WH" || iProcess->first=="ZH" ) continue;
      hist = getACopy(map_nom[channel][iProcess->first],color);
      h_nom->Add(hist);
      submap_nom[iProcess->first] = hist;
      hist = getACopy(map_b[channel][iProcess->first],color);
      h_b->Add(hist);
      submap_b[iProcess->first] = hist;
    }
    TH1* h_data = 0;
    if (channel == "n0of") h_data = histo_Data_n0of;
    if (channel == "n0sf") h_data = histo_Data_n0sf;
    if (channel == "n1of") h_data = histo_Data_n1of;
    if (channel == "n1sf") h_data = histo_Data_n1sf;
    assert(h_data);
    double maxY = std::max(std::max(h_nom->GetMaximum(),h_data->GetMaximum()),
			   std::max(h_b->GetMaximum(),h_sb->GetMaximum()));
    maxY *= 1.2;
    c3->cd(1);
    h_nom->SetMaximum(maxY);
    h_nom->Draw("hist");
    DrawLegends(submap_nom,h_data);
    h_data->Draw("e same");
    c3->cd(2);
    h_b->SetMaximum(maxY);
    h_b->Draw("hist");
    DrawLegends(submap_b,h_data);
    h_data->Draw("e same");
    c3->cd(3);
    h_sb->SetMaximum(maxY);
    h_sb->Draw("hist");
    DrawLegends(submap_sb,h_data,true);
    h_data->Draw("e same");
    c3->Print(Form("bdt_%u_%s.pdf", mass, channel.Data()));
  }
}
