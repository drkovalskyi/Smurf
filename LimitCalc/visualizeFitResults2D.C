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
#include "TH2D.h"
#include <algorithm>
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "tdrstyle.C"

#ifndef __CINT__
// typedef std::map<TString,TGraphAsymmErrors*> SubMap;
typedef std::map<TString,TH1*> SubMap;
typedef std::map<TString,SubMap > Map;
const float markersize = 1.3;

TH1* makeHistogram(const TGraphAsymmErrors* graph){
  TH1D* hist = new TH1D(Form("h_%s",graph->GetName()),Form("%s",graph->GetName()),126,-1,1);
  hist->SetDirectory(0);
  for (int i=0; i<graph->GetN();++i){
    Int_t bin = hist->FindBin(graph->GetX()[i]);
    hist->SetBinContent(bin,graph->GetY()[i]);
  }
  hist->SetStats(0);
  cout << graph->GetName() << " : " << hist->Integral()  << endl;  //FIXME  
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

// Zoom
TH2* Zoom(TH2* h2) {

    float mtbins[7]    =   {60,70,80,90,100,110,120};
    float mllbins[6]   =   {12,30,45,60,75,100};

    TH2D* h2zoom = new TH2D(Form("%s_zoom", h2->GetName()), Form("%s_zoom", h2->GetTitle()), 6, mtbins, 5, mllbins);

    for(unsigned int ibinX=1; ibinX <= 6; ++ibinX) {
        for(unsigned int ibinY=1; ibinY <= 5; ++ibinY) {
            h2zoom->SetBinContent(    ibinX, ibinY, h2->GetBinContent(ibinX, ibinY) );
            h2zoom->SetBinError(  ibinX, ibinY, h2->GetBinError(ibinX, ibinY) );
        }
    }
    h2zoom->SetStats(0);
    h2zoom->SetXTitle("M_{T} (GeV)");
    h2zoom->SetYTitle("M_{ll} (GeV)");
    h2zoom->SetTitleOffset(1.5,"Y");;
    h2zoom->SetMinimum(0);
    
    return h2zoom;
}

// rolling 
TH2* Roll1DTo2D(TH1 *h1, bool doZoom=true) {

    // style
    gStyle->SetPaintTextFormat(".1f");
    gStyle->SetMarkerSize(2.5);

    // palette
    //loadColorPalette();

    unsigned int nbins  = h1->GetXaxis()->GetNbins();
    
    float mtbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
    float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};
    TH2D* h2 = new TH2D(h1->GetName(), h1->GetTitle(), 14, mtbins, 9, mllbins);
    for(unsigned int ibin=1; ibin <= nbins; ++ibin) {

        int ibinX;
        int ibinY;
        if(ibin<=9*6) {
            ibinX = (int) (ibin-1)/9+1;
            ibinY = (int) (ibin-1)%9+1;
        } else {
            ibinX = (int) (ibin-1)/9+1;
            ibinY = (int) (ibin-1)%9+1;
        }

        h2->SetBinContent(  ibinX, ibinY, h1->GetBinContent(ibin) );
        h2->SetBinError(    ibinX, ibinY, h1->GetBinError(ibin) );
    }
    h2->SetStats(0); 
    h2->SetXTitle("M_{T} (GeV)");
    h2->SetYTitle("M_{ll} (GeV)");
    h2->SetTitleOffset(1.5,"Y");;
    h2->SetMinimum(0);
   
    // do zoom or not
    TH2D *h2zoom = (TH2D*)Zoom(h2);

    if(doZoom) return h2zoom;
    else return h2;
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

void makePlots(TString channel, const char* prefix, const char* rootfile, const char* prefix_output=0, int njets=0, int MASS=125)
{

  if (prefix_output==0) prefix_output = prefix;
  unsigned int mass = 0;
  std::string file_sb  = Form("%s_fittedShape_floatMu.root",prefix);
  std::string file_b   = Form("%s_fittedShape_mu0.root",prefix);
  std::string file_nom = Form("%s_nominalShape.root",prefix);
  
  
  gROOT->SetStyle("Plain"); 
  //gStyle->SetTitleBorderSize(1); 
  gStyle->SetPaintTextFormat("6.1f");
  cout << ".... S+B fit ...." << endl;
  Map map_sb(loadPlots(file_sb));
  cout << ".... B fit ...." << endl;
  Map map_b(loadPlots(file_b));
  assert(map_sb.size()==map_b.size());
  Map map_nom(loadPlots(file_nom));
  assert(map_sb.size()==map_nom.size());

  TFile* f = TFile::Open(rootfile);
  assert(f);
  TH1* histo_Data = (TH1*)f->Get("histo_Data");
  histo_Data->SetMarkerStyle(20);
  histo_Data->SetDirectory(0);
  
  std::map<TString,Color_t> color_map;
  color_map["Top"] = kYellow;
  color_map["VV"] = kAzure-2;
  color_map["WH"] = kRed+1;
  color_map["Wgamma"] = kBlue;
  color_map["Wg3l"] = kBlue;
  color_map["WjetsE"] = kGray+1;
  color_map["WjetsM"] = kGray+1;
  color_map["ZH"] = kRed+1;
  color_map["Zjets"] = kGreen+2;
  color_map["Ztt"] = kGreen+2;
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

  //for(Map::iterator itr = map_sb.begin(); itr!=map_sb.end(); ++itr,++iChannel){
  // if (nProcesses==0) nProcesses = itr->second.size();
  // assert(nProcesses==itr->second.size());
  // TString channel = itr->first;
  // printf("Channel %s, nProcesses: %u\n",channel.Data(),nProcesses);
  double maxValue(0);
  double maxInt(0);
  std::vector<StringPair> allProcesses;
  std::vector<StringPair> bkgProcesses;
  SubMap smap(map_sb[channel]);
  unsigned int nProcesses = smap.size();
  TCanvas* c = new TCanvas("c","",900,900);
  unsigned int sidex = int(sqrt(nProcesses));
  unsigned int sidey = sidex;
  if (sidex*sidey < nProcesses) sidex++;
  if (sidex*sidey < nProcesses) sidey++;
  c->Divide(sidex,sidey); 
  for(SubMap::const_iterator itr_p = smap.begin(); itr_p!=smap.end(); ++itr_p){
    double value = itr_p->second->GetMaximum();
    double integral = itr_p->second->Integral();
    if (value>maxValue) maxValue = value;
    if (integral>maxInt) maxInt = integral;
  }
  // std::cout << "maxInt: " << maxInt << std::endl;
  // const double cutoff = 0.10; // ignore contributions that are less than the threshold, which is a fraction with respect to the biggest contribution
  for(SubMap::const_iterator itr_p = smap.begin(); itr_p!=smap.end(); ++itr_p){
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
  for(SubMap::iterator itr_p = smap.begin(); itr_p!=smap.end(); ++itr_p){
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
  c->Print(Form("%s-allshapes.pdf",prefix_output));
  TCanvas* c2 = new TCanvas("c2","",1600,400);
  c2->Divide(4,1);
  for ( unsigned int j=0; j<4; ++j ){
    c2->cd(j+1);
    TString process = bkgProcesses.at(j).first;
    map_sb[channel][process]->Draw();
    map_b[channel][process]->Draw("same");
    map_nom[channel][process]->Draw("same");
    DrawLegends2(map_sb[channel][process],map_b[channel][process],map_nom[channel][process]);
  }
  c2->Print(Form("%s-mainbkg.pdf", prefix_output)); 
  TCanvas* c3 = new TCanvas("c3","",1500,500);
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
  double maxY = std::max(std::max(h_nom->GetMaximum(),histo_Data->GetMaximum()),
			   std::max(h_b->GetMaximum(),h_sb->GetMaximum()));
  maxY *= 1.2;
  c3->cd(1);
  h_nom->SetMaximum(maxY);
  h_nom->Draw("hist");
  DrawLegends(submap_nom,histo_Data);
  histo_Data->Draw("e same");
  c3->cd(2);
  h_b->SetMaximum(maxY);
  h_b->Draw("hist");
  DrawLegends(submap_b,histo_Data);
  histo_Data->Draw("e same");
  c3->cd(3);
  h_sb->SetMaximum(maxY);
  h_sb->Draw("hist");
  DrawLegends(submap_sb,histo_Data,true);
  histo_Data->Draw("e same");
  c3->Print(Form("%s-comp.pdf", prefix_output)); 

  //
  // 2D 
  // 
  // --- before fit
  TH1* hbkg = (TH1*)histo_Data->Clone(); 	hbkg->Reset();
  hbkg->SetTitle("Bkg before fit");
  hbkg->SetName("hbkg");
  TH1* histo_Top    = (TH1*)f->Get("histo_Top"); 		hbkg->Add(histo_Top);
  TH1* histo_VV     = (TH1*)f->Get("histo_VV"); 		hbkg->Add(histo_VV);
  TH1* histo_Wgamma = (TH1*)f->Get("histo_Wgamma");     hbkg->Add(histo_Wgamma);
  TH1* histo_Wg3l   = (TH1*)f->Get("histo_Wg3l");       hbkg->Add(histo_Wg3l);
  TH1* histo_WjetsE = (TH1*)f->Get("histo_WjetsE"); 	hbkg->Add(histo_WjetsE);
  TH1* histo_WjetsM = (TH1*)f->Get("histo_WjetsM"); 	hbkg->Add(histo_WjetsM);
  TH1* histo_Zjets  = (TH1*)f->Get("histo_Zjets"); 	    hbkg->Add(histo_Zjets);
  TH1* histo_ggWW   = (TH1*)f->Get("histo_ggWW"); 	    hbkg->Add(histo_ggWW);
  TH1* histo_qqWW   = (TH1*)f->Get("histo_qqWW");		hbkg->Add(histo_qqWW);
  TH1* histo_Ztt    = (TH1*)f->Get("histo_Ztt"); 	    hbkg->Add(histo_Ztt);
  
  TH1* hsig = (TH1*)histo_Data->Clone(); 	hsig->Reset();
  hsig->SetTitle("Signal before fit");
  hsig->SetName("hsig");
  TH1* histo_ggH = (TH1*)f->Get("histo_ggH"); 	hsig->Add(histo_ggH);
  TH1* histo_qqH = (TH1*)f->Get("histo_qqH"); 	hsig->Add(histo_qqH);
  TH1* histo_WH  = (TH1*)f->Get("histo_WH"); 	hsig->Add(histo_WH);
  TH1* histo_ZH  = (TH1*)f->Get("histo_ZH"); 	hsig->Add(histo_ZH);

  // --- after fit
  TH1* h1_s 		= (TH1*)f->Get("histo_ggH")->Clone();	h1_s->Reset();	
  TH1* h1_b 		= (TH1*)f->Get("histo_ggH")->Clone();	h1_b->Reset(); 
  TH1* h1_b_only 	= (TH1*)f->Get("histo_ggH")->Clone();	h1_b_only->Reset(); 
  for (std::vector<StringPair>::const_iterator iProcess = allProcesses.begin();
       iProcess != allProcesses.end(); ++iProcess){
    Color_t color = color_map[iProcess->first];
    TH1* hist = getACopy(map_sb[channel][iProcess->first],color);   

	// ----------------------------------------------------------------------
	//		Add systematics 
	// ----------------------------------------------------------------------
	float syst = 0.;
	float WWnorm[2] 	= { 0.0483, 0.09037	 };         // float WWnorm[2] 	= { 0.06153, 0.1281	 }; // HCP
	float Topnorm[2] 	= { 0.192, 0.03155 };           // float Topnorm[2] 	= { 0.20772, 0.03265 };  // HCP
	// instrumental : CMS_eff_m CMS_eff_e  CMS_scale_m CMS_scale_e CMS_hww_met_resolution  CMS_scale_j 
	// WW  
 	if ( iProcess->first=="qqWW" || iProcess->first=="ggWW") { 
		syst = TMath::Sqrt(syst*syst + WWnorm[njets]*WWnorm[njets]); 	
		syst = TMath::Sqrt(syst*syst + 0.03*0.03 + 0.04*0.04 + 0.015*0.015 + 0.02*0.02 + 0.02*0.02 + 0.02*0.02 ); // instrumental : ~ 6 %   	
		syst = TMath::Sqrt(syst*syst + 0.06*0.06); // QCDscale_WW_EXTRAP  	
		syst = TMath::Sqrt(syst*syst + 0.04*0.04); // pdf_gg or pdf_qqbar  	
	}
	// Top  
 	if ( iProcess->first=="Top") {
		syst = TMath::Sqrt(syst*syst + Topnorm[njets]*WWnorm[njets]); 	
		syst = TMath::Sqrt(syst*syst + 0.015*0.015 + 0.02*0.02 + 0.02*0.02 + 0.02*0.02 ); // instrumental : ~ 4 %  	
	}
	// Wjets 
 	if ( iProcess->first=="WjetsE" || iProcess->first=="WjetsM" ) {
		syst = TMath::Sqrt(syst*syst + 0.36*0.36); 	
	}
	// Wgamma
 	if ( iProcess->first=="Wgamma" || iProcess->first=="Wg3l") {
		syst = TMath::Sqrt(syst*syst + 0.4*0.4); 	
	}

	cout << "Syst: " << iProcess->first << " : " <<  syst << endl;

	// add stat uncert from the template before fit
	TH1* temp = (TH1*)f->Get(Form("histo_%s", (iProcess->first).Data()));
    unsigned int nbins  = hist->GetXaxis()->GetNbins();
    for(unsigned int ibin=1; ibin <= nbins; ++ibin) { 
			float systbin = syst * hist->GetBinContent(ibin);	

			if (iProcess->first=="qqWW") {  
				TH1* tempww = (TH1*)f->Get(Form("histo_%s_CMS_hww_MVAWWBoundingUp", (iProcess->first).Data()));   
				float diffbin = TMath::Abs(tempww->GetBinContent(ibin) - temp->GetBinContent(ibin)); 
				systbin = TMath::Sqrt(systbin*systbin + diffbin*diffbin);
				TH1* tempwwmcnlo = (TH1*)f->Get(Form("histo_%s_CMS_hww_MVAWWNLOBoundingUp", (iProcess->first).Data())); 
				diffbin = TMath::Abs(tempwwmcnlo->GetBinContent(ibin) - temp->GetBinContent(ibin)); 
				systbin = TMath::Sqrt(systbin*systbin + diffbin*diffbin);
			}

			if (iProcess->first=="Top") {
				TH1* temptop = (TH1*)f->Get(Form("histo_%s_CMS_hww_MVATopBoundingUp", (iProcess->first).Data())); 
				float diffbin = TMath::Abs(temptop->GetBinContent(ibin) - temp->GetBinContent(ibin)); 
				systbin = TMath::Sqrt(systbin*systbin + diffbin*diffbin);
			}

        	hist->SetBinError( ibin, TMath::Sqrt(systbin*systbin + temp->GetBinError(ibin)*temp->GetBinError(ibin)));
			//hist->SetBinError( ibin, TMath::Sqrt(systbin*systbin) );
	}
	// ----------------------------------------------------------------------

	if ( iProcess->first=="ggH" || iProcess->first=="qqH" || iProcess->first=="WH" || iProcess->first=="ZH" ) {
		h1_s->Add(hist);  
	}
	else {
		h1_b->Add(hist);  
    	hist = getACopy(map_b[channel][iProcess->first],color);
    	h1_b_only->Add(hist); 
	}
  }

  TFile *fout = new TFile(Form("%s_2D_postfit.root", prefix_output), "RECREATE"); // needs some editing to use this 
  //TFile *fout = new TFile("2D_postfit_125_of_0j.root", "RECREATE");

  TCanvas* c2d = new TCanvas("c2d","",2000,800);
  c2d->Divide(4,2);
  // ------ background before fit 
  c2d->cd(1);
  TH2* h2bkg = Roll1DTo2D(hbkg);  
  h2bkg->SetTitle("Background before fit");
  h2bkg->SetName("h2bkg");
  h2bkg->SetMarkerSize(markersize);
  h2bkg->Draw("colz text e");
  // ------ signal before fit 
  c2d->cd(2);
  TH2* h2sig = Roll1DTo2D(hsig);  
  h2sig->SetTitle("Signal before fit");
  h2sig->SetName("h2sig");
  h2sig->SetMarkerSize(markersize);
  h2sig->Draw("colz text e");
  // ------ Data
  c2d->cd(3);
  TH2* h2data = Roll1DTo2D(histo_Data);
  h2data->SetTitle("Data");
  h2data->SetName("h2data");
  h2data->SetMarkerSize(markersize);
  h2data->Draw("colz text e"); 
  // ------ Background only fit 
  c2d->cd(4);
//  TH2* h2_b_only = Roll1DTo2D(h1_b_only);
//  h2_b_only->SetTitle("Background after B-only fit");
//  h2_b_only->SetName("h2_b_only");
//  h2_b_only->SetMarkerSize(markersize);
//  h2_b_only->Draw("colz text e");
  // ------ background after fit 
  c2d->cd(5);
  TH2* h2_b = Roll1DTo2D(h1_b); 
  h2_b->SetTitle("Background after S+B fit");
  h2_b->SetName("h2_b");
  h2_b->SetMarkerSize(markersize);
  h2_b->Draw("colz text e");
  // ------ signal after fit 
  c2d->cd(6);
  TH2* h2_s = Roll1DTo2D(h1_s);  
  h2_s->SetTitle("Signal after S+B fit");
  h2_s->SetName("h2_s");
  h2_s->SetMarkerSize(markersize);
  h2_s->Draw("colz text e"); 
  // ------ Data - background fit 
  c2d->cd(7);
  TH2* h2_dataminusb = Roll1DTo2D(histo_Data);
  h2_dataminusb->SetTitle("Data - Background after S+B fit");
  h2_dataminusb->SetName("h2_dataminusb");
  h2_dataminusb->SetMarkerSize(markersize);
  h2_dataminusb->Add(h2_b, -1);
  h2_dataminusb->Draw("colz text e");
  // ------ Data - background  - sig after fit (S+B) 
  c2d->cd(8);
  TH2* h2_dataminusbminussig = Roll1DTo2D(histo_Data);
  h2_dataminusbminussig->SetTitle("Data - Background -Signal after S+B fit");
  h2_dataminusbminussig->SetName("h2_dataminusbminussig");
  h2_dataminusbminussig->SetMarkerSize(markersize);
  h2_dataminusbminussig->Add(h2_b, -1);
  h2_dataminusbminussig->Add(h2_s, -1);
  h2_dataminusbminussig->Draw("colz text e");
  
  c2d->Print(Form("%s-2d.pdf", prefix_output)); 
 
  // for PAS 
  TLatex* Lumi = new TLatex(90, 106, "CMS preliminary L = 19.5 fb^{-1} (8TeV)");
  Lumi ->SetTextAlign(22);
  Lumi ->SetTextFont(42);
  Lumi ->SetTextSize(.05);
  //Lumi ->SetNDC(1);
  TLatex* Lumi_zoom = new TLatex(120, 107, "CMS preliminary L = 19.5 fb^{-1} (8TeV)");
  Lumi_zoom->SetTextAlign(22);
  Lumi_zoom->SetTextFont(42);
  Lumi_zoom->SetTextSize(.06);
  //Lumi ->SetNDC(1);

  TCanvas* c2sig = new TCanvas("c2sig","",650,500);
  c2sig->cd(1);
  c2sig->SetTopMargin(0.15);
  //c2sig->SetRightMargin(0.15);
  //h2sig->SetTitle(Form("M_{H} = %i GeV", MASS))  
  TLatex* Title_sig = new TLatex(90, 113, Form("M_{H} = %i GeV", MASS));
  Title_sig->SetTextAlign(22);
  Title_sig->SetTextFont(42);
  Title_sig->SetTextSize(.06);;
  h2sig->SetTitle(""); 
  h2sig->SetTitleOffset(1.2,"Y");
  h2sig->SetXTitle("M_{T} (GeV)");
  h2sig->SetYTitle("M_{ll} (GeV)");
  h2sig->Draw("colz");
  Title_sig->Draw("same");
  Lumi ->Draw("same");
//  c2sig->Print(Form("plotPAS/2d_prefit_%ij_%i_sig_PAS.pdf", njets, MASS)); 
//  c2sig->Print(Form("plotPAS/2d_prefit_%ij_%i_sig_PAS.png", njets, MASS)); 
//  c2sig->Print(Form("plotPAS/2d_prefit_%ij_%i_sig_PAS.eps", njets, MASS)); 
  c2sig->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_sig_PAS_v1.pdf", njets, MASS)); 
  c2sig->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_sig_PAS_v1.png", njets, MASS)); 
  c2sig->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_sig_PAS_v1.eps", njets, MASS)); 
  c2sig->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_sig_PAS_v1.root", njets, MASS)); 
  
  TCanvas* c2data = new TCanvas("c2data","",650,500);
  c2data->cd(1); 
  c2data->SetTopMargin(0.15);
  //c2data->SetRightMargin(0.15);
  //h2data->SetTitle("Data"); 
  TLatex* Title_data = new TLatex(90, 113, "Data");
  Title_data->SetTextAlign(22);
  Title_data->SetTextFont(42);
  Title_data->SetTextSize(.06);
  h2data->SetTitle("");
  h2data->SetTitleOffset(1.2,"Y");
  h2data->SetXTitle("M_{T} (GeV)");
  h2data->SetYTitle("M_{ll} (GeV)");
  h2data->Draw("colz");
  Title_data->Draw("same");
  Lumi ->Draw("same");
//  c2data->Print(Form("plotPAS/2d_prefit_%ij_%i_data_PAS.pdf", njets, MASS)); 
//  c2data->Print(Form("plotPAS/2d_prefit_%ij_%i_data_PAS.png", njets, MASS)); 
//  c2data->Print(Form("plotPAS/2d_prefit_%ij_%i_data_PAS.eps", njets, MASS)); 
  c2data->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_data_PAS_v1.pdf", njets, MASS)); 
  c2data->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_data_PAS_v1.png", njets, MASS)); 
  c2data->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_data_PAS_v1.eps", njets, MASS)); 
  c2data->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_data_PAS_v1.root", njets, MASS)); 
  
  TCanvas* c2bkg = new TCanvas("c2bkg","",650,500);
  c2bkg->cd(1); 
  c2bkg->SetTopMargin(0.15);
  //c2bkg->SetRightMargin(0.15);
  //h2bkg->SetTitle("Background"); 
  TLatex* Title_bkg = new TLatex(90, 113, "Background");
  Title_bkg->SetTextAlign(22);
  Title_bkg->SetTextFont(42);
  Title_bkg->SetTextSize(.06);
  h2bkg->SetTitle("");
  h2bkg->SetTitleOffset(1.2,"Y");
  h2bkg->SetXTitle("M_{T} (GeV)");
  h2bkg->SetYTitle("M_{ll} (GeV)");
  h2bkg->Draw("colz");
  Title_bkg->Draw("same");
  Lumi ->Draw("same");
//  c2bkg->Print(Form("plotPAS/2d_prefit_%ij_%i_bkg_PAS.pdf", njets, MASS)); 
//  c2bkg->Print(Form("plotPAS/2d_prefit_%ij_%i_bkg_PAS.png", njets, MASS)); 
//  c2bkg->Print(Form("plotPAS/2d_prefit_%ij_%i_bkg_PAS.eps", njets, MASS)); 
  c2bkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_bkg_PAS_v1.pdf", njets, MASS)); 
  c2bkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_bkg_PAS_v1.png", njets, MASS)); 
  c2bkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_bkg_PAS_v1.eps", njets, MASS)); 
  c2bkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_bkg_PAS_v1.root", njets, MASS)); 
  
  TCanvas* c2sig_post = new TCanvas("c2sig_post","",650,500);
  c2sig_post->cd(1); 
  //c2sig_post->SetTopMargin(0.15);
  h2_s->SetMarkerSize(markersize*2);
  //h2_s->SetTitle("Signal post-fit");
  h2_s->SetTitle("");
  h2_s->SetTitleOffset(1.2,"Y");
  h2_s->SetXTitle("M_{T} (GeV)");
  h2_s->SetYTitle("M_{ll} (GeV)");
  h2_s->Draw("colz");
  Lumi ->Draw("same");
//  c2sig_post->Print(Form("plotPAS/2d_postfit_%ij_%i_sig_PAS.pdf", njets, MASS)); 
//  c2sig_post->Print(Form("plotPAS/2d_postfit_%ij_%i_sig_PAS.png", njets, MASS)); 
//  c2sig_post->Print(Form("plotPAS/2d_postfit_%ij_%i_sig_PAS.eps", njets, MASS)); 
  
  TCanvas* c2bkg_post = new TCanvas("c2bkg_post","",650,500);
  c2bkg_post->cd(1); 
  c2bkg_post->SetTopMargin(0.15);
  //c2bkg->SetRightMargin(0.15);
  h2_b->SetTitle("Background");
  h2_b->SetTitleOffset(1.2,"Y");
  h2_b->SetXTitle("M_{T} (GeV)");
  h2_b->SetYTitle("M_{ll} (GeV)");
  h2_b->Draw("colz");
  Lumi ->Draw("same");
//  c2bkg_post->Print(Form("plotPAS/2d_postfit_%ij_%i_bkg_PAS.pdf", njets, MASS)); 
//  c2bkg_post->Print(Form("plotPAS/2d_postfit_%ij_%i_bkg_PAS.png", njets, MASS)); 
//  c2bkg_post->Print(Form("plotPAS/2d_postfit_%ij_%i_bkg_PAS.eps", njets, MASS)); 
  
  TCanvas* c2dataminusbkg = new TCanvas("c2dataminusbkg","",600,500);
  c2dataminusbkg->cd(1); 
  h2_dataminusb->SetMarkerSize(markersize*1.5);
  c2dataminusbkg->SetTopMargin(0.15);
  //c2dataminusbkg->SetRightMargin(0.15);
  //h2_dataminusb->SetTitle("");  
  //gStyle->SetTitleFont(42);
  //gStyle->SetTitleFontSize(0.08);
  //gStyle->SetTitleX(0.2);
  //gStyle->SetTitleBorderSize(0.);
  //h2_dataminusb->SetTitle("Data - Background");
  //h2_dataminusb->SetTitleOffset(1.2,"Y");
  TLatex* Title_dataminusb = new TLatex(90, 113, "Data - Background");
  Title_dataminusb->SetTextAlign(22);
  Title_dataminusb->SetTextFont(42);
  Title_dataminusb->SetTextSize(.06);
  h2_dataminusb->SetXTitle("M_{T} (GeV)");
  h2_dataminusb->SetYTitle("M_{ll} (GeV)");
  h2_dataminusb->SetTitleOffset(1.2,"Y");
  h2_dataminusb->Draw("colz text e");
  h2_dataminusb->SetTitle("");
  Title_dataminusb->Draw("same");
  Lumi->Draw("same");
  //c2dataminusbkg->Print(Form("plotPAS/2d_postfit_%ij_%i_dataminusbkg_PAS.pdf", njets, MASS)); 
  //c2dataminusbkg->Print(Form("plotPAS/2d_postfit_%ij_%i_dataminusbkg_PAS.png", njets, MASS)); 
  //c2dataminusbkg->Print(Form("plotPAS/2d_postfit_%ij_%i_dataminusbkg_PAS.eps", njets, MASS)); 
  //c2dataminusbkg->Print(Form("plotPAS/2d_postfit_%ij_%i_dataminusbkg_PAS.root", njets, MASS)); 
  c2dataminusbkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_postfit_%ij_%i_dataminusbkg_PAS_v1.pdf", njets, MASS)); 
  c2dataminusbkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_postfit_%ij_%i_dataminusbkg_PAS_v1.png", njets, MASS)); 
  c2dataminusbkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_postfit_%ij_%i_dataminusbkg_PAS_v1.eps", njets, MASS)); 
  c2dataminusbkg->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_postfit_%ij_%i_dataminusbkg_PAS_v1.root", njets, MASS)); 


  // spin 2 2D plots 
  TCanvas* c2spin2 = new TCanvas("c2spin2","",600,500);
  c2spin2->cd(1); 
  c2spin2->SetTopMargin(0.15);
  
  TFile *File_spin2        = TFile::Open(Form("cards/hwwjcp_19p5fb/125/xwwof_%ij.input_8TeV.root",njets), "READ");
  TH1F *h1_ggH_spin2    = (TH1F*) File_spin2->Get("histo_ggH"); 
  TH2* h2_ggH_spin2     = Roll1DTo2D(h1_ggH_spin2);   
  //h2_ggH_spin2->SetTitle("2^{+}_{min} (125 GeV)"); 
  TLatex* Title_spin2 = new TLatex(90, 113, "2^{+}_{min} (125 GeV)");
  Title_spin2->SetTextAlign(22);
  Title_spin2->SetTextFont(42);
  Title_spin2->SetTextSize(.06);
  h2_ggH_spin2->SetTitle("");
  h2_ggH_spin2->SetTitleOffset(1.2,"Y");
  h2_ggH_spin2->SetName("h2sig");
  h2_ggH_spin2->SetMarkerSize(markersize);
  h2_ggH_spin2->Draw("colz");
  Title_spin2 ->Draw("same");
  Lumi ->Draw("same");
//  c2spin2->Print(Form("plotPAS/2d_prefit_%ij_%i_spin2_PAS.pdf", njets, MASS)); 
//  c2spin2->Print(Form("plotPAS/2d_prefit_%ij_%i_spin2_PAS.png", njets, MASS)); 
//  c2spin2->Print(Form("plotPAS/2d_prefit_%ij_%i_spin2_PAS.eps", njets, MASS)); 
  c2spin2->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_spin2_PAS_v1.pdf", njets, MASS)); 
  c2spin2->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_spin2_PAS_v1.png", njets, MASS)); 
  c2spin2->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_spin2_PAS_v1.eps", njets, MASS)); 
  c2spin2->Print(Form("/afs/cern.ch/user/j/jaehyeok/www/HWW/Moriond2013/PAS/plots/2d_prefit_%ij_%i_spin2_PAS_v1.root", njets, MASS)); 

  fout->Write(); 

}
