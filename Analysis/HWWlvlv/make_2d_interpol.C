#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <TH1D.h>
#include <TGraphErrors.h>

Double_t funPOL(Double_t* x, Double_t* par);

void make_2d_interpol(int Njet = 0, int ECM = 8){
  char finalStateName[2];
  sprintf(finalStateName,"of");

  const int numberOfBins = 126;
  const int nmass = 14;
  const double mH[nmass]   = {110,115,120,125,130,135,140,145,150,160,170,180,190,200};
  const double mHE[nmass]  = {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
  double effmea[nmass],effmeaE[nmass];

  TCanvas *canvas[4][numberOfBins];
  TGraphErrors* gr[4][numberOfBins];
  TH1D* histoI[4][nmass];
  TFile *_fileInput[nmass];
  double eff[4][nmass][numberOfBins],effE[4][nmass][numberOfBins];
  
  for(int i=0; i<nmass; i++){
    _fileInput[i] = TFile::Open(Form("/data/smurf/ceballos/inputLimits/ana_EPS/%d/hww%s_%dj.input_%dTeV.root",(int)mH[i],finalStateName,Njet,ECM));
    _fileInput[i]->cd();
    histoI[0][i] = (TH1D*) _fileInput[i]->Get("histo_ggH");
    histoI[1][i] = (TH1D*) _fileInput[i]->Get("histo_qqH");
    histoI[2][i] = (TH1D*) _fileInput[i]->Get("histo_WH");
    histoI[3][i] = (TH1D*) _fileInput[i]->Get("histo_ZH");
    histoI[0][i]->SetDirectory(0);
    histoI[1][i]->SetDirectory(0);
    histoI[2][i]->SetDirectory(0);
    histoI[3][i]->SetDirectory(0);
    histoI[0][i]->Scale(1.0/histoI[0][i]->GetSumOfWeights());
    histoI[1][i]->Scale(1.0/histoI[1][i]->GetSumOfWeights());
    histoI[2][i]->Scale(1.0/histoI[2][i]->GetSumOfWeights());
    histoI[3][i]->Scale(1.0/histoI[3][i]->GetSumOfWeights());
    for(Int_t nb=1;nb<=numberOfBins;nb++){
      eff[0][i][nb-1] = histoI[0][i]->GetBinContent(nb); effE[0][i][nb-1] = histoI[0][i]->GetBinError(nb);
      eff[1][i][nb-1] = histoI[1][i]->GetBinContent(nb); effE[1][i][nb-1] = histoI[1][i]->GetBinError(nb);
      eff[2][i][nb-1] = histoI[2][i]->GetBinContent(nb); effE[2][i][nb-1] = histoI[2][i]->GetBinError(nb);
      eff[3][i][nb-1] = histoI[3][i]->GetBinContent(nb); effE[3][i][nb-1] = histoI[3][i]->GetBinError(nb);
    }
    _fileInput[i]->Close();
  }

  if(numberOfBins != histoI[0][0]->GetNbinsX()) assert(0);

  for(int nchan=0;nchan<4;nchan++){
    for(int nb=1;nb<=numberOfBins;nb++){
      for(int i=0; i<nmass; i++){
        effmea[i]  = eff[nchan][i][nb];
        effmeaE[i] = effE[nchan][i][nb];
      }
      TF1 *func = new TF1("func",funPOL,mH[0],mH[nmass],4);
      canvas[nchan][nb] = new TCanvas(Form("c_%d_%d",nchan,nb),Form("c_%d_%d",nchan,nb),500,500);
      canvas[nchan][nb]->cd();
      gr[nchan][nb] = new TGraphErrors(nmass, mH, effmea, mHE, effmeaE);
      gr[nchan][nb]->SetTitle(Form("c_%d_%d",nchan,nb));
      gr[nchan][nb]->SetMarkerStyle(20);
      gr[nchan][nb]->SetMarkerColor(2);
      gr[nchan][nb]->GetXaxis()->SetTitle("m_H [GeV}");
      gr[nchan][nb]->GetXaxis()->SetTitleOffset(0.7);
      gr[nchan][nb]->GetXaxis()->SetTitleSize(0.065);
      gr[nchan][nb]->GetXaxis()->SetLabelSize(0.065);
      gr[nchan][nb]->GetYaxis()->SetTitle("Fraction");
      gr[nchan][nb]->GetYaxis()->SetTitleOffset(1.3);
      gr[nchan][nb]->GetYaxis()->SetTitleSize(0.05);
      gr[nchan][nb]->GetYaxis()->SetLabelSize(0.05);
      gr[nchan][nb]->GetYaxis()->CenterTitle(kTRUE);
      gr[nchan][nb]->Fit(func,"rome","e",mH[0],mH[nmass-1]);
      gr[nchan][nb]->Draw("A*");
      canvas[nchan][nb]->Update();
      canvas[nchan][nb]->SaveAs(Form("plots/make_2d_interpol_%d_%d.png",nchan,nb));
    }
  }
}
Double_t funPOL(Double_t* x, Double_t* par) {

  Double_t a = par[0];
  Double_t b = par[1];
  Double_t c = par[2];
  Double_t d = par[3];

  return a + b*x[0] + c*x[0]*x[0] + d*x[0]*x[0]*x[0];
}
