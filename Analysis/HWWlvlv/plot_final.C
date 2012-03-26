#include "TStyle.h"
void plot_final(int nsel = 0, int ReBin = 1, char XTitle[300] = "MVA Output", char units[300] = "", 
               char plotName[300] = "histo_tmva_ntuples_120train_1jets_hww120_chan4.root", char outputName[300] = "njets",
               bool isLogY = false, char MassH[300] = "160") {
   
  gROOT->SetStyle("Plain");
  //setTDRStyle();
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  TCanvas* cv = new TCanvas("cv", "cv", 700, 700);
  cv->Divide(1,3);

  TFile *_file0;
  
  char myRootFile[300];
  sprintf(myRootFile,"%s",plotName);
  cout << myRootFile<<endl;
  TFile *_file0 = TFile::Open(myRootFile);

  char YTitle[300];
  double binw1 = histoB->GetBinWidth(1);
  int	 binw2 = histoB->GetBinWidth(1);
  if(binw1 == binw2	) sprintf(YTitle,"Events",binw1,units);
  else  		  sprintf(YTitle,"Events / %4.2f %s",binw1,units);
  histoB->SetLabelSize(0.05, "X");
  histoB->SetTitleSize(0.06, "X");
  histoB->SetLabelSize(0.05, "Y");
  histoB->SetTitleSize(0.06, "Y");
  histoB->GetXaxis()->SetTitleFont(132);
  histoB->GetXaxis()->SetLabelFont(132);
  histoB->GetYaxis()->SetTitleFont(132); 
  histoB->GetYaxis()->SetLabelFont(132); 
  histoB->GetXaxis()->SetTitleOffset(0.9);
  histoB->GetYaxis()->SetTitleOffset(1.3);
  histoB->Rebin(ReBin);

  histoS->SetLabelSize(0.05, "X");
  histoS->SetTitleSize(0.06, "X");
  histoS->SetLabelSize(0.05, "Y");
  histoS->SetTitleSize(0.06, "Y");
  histoS->GetXaxis()->SetTitleFont(132);
  histoS->GetXaxis()->SetLabelFont(132);
  histoS->GetYaxis()->SetTitleFont(132); 
  histoS->GetYaxis()->SetLabelFont(132); 
  histoS->GetXaxis()->SetTitleOffset(0.9);
  histoS->GetYaxis()->SetTitleOffset(1.3);
  histoS->Scale(50.0);
  histoS->Rebin(ReBin);

  histoD->Rebin(ReBin);

  TH1F* histoB_empty =  histoB->Clone();
  histoB_empty->Scale(0.0);
  double max = histoB->GetMaximum();
  if(histoD->GetMaximum() > max) max = histoD->GetMaximum();
  if(histoD->GetMaximum() > max) max = histoD->GetMaximum();
  if(histoS->GetMaximum() > max) max = histoS->GetMaximum();
  if(histoS->GetMaximum() > max) max = histoS->GetMaximum();
  histoB->SetMaximum(max*1.3);
  histoB->SetMinimum(0.001);
  histoS->SetMaximum(max*1.3);
  histoS->SetMinimum(0.001);
  histoB_empty->SetMaximum(max*1.3);
  histoB_empty->SetMinimum(0.001);
  if(isLogY == true){
    histoB->SetMinimum(0.01);
    histoS->SetMinimum(0.01);
    histoB_empty->SetMinimum(0.01);
    histoB_empty->SetMaximum(max*3.8);
  }
  
  cv->cd(1);
  histoB_empty->Draw("hist");
  histoB->Draw("hist,same");
  if	 (nsel == 2 ||  nsel == 3 ||  nsel == 4){
    //histoB_empty->Draw("hist");
    histoS->Draw("hist,same");
  }
  //histoB->Draw("hist,same");
  histoD->Draw("e,same");

  labelcms  = new TPaveText(0.28,0.90,0.28,0.90,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.05);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS, #sqrt{s} = 7 TeV, L_{ int} = 1.1 fb^{-1}");
  labelcms->SetBorderSize(0);
  labelcms->SetTextFont(132);
  labelcms->SetLineWidth(2);
  labelcms->Draw();

  TLegend* leg;
  leg = new TLegend(0.70,0.70,0.85,0.85); 
  leg ->SetFillStyle(0);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(0);
  leg ->SetTextSize(0.04);
  leg ->SetTextFont(132);
  leg ->SetLineWidth(1);
  
  if(nsel < 31){
    leg ->AddEntry(histoD,"data","P"); 
    if(nsel == 2 || nsel == 3 || nsel == 4){
     char hwwTitle[300];
     sprintf(hwwTitle,"H(%s) #rightarrow WW",MassH);      
      leg ->AddEntry(histoS,hwwTitle,"F");
    }
    leg ->AddEntry(histoB,"background","F");
  }

  leg ->Draw();

  double S = 0.0;
  double B = 0.0;
  double fsig = 0.30;
  cv->cd(2);
  TH1D *hDSignif1 = new TH1D("hDSignif1","hDSignif1",histoS->GetNbinsX(),
					    -0.5,histoS->GetNbinsX()-0.5);
  for(int i=1; i<=histoS->GetNbinsX(); i++){
    S = S + histoS->GetBinContent(i);
    B = B + histoB->GetBinContent(i);
    if(B>0){
      double var = TMath::Sqrt(B+fsig*fsig*B*B);
      //var = TMath::Sqrt(S+B);
      hDSignif1->SetBinContent(i,S/var);
    }
    else {
      hDSignif1->SetBinContent(i,0.0);
    }
  }
  hDSignif1->SetDirectory(0);
  hDSignif1->Draw();

  cv->cd(3);
  S = 0.0;
  B = 0.0;
  TH1D *hDSignif2 = new TH1D("hDSignif2","hDSignif2",histoS->GetNbinsX(),
					    -0.5,histoS->GetNbinsX()-0.5);
  for(int i=histoS->GetNbinsX(); i>=1; i--){
    S = S + histoS->GetBinContent(i);
    B = B + histoB->GetBinContent(i);
    if(B>0){
      double var = TMath::Sqrt(B+fsig*fsig*B*B);
      //var = TMath::Sqrt(S+B);
      hDSignif2->SetBinContent(i,S/var);
    }
    else {
      hDSignif2->SetBinContent(i,0.0);
    }
  }
  hDSignif2->SetDirectory(0);
  hDSignif2->Draw();
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600);
  tdrStyle->SetCanvasDefW(600);
  tdrStyle->SetCanvasDefX(0);
  tdrStyle->SetCanvasDefY(0);

  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);

  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  tdrStyle->SetMarkerSize(2);

  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  tdrStyle->SetOptDate(0);

  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); 
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.027);

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(20, "XYZ");
  //tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXSize(0.02); 
  tdrStyle->SetTitleYSize(0.02);
  tdrStyle->SetTitleXOffset(1.5);
  tdrStyle->SetTitleYOffset(1.7);
  // tdrStyle->SetTitleOffset(1.1, "XYZ"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(18, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->cd();
}
