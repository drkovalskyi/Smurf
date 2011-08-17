#include "TStyle.h"
void plot_syst(int nsel = 0, int ReBin = 10, char XTitle[300] = "MVA Output", char units[300] = "", 
               char plotName[300] = "histo_tmva_ntuples_130train_0jets_hww130_chan5.root", char outputName[300] = "njets",
               bool isLogY = false, char MassH[300] = "160") {
   
  gROOT->SetStyle("Plain");
  setTDRStyle();
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);

  TFile *_file0;
  
  char myRootFile[300];
  sprintf(myRootFile,"%s",plotName);
  cout << myRootFile<<endl;
  TFile *_file0 = TFile::Open(myRootFile);

  TCanvas* c1 = new TCanvas();
  if(isLogY == true) c1->SetLogy();

  char YTitle[300];
  double binw1 = histoB->GetBinWidth(1);
  int	 binw2 = histoB->GetBinWidth(1);
  if(binw1 == binw2	) sprintf(YTitle,"Events",binw1,units);
  else  		  sprintf(YTitle,"Events / %4.2f %s",binw1,units);

  printf("kRed(%d) == lepton eff +\n",kBlack+0+1);
  printf("kGreen(%d) == lepton eff -\n",kBlack+1+1);
  printf("kBlue(%d) == lepton smear +\n",kBlack+2+1);
  printf("kYellow(%d) == lepton smear -\n",kBlack+3+1);
  printf("kPink(%d) == met smear\n",kBlack+4+1);
  
  TH1D *sigMVA[5][6],*sigMVASyst[5][6],*bgdMVADecays[5][8],*bgdMVADecaysSyst[5][8];
  for(int i=0; i<6; i++){
    sigMVA[4][i] = (TH1D*)_file0->Get(Form("sigMVA_4_%d",i));
    sigMVA[4][i]->Rebin(ReBin);
    sigMVA[4][i]->SetMarkerColor(kBlack);
    sigMVA[4][i]->SetMarkerStyle(20);
    sigMVA[4][i]->SetLineWidth(2);
    sigMVA[4][i]->SetMarkerStyle(kFullDotLarge);
    sigMVA[4][i]->SetTitle("");
    sigMVA[4][i]->SetLabelSize(0.05, "X");
    sigMVA[4][i]->SetTitleSize(0.06, "X");
    sigMVA[4][i]->SetLabelSize(0.05, "Y");
    sigMVA[4][i]->SetTitleSize(0.06, "Y");
    sigMVA[4][i]->GetXaxis()->SetTitleFont(132);
    sigMVA[4][i]->GetXaxis()->SetLabelFont(132);
    sigMVA[4][i]->GetYaxis()->SetTitleFont(132); 
    sigMVA[4][i]->GetYaxis()->SetLabelFont(132); 
    sigMVA[4][i]->GetXaxis()->SetTitleOffset(0.9);
    sigMVA[4][i]->GetYaxis()->SetTitleOffset(1.3);
    sigMVA[4][i]->SetYTitle(YTitle);
    sigMVA[4][i]->SetXTitle(XTitle);
    for(int j=0; j<5; j++){
      sigMVASyst[j][i] = (TH1D*)_file0->Get(Form("sigMVASyst_%d_%d",j,i));
      sigMVASyst[j][i]->Rebin(ReBin);
      sigMVASyst[j][i]->SetMarkerColor(kBlack+j+1);
      sigMVASyst[j][i]->SetMarkerStyle(20);
      sigMVASyst[j][i]->SetLineWidth(2);
      sigMVASyst[j][i]->SetMarkerStyle(kFullDotLarge);
      sigMVASyst[j][i]->SetTitle("");
      sigMVASyst[j][i]->SetLabelSize(0.05, "X");
      sigMVASyst[j][i]->SetTitleSize(0.06, "X");
      sigMVASyst[j][i]->SetLabelSize(0.05, "Y");
      sigMVASyst[j][i]->SetTitleSize(0.06, "Y");
      sigMVASyst[j][i]->GetXaxis()->SetTitleFont(132);
      sigMVASyst[j][i]->GetXaxis()->SetLabelFont(132);
      sigMVASyst[j][i]->GetYaxis()->SetTitleFont(132); 
      sigMVASyst[j][i]->GetYaxis()->SetLabelFont(132); 
      sigMVASyst[j][i]->GetXaxis()->SetTitleOffset(0.9);
      sigMVASyst[j][i]->GetYaxis()->SetTitleOffset(1.3);
      sigMVASyst[j][i]->SetYTitle(YTitle);
      sigMVASyst[j][i]->SetXTitle(XTitle);
    }
  }
  for(int i=0; i<8; i++) {
    bgdMVADecays[4][i] = (TH1D*)_file0->Get(Form("bgdMVADecays_4_%d",i));
    bgdMVADecays[4][i]->Rebin(ReBin);
    bgdMVADecays[4][i]->SetMarkerColor(kBlack);
    bgdMVADecays[4][i]->SetMarkerStyle(20);
    bgdMVADecays[4][i]->SetLineWidth(2);
    bgdMVADecays[4][i]->SetMarkerStyle(kFullDotLarge);
    bgdMVADecays[4][i]->SetTitle("");
    bgdMVADecays[4][i]->SetLabelSize(0.05, "X");
    bgdMVADecays[4][i]->SetTitleSize(0.06, "X");
    bgdMVADecays[4][i]->SetLabelSize(0.05, "Y");
    bgdMVADecays[4][i]->SetTitleSize(0.06, "Y");
    bgdMVADecays[4][i]->GetXaxis()->SetTitleFont(132);
    bgdMVADecays[4][i]->GetXaxis()->SetLabelFont(132);
    bgdMVADecays[4][i]->GetYaxis()->SetTitleFont(132); 
    bgdMVADecays[4][i]->GetYaxis()->SetLabelFont(132); 
    bgdMVADecays[4][i]->GetXaxis()->SetTitleOffset(0.9);
    bgdMVADecays[4][i]->GetYaxis()->SetTitleOffset(1.3);
    bgdMVADecays[4][i]->SetYTitle(YTitle);
    bgdMVADecays[4][i]->SetXTitle(XTitle);
    for(int j=0; j<5; j++) {
      bgdMVADecaysSyst[j][i] = (TH1D*)_file0->Get(Form("bgdMVADecaysSyst_%d_%d",j,i));
      bgdMVADecaysSyst[j][i]->Rebin(ReBin);
      bgdMVADecaysSyst[j][i]->SetMarkerColor(kBlack+j+1);
      bgdMVADecaysSyst[j][i]->SetMarkerStyle(20);
      bgdMVADecaysSyst[j][i]->SetLineWidth(2);
      bgdMVADecaysSyst[j][i]->SetMarkerStyle(kFullDotLarge);
      bgdMVADecaysSyst[j][i]->SetTitle("");
      bgdMVADecaysSyst[j][i]->SetLabelSize(0.05, "X");
      bgdMVADecaysSyst[j][i]->SetTitleSize(0.06, "X");
      bgdMVADecaysSyst[j][i]->SetLabelSize(0.05, "Y");
      bgdMVADecaysSyst[j][i]->SetTitleSize(0.06, "Y");
      bgdMVADecaysSyst[j][i]->GetXaxis()->SetTitleFont(132);
      bgdMVADecaysSyst[j][i]->GetXaxis()->SetLabelFont(132);
      bgdMVADecaysSyst[j][i]->GetYaxis()->SetTitleFont(132); 
      bgdMVADecaysSyst[j][i]->GetYaxis()->SetLabelFont(132); 
      bgdMVADecaysSyst[j][i]->GetXaxis()->SetTitleOffset(0.9);
      bgdMVADecaysSyst[j][i]->GetYaxis()->SetTitleOffset(1.3);
      bgdMVADecaysSyst[j][i]->SetYTitle(YTitle);
      bgdMVADecaysSyst[j][i]->SetXTitle(XTitle);
    }
  }

  if      (nsel <10){
    sigMVA[4][nsel]->SetMaximum(sigMVA[4][nsel]->GetMaximum()*1.3);
    sigMVA[4][nsel]->SetMinimum(0.001);
    sigMVA[4][nsel]->Draw("e,hist");
    sigMVASyst[0][nsel]->Draw("e,same,hist");
    sigMVASyst[1][nsel]->Draw("e,same,hist");
    sigMVASyst[2][nsel]->Draw("e,same,hist");
    sigMVASyst[3][nsel]->Draw("e,same,hist");
    sigMVASyst[4][nsel]->Draw("e,same,hist");
  }
  else if(nsel <20){
    sigMVASyst[0][nsel-10]->SetMaximum(1.3);
    sigMVASyst[0][nsel-10]->SetMinimum(0.7);
    sigMVASyst[0][nsel-10]->Divide(sigMVA[4][nsel-10]);
    sigMVASyst[1][nsel-10]->Divide(sigMVA[4][nsel-10]);
    sigMVASyst[2][nsel-10]->Divide(sigMVA[4][nsel-10]);
    sigMVASyst[3][nsel-10]->Divide(sigMVA[4][nsel-10]);
    sigMVASyst[4][nsel-10]->Divide(sigMVA[4][nsel-10]);
    sigMVASyst[0][nsel-10]->Draw("e     ,hist");
    sigMVASyst[1][nsel-10]->Draw("e,same,hist");
    sigMVASyst[2][nsel-10]->Draw("e,same,hist");
    sigMVASyst[3][nsel-10]->Draw("e,same,hist");
    sigMVASyst[4][nsel-10]->Draw("e,same,hist");
  }
  else if(nsel <30){
    bgdMVADecays[4][nsel-20]->SetMaximum(bgdMVADecays[4][nsel-20]->GetMaximum()*1.3);
    bgdMVADecays[4][nsel-20]->SetMinimum(0.001);
    bgdMVADecays[4][nsel-20]->Draw("e,hist");
    bgdMVADecaysSyst[0][nsel-20]->Draw("e,same,hist");
    bgdMVADecaysSyst[1][nsel-20]->Draw("e,same,hist");
    bgdMVADecaysSyst[2][nsel-20]->Draw("e,same,hist");
    bgdMVADecaysSyst[3][nsel-20]->Draw("e,same,hist");
    bgdMVADecaysSyst[4][nsel-20]->Draw("e,same,hist");
  }
  else if(nsel <40){
    bgdMVADecaysSyst[0][nsel-30]->Divide(bgdMVADecays[4][nsel-30]);
    bgdMVADecaysSyst[1][nsel-30]->Divide(bgdMVADecays[4][nsel-30]);
    bgdMVADecaysSyst[2][nsel-30]->Divide(bgdMVADecays[4][nsel-30]);
    bgdMVADecaysSyst[3][nsel-30]->Divide(bgdMVADecays[4][nsel-30]);
    bgdMVADecaysSyst[4][nsel-30]->Divide(bgdMVADecays[4][nsel-30]);
    bgdMVADecaysSyst[0][nsel-30]->SetMaximum(1.3);
    bgdMVADecaysSyst[0][nsel-30]->SetMinimum(0.7);
    bgdMVADecaysSyst[0][nsel-30]->Draw("e     ,hist");
    bgdMVADecaysSyst[1][nsel-30]->Draw("e,same,hist");
    bgdMVADecaysSyst[2][nsel-30]->Draw("e,same,hist");
    bgdMVADecaysSyst[3][nsel-30]->Draw("e,same,hist");
    bgdMVADecaysSyst[4][nsel-30]->Draw("e,same,hist");
  }

  char myOutputFile[300];
  //sprintf(myOutputFile,"plots/%s.eps",outputName);
  //c1->SaveAs(myOutputFile);
  //sprintf(myOutputFile,"plots/%s.png",outputName);
  //c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"plots/%s.pdf",outputName);
  c1->SaveAs(myOutputFile);

  return;
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
