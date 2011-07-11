// R. Gonzalez & G. Gomez-Ceballos
// macro used for HWW
#include "TStyle.h"
// nsel = 0 --> only backgrounds
//        1 --> only for ww cut evolution
//        2 --> signal on top of background
//        3 --> signal separated from background
//        4 --> signal on top of background and showing cut evolution
void plot_data(int nsel = 0, int ReBin = 10, char XTitle[300] = "N_{jets}", char units[300] = "", 
                char plotName[300] = "histo_tmva_ntuples_160train_0jets_hww160_chan4", char outputName[300] = "njets",
                bool isLogY = false, char MassH[300] = "160") {
   
  gROOT->SetStyle("Plain");
  setTDRStyle();
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);

  TFile *_file0;
  
  char myRootFile[300];
  sprintf(myRootFile,"%s.root",plotName);
  cout << myRootFile<<endl;
  TFile *_file0 = TFile::Open(myRootFile);

  TCanvas* c1 = new TCanvas();
  if(isLogY == true) c1->SetLogy();

  if(nsel < 40){   
    histoS->Rebin(ReBin);
    histoB->Rebin(ReBin);
    histoD->Rebin(ReBin);

    histoB->SetFillColor(kBlue);
    histoB->SetFillStyle(1001);
    histoB->SetLineColor(kBlue);
    histoB->SetLineStyle(0);
    histoB->SetLineWidth(0);
    histoB->SetTitle("");
   
    histoD->Sumw2();
    histoD->SetMarkerColor(kBlack);
    histoD->SetMarkerStyle(20);
    histoD->SetLineWidth(2);
    histoD->SetMarkerStyle(kFullDotLarge);
    histoD->SetTitle("");

    if(nsel == 2 || nsel == 4){
      histoS->Sumw2();
      histoS->SetFillColor(kWhite);
      histoS->SetFillStyle(1001);
      histoS->SetLineStyle(0);
      histoS->SetLineWidth(3);
      histoS->SetTitle("");
    }
    else {
      TH1F* histo_B =  histoB->Clone();
      histoS->Add(histo_B);
      histoS->Sumw2();
      histoS->SetFillColor(kBlack);
      histoS->SetFillStyle(1);
      histoS->SetLineStyle(1);
      histoS->SetLineWidth(4);
      histoS->SetTitle("");
    }

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

    char YTitle[300];
    double binw1 = histoB->GetBinWidth(1);
    int    binw2 = histoB->GetBinWidth(1);
    if     (nsel == 2 || nsel == 4) sprintf(YTitle,"Events",binw1,units);
    else if(binw1 == binw2        ) sprintf(YTitle,"Events",binw1,units);
    else                            sprintf(YTitle,"Events / %4.2f %s",binw1,units);

    histoS->SetYTitle(YTitle);
    histoB->SetYTitle(YTitle);
    histoS->SetXTitle(XTitle);
    histoB->SetXTitle(XTitle);
     
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
    if     (nsel == 2 ||  nsel == 3 ||  nsel == 4){
      //histoB_empty->Draw("hist");
      histoS->Draw("hist,same");
    }        
    //histoB->Draw("hist,same");
    histoD->Draw("e,same");
  } else {
    printf("Wrong option: %d\n",nsel);
    return;
  }

  labelcms  = new TPaveText(0.28,0.90,0.28,0.90,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.05);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS, #sqrt{s} = 7 TeV, L_{ int} = 188 pb^{-1}");
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
  char myOutputFile[300];

  //sprintf(myOutputFile,"plots/%s.eps",outputName);
  //c1->SaveAs(myOutputFile);
  //sprintf(myOutputFile,"plots/%s.png",outputName);
  //c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"plots/%s.pdf",outputName);
  c1->SaveAs(myOutputFile);
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
