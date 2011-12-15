#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TLine.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "tdrstyle.C"

// Can put multiple filenames in LimitFiles and labels in LimitNames separated by comma
// A latex table with additional column for observed from an additional file is made 
//   in the same style as makeTable.sh
// Ratio: 0 = don't take ratio, 1 = simple ratio, 2 = ratio in units of sigma

void PlotLimit(string LimitFiles = "ana_v7_Full2011_limits_nj_shape.txt",
	       string filePrefix = "comb",
	       string LimTitle   = "H #rightarrow WW (BDT based)",
	       bool   DoObsLim   = 1,
	       bool   DoExpLim   = 1,
	       int    Ratio      = 0,
	       string lumi       = "4.6 fb^{-1}",
	       string LimitNames = "",
	       string LimTable   = "",
	       string TabLegend  = "",
	       string TabChannel = "",
	       string TabLabel   = "")
{
  setTDRStyle();

  // Get input files
  stringstream fnames(LimitFiles);
  string fname;
  vector<string> LimitFile;    
  while ( getline(fnames, fname, ',') ) {
    printf("Using files %s\n", fname.c_str());
    LimitFile.push_back(fname);
  } 
  stringstream flabels(LimitNames);
  string flabel;
  vector<string> LimitName;    
  while ( getline(flabels, flabel, ',') ) {
    printf("Using file corresponding to %s\n", flabel.c_str());
    LimitName.push_back(flabel);
  } 
  if (LimitFile.size()>4) { printf("Not setup to do %d-d comparison\n", (Int_t)LimitFile.size()); return; }

  vector<float> vMass           ;
  vector<float> vObsLimit       ; 
  vector<float> vMeanExpLimit   ; 
  vector<float> vMedianExpLimit ; 
  vector<float> vExpLim68Down   ; 
  vector<float> vExpLim68Up     ; 
  vector<float> vExpLim95Down   ; 
  vector<float> vExpLim95Up     ;

  float Mass           ;
  float ObsLimit       ; 
  float MeanExpLimit   ; 
  float MedianExpLimit ; 
  float ExpLim68Down   ; 
  float ExpLim68Up     ; 
  float ExpLim95Down   ; 
  float ExpLim95Up     ;

  // Read in the nominal limits
  ifstream indata;
  indata.open(LimitFile[0].c_str());
  if(!indata) { // file couldn't be opened
    cerr << "Error: file could not be opened" << endl;
    return;
  }
  while ( indata >> Mass >> ObsLimit >> MeanExpLimit >> MedianExpLimit >> ExpLim95Down >> ExpLim68Down >> ExpLim68Up >> ExpLim95Up )  {
    cout << Mass << " " << ObsLimit  << " " << MedianExpLimit <<" "<< ExpLim68Down <<" "<< ExpLim68Up <<" "<< ExpLim95Down <<" "<< ExpLim95Up << endl;

    vMass           .push_back(Mass           );
    vObsLimit       .push_back(ObsLimit       ); 
    vMeanExpLimit   .push_back(MeanExpLimit   ); 
    vMedianExpLimit .push_back(MedianExpLimit ); 
    vExpLim68Down   .push_back(ExpLim68Down   ); 
    vExpLim68Up     .push_back(ExpLim68Up     ); 
    vExpLim95Down   .push_back(ExpLim95Down   ); 
    vExpLim95Up     .push_back(ExpLim95Up     );
  }

  TCanvas* cLimit = new TCanvas("cLimit", "cLimit");

  cLimit->SetLeftMargin  (1.30 * cLimit->GetLeftMargin());
  cLimit->SetRightMargin (2.10 * cLimit->GetRightMargin());
  cLimit->SetBottomMargin(1.35 * cLimit->GetBottomMargin());

  float x1 = vMass.at(0) - 5.;
  float x2 = vMass.at(vMass.size()-1) + 5.; 

  // Expected Limit
  TGraph* ExpLim = NULL ;
  TGraphAsymmErrors* ExpBand68 = NULL ;
  TGraphAsymmErrors* ExpBand95 = NULL ;
  float min = 999999., max = 0;
  if ( DoExpLim ) {
    float x[100];
    float ex[100];
    float y[100];
    float yu68[100];
    float yd68[100];
    float yu95[100];
    float yd95[100]; 
    for ( int i = 0 ; i < (signed) vMass.size() ; ++i ) {
      x[i] = vMass.at(i) ; ex[i] = 0 ; 
      y[i] = vMedianExpLimit.at(i);
      if (Ratio==1) y[i] /= vMedianExpLimit.at(i);  
      if (Ratio==2) y[i] -= vMedianExpLimit.at(i);  
      if(y[i]    > max) max = y[i]   ; if(y[i]    < min) min = y[i]   ;
      yu68[i] = vExpLim68Up.at(i) -y[i];
      if (Ratio>=1) yu68[i] /= vMedianExpLimit.at(i); 
      if(yu68[i] > max) max = yu68[i]; if(yu68[i] < min) min = yu68[i];
      yd68[i] = y[i] - vExpLim68Down.at(i);  
      if (Ratio>=1) yd68[i] /= vMedianExpLimit.at(i); 
      if(yd68[i] > max) max = yd68[i]; if(yd68[i] < min) min = yd68[i];
      yu95[i] = vExpLim95Up.at(i) -y[i];   
      if (Ratio>=1) yu95[i] /= vMedianExpLimit.at(i); 
      if(yu95[i] > max) max = yu95[i]; if(yu95[i] < min) min = yu95[i];
      yd95[i] = y[i] - vExpLim95Down.at(i);  
      if (Ratio>=1) yd95[i] /= vMedianExpLimit.at(i); 
      if(yd95[i] > max) max = yd95[i]; if(yd95[i] < min) min = yd95[i];
      if (Ratio==1) {
        yu68[i] = vExpLim68Up.at(i) / vMedianExpLimit.at(i) - 1;  
        yd68[i] = 1 - vExpLim68Down.at(i) / vMedianExpLimit.at(i);  
        yu95[i] = vExpLim95Up.at(i) / vMedianExpLimit.at(i) - 1;  
        yd95[i] = 1 - vExpLim95Down.at(i) / vMedianExpLimit.at(i);  
      }
      if (Ratio>=2) {
        yu68[i]=1; yd68[i]=1;
        yu95[i]=2; yd95[i]=2;
      }

    }
    ExpBand95 = new TGraphAsymmErrors((signed) vMass.size(),x,y,ex,ex,yd95,yu95);
    ExpBand95->SetFillColor(90); 
    ExpBand95->GetYaxis()->SetRangeUser(0.,50);
    ExpBand95->GetXaxis()->SetRangeUser(x1,x2);
    ExpBand95->GetXaxis()->SetTitle("Higgs mass [GeV/c^{2}]");
    ExpBand95->GetYaxis()->SetTitle("95% CL limit on #sigma/#sigma_{SM}");
    if (Ratio==1) ExpBand95->GetYaxis()->SetTitle("ratio to expected ");
    if (Ratio==2) ExpBand95->GetYaxis()->SetTitle("#Delta(observed, expected) / 1#sigma");
    ExpBand95->Draw("A3");

    ExpBand95->GetXaxis()->SetTitleOffset(1.3);
    ExpBand95->GetYaxis()->SetTitleOffset(1.7);
    ExpBand95->GetXaxis()->SetNdivisions(510);
    ExpBand95->GetYaxis()->SetNdivisions(505);

    ExpBand95->GetYaxis()->SetRangeUser(0.,50);
    ExpBand95->Draw("A3");
    ExpBand68 = new TGraphAsymmErrors((signed) vMass.size(),x,y,ex,ex,yd68,yu68);
    ExpBand68->SetFillColor(211); 
    ExpBand68->Draw("3");

    ExpLim = new TGraph((signed) vMass.size(),x,y);    
    ExpLim->SetLineWidth(2);
    ExpLim->SetLineStyle(2);
    ExpLim->Draw("l");
  }

  // Observed Limit
  TGraph* ObsLim = NULL ;
  TGraph *ObsLimComp[4];
  TGraph *ExpLimComp[4];
  float xcomp[4][100];
  float ycomp[4][100];    
  float ycompe[4][100];    
  if ( DoObsLim ) {
    float x[100];
    float y[100];    
    for ( int i = 0 ; i < (signed) vMass.size() ; ++i ) { 
      x[i] = vMass.at(i) ; 
      y[i] = vObsLimit.at(i);
      if (Ratio>=1) y[i] = vObsLimit.at(i) / vMedianExpLimit.at(i); 
      if (Ratio>=2) {
        float onesigma = 
        ( fabs(log(vExpLim68Down.at(i)) - log(vMedianExpLimit.at(i)))
        + fabs(log(vExpLim68Up.at(i))   - log(vMedianExpLimit.at(i)))
        + fabs(log(vExpLim95Down.at(i)) - log(vMedianExpLimit.at(i)))
        + fabs(log(vExpLim95Up.at(i))   - log(vMedianExpLimit.at(i))) ) / 6;
        y[i] = (log(vObsLimit.at(i))-log(vMedianExpLimit.at(i))) / onesigma;
      }
      if(y[i] > max) max = y[i]; if(y[i] < min) min = y[i];
    }
    int colobs = kBlack; if (Ratio) colobs = kBlue+1;
    ObsLim = new TGraph((signed) vMass.size(),x,y);
    ObsLim->SetMarkerColor(colobs);
    ObsLim->SetLineWidth(2);
    ObsLim->SetLineColor(colobs);
    //ObsLim->SetLineStyle(2);
    ObsLim->SetMarkerStyle(kFullDotLarge);
    ObsLim->SetMarkerSize(1.0);
    if   (DoExpLim) ObsLim->Draw("lp");
    else {
      ObsLim->GetYaxis()->SetRangeUser(0.,10); 
      ObsLim->GetXaxis()->SetRangeUser(x1,x2);
      ObsLim->GetXaxis()->SetTitle("Higgs mass [GeV/c^{2}]");
      ObsLim->GetYaxis()->SetTitle("95% CL Limit on #sigma/#sigma_{SM}");
      ObsLim->Draw("alp");
    }

    // Add comparison of observed limits by other methods
    if (LimitFile.size()>1) { 

      int cols[3] = { kAzure, kMagenta+2, kGreen+3 };
      if (Ratio) { cols[0] = kRed+1; cols[1] = kGreen+3; cols[2] = kMagenta+2; }
      int maks[3] = { 24, 25, 26 };
      for ( int ifile = 1 ; ifile < (signed) LimitFile.size() ; ++ifile ) {

        ifstream indatacomp;
        indatacomp.open(LimitFile[ifile].c_str());
        if(!indatacomp) { // file couldn't be opened
          cerr << "Error: file could not be opened" << endl;
          return;
        }
        int npt(0);
        while ( indatacomp >> Mass >> ObsLimit >> MeanExpLimit >> MedianExpLimit >> ExpLim95Down >> ExpLim68Down >> ExpLim68Up >> ExpLim95Up )  {
          cout << Mass << " " << ObsLimit  << " " << MedianExpLimit <<" "<< ExpLim68Down <<" "<< ExpLim68Up <<" "<< ExpLim95Down <<" "<< ExpLim95Up << endl;

          xcomp[ifile][npt] = Mass           ;
          ycomp[ifile][npt] = ObsLimit;
          ycompe[ifile][npt] = MedianExpLimit;
          if (Ratio==1) {
            ycomp[ifile][npt] = ObsLimit / MedianExpLimit;
            ycompe[ifile][npt] = MedianExpLimit / MedianExpLimit;            
          }
          if (Ratio==2) {
            float onesigma = 
            ( fabs(log(ExpLim68Down) - log(MedianExpLimit))
            + fabs(log(ExpLim68Up)   - log(MedianExpLimit))
            + fabs(log(ExpLim95Down) - log(MedianExpLimit))
            + fabs(log(ExpLim95Up)   - log(MedianExpLimit)) ) / 6.;
            ycomp[ifile][npt] = (log(ObsLimit) - log(MedianExpLimit)) / onesigma; 
            ycompe[ifile][npt] = MedianExpLimit - MedianExpLimit       ;                 
          }
          ++npt;
        }
        if (npt != (signed)vMass.size()) 
          printf("*** DIFFERENT NUMBER OF MASS POINTS!! %d %d*** \n", (Int_t)vMass.size(), npt);
        ObsLimComp[ifile-1] = new TGraph(npt,xcomp[ifile],ycomp[ifile]);
        ObsLimComp[ifile-1]->SetMarkerColor(cols[ifile-1]);
        ObsLimComp[ifile-1]->SetMarkerStyle(maks[ifile-1]);
        ObsLimComp[ifile-1]->SetLineWidth(2);
        ObsLimComp[ifile-1]->SetLineStyle(1);
        ObsLimComp[ifile-1]->SetLineColor(cols[ifile-1]);
        ObsLimComp[ifile-1]->Draw("lp");
        ExpLimComp[ifile-1] = new TGraph(npt,xcomp[ifile],ycompe[ifile]);
        ExpLimComp[ifile-1]->SetLineWidth(2);
        ExpLimComp[ifile-1]->SetLineStyle(ifile+2);
        ExpLimComp[ifile-1]->SetLineColor(cols[ifile-1]);
        //ExpLimComp[ifile-1]->Draw("l");

        printf("made graph for %d\n", ifile-1);

      }
    }    
  }

  TText* CMS = new TText(.92,.87,"CMS preliminary");
  CMS ->SetTextAlign(32);
  CMS ->SetTextFont(42);
  CMS ->SetTextSize(.035);
  CMS ->SetNDC(1);
  CMS ->Draw("same");

  TLatex* title = new TLatex(.92,.82,LimTitle.c_str());
  title->SetTextAlign(32);
  title->SetTextFont(42);
  title->SetTextSize(.035);
  title->SetNDC(1);
  title->Draw("same");

  TLatex* Lumi = new TLatex(.92,.77, TString("L = "+lumi).Data());
  Lumi ->SetTextAlign(32);
  Lumi ->SetTextFont(42);
  Lumi ->SetTextSize(.035);
  Lumi ->SetNDC(1);
  Lumi ->Draw("same");

  TLegend* leg = NULL ;

  float legymin = 0.72 - 0.04*(LimitFile.size()-1);
  leg = new TLegend(0.23, legymin, 0.60, 0.87, "");
  leg->SetTextFont(42);
  leg->SetTextSize(.035);

  ExpBand68->SetLineColor(10);
  ExpBand95->SetLineColor(10);

  TString algo(""); if (LimitFile.size()>1) algo = " ("+LimitName[0]+")";
  if (DoExpLim) leg->AddEntry(ExpLim,   " median expected"+algo,"l");
  if (DoExpLim) leg->AddEntry(ExpBand68," expected #pm 1#sigma","f");
  if (DoExpLim) leg->AddEntry(ExpBand95," expected #pm 2#sigma","f");
  if (DoObsLim) leg->AddEntry(ObsLim," observed"+algo,"lp");
  for (int i=1; i<signed(LimitFile.size()); ++i) 
    leg->AddEntry(ObsLimComp[i-1],TString(" observed ("+LimitName[i]+")").Data(),"lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  leg->Draw("same");

  float maxMass = 300;
  max=10;
  // hard-code max for different channels
  TString name(filePrefix);
  if      (name.Contains("comb")) {
    max = 10;
  } 
  else if (name.Contains("wh3l")) {
    max = 20; maxMass = 200; x1 = 120;
  }
  else if (name.Contains("0j")) {
    if      (name.Contains("of")) max = 14;
    else if (name.Contains("sf")) max = 25;
    else                          max = 12;
  } else if (name.Contains("1j")) {
    if      (name.Contains("of")) max = 30;
    else if (name.Contains("sf")) max = 65;
    else                          max = 30;
  } else if (name.Contains("2j")) {
    max = 32;
  } else {
    printf("*** CHANNEL NOT FOUND ***\n");
  }   

  if (Ratio)    { max = 5; min = 0; }
  if (Ratio>=2) { min = -1; }

  vector<string> extensions;
  extensions.push_back(".png");
  extensions.push_back(".pdf");
  extensions.push_back(".eps");

  ExpBand95->GetXaxis()->SetRangeUser(x1,x2);
  //ExpBand95->GetYaxis()->SetRangeUser(min-0.2,max);
  ExpBand95->GetYaxis()->SetRangeUser(0.1,max);

  cLimit->GetFrame()->DrawClone();
  cLimit->Update();

  TLine *l = new TLine(cLimit->GetUxmin(),1,cLimit->GetUxmax(),1);
  l->SetLineWidth(2);
  l->SetLineColor(kRed+1);
  if (!Ratio) l->Draw("same");

  for(size_t i=0;i<extensions.size();++i) cLimit->Print( ("plots/"+filePrefix+"_lin"+extensions[i]).c_str() );

  ExpBand95->GetXaxis()->SetRangeUser(x1,maxMass);
  //ExpBand95->GetYaxis()->SetRangeUser(min-0.2,max);
  ExpBand95->GetYaxis()->SetRangeUser(0.1,max);
  cLimit->Modified();
  cLimit->Update();
  l->SetX1(cLimit->GetUxmin());
  l->SetX2(cLimit->GetUxmax());

  for(size_t i=0;i<extensions.size();++i) cLimit->Print( ("plots/"+filePrefix+"_zoom_lin"+extensions[i]).c_str() );

  if (!Ratio) {

    cLimit->SetLogy();

    ExpBand95->GetXaxis()->SetRangeUser(x1,x2);
    //ExpBand95->GetYaxis()->SetRangeUser(min,max*10);
    ExpBand95->GetYaxis()->SetRangeUser(0.1,max*10);
    cLimit->Modified();
    cLimit->Update();
    l->SetX1(cLimit->GetUxmin());
    l->SetX2(cLimit->GetUxmax());
    for(size_t i=0;i<extensions.size();++i) cLimit->Print( ("plots/"+filePrefix+"_log"+extensions[i]).c_str() );

    ExpBand95->GetXaxis()->SetRangeUser(x1,maxMass);
    //ExpBand95->GetYaxis()->SetRangeUser(min,max*10);
    ExpBand95->GetYaxis()->SetRangeUser(0.1,max*10);
    cLimit->Modified();
    cLimit->Update();
    l->SetX1(cLimit->GetUxmin());
    l->SetX2(cLimit->GetUxmax());
    for(size_t i=0;i<extensions.size();++i) cLimit->Print( ("plots/"+filePrefix+"_zoom_log"+extensions[i]).c_str() );

  }
  
  // Make latex table with nominal result + one column of observed from 2nd method                                                                                                                                                         
  if (LimitFile.size()>1) {

    printf("Writing table to %s\n", LimTable.c_str());
    ofstream fout; fout.open(LimTable.c_str());
    fout<<"\\begin{table}[h!]\\begin{center}"<<endl;
    fout<<Form("\\caption{Expected and observed %s upper limits and uncertainty band %s.}",
      TabLegend.c_str(), TabChannel.c_str())<<endl;
    fout<<Form(" \\label{table:%s}", TabLabel.c_str())<<endl;
    fout<<" \\small{"<<endl;
    fout<<"    \\begin{tabular}{c|c c c c c}"<<endl;
    fout<<"    \\hline"<<endl;
    fout<<"    $m_{\\mathrm{H}}$~[\\GeV] &"<<endl;
    fout<<"    observed &"<<endl;
    fout<<"    observed &"<<endl;
    fout<<"    median expected &"<<endl;
    fout<<"    expected range &"<<endl;
    fout<<"    expected range \\\\ "<<endl;
    fout<<"    &(CLs)&"<<endl;
    fout<<"    (Bayesian)&"<<endl;
    fout<<"    (CLs) &"<<endl;
    fout<<"    for 68\\% &"<<endl;
    fout<<"    for 95\\% \\\\ "<<endl;
    fout<<"    \\hline"<<endl;
    for (int im=0; im<(signed)vMass.size(); ++im) {
      fout<<Form("%.0f \t&\t %6.1f \t&\t %6.1f \t&\t %6.1f \t&\t[ %6.1f , %6.1f ]\t&\t[ %6.1f , %6.1f ] \\\\ \n",
        vMass[im], vObsLimit[im], ycomp[1][im], vMedianExpLimit[im],
        vExpLim68Down[im], vExpLim68Up[im], vExpLim95Down[im], vExpLim95Up[im]);
    }
    fout<<"    \\hline"<<endl;
    fout<<"    \\end{tabular}"<<endl;
    fout<<"  }"<<endl;
    fout<<"  \\end{center}"<<endl;
    fout<<"\\end{table}"<<endl;
    fout.close();

  }

  return;
}

