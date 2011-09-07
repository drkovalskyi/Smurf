#include <iostream>
#include <fstream>
#include "TString.h"

void shapeStat(int mH = 140, int njet = 0, const char *flavor = "of", TString inputdir = "ana_v6_1500pb_LP_FINAL/")
{ 
  inputdir = inputdir + Form("%i/",mH);
  cout << "Accessing directory " << inputdir  << endl;
  TString shapecard = inputdir + Form("hww%s_%ij_shape.txt", flavor, njet);
  TString newcardname = inputdir + Form("shapeStat/hww%s_%ij_shape.txt", flavor, njet);
  TString shapeHistName = inputdir + Form("/shapeStat/hww%s_%ij.input.root", flavor, njet);
  
  writeNewHist(njet, flavor, shapeHistName);
  writeNewCard(mH, shapecard, newcardname);
  
}

void writeNewCard(int mH, TString shapecard, TString newcardname) 
{
  ifstream inputfile(shapecard);
  ofstream newcard;
  
  std::cout << "Writing " << newcardname << std::endl;
  newcard.open(newcardname);
  string line;
  if (inputfile.is_open())
    {
      while ( inputfile.good() )
	{
	  getline (inputfile,line);
	  if (TString(line).Contains("histo_$PROCESS", TString::kExact)  && TString(line).Contains("shapes", TString::kExact)) 
	    newcard << TString(line).ReplaceAll(Form("/%i/", mH), Form("/%i/shapeStat/", mH)) << " histo_$PROCESS_$SYSTEMATIC" << endl;
	  else if (TString(line).Contains("histo_Data", TString::kExact)  && TString(line).Contains("shapes", TString::kExact)) 
	    newcard << TString(line).ReplaceAll(Form("/%i%/", mH), Form("/%i%/shapeStat/", mH)) << endl;
	  else if (TString(line).Contains("stat", TString::kExact)) {
	    newcard << TString(line).ReplaceAll("lnN", "shapeStat") << endl;
	  }
	  else 
	    newcard << line << endl;
	}
      inputfile.close();
      newcard.close();
    }
  else 
    cout << "Unable to open " << shapecard << endl;
}

void writeNewHist(const int njet, const char *flavor, TString shapeHistName) 
{
  TFile *file = TFile::Open(shapeHistName, "UPDATE");
  if ( file == 0x0 ) {
    std::cout << "Error: " << shapeHistName << " doesn't exist! exitting the script.\n";
    return;
  }
  std::cout << "Writing histogram " << shapeHistName << "\n";
  
  // ggH 
  TH1F *histo_ggH_Up = (TH1F*)file->Get("histo_ggH") ->Clone(Form("histo_ggH_CMS_hww%s_stat_%ij_ggH_bin1Up", flavor, njet));
  if (histo_ggH_Up->GetEntries() > 0 )  modifyHist( histo_ggH_Up, "Up");
  TH1F *histo_ggH_Down = (TH1F*)file->Get("histo_ggH") ->Clone(Form("histo_ggH_CMS_hww%s_stat_%ij_ggH_bin1Down", flavor, njet));
  if (histo_ggH_Down->GetEntries() > 0 )  modifyHist( histo_ggH_Down, "Down");
  histo_ggH_Up->Write();
  histo_ggH_Down->Write();

  // qqH 
  TH1F *histo_qqH_Up = (TH1F*)file->Get("histo_qqH") ->Clone(Form("histo_qqH_CMS_hww%s_stat_%ij_qqH_bin1Up", flavor, njet));
  if (histo_qqH_Up->GetEntries() > 0 )  modifyHist( histo_qqH_Up, "Up");
  TH1F *histo_qqH_Down = (TH1F*)file->Get("histo_qqH") ->Clone(Form("histo_qqH_CMS_hww%s_stat_%ij_qqH_bin1Down", flavor, njet));
  if (histo_qqH_Down->GetEntries() > 0 )  modifyHist( histo_qqH_Down, "Down");
  histo_qqH_Up->Write();
  histo_qqH_Down->Write();

  // WH 
  TH1F *histo_WH_Up = (TH1F*)file->Get("histo_WH") ->Clone(Form("histo_WH_CMS_hww%s_stat_%ij_WH_bin1Up", flavor, njet));
  if (histo_WH_Up->GetEntries() > 0 )  modifyHist( histo_WH_Up, "Up");
  TH1F *histo_WH_Down = (TH1F*)file->Get("histo_WH") ->Clone(Form("histo_WH_CMS_hww%s_stat_%ij_WH_bin1Down", flavor, njet));
  if (histo_WH_Down->GetEntries() > 0 )  modifyHist( histo_WH_Down, "Down");
  histo_WH_Up->Write();
  histo_WH_Down->Write();
  
  // ZH 
  TH1F *histo_ZH_Up = (TH1F*)file->Get("histo_ZH") ->Clone(Form("histo_ZH_CMS_hww%s_stat_%ij_ZH_bin1Up", flavor, njet));
  if (histo_ZH_Up->GetEntries() > 0 )  modifyHist( histo_ZH_Up, "Up");
  TH1F *histo_ZH_Down = (TH1F*)file->Get("histo_ZH") ->Clone(Form("histo_ZH_CMS_hww%s_stat_%ij_ZH_bin1Down", flavor, njet));
  if (histo_ZH_Down->GetEntries() > 0 )  modifyHist( histo_ZH_Down, "Down");
  histo_ZH_Up->Write();
  histo_ZH_Down->Write();
  
  // qqWW 
  TH1F *histo_qqWW_Up = (TH1F*)file->Get("histo_qqWW") ->Clone(Form("histo_qqWW_CMS_hww%s_stat_%ij_WW_bin1Up", flavor, njet));
  if (histo_qqWW_Up->GetEntries() > 0 )  modifyHist( histo_qqWW_Up, "Up");
  TH1F *histo_qqWW_Down = (TH1F*)file->Get("histo_qqWW") ->Clone(Form("histo_qqWW_CMS_hww%s_stat_%ij_WW_bin1Down", flavor, njet));
  if (histo_qqWW_Down->GetEntries() > 0 )  modifyHist( histo_qqWW_Down, "Down");
  histo_qqWW_Up->Write();
  histo_qqWW_Down->Write();
  
  // ggWW 
  TH1F *histo_ggWW_Up = (TH1F*)file->Get("histo_ggWW") ->Clone(Form("histo_ggWW_CMS_hww%s_stat_%ij_ggWW_bin1Up", flavor, njet));
  if (histo_ggWW_Up->GetEntries() > 0 )  modifyHist( histo_ggWW_Up, "Up");
  TH1F *histo_ggWW_Down = (TH1F*)file->Get("histo_ggWW") ->Clone(Form("histo_ggWW_CMS_hww%s_stat_%ij_ggWW_bin1Down", flavor, njet));
  if (histo_ggWW_Down->GetEntries() > 0 )  modifyHist( histo_ggWW_Down, "Down");
  histo_ggWW_Up->Write();
  histo_ggWW_Down->Write();

  // VV 
  TH1F *histo_VV_Up = (TH1F*)file->Get("histo_VV") ->Clone(Form("histo_VV_CMS_hww%s_stat_%ij_VV_bin1Up", flavor, njet));
  if (histo_VV_Up->GetEntries() > 0 )  modifyHist( histo_VV_Up, "Up");
  TH1F *histo_VV_Down = (TH1F*)file->Get("histo_VV") ->Clone(Form("histo_VV_CMS_hww%s_stat_%ij_VV_bin1Down", flavor, njet));
  if (histo_VV_Down->GetEntries() > 0 )  modifyHist( histo_VV_Down, "Down");
  histo_VV_Up->Write();
  histo_VV_Down->Write();

  // Top 
  TH1F *histo_Top_Up = (TH1F*)file->Get("histo_Top") ->Clone(Form("histo_Top_CMS_hww%s_stat_%ij_ttbar_bin1Up", flavor, njet));
  if (histo_Top_Up->GetEntries() > 0 )  modifyHist( histo_Top_Up, "Up");
  TH1F *histo_Top_Down = (TH1F*)file->Get("histo_Top") ->Clone(Form("histo_Top_CMS_hww%s_stat_%ij_ttbar_bin1Down", flavor, njet));
  if (histo_Top_Down->GetEntries() > 0 )  modifyHist( histo_Top_Down, "Down");
  histo_Top_Up->Write();
  histo_Top_Down->Write();

  // Zjets 
  TH1F *histo_Zjets_Up = (TH1F*)file->Get("histo_Zjets") ->Clone(Form("histo_Zjets_CMS_hww%s_stat_%ij_Z_bin1Up", flavor, njet));
  if (histo_Zjets_Up->GetEntries() > 0 )  modifyHist( histo_Zjets_Up, "Up");
  TH1F *histo_Zjets_Down = (TH1F*)file->Get("histo_Zjets") ->Clone(Form("histo_Zjets_CMS_hww%s_stat_%ij_Z_bin1Down", flavor, njet));
  if (histo_Zjets_Down->GetEntries() > 0 )  modifyHist( histo_Zjets_Down, "Down");
  histo_Zjets_Up->Write();
  histo_Zjets_Down->Write();

  // Wjets 
  TH1F *histo_Wjets_Up = (TH1F*)file->Get("histo_Wjets") ->Clone(Form("histo_Wjets_CMS_hww%s_stat_%ij_Wjets_bin1Up", flavor, njet));
  if (histo_Wjets_Up->GetEntries() > 0 )  modifyHist( histo_Wjets_Up, "Up");
  TH1F *histo_Wjets_Down = (TH1F*)file->Get("histo_Wjets") ->Clone(Form("histo_Wjets_CMS_hww%s_stat_%ij_Wjets_bin1Down", flavor, njet));
  if (histo_Wjets_Down->GetEntries() > 0 )  modifyHist( histo_Wjets_Down, "Down");
  histo_Wjets_Up->Write();
  histo_Wjets_Down->Write();

  // Wgamma 
  TH1F *histo_Wgamma_Up = (TH1F*)file->Get("histo_Wgamma") ->Clone(Form("histo_Wgamma_CMS_hww%s_stat_%ij_Wgamma_bin1Up", flavor, njet));
  if (histo_Wgamma_Up->GetEntries() > 0 )  modifyHist( histo_Wgamma_Up, "Up");
  TH1F *histo_Wgamma_Down = (TH1F*)file->Get("histo_Wgamma") ->Clone(Form("histo_Wgamma_CMS_hww%s_stat_%ij_Wgamma_bin1Down", flavor, njet));
  if (histo_Wgamma_Down->GetEntries() > 0 )  modifyHist( histo_Wgamma_Down, "Down");
  histo_Wgamma_Up->Write();
  histo_Wgamma_Down->Write();
  
  // Ztt 
  TH1F *histo_Ztt_Up = (TH1F*)file->Get("histo_Ztt") ->Clone(Form("histo_Ztt_CMS_hww%s_stat_%ij_Ztt_bin1Up", flavor, njet));
  if (histo_Ztt_Up->GetEntries() > 0 )  modifyHist( histo_Ztt_Up, "Up");
  TH1F *histo_Ztt_Down = (TH1F*)file->Get("histo_Ztt") ->Clone(Form("histo_Ztt_CMS_hww%s_stat_%ij_Ztt_bin1Down", flavor, njet));
  if (histo_Ztt_Down->GetEntries() > 0 )  modifyHist( histo_Ztt_Down, "Down");
  histo_Ztt_Up->Write();
  histo_Ztt_Down->Write();

  file->Close();

}

void modifyHist( TH1F*& histo, TString option)
{
  if (histo == 0x0 || histo->GetEntries() == 0) return;
  for ( int i = 0 ; i < histo->GetNbinsX(); i++) {
    double bincontent = histo->GetBinContent(i);
    double binerr = histo->GetBinError(i);
    double newbincontent (bincontent), newbinerr(binerr);

    if ( bincontent <= 0.0) {
      float histweight = histo->Integral(0,1000)/float(histo->GetEntries());
      bincontent = 1.*histweight;
      binerr = 1.*histweight;
    }
    if ( option == "Up")  newbincontent = bincontent + binerr;
    if ( option == "Down") newbincontent = bincontent - binerr;
    
    if ( newbincontent <= 0.0) newbincontent = 0.001;
    histo->SetBinContent(i, newbincontent);
  }
}
