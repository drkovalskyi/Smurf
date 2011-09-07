#include <iostream>
#include <fstream>

void splitcard(int mH = 140, int njet = 0, const char *flavor = "of", TString inputdir = "EPS2011_StatTest/")
{ 
  inputdir = inputdir + Form("%i/",mH);
  
  cout << inputdir << endl;
  TString shapecard = inputdir + Form("hww%s_%ij_shape.txt", flavor, njet);
  TString shapehist = inputdir + Form("hww%s_%ij.input.root", flavor, njet);
  
  TFile *file = TFile::Open(shapehist);
  gROOT->cd();
  
  TH1F *histo_Data = (TH1F*)file->Get("histo_Data");
  
  TH1F *histo_ZH = (TH1F*)file->Get("histo_ZH");
  TH1F *histo_WH = (TH1F*)file->Get("histo_WH");
  TH1F *histo_qqH = (TH1F*)file->Get("histo_qqH");
  TH1F *histo_ggH = (TH1F*)file->Get("histo_ggH");
  TH1F *histo_qqWW = (TH1F*)file->Get("histo_qqWW");
  TH1F *histo_ggWW = (TH1F*)file->Get("histo_ggWW");
  TH1F *histo_Wjets = (TH1F*)file->Get("histo_Wjets");
  TH1F *histo_Wgamma = (TH1F*)file->Get("histo_Wgamma");
  TH1F *histo_Zjets = (TH1F*)file->Get("histo_Zjets");
  TH1F *histo_VV = (TH1F*)file->Get("histo_VV");
  TH1F *histo_Top = (TH1F*)file->Get("histo_Top");
  TH1F *histo_Ztt = (TH1F*)file->Get("histo_Ztt");

  histcheck(histo_Wjets);
  float nZH(0.0), nWH(0.0), nqqH(0.0), nggH(0.0), nqqWW(0.0), nggWW(0.0), nWjets(0.0), nWgamma(0.0), nZjets(0.0), nVV(0.0), nTop(0.0), nZtt(0.0);
  float sigmanZH(0.0), sigmanWH(0.0), sigmanqqH(0.0), sigmanggH(0.0), sigmanqqWW(0.0), sigmanggWW(0.0), sigmanWjets(0.0), sigmanWgamma(0.0), sigmanZjets(0.0), sigmanVV(0.0), sigmanTop(0.0), sigmanZtt(0.0);
  
  for (int i = 1; i < histo_Data->GetNbinsX()+1; i++ ) {
    
    getbin(histo_ZH, i, nZH, sigmanZH);
    getbin(histo_WH, i, nWH, sigmanWH);
    getbin(histo_qqH, i, nqqH, sigmanqqH);
    getbin(histo_ggH, i, nggH, sigmanggH);
    getbin(histo_qqWW, i, nqqWW, sigmanqqWW);
    getbin(histo_ggWW, i, nggWW, sigmanggWW);
    getbin(histo_Wjets, i, nWjets, sigmanWjets);
    getbin(histo_Wgamma, i, nWgamma, sigmanWgamma);
    getbin(histo_Zjets, i, nZjets, sigmanZjets);
    getbin(histo_VV, i, nVV, sigmanVV);
    getbin(histo_Top, i, nTop, sigmanTop);
    getbin(histo_Ztt, i, nZtt, sigmanZtt);
    
    
    ifstream inputfile(shapecard);
    ofstream newcard;
    TString newcardname = inputdir + Form("splitcards/hww%s_%ij_shape_Bin%i.txt", flavor, njet,i);
    cout << newcardname << endl;
    newcard.open(newcardname);
    string line;
    if (inputfile.is_open())
      {
	while ( inputfile.good() )
	  {
	    getline (inputfile,line);
	    if (TString(line).Contains("Observation", TString::kExact)) 
	      newcard << "Observation " << histo_Data->GetBinContent(i) << "\n";
	    else if (TString(line).Contains("shapes", TString::kExact)) 
	      newcard << "";
	    else if (TString(line).Contains("rate", TString::kExact))
	      newcard << Form("rate  %.3f   %.3f   %.3f   %.3f  %.3f  %.3f   %.3f   %.3f   %.3f  %.3f   %.3f   %.3f\n", 
			      nZH, nWH, nqqH, nggH, nqqWW, nggWW, nVV, nTop, nZjets, nWjets, nWgamma, nZtt);
	
	    else if (TString(line).Contains("stat", TString::kExact)) {
	      if (TString(line).Contains("ZH", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_ZH_bin%i     lnN %4.3f   -     -     -     -     -     -     -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanZH);
	      else if (TString(line).Contains("WH", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_WH_bin%i     lnN   -   %4.3f   -     -     -     -     -     -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanWH);
	      else if (TString(line).Contains("qqH_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_qqH_bin%i    lnN   -     -   %4.3f   -     -     -     -     -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanqqH);
	      else if (TString(line).Contains("ggH_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_ggH_bin%i    lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanggH);
	      else if (TString(line).Contains("WW_", TString::kExact) && !(TString(line).Contains("ggWW", TString::kExact)))
		newcard << Form("CMS_hww%s_stat_%ij_WW_bin%i     lnN   -     -     -     -   %4.3f   -     -     -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanqqWW);
	      else if (TString(line).Contains("ggWW_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_ggWW_bin%i   lnN   -     -     -     -     -   %4.3f   -     -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanggWW);
	      else if (TString(line).Contains("VV_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_VV_bin%i     lnN   -     -     -     -     -     -   %4.3f   -     -     -     -     -  \n", flavor, njet, i, 1.0+sigmanVV);
	      else if (TString(line).Contains("ttbar_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_ttbar_bin%i  lnN   -     -     -     -     -     -     -   %4.3f   -     -     -     -  \n", flavor, njet, i, 1.0+sigmanTop);
	      else if (TString(line).Contains("Z_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_Z_bin%i      lnN   -     -     -     -     -     -     -     -   %4.3f   -     -     -  \n", flavor, njet, i, 1.0+sigmanZjets);
	      else if (TString(line).Contains("Wjets_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_Wjets_bin%i  lnN   -     -     -     -     -     -     -     -     -   %4.3f   -     -  \n", flavor, njet, i, 1.0+sigmanWjets);
	      else if (TString(line).Contains("Wgamma_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_Wgamma_bin%i lnN   -     -     -     -     -     -     -     -     -     -   %4.3f   -  \n", flavor, njet, i, 1.0+sigmanWgamma);
	      else if (TString(line).Contains("Ztt_", TString::kExact))
		newcard << Form("CMS_hww%s_stat_%ij_Ztt_bin%i    lnN   -     -     -     -     -     -     -     -     -     -     -   %4.3f \n", flavor, njet, i, 1.0+sigmanZtt);
	      else 
		newcard << line << endl;
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
  
}

// simple sanity check that the histograms are filled correctly
void histcheck(TH1F*& hist) {
  
  if(hist == 0x0 ) return;
  
  float yield;
  float yielderr;
  
  for (int i = 1; i < hist->GetNbinsX() + 1; i++) {
  
    yield += hist->GetBinContent(i);
    yielderr += pow(hist->GetBinError(i),2);
  }

  //cout << sqrt(yielderr)/yield <<endl;

}

void getbin(TH1F* & hist, int& bin, float & bincontent, float & binerror)
{
  
  if(hist == 0x0) {
    cout << "Histogram is not found.... Bins and Bin errors are set to 0...";
    return;
  }
  
  bincontent = hist->GetBinContent(bin);
  
  // ignore the bin-statistical error....
  /*
  float tot (0.0), toterr (0.0);
  for (int i = 1; i < hist->GetNbinsX()+1; i++) {
    tot += hist->GetBinContent(i);
    toterr += pow(hist->GetBinError(i),2);
  }
  binerror = sqrt(toterr)/tot;
  */

  // realistic uncertainties
  
  if (bincontent <= 0) {
    // flucturate the original event count by 1...
    if ( hist->GetEntries() == 0.) {
      bincontent = 0.0;
      binerror = 0.0;
    }
    else {
      float histweight = hist->Integral(0,1000)/float(hist->GetEntries());
      bincontent = 1.*histweight;
      binerror = 1.0;
    }
  }
  else 
    binerror = hist->GetBinError(bin) / bincontent;
  
}
