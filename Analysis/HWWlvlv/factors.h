#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"

void atributes(TH1D *histo, Char_t xtitle[]="", Int_t COLOR = 1, Char_t ytitle[]="Fraction");
double scpFast(double sig, double bkg, double sigma_b, double delta_b = 0.0);
double scaleFactor(double pt1, double eta1, double pt2, double eta2, int type, int nsel);
double enhancementFactor(double mass, bool isBRFactor = false);
double fakeRate(double pt, double eta, TH2D *fhDFRMu, TH2D *fhDFREl, int fm, int fe);
double leptonEfficiency(double pt, double eta, TH2D *fhDEffMu, TH2D *fhDEffEl, int lid, int syst = 0);
double hzz2l_cuts(double mass, int opt);

void atributes(TH1D *histo, Char_t xtitle[], Int_t COLOR, Char_t ytitle[]){
  //InitHist(histo,xtitle,ytitle);

  histo->ResetAttLine();
  histo->ResetAttFill();
  histo->ResetAttMarker();
  histo->GetYaxis()->SetNdivisions(505);
  histo->GetXaxis()->SetNdivisions(505);
  histo->SetTitle("");
  //histo->SetMarkerColor(COLOR);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.8);
  histo->SetLineWidth(4);
  histo->SetLineColor(COLOR);
  histo->SetMarkerStyle(kFullDotLarge);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->SetTitleOffset(1.);
  //histo->GetXaxis()->SetTitleSize(0.036);
  //double binw = histo->GetBinWidth(1);
  //char yTitle[300];
  //sprintf(yTitle,"Number of Entries / %3.1f MeV",binw*1000.0);
  //sprintf(yTitle,"Number of Entries / %1.2f cm",binw);
  //sprintf(yTitle,"Number of Entries");
  //sprintf(yTitle,"Ratio");
  //sprintf(yTitle,"(P-N)/(P+N)");
  //sprintf(yTitle,"Fraction");
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->GetYaxis()->SetTitle(ytitle);
  //histo->GetYaxis()->SetTitleSize(0.05);
  //histo->GetYaxis()->SetLabelSize(0.02);
  histo->GetYaxis()->CenterTitle(kTRUE);
}

double scpFast(double sig, double bkg, double sigma_b, double delta_b)
{
double fac2  = sqrt(bkg+delta_b)/sqrt(bkg+sigma_b*sigma_b+delta_b);
double Sc12_sys = 2*(sqrt(sig+bkg)-sqrt(bkg+delta_b))*fac2;
return Sc12_sys;
//double Sc12_sys = sig/sqrt(bkg+sigma_b*sigma_b+0.0*delta_b);
//return Sc12_sys;
}

double scaleFactor(double pt1, double eta1, double pt2, double eta2, int type, int nsel){
  // type == 10/11/12/13 ->mm/ee/em/me
  // hardcoded, it's not used much
  int syst = 0;
  double scaleE[2] = {0.969, 0.992};
  if(syst == 1){
    scaleE[0] = scaleE[0] - 0.019;
    scaleE[1] = scaleE[1] - 0.026;
  };
  double scaleM[2] = {1.000, 1.000};
  if(syst == 1){
    scaleM[0] = scaleM[0] + 0.002;
    scaleM[1] = scaleM[1] + 0.002;
  };
  double weight = 1.0;
  
  if(nsel == 0){ // electron scale factor
    if     (type == 11){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
      if(fabs(eta2) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
    }
    else if(type == 12){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
    }
    else if(type == 13){ 
      if(fabs(eta2) < 1.479) weight = weight * scaleE[0];
      else                   weight = weight * scaleE[1];
    }
  }
  else if(nsel == 1){ // muon scale factor
    if     (type == 10){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
      if(fabs(eta2) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
    }
    else if(type == 12){ 
      if(fabs(eta2) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
    }
    else if(type == 13){ 
      if(fabs(eta1) < 1.479) weight = weight * scaleM[0];
      else                   weight = weight * scaleM[1];
    }
  }
  else if(nsel == 2){ // luminosity scale factor
    if     (type == 0) weight = weight * 0.00000;
    else if(type == 1) weight = weight * 0.23111;
    else if(type == 2) weight = weight * 0.82486;
    else if(type == 3) weight = weight * 1.47045;
    else if(type == 4) weight = weight * 1.83450;
    else if(type == 5) weight = weight * 1.79663;
    else if(type == 6) weight = weight * 1.46496;
    else if(type == 7) weight = weight * 1.05682;
    else if(type == 8) weight = weight * 0.70823;
    else if(type == 9) weight = weight * 0.47386;
    else if(type ==10) weight = weight * 0.32382;
    else if(type ==11) weight = weight * 0.22383;
    else if(type ==12) weight = weight * 0.17413;
    else if(type ==13) weight = weight * 0.10930;
    else if(type ==14) weight = weight * 0.09563;
    else if(type ==15) weight = weight * 0.08367;
    else if(type ==16) weight = weight * 0.05418;
    else if(type ==17) weight = weight * 0.04891	;
    else if(type ==18) weight = weight * 0.03515	;
    else if(type >=19) weight = weight * 0.01000;
  }
  else if(nsel == 3){ // momentum scale factor
    if     (type == 10){ 
      if(pt1*0.99 < 20) weight = 0.0;
      if(pt2*0.99 < 20) weight = 0.0;
    }
    else if(type == 11){ 
      if     (fabs(eta1) <  1.479 && pt1*0.98 < 20) weight = 0.0;
      else if(fabs(eta1) >= 1.479 && pt1*0.96 < 20) weight = 0.0;
      if     (fabs(eta2) <  1.479 && pt2*0.98 < 20) weight = 0.0;
      else if(fabs(eta2) >= 1.479 && pt2*0.96 < 20) weight = 0.0;
    }
    else if(type == 12){ 
      if     (fabs(eta1) <  1.479 && pt1*0.98 < 20) weight = 0.0;
      else if(fabs(eta1) >= 1.479 && pt1*0.96 < 20) weight = 0.0;
      if(pt2*0.99 < 20) weight = 0.0;
    }
    else if(type == 13){ 
      if(pt1*0.99 < 20) weight = 0.0;
      if     (fabs(eta2) <  1.479 && pt2*0.98 < 20) weight = 0.0;
      else if(fabs(eta2) >= 1.479 && pt2*0.96 < 20) weight = 0.0;
    }
  }
  
  return weight;
}

double enhancementFactor(double mass, bool isBRFactor){

if(isBRFactor == false){
  if     (mass==110.000) return 9.203 ;
  else if(mass==115.000) return 9.170 ;
  else if(mass==120.000) return 9.129 ;
  else if(mass==130.000) return 9.053 ;
  else if(mass==140.000) return 8.966 ;
  else if(mass==150.000) return 8.867 ;
  else if(mass==160.000) return 8.768 ;
  else if(mass==170.000) return 8.675 ;
  else if(mass==180.000) return 8.562 ;
  else if(mass==190.000) return 8.454 ;
  else if(mass==200.000) return 8.337 ;
  else if(mass==210.000) return 8.218 ;
  else if(mass==220.000) return 8.096 ;
  else if(mass==230.000) return 7.962 ;
  else if(mass==250.000) return 7.685 ;
  else if(mass==300.000) return 6.824 ;
  else if(mass==350.000) return 5.304 ;
  else if(mass==400.000) return 4.498 ;
  else if(mass==450.000) return 4.3895;
  else if(mass==500.000) return 4.450 ;
  else if(mass==550.000) return 4.6065;
  else if(mass==600.000) return 4.837 ;
}
else {
  if     (mass==110.000) return 0.574;
  else if(mass==115.000) return 0.569;
  else if(mass==120.000) return 0.570;
  else if(mass==130.000) return 0.600;
  else if(mass==140.000) return 0.672;
  else if(mass==150.000) return 0.781;
  else if(mass==160.000) return 0.931;
  else if(mass==170.000) return 0.981;
  else if(mass==180.000) return 0.986;
  else if(mass==190.000) return 0.988;
  else if(mass==200.000) return 0.988;
  else if(mass==210.000) return 0.989;
  else if(mass==220.000) return 0.990;
  else if(mass==230.000) return 0.991;
  else if(mass==250.000) return 0.991;
  else if(mass==300.000) return 0.993;
  else if(mass==350.000) return 0.992;
  else if(mass==400.000) return 0.989;
  else if(mass==450.000) return 0.989;
  else if(mass==500.000) return 0.991;
  else if(mass==550.000) return 0.993;
  else if(mass==600.000) return 0.996;
}

return 1.0;

}

double fakeRate(double pt, double eta, TH2D *fhDFRMu, TH2D *fhDFREl, int fm, int fe){
  // fm == apply muon fake rate, fe == apply electron fake rate
  if(fm == 0 && fe == 0) return 1.0;
  double mypt   = TMath::Min(pt,34.999);
  double myeta  = TMath::Min(fabs(eta),2.4999);
  double prob = 1.0;
  if     (fm == 1){
    Int_t ptbin = fhDFRMu->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDFRMu->GetYaxis()->FindBin(myeta);        
    prob = fhDFRMu->GetBinContent(ptbin,etabin);
    //Int_t ptbin = fhDFRMu->GetYaxis()->FindBin(mypt);
    //Int_t etabin = fhDFRMu->GetXaxis()->FindBin(myeta);  
    //prob = fhDFRMu->GetBinContent(etabin,ptbin);
  }
  else if(fe == 1){
    Int_t ptbin = fhDFREl->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDFREl->GetYaxis()->FindBin(myeta);        
    prob = fhDFREl->GetBinContent(ptbin,etabin);
    //Int_t ptbin = fhDFREl->GetYaxis()->FindBin(mypt);
    //Int_t etabin = fhDFREl->GetXaxis()->FindBin(myeta);  
    //prob = fhDFREl->GetBinContent(etabin,ptbin);
  }
  return prob/(1-prob);
}

double leptonEfficiency(double pt, double eta, TH2D *fhDEffMu, TH2D *fhDEffEl, int lid, int syst){
  // lid == 13 (muon), 11 (electron)
  double mypt   = TMath::Min(pt,49.999);
  double myeta  = TMath::Min(fabs(eta),2.4999);
  double prob = 1.0;
  if     (TMath::Abs(lid) == 13){
    Int_t ptbin = fhDEffMu->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDEffMu->GetYaxis()->FindBin(myeta);	 
    prob = fhDEffMu->GetBinContent(ptbin,etabin);
    if     (syst > 0) prob += fhDEffMu->GetBinError(ptbin,etabin);
    else if(syst < 0) prob -= fhDEffMu->GetBinError(ptbin,etabin);
  }
  else if(TMath::Abs(lid) == 11){
    Int_t ptbin = fhDEffEl->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDEffEl->GetYaxis()->FindBin(myeta);	 
    prob = fhDEffEl->GetBinContent(ptbin,etabin);
    if     (syst > 0) prob += fhDEffEl->GetBinError(ptbin,etabin);
    else if(syst < 0) prob -= fhDEffEl->GetBinError(ptbin,etabin);
  }
  return prob;
}

double hzz2l_cuts(double mass, int opt){    
double x        = mass/50. - 2.0;
double cutValue = 0.0;

if     (opt == 0){ // min(met)
 cutValue = 10.933 * x + 41.16;
}
else if(opt == 1){ // max(met)
 cutValue = 79.952 * x + 234.67;
}
else if(opt == 2){ // min(mthiggs)
 cutValue = -1.6174 * x*x + 46.683 * x + 105.42;
}
else if(opt == 3){ // max(mthiggs)
 if(mass <= 450){
   cutValue =  5.2857 * x*x + 11.071 * x + 184.43;
 }
 else {
   cutValue = 7000.0;
 }
}
return cutValue;
}
