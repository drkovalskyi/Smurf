#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double Unroll2VarTo1Var(double varA, double varB, int binsA, int binsB, bool verbose = false);
TH1F* Unroll2DTo1D(TH2F* h2, const char* hname);
TH1D * SmurfRebin(const TH1D *old, const unsigned int rebin);
double CalcGammaMRstar(LorentzVector ja, LorentzVector jb);
double CalcMTR(LorentzVector ja, LorentzVector jb, double metValue, double metPhi);
Double_t CalcMR(LorentzVector ja, LorentzVector jb);
double DeltaPhi(double phi1, double phi2);
void atributes(TH1D *histo, Char_t xtitle[]="", Int_t COLOR = 1, Char_t ytitle[]="Fraction");
double scpFast(double sig, double bkg, double sigma_b, double delta_b = 0.0);
double scaleFactor(double pt1, double eta1, double pt2, double eta2, int type, int nsel);
double enhancementFactor(double mass, int type = 0);
double fakeRate(double pt, double eta, TH2D *fhDFRMu, TH2D *fhDFREl, int fm, int fe);
double leptonEfficiency(double pt, double eta, TH2D *fhDEffMu, TH2D *fhDEffEl, int lid, int syst = 0);
double hzz2l_cuts(double mass, int opt);
double nVtxScaleFactor(TH1D *fhDNvtx, int nvtx);
double nPUScaleFactor2011(TH1D *fhDPU, int npu);
double nPUScaleFactor2012(TH1D *fhDPU, int npu);
double mt_atlas(LorentzVector dilep, double met, double metPhi);
double poorManMetSyst(LorentzVector l1, LorentzVector l2, LorentzVector l3, 
                      int lid1, int lid2, int lid3, double met, double metPhi,
                      LorentzVector j1, LorentzVector j2, LorentzVector j3, LorentzVector j4, int nsyst);

double Unroll2VarTo1Var(double varA, double varB, int binsA, int binsB, bool verbose){
  // variables are supposed to be in the [0,1] range
  if(varA < 0 || varA > 1 || varB < 0 || varB > 1) return -1.0;
  if(varA < 0 || varA > 1 || varB < 0 || varB > 1) assert(0);

  double newVarA = (int)(varA*binsA);
  double newVarB = (int)(varB*binsB);
  double finalVar = binsB*newVarA+newVarB;
  double binRange = binsA*binsB-1;
  if(binsB > binsA) {finalVar = binsA*newVarB+newVarA;}
  if(verbose == true) std::cout << newVarA  << " " << varA     << " " <<  binsA     	   << " " <<
                                   newVarB  << " " << varB     << " " <<  binsB      	   << " " <<
			           finalVar << " " << binRange << " " << finalVar/binRange << std::endl;

  // output is forced to be in the range [-1,1]
  return ((finalVar/binRange)-0.5)*2.0;
}

TH1F* Unroll2DTo1D(TH2F* h2, const char* hname) {
  
  unsigned int nbinsX = h2->GetXaxis()->GetNbins();
  unsigned int nbinsY = h2->GetYaxis()->GetNbins();
  unsigned int nbins  = nbinsX*nbinsY;

  TH1F* h1_unroll = new TH1F(hname, hname, nbins, 0.5, 0.5+nbins);

  for(unsigned int x=1; x <= nbinsX; ++x) {
    for(unsigned int y=1; y <= nbinsY; ++y) {
      h1_unroll->SetBinContent( (x-1)*nbinsY + y, h2->GetBinContent(x, y) );
      h1_unroll->SetBinError( (x-1)*nbinsY + y, h2->GetBinError(x, y) );
    }
  }

  // storing entries is needed to calculate stat up bounding for empty bins 
  h1_unroll->SetEntries(h2->GetEntries()); ;

  return h1_unroll;
}


TH1D * SmurfRebin(const TH1D *old, const unsigned int rebin)
{
  TH1D *h1_new = (TH1D*)old->Clone(TString(old->GetName()) + "_rebinned");
  TH1D *h1_old_tmp = (TH1D*)old->Clone("old_tmp");
  h1_old_tmp->Rebin(rebin);

  for(int i=1; i<=h1_new->GetNbinsX(); ++i){
    int bin = h1_old_tmp->FindBin(h1_new->GetBinCenter(i));
    h1_new->SetBinContent(i, h1_old_tmp->GetBinContent(bin)/rebin);
    h1_new->SetBinError  (i, h1_old_tmp->GetBinError  (bin)/rebin);
  }

  delete h1_old_tmp;
  return h1_new;
}

double CalcGammaMRstar(LorentzVector ja, LorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  LorentzVector jaT(ja.Px(),ja.Py(),0.0,ja.Pt());
  LorentzVector jbT(jb.Px(),jb.Py(),0.0,jb.Pt());
  double ATBT = (jaT+jbT).mag2();

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  //gamma times MRstar
  temp *= mygamma;

  return temp;
}

double CalcMTR(LorentzVector ja, LorentzVector jb, double metValue, double metPhi){

  Double_t Pt0    = ja.Pt();
  Double_t Pt1    = jb.Pt();
  Double_t etmiss = metValue;

  Double_t Px0    = ja.Px();
  Double_t Px1    = jb.Px();
  Double_t metx   = metValue*cos(metPhi);
  Double_t Py0    = ja.Py();
  Double_t Py1    = jb.Py();
  Double_t mety   = metValue*sin(metPhi);

  return TMath::Sqrt(TMath::Abs(0.5*etmiss*(Pt0 + Pt1) - 0.5*(metx*(Px0 + Px1) + mety*(Py0 + Py1))));

}

Double_t CalcMR(LorentzVector ja, LorentzVector jb){

  Double_t E0  = ja.E();
  Double_t E1  = jb.E();
  Double_t Pz0 = ja.Pz();
  Double_t Pz1 = jb.Pz();

  Double_t den =  TMath::Power(Pz0-Pz1, 2) - TMath::Power(E0-E1,2);
  if(den <= 0) return -100.0;

  return 2.0*TMath::Sqrt(TMath::Power(E0*Pz1 - E1*Pz0, 2)/den);
}

double DeltaPhi(double phi1, double phi2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].
  double dphi = TMath::Abs(phi1-phi2);
  while (dphi>TMath::Pi())
    dphi = TMath::Abs(dphi - TMath::TwoPi());
  return(dphi);
}
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

double nVtxScaleFactor(TH1D *fhDNvtx, int nvtx){
  double mynvtx = TMath::Min((double)nvtx,19.499);
  Int_t nvtxbin = fhDNvtx->GetXaxis()->FindBin(mynvtx);
  return fhDNvtx->GetBinContent(nvtxbin);
}

double nPUScaleFactor2011(TH1D *fhDPU, int npu){
  double mynpu = TMath::Min((double)npu,49.499);
  Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
  return fhDPU->GetBinContent(npuxbin);
}

double nPUScaleFactor2012(TH1D *fhDPU, int npu){
  double mynpu = TMath::Min((double)npu,49.499);
  Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
  return fhDPU->GetBinContent(npuxbin);
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
    if     (type == 0) weight = weight * 0	;
    else if(type == 1) weight = weight * 0.14306;
    else if(type == 2) weight = weight * 0.47351;
    else if(type == 3) weight = weight * 0.65951;
    else if(type == 4) weight = weight * 0.87082;
    else if(type == 5) weight = weight * 0.97577;
    else if(type == 6) weight = weight * 1.02573;
    else if(type == 7) weight = weight * 1.01917;
    else if(type == 8) weight = weight * 1.09291;
    else if(type == 9) weight = weight * 1.17922;
    else if(type ==10) weight = weight * 1.31628;
    else if(type ==11) weight = weight * 1.56301;
    else if(type ==12) weight = weight * 1.86699;
    else if(type ==13) weight = weight * 2.54774;
    else if(type ==14) weight = weight * 3.13798;
    else if(type ==15) weight = weight * 4.33400;
    else if(type ==16) weight = weight * 5.11852;
    else if(type ==17) weight = weight * 8.05315;
    else if(type ==18) weight = weight * 10.1783;
    else if(type ==19) weight = weight * 12.7989;
    else if(type ==20) weight = weight * 18.2182;
    else if(type ==21) weight = weight * 33.3586;
    else if(type ==22) weight = weight * 43.8323;
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

double enhancementFactor(double mass, int type){

if(type == 0){ // 4th-generation gg->H enhancement factor
  if     (mass==110.000) return 10.0302;
  else if(mass==111.000) return 10.0411;
  else if(mass==112.000) return 10.047 ;
  else if(mass==113.000) return 10.0479;
  else if(mass==114.000) return 10.0488;
  else if(mass==115.000) return 10.0386;
  else if(mass==116.000) return 10.0224;
  else if(mass==117.000) return 10.0057;
  else if(mass==118.000) return 9.98257;
  else if(mass==119.000) return 9.95272;
  else if(mass==120.000) return 9.92183;
  else if(mass==121.000) return 9.92054;
  else if(mass==122.000) return 9.92537;
  else if(mass==123.000) return 9.91783;
  else if(mass==124.000) return 9.91003;
  else if(mass==125.000) return 9.89549;
  else if(mass==126.000) return 9.88048;
  else if(mass==127.000) return 9.8583 ;
  else if(mass==128.000) return 9.83539;
  else if(mass==129.000) return 9.80488;
  else if(mass==130.000) return 9.77337;
  else if(mass==131.000) return 9.77698;
  else if(mass==132.000) return 9.77356;
  else if(mass==133.000) return 9.77003;
  else if(mass==134.000) return 9.75904;
  else if(mass==135.000) return 9.74771;
  else if(mass==136.000) return 9.73602;
  else if(mass==137.000) return 9.72397;
  else if(mass==138.000) return 9.696  ;
  else if(mass==139.000) return 9.67506;
  else if(mass==140.000) return 9.64551;
  else if(mass==141.000) return 9.62343;
  else if(mass==142.000) return 9.59253;
  else if(mass==143.000) return 9.56897;
  else if(mass==144.000) return 9.52797;
  else if(mass==145.000) return 9.49423;
  else if(mass==146.000) return 9.48201;
  else if(mass==147.000) return 9.4781 ;
  else if(mass==148.000) return 9.47407;
  else if(mass==149.000) return 9.46103;
  else if(mass==150.000) return 9.44762;
  else if(mass==151.000) return 9.42085;
  else if(mass==152.000) return 9.41234;
  else if(mass==153.000) return 9.39424;
  else if(mass==154.000) return 9.37185;
  else if(mass==155.000) return 9.36192;
  else if(mass==156.000) return 9.34231;
  else if(mass==157.000) return 9.3231 ;
  else if(mass==158.000) return 9.31995;
  else if(mass==159.000) return 9.30081;
  else if(mass==160.000) return 9.30617;
  else if(mass==162.000) return 9.36146;
  else if(mass==164.000) return 9.42705;
  else if(mass==166.000) return 9.46276;
  else if(mass==168.000) return 9.45674;
  else if(mass==170.000) return 9.44495;
  else if(mass==172.000) return 9.43446;
  else if(mass==174.000) return 9.41176;
  else if(mass==176.000) return 9.40664;
  else if(mass==178.000) return 9.38899;
  else if(mass==180.000) return 9.36341;
  else if(mass==182.000) return 9.36975;
  else if(mass==184.000) return 9.37451;
  else if(mass==186.000) return 9.38506;
  else if(mass==188.000) return 9.3719 ;
  else if(mass==190.000) return 9.36228;
  else if(mass==192.000) return 9.36185;
  else if(mass==194.000) return 9.34473;
  else if(mass==196.000) return 9.32386;
  else if(mass==198.000) return 9.29903;
  else if(mass==200.000) return 9.25891;
  else if(mass==202.000) return 9.24455;
  else if(mass==204.000) return 9.2222 ;
  else if(mass==206.000) return 9.18765;
  else if(mass==208.000) return 9.14973;
  else if(mass==210.000) return 9.10438;
  else if(mass==212.000) return 9.07991;
  else if(mass==214.000) return 9.05045;
  else if(mass==216.000) return 9.00763;
  else if(mass==218.000) return 8.96315;
  else if(mass==220.000) return 8.90858;
  else if(mass==222.000) return 8.88995;
  else if(mass==224.000) return 8.86628;
  else if(mass==226.000) return 8.83296;
  else if(mass==228.000) return 8.79397;
  else if(mass==230.000) return 8.75128;
  else if(mass==232.000) return 8.72623;
  else if(mass==234.000) return 8.69796;
  else if(mass==236.000) return 8.65929;
  else if(mass==238.000) return 8.61927;
  else if(mass==240.000) return 8.57302;
  else if(mass==242.000) return 8.54385;
  else if(mass==244.000) return 8.50634;
  else if(mass==246.000) return 8.46514;
  else if(mass==248.000) return 8.41761;
  else if(mass==250.000) return 8.36353;
  else if(mass==252.000) return 8.33487;
  else if(mass==254.000) return 8.30012;
  else if(mass==256.000) return 8.26169;
  else if(mass==258.000) return 8.21681;
  else if(mass==260.000) return 8.17057;
  else if(mass==262.000) return 8.13738;
  else if(mass==264.000) return 8.10322;
  else if(mass==266.000) return 8.05707;
  else if(mass==268.000) return 8.0124 ;
  else if(mass==270.000) return 7.96089;
  else if(mass==272.000) return 7.92079;
  else if(mass==274.000) return 7.87683;
  else if(mass==276.000) return 7.82609;
  else if(mass==278.000) return 7.77126;
  else if(mass==280.000) return 7.71513;
  else if(mass==282.000) return 7.68018;
  else if(mass==284.000) return 7.64147;
  else if(mass==286.000) return 7.59892;
  else if(mass==288.000) return 7.55245;
  else if(mass==290.000) return 7.50196;
  else if(mass==295.000) return 7.39919;
  else if(mass==300.000) return 7.26672;
  else if(mass==305.000) return 7.15492;
  else if(mass==310.000) return 7.01981;
  else if(mass==315.000) return 6.90487;
  else if(mass==320.000) return 6.76457;
  else if(mass==325.000) return 6.61864;
  else if(mass==330.000) return 6.44283;
  else if(mass==335.000) return 6.28702;
  else if(mass==340.000) return 6.14195;
  else if(mass==345.000) return 5.88764;
  else if(mass==350.000) return 5.68083;
  else if(mass==360.000) return 5.33672;
  else if(mass==370.000) return 5.12601;
  else if(mass==380.000) return 4.94263;
  else if(mass==390.000) return 4.81928;
  else if(mass==400.000) return 4.71949;
  else if(mass==420.000) return 4.57859;
  else if(mass==440.000) return 4.50067;
  else if(mass==450.000) return 4.45017;
  else if(mass==460.000) return 4.43816;
  else if(mass==480.000) return 4.38596;
  else if(mass==500.000) return 4.35756;
  else if(mass==520.000) return 4.45904;
  else if(mass==540.000) return 4.40678;
  else if(mass==550.000) return 4.28314;
  else if(mass==560.000) return 4.37644;
  else if(mass==580.000) return 4.42191;
  else if(mass==600.000) return 4.28528;
  else assert(0);
}
else if(type == 1){ // 4th-generation BR(H->WW) enhancement factor
  if     (mass==110.000) return 0.149378;
  else if(mass==111.000) return 0.160256;
  else if(mass==112.000) return 0.167208;
  else if(mass==113.000) return 0.170996;
  else if(mass==114.000) return 0.17268 ;
  else if(mass==115.000) return 0.172434;
  else if(mass==116.000) return 0.171162;
  else if(mass==117.000) return 0.168692;
  else if(mass==118.000) return 0.166102;
  else if(mass==119.000) return 0.162692;
  else if(mass==120.000) return 0.158741;
  else if(mass==121.000) return 0.167949;
  else if(mass==122.000) return 0.174706;
  else if(mass==123.000) return 0.179459;
  else if(mass==124.000) return 0.1835	;
  else if(mass==125.000) return 0.186111;
  else if(mass==126.000) return 0.187554;
  else if(mass==127.000) return 0.1888	;
  else if(mass==128.000) return 0.189179;
  else if(mass==129.000) return 0.18951 ;
  else if(mass==130.000) return 0.18918 ;
  else if(mass==131.000) return 0.200093;
  else if(mass==132.000) return 0.209796;
  else if(mass==133.000) return 0.217879;
  else if(mass==134.000) return 0.225117;
  else if(mass==135.000) return 0.231638;
  else if(mass==136.000) return 0.237541;
  else if(mass==137.000) return 0.242912;
  else if(mass==138.000) return 0.247819;
  else if(mass==139.000) return 0.252319;
  else if(mass==140.000) return 0.256461;
  else if(mass==141.000) return 0.268451;
  else if(mass==142.000) return 0.279558;
  else if(mass==143.000) return 0.289876;
  else if(mass==144.000) return 0.299485;
  else if(mass==145.000) return 0.30897 ;
  else if(mass==146.000) return 0.32818 ;
  else if(mass==147.000) return 0.34625 ;
  else if(mass==148.000) return 0.362727;
  else if(mass==149.000) return 0.378792;
  else if(mass==150.000) return 0.393983;
  else if(mass==151.000) return 0.414226;
  else if(mass==152.000) return 0.436141;
  else if(mass==153.000) return 0.460317;
  else if(mass==154.000) return 0.487742;
  else if(mass==155.000) return 0.518239;
  else if(mass==156.000) return 0.553922;
  else if(mass==157.000) return 0.589499;
  else if(mass==158.000) return 0.637631;
  else if(mass==159.000) return 0.688136;
  else if(mass==160.000) return 0.748899;
  else if(mass==162.000) return 0.836691;
  else if(mass==164.000) return 0.878788;
  else if(mass==166.000) return 0.900208;
  else if(mass==168.000) return 0.9139	;
  else if(mass==170.000) return 0.923237;
  else if(mass==172.000) return 0.929387;
  else if(mass==174.000) return 0.934375;
  else if(mass==176.000) return 0.938285;
  else if(mass==178.000) return 0.942978;
  else if(mass==180.000) return 0.946352;
  else if(mass==182.000) return 0.952381;
  else if(mass==184.000) return 0.958237;
  else if(mass==186.000) return 0.96256 ;
  else if(mass==188.000) return 0.963885;
  else if(mass==190.000) return 0.966921;
  else if(mass==192.000) return 0.970466;
  else if(mass==194.000) return 0.970302;
  else if(mass==196.000) return 0.969456;
  else if(mass==198.000) return 0.970509;
  else if(mass==200.000) return 0.968961;
  else if(mass==202.000) return 0.972011;
  else if(mass==204.000) return 0.97377 ;
  else if(mass==206.000) return 0.975549;
  else if(mass==208.000) return 0.976	;
  else if(mass==210.000) return 0.975104;
  else if(mass==212.000) return 0.975867;
  else if(mass==214.000) return 0.976634;
  else if(mass==216.000) return 0.977406;
  else if(mass==218.000) return 0.978182;
  else if(mass==220.000) return 0.977591;
  else if(mass==222.000) return 0.979494;
  else if(mass==224.000) return 0.980028;
  else if(mass==226.000) return 0.980563;
  else if(mass==228.000) return 0.9811	;
  else if(mass==230.000) return 0.981638;
  else if(mass==232.000) return 0.982178;
  else if(mass==234.000) return 0.98272 ;
  else if(mass==236.000) return 0.98187 ;
  else if(mass==238.000) return 0.982411;
  else if(mass==240.000) return 0.982955;
  else if(mass==242.000) return 0.983239;
  else if(mass==244.000) return 0.984922;
  else if(mass==246.000) return 0.98661 ;
  else if(mass==248.000) return 0.986895;
  else if(mass==250.000) return 0.988588;
  else if(mass==252.000) return 0.988017;
  else if(mass==254.000) return 0.988857;
  else if(mass==256.000) return 0.988286;
  else if(mass==258.000) return 0.989127;
  else if(mass==260.000) return 0.988555;
  else if(mass==262.000) return 0.989112;
  else if(mass==264.000) return 0.988252;
  else if(mass==266.000) return 0.987393;
  else if(mass==268.000) return 0.987948;
  else if(mass==270.000) return 0.987088;
  else if(mass==272.000) return 0.988506;
  else if(mass==274.000) return 0.988506;
  else if(mass==276.000) return 0.988506;
  else if(mass==278.000) return 0.989928;
  else if(mass==280.000) return 0.989928;
  else if(mass==282.000) return 0.989353;
  else if(mass==284.000) return 0.990202;
  else if(mass==286.000) return 0.989625;
  else if(mass==288.000) return 0.989049;
  else if(mass==290.000) return 0.989899;
  else if(mass==295.000) return 0.989884;
  else if(mass==300.000) return 0.988439;
  else if(mass==305.000) return 0.986975;
  else if(mass==310.000) return 0.985507;
  else if(mass==315.000) return 0.984783;
  else if(mass==320.000) return 0.985486;
  else if(mass==325.000) return 0.984761;
  else if(mass==330.000) return 0.985465;
  else if(mass==335.000) return 0.985444;
  else if(mass==340.000) return 0.983988;
  else if(mass==345.000) return 0.976642;
  else if(mass==350.000) return 0.933432;
  else if(mass==360.000) return 0.827957;
  else if(mass==370.000) return 0.753185;
  else if(mass==380.000) return 0.702791;
  else if(mass==390.000) return 0.666667;
  else if(mass==400.000) return 0.642612;
  else if(mass==420.000) return 0.613475;
  else if(mass==440.000) return 0.597473;
  else if(mass==450.000) return 0.592928;
  else if(mass==460.000) return 0.590164;
  else if(mass==480.000) return 0.587912;
  else if(mass==500.000) return 0.589744;
  else if(mass==520.000) return 0.595247;
  else if(mass==540.000) return 0.599636;
  else if(mass==550.000) return 0.601272;
  else if(mass==560.000) return 0.605072;
  else if(mass==580.000) return 0.612613;
  else if(mass==600.000) return 0.620072;
  else assert(0);
}
else if(type == 2){ // fermiophobic BR(H->WW) enhancement factor
  if     (mass==90 ) return 214.832536;
  else if(mass==95 ) return 127.754237;
  else if(mass==100) return  66.486486;
  else if(mass==105) return  33.621399;
  else if(mass==110) return  17.717842;
  else if(mass==110) return  16.686160;
  else if(mass==111) return  15.714286;
  else if(mass==111) return  14.810345;
  else if(mass==112) return  13.977273;
  else if(mass==112) return  13.180428;
  else if(mass==113) return  12.467532;
  else if(mass==113) return  11.784741;
  else if(mass==114) return  11.159794;
  else if(mass==114) return  10.548112;
  else if(mass==115) return  10.000000;
  else if(mass==115) return   9.496718;
  else if(mass==116) return   9.004149;
  else if(mass==116) return   8.519608;
  else if(mass==117) return   8.121495;
  else if(mass==117) return   7.758929;
  else if(mass==118) return   7.372881;
  else if(mass==118) return   7.016129;
  else if(mass==119) return   6.692308;
  else if(mass==119) return   6.397059;
  else if(mass==120) return   6.083916;
  else if(mass==120) return   5.838926;
  else if(mass==121) return   5.576923;
  else if(mass==121) return   5.337423;
  else if(mass==122) return   5.117647;
  else if(mass==122) return   4.915254;
  else if(mass==123) return   4.702703;
  else if(mass==123) return   4.526042;
  else if(mass==124) return   4.345000;
  else if(mass==124) return   4.177885;
  else if(mass==125) return   4.023148;
  else if(mass==125) return   3.879464;
  else if(mass==126) return   3.729614;
  else if(mass==126) return   3.601660;
  else if(mass==127) return   3.472000;
  else if(mass==127) return   3.351351;
  else if(mass==128) return   3.238806;
  else if(mass==128) return   3.133574;
  else if(mass==129) return   3.034965;
  else if(mass==129) return   2.938983;
  else if(mass==130) return   2.842623;
  else if(mass==130) return   2.761146;
  else if(mass==131) return   2.675926;
  else if(mass==131) return   2.603604;
  else if(mass==132) return   2.527697;
  else if(mass==132) return   2.456091;
  else if(mass==133) return   2.388430;
  else if(mass==133) return   2.324397;
  else if(mass==134) return   2.263708;
  else if(mass==134) return   2.211735;
  else if(mass==135) return   2.151365;
  else if(mass==135) return   2.099274;
  else if(mass==136) return   2.049645;
  else if(mass==136) return   2.002309;
  else if(mass==137) return   1.957111;
  else if(mass==137) return   1.913907;
  else if(mass==138) return   1.874730;
  else if(mass==138) return   1.835095;
  else if(mass==139) return   1.797101;
  else if(mass==139) return   1.760649;
  else if(mass==140) return   1.727634;
  else if(mass==141) return   1.663480;
  else if(mass==142) return   1.604052;
  else if(mass==143) return   1.548845;
  else if(mass==144) return   1.497427;
  else if(mass==145) return   1.453488;
  else if(mass==146) return   1.412238;
  else if(mass==147) return   1.373437;
  else if(mass==148) return   1.334848;
  else if(mass==149) return   1.301915;
  else if(mass==150) return   1.270774;
  else if(mass==151) return   1.241283;
  else if(mass==152) return   1.214674;
  else if(mass==153) return   1.189153;
  else if(mass==154) return   1.166452;
  else if(mass==155) return   1.143396;
  else if(mass==156) return   1.122549;
  else if(mass==157) return   1.102625;
  else if(mass==158) return   1.082462;
  else if(mass==159) return   1.064407;
  else if(mass==160) return   1.048458;
  else if(mass==162) return   1.026511;
  else if(mass==164) return   1.017764;
  else if(mass==166) return   1.014553;
  else if(mass==168) return   1.012448;
  else if(mass==170) return   1.011411;
  else if(mass==172) return   1.009346;
  else if(mass==174) return   1.009375;
  else if(mass==176) return   1.007322;
  else if(mass==178) return   1.007392;
  else if(mass==180) return   1.006438;
  else if(mass==182) return   1.005537;
  else if(mass==184) return   1.004640;
  else if(mass==186) return   1.003623;
  else if(mass==188) return   1.003736;
  else if(mass==190) return   1.002545;
  else if(mass==192) return   1.002591;
  else if(mass==194) return   1.002628;
  else if(mass==196) return   1.002656;
  else if(mass==198) return   1.002681;
  else if(mass==200) return   1.001350;
  else if(mass==202) return   1.001359;
  else if(mass==204) return   1.001366;
  else if(mass==206) return   1.002747;
  else if(mass==208) return   1.002759;
  else if(mass==210) return   1.001383;
  else if(mass==212) return   1.001387;
  else if(mass==214) return   1.001391;
  else if(mass==216) return   1.001395;
  else if(mass==218) return   1.002797;
  else if(mass==220) return   1.001401;
  else if(mass==222) return   1.002809;
  else if(mass==224) return   1.002813;
  else if(mass==226) return   1.002817;
  else if(mass==228) return   1.001410;
  else if(mass==230) return   1.001412;
  else if(mass==232) return   1.002829;
  else if(mass==234) return   1.002833;
  else if(mass==236) return   1.001416;
  else if(mass==238) return   1.001418;
  else if(mass==240) return   1.002841;
  else if(mass==242) return   1.001420;
  else if(mass==244) return   1.001422;
  else if(mass==246) return   1.002849;
  else if(mass==248) return   1.001425;
  else if(mass==250) return   1.001427;
  else if(mass> 250) return   1.000000;
  else assert(0);
}
else if(type == 3){ // 4th-generation BR(H->tautau) enhancement factor
  if     (mass==110.000) return 0.706983;
  else if(mass==111.000) return 0.710792;
  else if(mass==112.000) return 0.7146	;
  else if(mass==113.000) return 0.718409;
  else if(mass==114.000) return 0.722218;
  else if(mass==115.000) return 0.726027;
  else if(mass==116.000) return 0.729835;
  else if(mass==117.000) return 0.733644;
  else if(mass==118.000) return 0.737453;
  else if(mass==119.000) return 0.741261;
  else if(mass==120.000) return 0.74507 ;
  else if(mass==121.000) return 0.758519;
  else if(mass==122.000) return 0.771968;
  else if(mass==123.000) return 0.785418;
  else if(mass==124.000) return 0.798867;
  else if(mass==125.000) return 0.812316;
  else if(mass==126.000) return 0.825765;
  else if(mass==127.000) return 0.839214;
  else if(mass==128.000) return 0.852664;
  else if(mass==129.000) return 0.866113;
  else if(mass==130.000) return 0.879562;
  else if(mass==131.000) return 0.909402;
  else if(mass==132.000) return 0.939243;
  else if(mass==133.000) return 0.969083;
  else if(mass==134.000) return 0.998924;
  else if(mass==135.000) return 1.02876 ;
  else if(mass==136.000) return 1.0586	;
  else if(mass==137.000) return 1.08844 ;
  else if(mass==138.000) return 1.11829 ;
  else if(mass==139.000) return 1.14813 ;
  else if(mass==140.000) return 1.17797 ;
  else if(mass==141.000) return 1.22973 ;
  else if(mass==142.000) return 1.28149 ;
  else if(mass==143.000) return 1.33326 ;
  else if(mass==144.000) return 1.38502 ;
  else if(mass==145.000) return 1.43678 ;
  else if(mass==146.000) return 1.50898 ;
  else if(mass==147.000) return 1.58117 ;
  else if(mass==148.000) return 1.65336 ;
  else if(mass==149.000) return 1.72556 ;
  else if(mass==150.000) return 1.79775 ;
  else if(mass==151.000) return 1.88957 ;
  else if(mass==152.000) return 1.98649 ;
  else if(mass==153.000) return 2.10526 ;
  else if(mass==154.000) return 2.21849 ;
  else if(mass==155.000) return 2.35238 ;
  else if(mass==156.000) return 2.5	;
  else if(mass==157.000) return 2.67606 ;
  else if(mass==158.000) return 2.83951 ;
  else if(mass==159.000) return 3.05019 ;
  else if(mass==160.000) return 3.25758 ;
  else if(mass==162.000) return 3.62114 ;
  else if(mass==164.000) return 3.84076 ;
  else if(mass==166.000) return 3.96774 ;
  else if(mass==168.000) return 4.04762 ;
  else if(mass==170.000) return 4.11317 ;
  else if(mass==172.000) return 4.15758 ;
  else if(mass==174.000) return 4.19441 ;
  else if(mass==176.000) return 4.23188 ;
  else if(mass==178.000) return 4.26101 ;
  else if(mass==180.000) return 4.29301 ;
  else if(mass==182.000) return 4.32584 ;
  else if(mass==184.000) return 4.32432 ;
  else if(mass==186.000) return 4.32494 ;
  else if(mass==188.000) return 4.34243 ;
  else if(mass==190.000) return 4.3617	;
  else if(mass==192.000) return 4.37162 ;
  else if(mass==194.000) return 4.38154 ;
  else if(mass==196.000) return 4.39422 ;
  else if(mass==198.000) return 4.40965 ;
  else if(mass==200.000) return 4.42509 ;
  else if(mass==202.000) return 4.43751 ;
  else if(mass==204.000) return 4.44992 ;
  else if(mass==206.000) return 4.46234 ;
  else if(mass==208.000) return 4.47476 ;
  else if(mass==210.000) return 4.48718 ;
  else if(mass==212.000) return 4.50097 ;
  else if(mass==214.000) return 4.51476 ;
  else if(mass==216.000) return 4.52854 ;
  else if(mass==218.000) return 4.54233 ;
  else if(mass==220.000) return 4.55612 ;
  else if(mass==222.000) return 4.56037 ;
  else if(mass==224.000) return 4.56463 ;
  else if(mass==226.000) return 4.56888 ;
  else if(mass==228.000) return 4.57313 ;
  else if(mass==230.000) return 4.57738 ;
  else if(mass==232.000) return 4.58742 ;
  else if(mass==234.000) return 4.59746 ;
  else if(mass==236.000) return 4.6075	;
  else if(mass==238.000) return 4.61754 ;
  else if(mass==240.000) return 4.62759 ;
  else if(mass==242.000) return 4.62963 ;
  else if(mass==244.000) return 4.63167 ;
  else if(mass==246.000) return 4.63371 ;
  else if(mass==248.000) return 4.63575 ;
  else if(mass==250.000) return 4.63779 ;
  else if(mass==252.000) return 4.64595 ;
  else if(mass==254.000) return 4.65411 ;
  else if(mass==256.000) return 4.66226 ;
  else if(mass==258.000) return 4.67042 ;
  else if(mass==260.000) return 4.67857 ;
  else if(mass==262.000) return 4.68086 ;
  else if(mass==264.000) return 4.68314 ;
  else if(mass==266.000) return 4.68543 ;
  else if(mass==268.000) return 4.68771 ;
  else if(mass==270.000) return 4.69	;
  else if(mass==272.000) return 4.69409 ;
  else if(mass==274.000) return 4.69819 ;
  else if(mass==276.000) return 4.70228 ;
  else if(mass==278.000) return 4.70637 ;
  else if(mass==280.000) return 4.71047 ;
  else if(mass==282.000) return 4.71275 ;
  else if(mass==284.000) return 4.71503 ;
  else if(mass==286.000) return 4.71731 ;
  else if(mass==288.000) return 4.7196	;
  else if(mass==290.000) return 4.72188 ;
  else if(mass==295.000) return 4.73151 ;
  else if(mass==300.000) return 4.74114 ;
  else if(mass==305.000) return 4.76578 ;
  else if(mass==310.000) return 4.79042 ;
  else if(mass==315.000) return 4.79717 ;
  else if(mass==320.000) return 4.80392 ;
  else if(mass==325.000) return 4.81759 ;
  else if(mass==330.000) return 4.83126 ;
  else if(mass==335.000) return 4.84288 ;
  else if(mass==340.000) return 4.86538 ;
  else if(mass==345.000) return 4.88048 ;
  else if(mass==350.000) return 4.64286 ;
  else if(mass==360.000) return 4.06619 ;
  else if(mass==370.000) return 3.65079 ;
  else if(mass==380.000) return 3.38235 ;
  else if(mass==390.000) return 3.1877	;
  else if(mass==400.000) return 3.03873 ;
  else if(mass==420.000) return 2.8642	;
  else if(mass==440.000) return 2.77359 ;
  else if(mass==450.000) return 2.75381 ;
  else if(mass==460.000) return 2.73404 ;
  else if(mass==480.000) return 2.71006 ;
  else if(mass==500.000) return 2.71895 ;
  else if(mass==520.000) return 2.74032 ;
  else if(mass==540.000) return 2.76168 ;
  else if(mass==550.000) return 2.77236 ;
  else if(mass==560.000) return 2.79044 ;
  else if(mass==580.000) return 2.82659 ;
  else if(mass==600.000) return 2.86274 ;
  else assert(0);
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
  if(pt<=10) return 1.0;
  double mypt   = TMath::Min(pt,49.999);
  double myeta  = TMath::Min(fabs(eta),2.4999);
  double prob = 1.0;
  if     (TMath::Abs(lid) == 13){
    Int_t ptbin = fhDEffMu->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDEffMu->GetYaxis()->FindBin(myeta);	 
    prob = fhDEffMu->GetBinContent(ptbin,etabin);
    if     (syst > 0) prob = prob + sqrt(fhDEffMu->GetBinError(ptbin,etabin)*fhDEffMu->GetBinError(ptbin,etabin) + 0.015*0.015);
    else if(syst < 0) prob = prob - sqrt(fhDEffMu->GetBinError(ptbin,etabin)*fhDEffMu->GetBinError(ptbin,etabin) + 0.015*0.015);
  }
  else if(TMath::Abs(lid) == 11){
    Int_t ptbin = fhDEffEl->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDEffEl->GetYaxis()->FindBin(myeta);	 
    prob = fhDEffEl->GetBinContent(ptbin,etabin);
    if     (syst > 0) prob = prob + sqrt(fhDEffEl->GetBinError(ptbin,etabin)*fhDEffEl->GetBinError(ptbin,etabin) + 0.02*0.02);
    else if(syst < 0) prob = prob - sqrt(fhDEffEl->GetBinError(ptbin,etabin)*fhDEffEl->GetBinError(ptbin,etabin) + 0.02*0.02);
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

double mt_atlas(LorentzVector dilep, double met, double metPhi){
  double deltaPhi = TMath::Abs(metPhi-dilep.Phi());
  while(deltaPhi>TMath::Pi()) deltaPhi = TMath::Abs(deltaPhi - 2*TMath::Pi());
  double aux = sqrt((dilep.M()*dilep.M()+dilep.Pt()*dilep.Pt())*met*met)-met*dilep.Pt()*cos(deltaPhi);
  
  double mt = dilep.M()*dilep.M()  + 2 * aux;
  if(mt >= 0) mt = sqrt(mt); else mt = 0;
  
  return mt;
}

double poorManMetSyst(LorentzVector l1, LorentzVector l2, LorentzVector l3, 
                      int lid1, int lid2, int lid3, double met, double metPhi,
                      LorentzVector j1, LorentzVector j2, LorentzVector j3, LorentzVector j4, int nsyst){
  double rec_x = met*cos(metPhi) + l1.Px() + l2.Px() + l3.Px();
  double rec_y = met*sin(metPhi) + l1.Py() + l2.Py() + l3.Py();
  if(j1.Pt() > 10) {
    rec_x = rec_x + j1.Px();
    rec_y = rec_y + j1.Py();
  }
  if(j2.Pt() > 10) {
    rec_x = rec_x + j2.Px();
    rec_y = rec_y + j2.Py();
  }
  if(j3.Pt() > 10) {
    rec_x = rec_x + j3.Px();
    rec_y = rec_y + j3.Py();
  }
  if(j4.Pt() > 10) {
    rec_x = rec_x + j4.Px();
    rec_y = rec_y + j4.Py();
  }
  double corr_x = met*cos(metPhi);
  double corr_y = met*sin(metPhi);
  if(nsyst > 0){
    corr_x = rec_x*1.1;
    if(abs(lid1) == 13) corr_x+= l1.Px()*1.002; if(abs(lid1) == 13) corr_y+= l1.Py()*1.002;
    if(abs(lid2) == 13) corr_x+= l2.Px()*1.002; if(abs(lid2) == 13) corr_y+= l2.Py()*1.002;
    if(abs(lid3) == 13) corr_x+= l3.Px()*1.002; if(abs(lid3) == 13) corr_y+= l3.Py()*1.002;
    if(abs(lid1) == 11) corr_x+= l1.Px()*1.006; if(abs(lid1) == 11) corr_y+= l1.Py()*1.006;
    if(abs(lid2) == 11) corr_x+= l2.Px()*1.006; if(abs(lid2) == 11) corr_y+= l2.Py()*1.006;
    if(abs(lid3) == 11) corr_x+= l3.Px()*1.006; if(abs(lid3) == 11) corr_y+= l3.Py()*1.006;
    if(j1.Pt() > 10)    corr_x+= j1.Px()*1.1;	if(j1.Pt() > 10)    corr_y+= j1.Py()*1.1;
    if(j2.Pt() > 10)    corr_x+= j2.Px()*1.1;	if(j2.Pt() > 10)    corr_y+= j2.Py()*1.1;
    if(j3.Pt() > 10)    corr_x+= j3.Px()*1.1;	if(j3.Pt() > 10)    corr_y+= j3.Py()*1.1;
    if(j4.Pt() > 10)    corr_x+= j4.Px()*1.1;	if(j4.Pt() > 10)    corr_y+= j4.Py()*1.1;
  } else if(nsyst < 0){
    corr_x = rec_x/1.1;
    if(abs(lid1) == 13) corr_x+= l1.Px()/1.002; if(abs(lid1) == 13) corr_y+= l1.Py()/1.002;
    if(abs(lid2) == 13) corr_x+= l2.Px()/1.002; if(abs(lid2) == 13) corr_y+= l2.Py()/1.002;
    if(abs(lid3) == 13) corr_x+= l3.Px()/1.002; if(abs(lid3) == 13) corr_y+= l3.Py()/1.002;
    if(abs(lid1) == 11) corr_x+= l1.Px()/1.006; if(abs(lid1) == 11) corr_y+= l1.Py()/1.006;
    if(abs(lid2) == 11) corr_x+= l2.Px()/1.006; if(abs(lid2) == 11) corr_y+= l2.Py()/1.006;
    if(abs(lid3) == 11) corr_x+= l3.Px()/1.006; if(abs(lid3) == 11) corr_y+= l3.Py()/1.006;
    if(j1.Pt() > 10)    corr_x+= j1.Px()/1.1;	if(j1.Pt() > 10)    corr_y+= j1.Py()/1.1;
    if(j2.Pt() > 10)    corr_x+= j2.Px()/1.1;	if(j2.Pt() > 10)    corr_y+= j2.Py()/1.1;
    if(j3.Pt() > 10)    corr_x+= j3.Px()/1.1;	if(j3.Pt() > 10)    corr_y+= j3.Py()/1.1;
    if(j4.Pt() > 10)    corr_x+= j4.Px()/1.1;	if(j4.Pt() > 10)    corr_y+= j4.Py()/1.1;
  } else return met;
  if(j4.Pt() > 10) return met;
  return sqrt(corr_x*corr_x+corr_y*corr_y);
}

