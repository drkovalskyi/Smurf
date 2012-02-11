#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double DileptonMassPreselectionCut( double mH ) {
  float dilmass_cut = 10000;   
  if     ( mH == 0) dilmass_cut = 14000;
  else if( mH >= 110 &&
           mH <  118 ) dilmass_cut =  70.0;
  else if( mH >= 118 &&
           mH <  125 ) dilmass_cut =  70.0;
  else if( mH >= 125 &&
           mH <= 130 ) dilmass_cut =  80.0;
  else if( mH >  130 &&
           mH <= 140 ) dilmass_cut =  90.0;
  else if( mH == 150 ) dilmass_cut = 100.0;
  else if( mH == 160 ) dilmass_cut = 100.0;
  else if( mH == 165 ) dilmass_cut = 100.0;
  else if( mH == 170 ) dilmass_cut = 100.0;
  else if( mH == 180 ) dilmass_cut = 110.0;
  else if( mH == 190 ) dilmass_cut = 120.0;
  else if( mH == 200 ) dilmass_cut = 130.0;
  else if( mH == 210 ) dilmass_cut = 140.0;
  else if( mH == 220 ) dilmass_cut = 150.0;
  else                 dilmass_cut = mH;
  
  return dilmass_cut;
}

int HiggsMassIndex ( double mH ) {

  int channel = 0;
  if     (mH == 0)   channel = 0;
  else if(mH >= 110 &&
          mH <  118) channel = 1;
  else if(mH == 118) channel = 2;
  else if(mH == 120) channel = 3;
  else if(mH == 122) channel = 4;
  else if(mH == 124) channel = 5;
  else if(mH == 126) channel = 6;
  else if(mH == 128) channel = 7;
  else if(mH == 130) channel = 8;
  else if(mH == 135) channel = 9;
  else if(mH == 140) channel = 10;
  else if(mH == 150) channel = 11;
  else if(mH == 160) channel = 12;
  else if(mH == 170) channel = 13;
  else if(mH == 180) channel = 14;
  else if(mH == 190) channel = 15;
  else if(mH == 200) channel = 16;
  else if(mH == 210) channel = 17;
  else if(mH == 220) channel = 18;
  else if(mH == 230) channel = 19;
  else if(mH == 250) channel = 20;
  else if(mH == 300) channel = 21;
  else if(mH == 350) channel = 22;
  else if(mH == 400) channel = 23;
  else if(mH == 450) channel = 24;
  else if(mH == 500) channel = 25;
  else if(mH == 550) channel = 26;
  else if(mH == 600) channel = 27;

  return channel;
}

double cutMassHigh ( double mH ) {
  double CutMassHigh[28]      = { 7000, 40, 40, 40, 41, 42, 43, 44, 45, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 28);
  return CutMassHigh[index];  
}

double cutPtMaxLow ( double mH ) {
  double CutPtMaxLow[28]      = { 20, 20, 20, 20, 21, 22, 23, 24, 25, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 28);
  return CutPtMaxLow[index];  
}

double cutPtMinLow ( double mH , Int_t finalstateType ) {
  double SFCutPtMinLow[28]      = { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double OFCutPtMinLow[28]      = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 12, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 28);
  if (finalstateType == 0 || finalstateType == 3) {
    return SFCutPtMinLow[index];
  } else {
    return OFCutPtMinLow[index];  
  }
}

double cutDeltaphiHigh ( double mH ) {
  double CutDeltaphiHigh[28]      = {180,115,115,115,110,105,100, 95, 90, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 28);
  return CutDeltaphiHigh[index];  
}

double cutMTLow ( double mH ) {
  double CutMTLow[28]      = { 0, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 28);
  return CutMTLow[index];  
}

double cutMTHigh ( double mH ) {
  double CutMTHigh[28]      = {70000,110,115,120,121,122,123,124,125,128,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 28);
  return CutMTHigh[index];  
}
