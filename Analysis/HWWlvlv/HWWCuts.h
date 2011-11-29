#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double DileptonMassPreselectionCut( double mH ) {
  float dilmass_cut = 10000;   
  if     ( mH >= 110 &&
           mH <  120 ) dilmass_cut =  70.0;
  else if( mH >= 120 &&
           mH <  125 ) dilmass_cut =  70.0;
  else if( mH >= 125 &&
           mH <  130 ) dilmass_cut =  75.0;
  else if( mH >= 130 &&
           mH <  140 ) dilmass_cut =  80.0;
  else if( mH == 140 ) dilmass_cut =  90.0;
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
          mH <  120) channel = 1;
  else if(mH >= 120 &&
          mH <  130) channel = 2;
  else if(mH >= 130 &&
          mH <  140) channel = 3;
  else if(mH == 140) channel = 4;
  else if(mH == 150) channel = 5;
  else if(mH == 160) channel = 6;
  else if(mH == 170) channel = 7;
  else if(mH == 180) channel = 8;
  else if(mH == 190) channel = 9;
  else if(mH == 200) channel = 10;
  else if(mH == 210) channel = 11;
  else if(mH == 220) channel = 12;
  else if(mH == 230) channel = 13;
  else if(mH == 250) channel = 14;
  else if(mH == 300) channel = 15;
  else if(mH == 350) channel = 16;
  else if(mH == 400) channel = 17;
  else if(mH == 450) channel = 18;
  else if(mH == 500) channel = 19;
  else if(mH == 550) channel = 20;
  else if(mH == 600) channel = 21;

  return channel;
}

double cutMassHigh ( double mH ) {
  double CutMassHigh[22]      = { 7000, 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 22);
  return CutMassHigh[index];  
}

double cutPtMaxLow ( double mH ) {
  double CutPtMaxLow[22]      = { 20, 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 22);
  return CutPtMaxLow[index];  
}

double cutPtMinLow ( double mH , Int_t finalstateType ) {
  double SFCutPtMinLow[22]      = { 15, 15, 15, 15, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double OFCutPtMinLow[22]      = { 10, 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 22);
  if (finalstateType == 0 || finalstateType == 3) {
    return SFCutPtMinLow[index]; 
  } else {
    return OFCutPtMinLow[index];  
  }
}

double cutDeltaphiHigh ( double mH ) {
  double CutDeltaphiHigh[22]      = {180, 115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 22);
  return CutDeltaphiHigh[index];  
}

double cutMTLow ( double mH ) {
  double CutMTLow[22]      = { 0, 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 22);
  return CutMTLow[index];  
}

double cutMTHigh ( double mH ) {
  double CutMTHigh[22]      = {70000, 110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  int index = HiggsMassIndex(mH);
  assert (index >= 0 && index < 22);
  return CutMTHigh[index];  
}

