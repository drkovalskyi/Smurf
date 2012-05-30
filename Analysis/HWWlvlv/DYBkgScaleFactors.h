Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.26144, 3.80073, 1.78984  };
  Double_t DYBkgScaleFactorHiggsSelection[3][20] = { 
    { 6.51897,6.52643,6.24007,6.42449,7.34506,8.44812,8.44812,6.80505,8.37395,8.67237,7.88954,10.1815,2.35414,1.05196,0.21233,0.162931,0.995529,1.20637,1.01923,1.46994},
    { 2.14023,2.44957,2.29237,2.80915,2.76834,2.56533,2.56533,4.14131,3.66603,3.45753,2.3403,5.74994,5.82765,4.34686,3.19114,3.47229,6.87839,5.48238,4.81245,2.62489},
    { 5.67183,6.41765,7.92906,8.50831,8.46921,8.47601,8.47601,6.94949,7.35015,7.59407,5.96643,5.96643,5.48774,4.53087,3.44055,4.06897,6.33791,6.50007,5.23597,5.89405} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 20 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.26536, 1.2509, 1.18681  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][20] = { 
    { 1.35666,1.37006,1.39372,1.44766,1.4931,1.55411,1.55411,1.7144,1.69444,1.69608,1.65323,1.66856,1.85834,3.15929,3.53033,3.89553,2.35041,1.92724,1.66837,1.67513},
    { 1.78533,1.76983,1.78812,1.79213,1.8187,1.83836,1.83836,1.80962,1.83254,1.84929,1.88959,1.73069,1.70109,1.61822,1.65932,1.58424,1.57378,1.55034,1.43284,2.23302},
    { 1.39199,1.3869,1.37928,1.32855,1.28916,1.32089,1.32089,1.3323,1.31092,1.30099,1.32535,1.32535,1.43096,1.71816,1.70908,1.78995,1.93977,1.88529,1.45397,1.59944} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 20 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

