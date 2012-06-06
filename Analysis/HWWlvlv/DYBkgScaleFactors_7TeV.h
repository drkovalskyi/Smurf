static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 2.71735, 3.56158, 7.86249  };
  Double_t DYBkgScaleFactorHiggsSelection[3][18] = { 
    { 3.4296,3.26764,3.47829,3.9614,4.44737,4.38557,3.46856,4.33784,4.6754,4.81333,3.31946,2.2153,2.55011,5.97335,20.5566,17.7712,0.0205386,0.172233},
    { 1.83522,2.49191,3.2039,2.46916,2.45932,2.28181,2.25866,2.12551,2.72826,2.41616,3.95606,3.47765,6.10595,7.17139,7.50276,7.1685,2.5476,4.29931},
    { 4.40323,4.29175,5.23227,5.6026,5.71493,4.99674,4.69874,5.08022,4.7121,5.03275,4.7184,5.00227,5.78293,8.33302,8.36391,7.55756,6.62105,6.83446} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

static Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.63116, 1.19547, 1.18884  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][18] = { 
    { 1.71448,1.7155,1.71668,1.5849,1.53941,1.46405,1.71491,1.67676,1.66811,1.66981,1.98007,2.95374,3.24026,3.60468,2.54213,2.14438,29.8815,9.9311},
    { 1.44405,1.43721,1.43336,1.50776,1.54845,1.51378,1.52481,1.5197,1.50964,1.50799,1.36325,1.4194,1.43002,1.42207,1.39809,1.35287,1.32499,1.32683},
    { 1.27245,1.26499,1.25959,1.29261,1.33199,1.33934,1.37158,1.40799,1.40593,1.40519,1.42584,1.56328,1.40966,1.28816,1.25695,1.26015,1.51675,1.46791} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

