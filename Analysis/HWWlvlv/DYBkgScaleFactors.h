Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.04789, 3.80056, 11.3794  };
  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { 
    { 3.71182,3.68507,4.32299,4.83513,3.37367,2.38566,1.99771,6.83772,22.7382,18.9507,0.162336,0.0437227},
    { 2.00234,3.44133,2.23393,2.54082,4.16117,3.66417,6.43091,7.56555,7.91577,7.56719,2.70016,4.55988},
    { 7.61378,11.3105,9.47995,11.4275,13.4377,14.1407,18.7867,19.1519,18.7112,16.5703,11.92,17.6918} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.66861, 1.25715, 1.38466  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { 
    { 1.80236,1.79617,1.77106,1.7624,2.08077,3.00986,3.56882,3.74332,2.83219,2.52016,5.17136,37.5824},
    { 1.57276,1.56505,1.60674,1.58084,1.50864,1.5491,1.64625,1.68053,1.59138,1.56197,1.45495,1.54134},
    { 1.80661,1.80089,1.70918,1.70557,1.82323,1.8564,1.96773,1.8899,1.86139,1.83612,1.86498,1.74341} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

