static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.41526, 4.28409, 2.00013  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 45.4944,63.3264,81.8409,85.0143,85.3732,93.0339,93.0339,93.6571,104.719,94.6107,81.7184,81.7184,37.9556,37.9556,18.8852,7.98769,4.40585,22.3879,20.6206,11.3141,11.3906},
    { 5.28728,8.9096,13.468,12.8478,12.7679,14.7039,14.7039,14.2863,14.2227,15.3428,14.341,14.341,11.7984,11.7984,8.81893,6.49877,7.52745,15.7283,16.6221,12.1183,7.97249},
    { 2.89198,3.36431,4.08532,5.09841,4.82074,4.17463,4.17463,4.02485,3.77453,3.96193,3.98473,3.98473,2.04409,2.04409,3.29151,5.25066,6.23518,6.15284,5.07825,6.95731,4.05701} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 21 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

static Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31535, 1.31933, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31826,1.31408,1.31136,1.31311,1.31627,1.32779,1.32779,1.32316,1.32393,1.32799,1.32952,1.32952,1.34217,1.34217,1.51267,1.88904,2.33764,1.38863,1.38075,1.35966,1.43876},
    { 1.43122,1.37226,1.34757,1.357,1.36434,1.35965,1.35965,1.37009,1.38551,1.37353,1.37238,1.37238,1.36408,1.36408,1.44914,1.51555,1.46809,1.35662,1.34468,1.3465,1.39283},
    { 1.46826,1.45563,1.44127,1.43507,1.44595,1.46294,1.46294,1.47979,1.50243,1.4939,1.49645,1.49645,1.57341,1.57341,1.8616,1.98306,1.59166,1.51026,1.47461,1.50291,1.60161} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 21 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.41526, 4.28409, 2.00013  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,73.3652},
    { 83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,51.579},
    { 19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,11.9131} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 21 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

static Double_t DYBkgScaleFactorBDTKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31535, 1.31933, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.32111},
    { 1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.32153},
    { 1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.54032} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 21 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

