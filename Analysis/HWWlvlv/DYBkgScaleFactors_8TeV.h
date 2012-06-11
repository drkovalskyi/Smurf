static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.48015, 3.83449, 1.78223  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 11.8114,12.0654,11.5912,12.8532,13.319,13.234,13.234,9.92246,12.114,12.2155,10.218,16.6501,4.42801,4.42801,3.73363,0.637033,0.629075,2.41865,2.24205,1.35329,1.57205},
    { 3.30428,3.95602,4.12211,3.47828,3.59562,3.5507,3.5507,6.24448,6.16396,5.84534,3.68661,11.8233,9.73669,9.73669,8.80798,8.05073,9.25977,12.7774,11.2528,8.26376,5.21201},
    { 9.94399,13.4021,16.4878,17.1825,18.5365,17.2436,17.2436,16.0188,15.4891,15.9359,13.9502,13.9502,11.159,11.159,8.88221,7.62407,8.13011,12.3709,12.2224,11.3578,11.3369} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.18667, 1.27736, 1.17733  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.26888,1.27956,1.29479,1.32825,1.38823,1.4599,1.4599,1.60026,1.58967,1.61678,1.61783,1.5238,1.82325,1.82325,2.8209,2.32651,2.12906,2.04043,1.7472,1.70848,1.72699},
    { 1.73523,1.71931,1.71998,1.80045,1.82682,1.84709,1.84709,1.75468,1.76759,1.7823,1.83586,1.47339,1.61279,1.61279,1.46184,1.51935,1.4813,1.52581,1.46656,1.27459,1.94158},
    { 1.3769,1.37028,1.36668,1.31614,1.26224,1.2748,1.2748,1.26801,1.27844,1.25445,1.26094,1.26094,1.61198,1.61198,1.89575,1.86637,2.01162,2.18551,2.09765,1.41947,1.44177} };
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

