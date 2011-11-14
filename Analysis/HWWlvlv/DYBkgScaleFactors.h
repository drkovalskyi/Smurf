Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.05333, 3.77023, 8.32348  };
  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { 
    { 3.62135,3.68817,4.60192,5.11402,3.53062,2.38553,2.78923,6.95108,22.814,21.144,0.164476,0.180441},
    { 1.96367,3.40266,2.23395,2.54084,4.1612,3.66422,6.43097,7.56567,7.91583,7.56728,2.70019,4.55992},
    { 4.89691,5.82205,5.85139,5.80057,5.59914,6.04599,6.50291,9.15083,9.03001,8.17839,6.80375,7.0233} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.62216, 1.195, 1.16668  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { 
    { 1.7141,1.71604,1.67435,1.66728,1.97881,2.94665,3.22312,3.46736,2.53005,2.06235,5.00285,9.83972},
    { 1.44345,1.4331,1.51945,1.50779,1.36303,1.41852,1.42939,1.42123,1.39771,1.35249,1.32435,1.32602},
    { 1.29867,1.28697,1.3183,1.31468,1.32316,1.38512,1.34464,1.31891,1.28118,1.28792,1.55915,1.51588} };
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

