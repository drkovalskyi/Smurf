Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.05333, 3.77023, 8.77368  };
  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { 
    { 3.62135,3.68817,4.60192,5.11402,3.53062,2.38553,2.78923,6.95108,22.814,21.144,0.164476,0.180441},
    { 1.96367,3.40266,2.23395,2.54084,4.1612,3.66422,6.43097,7.56567,7.91583,7.56728,2.70019,4.55992},
    { 7.33449,10.9111,8.93152,7.99434,5.85105,6.70623,7.58329,8.41707,7.79398,7.8258,6.18506,7.60171} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.62216, 1.195, 1.20837  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { 
    { 1.7141,1.71604,1.67435,1.66728,1.97881,2.94665,3.22312,3.46736,2.53005,2.06235,5.00285,9.83972},
    { 1.44345,1.4331,1.51945,1.50779,1.36303,1.41852,1.42939,1.42123,1.39771,1.35249,1.32435,1.32602},
    { 1.31339,1.3016,1.34244,1.33657,1.46566,1.47071,1.49686,1.44707,1.37993,1.36014,1.75929,1.46587} };
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

