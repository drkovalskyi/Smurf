Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.27351, 4.08024, 7.94197  };
  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { 
    { 10.3768,10.369,3.92054,3.05213,5.28014,1.00936,0.528444,0.0660013,30.2001,1.84261,0.150859,0.270647},
    { 2.12132,3.6175,2.25792,2.5705,3.89025,3.87764,5.28484,3.15614,5.03853,5.76586,2.84815,4.29011},
    { 3.30322,2.77334,2.63991,2.951,7.64305,8.74589,8.14114,8.77525,8.89862,8.82666,11.4185,10.7071} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.68136, 1.41036, 1.42692  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { 
    { 1.57883,1.57368,2.00993,2.20637,2.34774,2.81329,4.7402,10.3028,2.81371,10.2855,12.6231,10.0889},
    { 2.28799,2.28291,2.63335,2.62375,1.97502,1.79414,1.88022,2.39603,2.09786,1.89222,1.50605,1.58147},
    { 1.91925,1.84575,1.99833,1.98621,1.98113,1.77821,2.06667,2.17925,2.0841,1.94222,1.42452,1.40424} };
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

