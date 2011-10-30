Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.38833, 4.22787, 6.90608  };
  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { 
    { 10.3768,10.369,3.92054,3.45723,5.79417,1.00936,0.528444,0.0660013,35.3466,1.84261,0.150859,0.270647},
    { 2.12132,3.6175,2.25792,2.62923,3.95684,4.05744,5.55653,3.35359,5.24604,5.92567,2.99853,4.56334},
    { 5.17076,3.41434,2.73171,3.36498,9.36097,9.12765,11.0804,10.9941,8.84725,8.83281,7.28955,6.75973} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.67325, 1.40929, 1.24016  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { 
    { 1.57883,1.57368,2.00993,2.11315,2.3133,2.81329,4.7402,10.3028,2.68382,10.2278,12.6231,10.0889},
    { 2.28799,2.28291,2.63335,2.6232,1.97446,1.78875,1.87406,2.38816,2.09624,1.8913,1.50286,1.5778},
    { 1.62773,1.70795,1.86177,1.84678,2.02728,1.95265,2.22443,1.99373,1.64431,1.58545,1.28861,1.23734} };
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

