static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.05311, 3.84797, 2.18044  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 22.2751,30.1587,40.1426,42.2316,42.1997,46.5676,46.5676,46.5914,50.8588,48.7668,51.1144,51.1144,37.8227,37.8227,18.7232,7.90791,3.75931,21.9913,20.3873,11.2283,11.3839},
    { 3.40483,5.5976,7.92928,7.57119,7.44583,8.63733,8.63733,8.2066,8.49613,9.64566,10.1553,10.1553,11.7813,11.7813,8.80848,6.49507,7.52336,15.9668,16.8272,12.2098,8.07201},
    { 2.10888,2.27751,2.80275,3.51506,3.589,3.12521,3.12521,3.04691,2.85063,3.16642,3.1603,3.1603,2.04501,2.04501,3.29151,5.25066,6.23699,6.1551,5.08157,6.96001,4.0604} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31548, 1.31969, 1.3022  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.32123,1.31732,1.31368,1.31563,1.3199,1.32305,1.32305,1.32921,1.3321,1.33728,1.33619,1.33619,1.34208,1.34208,1.5122,1.89123,2.53281,1.39023,1.38163,1.35942,1.43755},
    { 1.42841,1.37442,1.35442,1.36508,1.37491,1.36909,1.36909,1.38357,1.39764,1.38722,1.38539,1.38539,1.36385,1.36385,1.44893,1.5149,1.46758,1.35531,1.34385,1.34571,1.39036},
    { 1.53082,1.52514,1.51095,1.53187,1.48772,1.48028,1.48028,1.49693,1.5191,1.50816,1.50855,1.50855,1.57324,1.57324,1.8616,1.98306,1.59159,1.51021,1.47453,1.50288,1.60155} };
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

