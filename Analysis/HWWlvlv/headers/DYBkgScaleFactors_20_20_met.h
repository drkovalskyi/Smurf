static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.10987, 3.84051, 2.17943  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 22.677,30.7523,40.7621,42.7698,42.8563,47.0214,47.0214,47.1316,51.5343,49.4254,51.7796,51.7796,37.9963,37.9963,18.9292,8.02787,4.4439,22.4183,20.6519,11.3286,11.4069},
    { 3.42373,5.62811,7.97237,7.60769,7.48024,8.66954,8.66954,8.23447,8.51834,9.67457,10.1855,10.1855,11.8286,11.8286,8.83761,6.51815,7.54846,15.7577,16.6525,12.1437,7.99441},
    { 2.10698,2.27561,2.80085,3.51295,3.58658,3.12261,3.12261,3.04691,2.85063,3.16642,3.1603,3.1603,2.04501,2.04501,3.29151,5.25066,6.237,6.15511,5.07983,6.9586,4.0585} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.3153, 1.31987, 1.30221  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.32084,1.31697,1.31347,1.31544,1.31957,1.32283,1.32283,1.32886,1.33167,1.33674,1.3357,1.3357,1.34206,1.34206,1.51169,1.88454,2.32592,1.38838,1.3805,1.35952,1.43841},
    { 1.42798,1.37416,1.35422,1.36488,1.37472,1.36892,1.36892,1.38343,1.39761,1.38713,1.3853,1.3853,1.36367,1.36367,1.44839,1.5142,1.46707,1.35636,1.34449,1.34628,1.3923},
    { 1.53097,1.52526,1.51104,1.53193,1.48781,1.48042,1.48042,1.49693,1.5191,1.50816,1.50855,1.50855,1.57324,1.57324,1.8616,1.98306,1.59159,1.5102,1.47456,1.50289,1.60157} };
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

