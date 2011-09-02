Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.02398, 2.81254, 5.08747  };
  Double_t DYBkgScaleFactorHiggsSelection[3][12] = { 
    { 5.06356,4.6418,3.06347,3.40455,1.56604,2.27749,3.82277,6.44,2.40173,1.38875,0.481823,2.37679},
    { 1.07078,1.84613,1.12533,1.04673,2.38312,2.28017,4.20695,1.96313,3.98251,4.85908,3.6442,3.59311},
    { 2.31884,3.3557,1.85614,3.09849,3.19556,2.59791,4.394,4.27537,5.16934,5.05894,5.28543,5.07749} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin][massIndex];
  }

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.6743, 1.55228, 1.46547  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][12] = { 
    { 3.91873,3.93674,3.33284,2.98309,6.20489,3.00385,3.75876,3.29036,4.19649,7.12789,8.10841,6.12231},
    { 3.29404,3.28525,4.86545,3.79078,2.94535,2.38163,2.19452,2.50524,2.70584,2.55319,1.56367,2.41627},
    { 3.17044,3.16377,3.94198,3.02316,5.22859,6.04353,5.98962,5.12945,3.04452,2.05939,1.30263,1.35387} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin][massIndex];
  }

