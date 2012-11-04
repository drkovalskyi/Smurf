static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 2.12934, 2.75084, 4.60605  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 7.02554,8.11659,9.45938,10.1396,10.5251,11.0867,11.0867,8.58461,8.11135,7.81623,8.94725,8.94725,7.06394,7.06394,6.60297,4.60429,1.87527,0.222625,0.438654,0.26298,0.220132},
    { 4.35128,6.12406,8.15462,8.11789,8.86022,9.56809,9.56809,9.9262,8.05405,9.78748,11.0918,11.0918,10.6304,10.6304,9.02557,9.94949,10.3637,18.3234,18.485,7.35354,5.88319},
    { 0.0903793,0.179345,0.449573,0.55084,0.710005,1.6697,1.6697,1.60525,1.50803,1.50169,1.50169,1.50169,1.49585,1.49585,9.62175e-06,9.60607e-06,2.9647e-05,7.93415e-05,9.91889e-05,0.000128277,8.83538e-05} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.38169, 1.32731, 1.32687  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.40592,1.403,1.40239,1.40735,1.41347,1.42546,1.42546,1.46236,1.63657,1.64711,1.63246,1.63246,1.50247,1.50247,2.08956,2.19291,2.35095,4.03838,2.88628,2.65613,8.07333},
    { 1.36941,1.35967,1.35299,1.35611,1.35787,1.36171,1.36171,1.36631,1.38393,1.37517,1.37043,1.37043,1.38578,1.38578,1.44265,1.43395,1.41344,1.36831,1.35837,1.37546,1.84467},
    { 2.45515,2.6203,2.52299,2.52967,2.54019,2.37039,2.37039,2.40132,2.44972,2.4507,2.4507,2.4507,2.56619,2.56619,2.08189,2.08352,1.65677,1.46593,1.43778,1.41111,1.45311} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 2.12934, 2.75084, 4.60605  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,24.4922,8.9589},
    { 44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,44.1642,24.563},
    { 3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,3.13575,0.987835} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.38169, 1.32731, 1.32687  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.37759,1.44566},
    { 1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.32513,1.34205},
    { 2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.20525,2.45612} };
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

