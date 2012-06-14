static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24068,1.23651,1.25029,1.18914,1.17201,1.17201,1.16936,1.16789,1.16309,1.16416},
    { 0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.948785,0.944046,0.970789,0.88353,0.88353,0.875504,0.866238,0.871906,0.848302} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

static Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.24055,1.18554,1.18554,1.18554,1.18554,1.18554,1.19982,1.18667},
    { 0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.951789,0.984764,0.984764,0.984764,0.984764,0.984764,1.01872,0.923297} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

static Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][18] = { 
    { 1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08764,1.08587,1.08415,1.08659,1.08659,1.08677,1.08711,1.0876,1.08787,},
    { 1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20947,1.21,1.21762,1.24976,1.24976,1.25202,1.25488,1.25205,1.2591} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

static Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][18] = { 
    { 1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08745,1.08479,1.08479,1.08479,1.08479,1.08479,1.0874,1.09316},
    { 1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.20855,1.21441,1.21441,1.21441,1.21441,1.21441,1.20896,1.24992} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

