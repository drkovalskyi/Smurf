static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.2071,1.20281,1.21643,1.15223,1.13461,1.13461,1.13193,1.13029,1.12542,1.12641},
    { 1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.00851,1.00424,1.03103,0.945868,0.945868,0.937885,0.928625,0.934293,0.910861} };
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
    { 1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.20697,1.14889,1.14889,1.14889,1.14889,1.14889,1.16149,1.14701},
    { 1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.01129,1.04438,1.04438,1.04438,1.04438,1.04438,1.07993,0.987097} };
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
    { 1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.09193,1.09012,1.0892,1.09189,1.09189,1.09209,1.09248,1.09302,1.09331,},
    { 1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.1914,1.1916,1.19968,1.22736,1.22736,1.22925,1.23169,1.22916,1.23506} };
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
    { 1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.0917,1.08981,1.08981,1.08981,1.08981,1.08981,1.09272,1.09894},
    { 1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19066,1.19712,1.19712,1.19712,1.19712,1.19712,1.1919,1.22783} };
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

