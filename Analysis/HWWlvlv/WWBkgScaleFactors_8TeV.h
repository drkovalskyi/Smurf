static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.08869,1.08869,1.08869,1.08869,1.08869,1.08869,1.08869,1.08869,1.08877,1.09389,1.092,1.092,1.07555,1.07555,1.07551,1.07131,1.06919,1.06503,1.06835},
    { 0.926326,0.926326,0.926326,0.926326,0.926326,0.926326,0.926326,0.926326,0.926337,0.931031,0.936225,0.936225,0.939683,0.939683,0.939251,0.935788,0.935659,0.930316,0.925624} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
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
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669,1.15669},
    { 1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232,1.04232} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
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
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][19] = { 
    { 1.05896,1.05896,1.05896,1.05896,1.05896,1.05896,1.05896,1.05896,1.05895,1.05817,1.05802,1.05802,1.05768,1.05768,1.0577,1.05805,1.05837,1.05885,1.05906,},
    { 1.10484,1.10484,1.10484,1.10484,1.10484,1.10484,1.10484,1.10484,1.10484,1.10251,1.10088,1.10088,1.10221,1.10221,1.10223,1.10236,1.10249,1.10305,1.1038} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
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
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][19] = { 
    { 1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524,1.0524},
    { 1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001,1.1001} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

