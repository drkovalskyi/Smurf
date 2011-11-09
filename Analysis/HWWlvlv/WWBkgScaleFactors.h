Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][9] = { 
    { 1.15835,1.15835,1.15836,1.16886,1.15096,1.15413,1.15559,1.14647,1.1453},
    { 1.20581,1.20581,1.20579,1.23919,1.25228,1.25299,1.2635,1.25401,1.25424} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 9 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][9] = { 
    { 1.18177,1.18177,1.15766,1.17771,1.15096,1.15413,1.15559,1.09295,1.09778},
    { 1.193,1.193,1.15257,1.24519,1.25228,1.25299,1.2635,1.18162,1.20722} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 9 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][9] = { 
    { 1.07949,1.07949,1.07949,1.07824,1.08213,1.08195,1.08207,1.08269,1.08311,},
    { 1.13279,1.13279,1.13281,1.12826,1.13181,1.13178,1.13068,1.13162,1.13189} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 9 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][9] = { 
    { 1.07048,1.07048,1.074,1.07474,1.08213,1.08195,1.08207,1.08974,1.09515},
    { 1.12169,1.12169,1.12921,1.12312,1.13181,1.13178,1.13068,1.14292,1.14766} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 9 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

