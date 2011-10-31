Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][9] = { 
    { 1.13488,1.13488,1.13489,1.14483,1.1287,1.13078,1.13449,1.1233,1.12048},
    { 1.11624,1.11624,1.11622,1.14734,1.15525,1.15594,1.16672,1.1549,1.15361} };
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
    { 1.14369,1.14369,1.12445,1.14199,1.1287,1.13078,1.13449,1.05936,1.04982},
    { 1.0356,1.0356,1.05916,1.12746,1.15525,1.15594,1.16672,1.04827,1.09567} };
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
    { 1.09863,1.09863,1.09863,1.09729,1.10118,1.10114,1.10131,1.10181,1.10244,},
    { 1.321,1.321,1.32101,1.31511,1.31997,1.31991,1.31843,1.32045,1.32129} };
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
    { 1.09208,1.09208,1.09374,1.0943,1.10118,1.10114,1.10131,1.10946,1.11641},
    { 1.32736,1.32736,1.32136,1.31375,1.31997,1.31991,1.31843,1.33307,1.32279} };
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

