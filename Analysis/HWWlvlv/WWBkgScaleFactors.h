Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][17] = { 
    { 1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37274,1.37796,1.36554,1.23495,1.15342,1.15399,1.15962,1.13284,1.12437},
    { 1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.23505,1.25054,1.20733,1.02627,1.02703,1.03477,1.02038,1.02353} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
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
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][17] = { 
    { 1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.37258,1.28811,1.24202,1.24202,1.24202,1.28863,1.24567},
    { 1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25734,1.25305,1.2153,1.2153,1.2153,1.19103,1.20752} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
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
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][17] = { 
    { 1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.1452,1.14602,1.14196,1.15468,1.1547,1.15468,1.15752,1.15891,},
    { 1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26876,1.26829,1.28896,1.35233,1.35232,1.35081,1.35167,1.35163} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
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
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][17] = { 
    { 1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.14516,1.18544,1.14029,1.14029,1.14029,1.14415,1.15783},
    { 1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.26503,1.36133,1.28469,1.28469,1.28469,1.29272,1.31696} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

