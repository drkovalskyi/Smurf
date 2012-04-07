Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12656,1.13385,1.12383,1.12692,1.12245,1.1134,1.11334},
    { 1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18537,1.18533,1.19354,1.20846,1.21743,1.21967,1.21038,1.19863,1.19779} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
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
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.12092,1.10437,1.10804},
    { 1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.18534,1.16807,1.2019} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
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
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07585,1.07533,1.07923,1.07908,1.07948,1.08028,1.08073,},
    { 1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12469,1.12288,1.12667,1.12662,1.12756,1.12862,1.12902} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
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
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.07629,1.08146,1.08631},
    { 1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.12539,1.13177,1.13746} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

