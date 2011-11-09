Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][9] = { 
    { 1.15822,1.15822,1.15823,1.16873,1.15083,1.154,1.15546,1.14633,1.14517},
    { 1.20737,1.20737,1.20735,1.24075,1.25382,1.25453,1.26504,1.25555,1.25578} };
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
    { 1.18165,1.18165,1.15754,1.17758,1.15083,1.154,1.15546,1.09281,1.09764},
    { 1.19459,1.19459,1.15421,1.24679,1.25382,1.25453,1.26504,1.18327,1.20902} };
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
    { 1.09709,1.09709,1.09709,1.0952,1.0988,1.09861,1.09876,1.09922,1.09971,},
    { 1.30117,1.30117,1.30118,1.29581,1.29991,1.29986,1.29871,1.3001,1.30067} };
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
    { 1.08974,1.08974,1.09185,1.09168,1.0988,1.09861,1.09876,1.10577,1.11101},
    { 1.29663,1.29663,1.29955,1.29129,1.29991,1.29986,1.29871,1.30366,1.29766} };
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

