Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[9] = {115,120,130,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][9] = { 
    { 1.13567,1.13567,1.13568,1.14562,1.12951,1.13159,1.1353,1.12411,1.1213},
    { 1.10537,1.10537,1.10535,1.13645,1.14413,1.14482,1.15558,1.14375,1.14244} };
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
    { 1.14444,1.14444,1.12523,1.14277,1.12951,1.13159,1.1353,1.06022,1.05071},
    { 1.02487,1.02487,1.04873,1.11676,1.14413,1.14482,1.15558,1.03767,1.0859} };
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
    { 1.09845,1.09845,1.09845,1.09711,1.10099,1.10095,1.10112,1.10162,1.10224,},
    { 1.3255,1.3255,1.32551,1.3195,1.32452,1.32445,1.32294,1.32502,1.32588} };
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
    { 1.09192,1.09192,1.09356,1.09412,1.10099,1.10095,1.10112,1.10924,1.11618},
    { 1.3321,1.3321,1.32576,1.31805,1.32452,1.32445,1.32294,1.33775,1.32672} };
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

