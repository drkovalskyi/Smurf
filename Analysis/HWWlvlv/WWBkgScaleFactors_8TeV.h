static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.12599,1.12599,1.12599,1.12599,1.12599,1.12599,1.12599,1.12599,1.12607,1.12929,1.12448,1.12448,1.09747,1.09747,1.09656,1.09259,1.09234,1.08977,1.09165},
    { 0.867757,0.867757,0.867757,0.867757,0.867757,0.867757,0.867757,0.867757,0.867768,0.878041,0.875875,0.875875,0.881663,0.881663,0.880409,0.881397,0.88224,0.880543,0.874115} };
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
    { 1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766,1.18766},
    { 1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779,1.03779} };
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
    { 1.06163,1.06163,1.06163,1.06163,1.06163,1.06163,1.06163,1.06163,1.06162,1.06091,1.06073,1.06073,1.06331,1.06331,1.06337,1.06373,1.06399,1.06448,1.06479,},
    { 1.13897,1.13897,1.13897,1.13897,1.13897,1.13897,1.13897,1.13897,1.13896,1.13552,1.13434,1.13434,1.13547,1.13547,1.13563,1.13468,1.13437,1.13455,1.136} };
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
    { 1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153,1.06153},
    { 1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281,1.1281} };
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

