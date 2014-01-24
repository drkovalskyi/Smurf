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
    { 1.06073,1.06073,1.06073,1.06073,1.06073,1.06073,1.06073,1.06073,1.06072,1.05999,1.05979,1.05979,1.06216,1.06216,1.06222,1.06256,1.06282,1.06329,1.06359,},
    { 1.14633,1.14633,1.14633,1.14633,1.14633,1.14633,1.14633,1.14633,1.14633,1.1429,1.14185,1.14185,1.1429,1.1429,1.14306,1.14216,1.14186,1.14208,1.14359} };
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
    { 1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771},
    { 1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879,1.12879} };
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

