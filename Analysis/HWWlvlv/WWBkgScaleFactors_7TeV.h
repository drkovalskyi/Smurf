static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10848,1.1085,1.11114,1.11846,1.10171,1.10476,1.10616,1.09746,1.09636},
    { 1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10143,1.10147,1.1134,1.13282,1.14253,1.14318,1.15305,1.14384,1.14383} };
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

static Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.10849,1.09168,1.1163},
    { 1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.10149,1.05089,1.07181} };
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

static Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07851,1.07763,1.08148,1.08129,1.08142,1.08203,1.08245,},
    { 1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.1418,1.1418,1.13996,1.13698,1.14138,1.14135,1.14017,1.14126,1.1416} };
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

static Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.07886,1.0834,1.08697},
    { 1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.14179,1.15087,1.15453} };
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

