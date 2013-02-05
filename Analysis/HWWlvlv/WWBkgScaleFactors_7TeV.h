static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.08532,1.08532,1.08532,1.08532,1.08532,1.08532,1.08532,1.08531,1.08531,1.08565,1.09198,1.09198,1.08916,1.08916,1.0921,1.09255,1.08818,1.0894,1.08468},
    { 1.0418,1.0418,1.0418,1.0418,1.0418,1.0418,1.0418,1.04183,1.04178,1.05214,1.06757,1.06757,1.07084,1.07084,1.0723,1.08088,1.07241,1.07189,1.05576} };
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
    { 1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098,1.08098},
    { 0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883,0.978883} };
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
    { 1.07409,1.07409,1.07409,1.07409,1.07409,1.07409,1.07409,1.07409,1.07409,1.07403,1.07343,1.07343,1.07536,1.07536,1.07516,1.0753,1.07564,1.07598,1.07669,},
    { 1.14436,1.14436,1.14436,1.14436,1.14436,1.14436,1.14436,1.14435,1.14436,1.14294,1.14062,1.14062,1.14667,1.14667,1.14661,1.14561,1.14675,1.14727,1.14943} };
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
    { 1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565,1.0565},
    { 1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833,1.13833} };
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

