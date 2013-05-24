static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07698,1.08327,1.08327,1.08105,1.08105,1.08416,1.08456,1.08015,1.08129,1.0765},
    { 1.01452,1.01452,1.01452,1.01452,1.01452,1.01452,1.01452,1.01455,1.01452,1.02488,1.04034,1.04034,1.04327,1.04327,1.04472,1.05323,1.04471,1.0441,1.02793} };
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
    { 1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196,1.07196},
    { 0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726,0.948726} };
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
    { 1.07515,1.07515,1.07515,1.07515,1.07515,1.07515,1.07515,1.07515,1.07515,1.07508,1.07447,1.07447,1.07632,1.07632,1.0761,1.07625,1.0766,1.07694,1.07768,},
    { 1.14949,1.14949,1.14949,1.14949,1.14949,1.14949,1.14949,1.14949,1.14949,1.14798,1.14554,1.14554,1.15188,1.15188,1.15181,1.15076,1.15199,1.15254,1.15486} };
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
    { 1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749,1.05749},
    { 1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448,1.1448} };
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

