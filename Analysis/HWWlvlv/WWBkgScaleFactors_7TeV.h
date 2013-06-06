static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07698,1.08327,1.08327,1.08106,1.08106,1.08416,1.08457,1.08015,1.0813,1.0765},
    { 1.0146,1.0146,1.0146,1.0146,1.0146,1.0146,1.0146,1.01463,1.0146,1.02496,1.04041,1.04041,1.04333,1.04333,1.04477,1.05329,1.04477,1.04415,1.02799} };
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
    { 0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793,0.948793} };
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
    { 1.14948,1.14948,1.14948,1.14948,1.14948,1.14948,1.14948,1.14948,1.14949,1.14797,1.14553,1.14553,1.15187,1.15187,1.15181,1.15076,1.15198,1.15254,1.15485} };
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

