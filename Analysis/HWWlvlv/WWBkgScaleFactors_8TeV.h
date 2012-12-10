static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07266,1.07266,1.07266,1.07266,1.07266,1.07266,1.07266,1.07266,1.07275,1.07733,1.0763,1.0763,1.05774,1.05774,1.05765,1.05303,1.05059,1.04656,1.04941},
    { 0.881697,0.881697,0.881697,0.881697,0.881697,0.881697,0.881697,0.881697,0.881707,0.888147,0.894656,0.894656,0.888421,0.888421,0.887976,0.887891,0.887103,0.881669,0.878226} };
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
    { 1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535,1.15535},
    { 1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289,1.02289} };
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
    { 1.06238,1.06238,1.06238,1.06238,1.06238,1.06238,1.06238,1.06238,1.06237,1.06153,1.06125,1.06125,1.06093,1.06093,1.06095,1.06134,1.06169,1.06219,1.06242,},
    { 1.12161,1.12161,1.12161,1.12161,1.12161,1.12161,1.12161,1.12161,1.1216,1.11851,1.11638,1.11638,1.11895,1.11895,1.11897,1.1187,1.11897,1.11965,1.12042} };
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
    { 1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525,1.05525},
    { 1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109,1.1109} };
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

