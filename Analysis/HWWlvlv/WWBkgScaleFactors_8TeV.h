static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.09869,1.09869,1.09869,1.09869,1.09869,1.09869,1.09869,1.09869,1.09877,1.10177,1.09761,1.09761,1.07667,1.07667,1.0766,1.07204,1.06984,1.06555,1.06849},
    { 0.9306,0.9306,0.9306,0.9306,0.9306,0.9306,0.9306,0.9306,0.930611,0.931113,0.931595,0.931595,0.926992,0.926992,0.926615,0.922805,0.922528,0.917023,0.91216} };
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
    { 1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407,1.19407},
    { 1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443,1.08443} };
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
    { 1.05864,1.05864,1.05864,1.05864,1.05864,1.05864,1.05864,1.05864,1.05863,1.05821,1.05833,1.05833,1.05839,1.05839,1.05841,1.05879,1.05912,1.05961,1.05986,},
    { 1.10278,1.10278,1.10278,1.10278,1.10278,1.10278,1.10278,1.10278,1.10278,1.10169,1.10121,1.10121,1.10393,1.10393,1.10395,1.10417,1.10432,1.10494,1.10572} };
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
    { 1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483,1.0483},
    { 1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037,1.09037} };
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

