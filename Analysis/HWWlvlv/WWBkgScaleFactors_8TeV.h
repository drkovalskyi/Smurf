static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.09277,1.09277,1.09277,1.09277,1.09277,1.09277,1.09277,1.09277,1.09285,1.09585,1.09173,1.09173,1.07096,1.07096,1.07088,1.06649,1.06422,1.05985,1.0629},
    { 0.935024,0.935024,0.935024,0.935024,0.935024,0.935024,0.935024,0.935025,0.935035,0.93555,0.936115,0.936115,0.932982,0.932982,0.932604,0.928786,0.928526,0.923005,0.918141} };
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
    { 1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184,1.19184},
    { 1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865,1.08865} };
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
    { 1.05265,1.05265,1.05265,1.05265,1.05265,1.05265,1.05265,1.05265,1.05264,1.05219,1.05226,1.05226,1.05168,1.05168,1.05169,1.05198,1.05225,1.05268,1.05287,},
    { 1.10268,1.10268,1.10268,1.10268,1.10268,1.10268,1.10268,1.10268,1.10268,1.10158,1.10109,1.10109,1.10362,1.10362,1.10363,1.10385,1.104,1.10461,1.10538} };
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
    { 1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503,1.04503},
    { 1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905,1.0905} };
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

