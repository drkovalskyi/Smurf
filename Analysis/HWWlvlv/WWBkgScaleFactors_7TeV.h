static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07625,1.07625,1.07625,1.07625,1.07625,1.07625,1.07625,1.07625,1.07625,1.07671,1.083,1.083,1.08083,1.08083,1.08393,1.08434,1.07992,1.08106,1.07627},
    { 1.01048,1.01048,1.01048,1.01048,1.01048,1.01048,1.01048,1.01051,1.01048,1.02085,1.03631,1.03631,1.03918,1.03918,1.04062,1.04913,1.0406,1.03997,1.0238} };
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
    { 1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164,1.07164},
    { 0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979,0.944979} };
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
    { 1.07519,1.07519,1.07519,1.07519,1.07519,1.07519,1.07519,1.07519,1.07519,1.07512,1.07451,1.07451,1.07635,1.07635,1.07613,1.07628,1.07663,1.07698,1.07771,},
    { 1.15012,1.15012,1.15012,1.15012,1.15012,1.15012,1.15012,1.15012,1.15012,1.1486,1.14613,1.14613,1.15251,1.15251,1.15244,1.15138,1.15262,1.15318,1.15551} };
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
    { 1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754,1.05754},
    { 1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538,1.14538} };
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

