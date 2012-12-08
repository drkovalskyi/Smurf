static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.08798,1.08798,1.08798,1.08798,1.08798,1.08798,1.08798,1.08797,1.08797,1.08843,1.09677,1.09677,1.09669,1.09669,1.09945,1.09988,1.09586,1.09712,1.09217},
    { 1.02305,1.02305,1.02305,1.02305,1.02305,1.02305,1.02305,1.02308,1.02303,1.03334,1.05166,1.05166,1.05656,1.05656,1.05802,1.06704,1.05852,1.05794,1.04177} };
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
    { 1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114,1.07114},
    { 0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944,0.952944} };
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
    { 1.07503,1.07503,1.07503,1.07503,1.07503,1.07503,1.07503,1.07503,1.07503,1.07494,1.07396,1.07396,1.07531,1.07531,1.07513,1.07528,1.07556,1.07589,1.07664,},
    { 1.15011,1.15011,1.15011,1.15011,1.15011,1.15011,1.15011,1.1501,1.15012,1.14854,1.14532,1.14532,1.15092,1.15092,1.15086,1.1497,1.1509,1.15144,1.15371} };
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
    { 1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053,1.06053},
    { 1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776,1.14776} };
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

