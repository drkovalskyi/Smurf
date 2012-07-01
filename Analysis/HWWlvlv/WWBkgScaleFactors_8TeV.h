static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20721,1.20847,1.22022,1.14966,1.10693,1.10693,1.10256,1.09776,1.09584,1.09797,1.0931},
    { 0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932208,0.935514,0.948936,0.986752,0.970156,0.970156,0.964753,0.96086,0.969458,0.956082,0.939551} };
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
    { 1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.20711,1.14341,1.14341,1.14341,1.14341,1.14341,1.15591,1.15261,1.16918},
    { 0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.932196,0.97225,0.97225,0.97225,0.97225,0.97225,1.04689,1.02926,1.00214} };
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
    { 1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08112,1.07958,1.0801,1.08355,1.08355,1.08383,1.08433,1.08467,1.08497,1.08587,},
    { 1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17253,1.16822,1.17979,1.1888,1.1888,1.18987,1.19081,1.18806,1.19073,1.1945} };
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
    { 1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08152,1.08123,1.08123,1.08123,1.08123,1.08123,1.08303,1.0864,1.08915},
    { 1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.17345,1.18421,1.18421,1.18421,1.18421,1.18421,1.16811,1.18423,1.1996} };
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

