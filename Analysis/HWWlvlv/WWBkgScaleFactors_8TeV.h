static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21332,1.21851,1.21558,1.09635,1.07368,1.07368,1.0742,1.0717,1.05218,1.05186},
    { 0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.912462,0.923497,0.918087,0.794405,0.794405,0.773776,0.778642,0.774324,0.755239} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
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
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.21318,1.1005,1.1005,1.1005,1.1005,1.1005,1.1491,1.14402},
    { 0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.935392,0.930273,0.930273,0.930273,0.930273,0.930273,0.939658,0.878026} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
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
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][18] = { 
    { 1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11825,1.11811,1.1189,1.12531,1.12531,1.12533,1.12591,1.12794,1.12843,},
    { 1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26948,1.26659,1.27957,1.33548,1.33548,1.34353,1.34249,1.34189,1.35035} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
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
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][18] = { 
    { 1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11823,1.11799,1.11799,1.11799,1.11799,1.11799,1.12072,1.12818},
    { 1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.26419,1.27592,1.27592,1.27592,1.27592,1.27592,1.27792,1.32142} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

