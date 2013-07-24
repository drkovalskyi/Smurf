static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07652,1.07698,1.08327,1.08327,1.08106,1.08106,1.08416,1.08457,1.08015,1.0813,1.0765},
    { 1.01464,1.01464,1.01464,1.01464,1.01464,1.01464,1.01464,1.01467,1.01463,1.02499,1.04045,1.04045,1.04337,1.04337,1.04481,1.05333,1.04481,1.04419,1.02802} };
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
    { 0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959,0.948959} };
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
    { 1.14947,1.14947,1.14947,1.14947,1.14947,1.14947,1.14947,1.14947,1.14948,1.14796,1.14552,1.14552,1.15186,1.15186,1.1518,1.15075,1.15197,1.15253,1.15484} };
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
    { 1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474,1.14474} };
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

