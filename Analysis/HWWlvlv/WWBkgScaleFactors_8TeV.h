static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.09622,1.09622,1.09622,1.09622,1.09622,1.09622,1.09622,1.09622,1.09631,1.09929,1.09519,1.09519,1.07449,1.07449,1.07442,1.07009,1.06788,1.06358,1.06668},
    { 0.927169,0.927169,0.927169,0.927169,0.927169,0.927169,0.927169,0.92717,0.92718,0.927669,0.92811,0.92811,0.924179,0.924179,0.923802,0.919986,0.919702,0.914186,0.909309} };
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
    { 1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047,1.19047},
    { 1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809,1.0809} };
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
    { 1.05888,1.05888,1.05888,1.05888,1.05888,1.05888,1.05888,1.05888,1.05888,1.05846,1.05857,1.05857,1.05861,1.05861,1.05863,1.05899,1.05932,1.05982,1.06004,},
    { 1.10342,1.10342,1.10342,1.10342,1.10342,1.10342,1.10342,1.10342,1.10341,1.10232,1.10185,1.10185,1.10446,1.10446,1.10448,1.1047,1.10486,1.10548,1.10627} };
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
    { 1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863,1.04863},
    { 1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097,1.09097} };
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

