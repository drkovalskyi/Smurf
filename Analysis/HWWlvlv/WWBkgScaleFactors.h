Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15748,1.16026,1.16797,1.15035,1.15352,1.15498,1.14585,1.14468},
    { 1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19942,1.19945,1.21211,1.23281,1.24654,1.24724,1.25774,1.24823,1.24844} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.15747,1.13758,1.16304},
    { 1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.19947,1.13291,1.13526} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07921,1.07831,1.08218,1.082,1.08213,1.08275,1.08317,},
    { 1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13366,1.13366,1.13193,1.12908,1.13254,1.13251,1.13141,1.13236,1.13264} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.07957,1.08439,1.08798},
    { 1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.13365,1.14454,1.15185} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

