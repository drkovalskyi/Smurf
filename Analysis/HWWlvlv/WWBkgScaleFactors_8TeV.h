static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.13341,1.13341,1.13341,1.13341,1.13341,1.13341,1.13341,1.13341,1.13349,1.13644,1.13071,1.13071,1.09516,1.09516,1.09414,1.08986,1.08937,1.08651,1.08801},
    { 0.855389,0.855389,0.855389,0.855389,0.855389,0.855389,0.855389,0.855389,0.8554,0.865574,0.862348,0.862348,0.855631,0.855631,0.854303,0.855105,0.855781,0.853858,0.847089} };
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
    { 1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373,1.22373},
    { 1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661,1.06661} };
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
    { 1.06299,1.06299,1.06299,1.06299,1.06299,1.06299,1.06299,1.06299,1.06298,1.06227,1.06214,1.06214,1.06502,1.06502,1.06509,1.06546,1.06574,1.06625,1.06659,},
    { 1.14344,1.14344,1.14344,1.14344,1.14344,1.14344,1.14344,1.14344,1.14344,1.13985,1.13872,1.13872,1.14099,1.14099,1.14117,1.1402,1.13989,1.14012,1.14172} };
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
    { 1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062,1.062},
    { 1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826,1.12826} };
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

