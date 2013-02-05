static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.09874,1.09874,1.09874,1.09874,1.09874,1.09874,1.09874,1.09874,1.09883,1.10182,1.09766,1.09766,1.07671,1.07671,1.07664,1.07208,1.06988,1.06559,1.06853},
    { 0.930978,0.930978,0.930978,0.930978,0.930978,0.930978,0.930978,0.930979,0.930989,0.931498,0.931984,0.931984,0.927381,0.927381,0.927005,0.923196,0.922916,0.917413,0.91255} };
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
    { 1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423,1.19423},
    { 1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779,1.08779} };
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
    { 1.05863,1.05863,1.05863,1.05863,1.05863,1.05863,1.05863,1.05863,1.05862,1.05821,1.05832,1.05832,1.05839,1.05839,1.05841,1.05879,1.05912,1.05961,1.05985,},
    { 1.10269,1.10269,1.10269,1.10269,1.10269,1.10269,1.10269,1.10269,1.10269,1.1016,1.10112,1.10112,1.10384,1.10384,1.10386,1.10408,1.10423,1.10485,1.10563} };
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
    { 1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829,1.04829},
    { 1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977,1.08977} };
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

