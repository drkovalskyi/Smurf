static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.09618,1.09917,1.09507,1.09507,1.07437,1.07437,1.0743,1.06997,1.06776,1.06346,1.06655},
    { 0.923524,0.923524,0.923524,0.923524,0.923524,0.923524,0.923524,0.923524,0.923535,0.924014,0.924474,0.924474,0.920466,0.920466,0.920087,0.916272,0.915991,0.910473,0.905592} };
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
    { 1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037,1.19037},
    { 1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785,1.0785} };
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
    { 1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.05889,1.05847,1.05858,1.05858,1.05862,1.05862,1.05864,1.059,1.05933,1.05983,1.06005,},
    { 1.1041,1.1041,1.1041,1.1041,1.1041,1.1041,1.1041,1.1041,1.1041,1.10301,1.10253,1.10253,1.10516,1.10516,1.10518,1.10541,1.10556,1.10619,1.10698} };
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
    { 1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864},
    { 1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144,1.09144} };
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

