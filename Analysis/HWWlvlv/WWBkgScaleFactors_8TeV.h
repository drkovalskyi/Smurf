static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.09618,1.09917,1.09507,1.09507,1.07438,1.07438,1.0743,1.06997,1.06776,1.06346,1.06656},
    { 0.923563,0.923563,0.923563,0.923563,0.923563,0.923563,0.923563,0.923564,0.923574,0.924055,0.924513,0.924513,0.920497,0.920497,0.920119,0.916305,0.916024,0.910507,0.905628} };
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
    { 1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038},
    { 1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845,1.07845} };
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
    { 1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.05889,1.05847,1.05859,1.05859,1.05862,1.05862,1.05864,1.059,1.05933,1.05983,1.06005,},
    { 1.10407,1.10407,1.10407,1.10407,1.10407,1.10407,1.10407,1.10407,1.10407,1.10298,1.1025,1.1025,1.10514,1.10514,1.10515,1.10538,1.10554,1.10616,1.10696} };
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
    { 1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143} };
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

