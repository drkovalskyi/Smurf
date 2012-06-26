static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20583,1.20708,1.21882,1.14855,1.1068,1.1068,1.10243,1.09763,1.09571,1.09784,1.09297},
    { 0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933211,0.93652,0.949952,0.988005,0.971388,0.971388,0.965986,0.962093,0.970692,0.957319,0.940795} };
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
    { 1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.20573,1.14233,1.14233,1.14233,1.14233,1.14233,1.15576,1.15242,1.16896},
    { 0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.933199,0.973487,0.973487,0.973487,0.973487,0.973487,1.04839,1.03083,1.00412} };
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
    { 1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08124,1.0797,1.0802,1.08356,1.08356,1.08385,1.08434,1.08469,1.08499,1.08588,},
    { 1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17233,1.16802,1.17955,1.18855,1.18855,1.18961,1.19055,1.1878,1.19047,1.19422} };
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
    { 1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08164,1.08134,1.08134,1.08134,1.08134,1.08134,1.08304,1.08642,1.08917},
    { 1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.17325,1.18396,1.18396,1.18396,1.18396,1.18396,1.16785,1.18394,1.19919} };
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

