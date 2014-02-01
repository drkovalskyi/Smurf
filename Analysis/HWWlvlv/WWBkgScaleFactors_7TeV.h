static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.06853,1.06853,1.06853,1.06853,1.06853,1.06853,1.06853,1.06852,1.06852,1.06899,1.07525,1.07525,1.07278,1.07278,1.07587,1.07625,1.0718,1.07288,1.06801},
    { 1.00652,1.00652,1.00652,1.00652,1.00652,1.00652,1.00652,1.00655,1.00652,1.01689,1.03234,1.03234,1.03506,1.03506,1.0365,1.045,1.03646,1.03582,1.01964} };
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
    { 1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618,1.06618},
    { 0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654,0.940654} };
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
    { 1.07293,1.07293,1.07293,1.07293,1.07293,1.07293,1.07293,1.07293,1.07293,1.07286,1.07224,1.07224,1.07395,1.07395,1.07373,1.07387,1.07419,1.07452,1.07521,},
    { 1.15083,1.15083,1.15083,1.15083,1.15083,1.15083,1.15083,1.15082,1.15083,1.14929,1.14681,1.14681,1.15327,1.15327,1.1532,1.15213,1.15338,1.15395,1.15631} };
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
    { 1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599,1.05599},
    { 1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632,1.14632} };
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

