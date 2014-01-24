static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.13846,1.13846,1.13846,1.13846,1.13846,1.13846,1.13846,1.13846,1.13855,1.1405,1.14811,1.14811,1.11186,1.11186,1.11053,1.1055,1.11212,1.11482,1.11054},
    { 0.891826,0.891826,0.891826,0.891826,0.891826,0.891826,0.891826,0.891826,0.891837,0.900526,0.90275,0.90275,0.922507,0.922507,0.917937,0.914456,0.917975,0.91507,0.906548} };
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
    { 1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186,1.19186},
    { 0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512,0.998512} };
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
    { 1.07445,1.07445,1.07445,1.07445,1.07445,1.07445,1.07445,1.07445,1.07444,1.074,1.07316,1.07316,1.07681,1.07681,1.07693,1.07735,1.07725,1.07747,1.07827,},
    { 1.1807,1.1807,1.1807,1.1807,1.1807,1.1807,1.1807,1.1807,1.1807,1.17839,1.17654,1.17654,1.17819,1.17819,1.17895,1.17908,1.17785,1.17843,1.18055} };
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
    { 1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545,1.06545},
    { 1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628,1.15628} };
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

