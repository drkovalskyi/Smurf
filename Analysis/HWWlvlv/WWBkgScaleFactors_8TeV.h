static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.09514,1.09514,1.09514,1.09514,1.09514,1.09514,1.09514,1.09514,1.09522,1.09821,1.09409,1.09409,1.07336,1.07336,1.07328,1.06895,1.06674,1.06243,1.0655},
    { 0.929027,0.929027,0.929027,0.929027,0.929027,0.929027,0.929027,0.929027,0.929037,0.929561,0.930164,0.930164,0.926821,0.926821,0.926393,0.922629,0.922407,0.916872,0.911937} };
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
    { 1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172,1.19172},
    { 1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833,1.07833} };
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
    { 1.05899,1.05899,1.05899,1.05899,1.05899,1.05899,1.05899,1.05899,1.05898,1.05856,1.05867,1.05867,1.05871,1.05871,1.05872,1.05909,1.05941,1.05991,1.06014,},
    { 1.10359,1.10359,1.10359,1.10359,1.10359,1.10359,1.10359,1.10359,1.10359,1.10251,1.10203,1.10203,1.10459,1.10459,1.10461,1.10483,1.10498,1.1056,1.1064} };
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
    { 1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866,1.04866},
    { 1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162,1.09162} };
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

