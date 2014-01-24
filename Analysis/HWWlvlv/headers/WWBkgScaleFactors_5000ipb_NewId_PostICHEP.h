static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18066,1.18153,1.19124,1.11921,1.06989,1.06989,1.06542,1.06056,1.0611,1.06373,1.05891},
    { 0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.90341,0.908656,0.910208,0.950039,0.929374,0.929374,0.923913,0.91985,0.928161,0.921467,0.904009} };
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
    { 1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.18052,1.1151,1.1151,1.1151,1.1151,1.1151,1.1253,1.12593,1.14293},
    { 0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.903411,0.945014,0.945014,0.945014,0.945014,0.945014,1.02566,0.977571,0.959391} };
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
    { 1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08188,1.08067,1.08212,1.0863,1.0863,1.0866,1.08721,1.08741,1.08771,1.08868,},
    { 1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17448,1.17366,1.18601,1.19741,1.19741,1.19862,1.1997,1.19701,1.19862,1.20276} };
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
    { 1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08222,1.08285,1.08285,1.08285,1.08285,1.08285,1.08499,1.08834,1.09119},
    { 1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.17589,1.18763,1.18763,1.18763,1.18763,1.18763,1.16966,1.19031,1.20469} };
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

