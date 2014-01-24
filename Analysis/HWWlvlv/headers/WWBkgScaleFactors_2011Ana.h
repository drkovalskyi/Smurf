Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][17] = { 
    { 1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15181,1.15635,1.15858,1.15858,1.05041,1.05113,1.05574,1.02532,1.03191},
    { 0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.761443,0.770051,0.770051,0.627558,0.628023,0.634624,0.592974,0.595757} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
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
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][17] = { 
    { 1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.15169,1.12053,1.15169,1.15169,1.15169,1.15596,1.10043},
    { 0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.784336,0.863107,0.784336,0.784336,0.784336,0.755881,0.622617} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
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
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][17] = { 
    { 1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15896,1.15893,1.15893,1.17727,1.17729,1.17728,1.18144,1.18144,},
    { 1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.45662,1.45611,1.45611,1.57953,1.57952,1.57516,1.6115,1.61148} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
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
  Int_t mHiggs[17] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][17] = { 
    { 1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.15889,1.21917,1.15889,1.15889,1.15889,1.16807,1.18422},
    { 1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.44542,1.52619,1.44542,1.44542,1.44542,1.45958,1.60842} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 17 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

