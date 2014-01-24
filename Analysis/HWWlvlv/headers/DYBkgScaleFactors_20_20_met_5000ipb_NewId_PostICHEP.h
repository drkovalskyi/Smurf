static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.52092, 3.66118, 1.9219  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 11.856,13.0655,13.2552,15.0269,15.9684,18.9024,18.9024,19.2378,19.4643,20.7095,20.6574,30.6428,13.8741,13.8741,9.82456,3.33938,1.56898,5.71158,7.41711,3.86353,3.47249},
    { 2.18955,3.04971,3.44186,2.87797,2.693,1.59214,1.59214,2.73409,1.8477,2.20915,1.65321,20.097,23.6927,23.6927,26.5715,23.3778,23.21,31.6519,29.6807,20.8509,15.312},
    { 1.73663,1.73593,1.73561,1.73479,1.73433,1.89207,2.21009,2.20945,2.36855,2.36212,2.52125,2.5204,2.35722,2.67133,2.828,2.98499,3.13926,3.44558,3.58937,2.46545,2.46103} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 21 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

static Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.29474, 1.15754, 1.22813  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.21243,1.21096,1.21371,1.27594,1.32523,1.38822,1.38822,1.39105,1.4493,1.44502,1.44532,1.61144,1.80729,1.80729,1.87229,1.64642,4.96014,2.32559,2.0224,1.71049,1.58394},
    { 2.10427,2.08432,2.08013,2.10472,2.1264,2.25221,2.25221,1.6939,1.95968,1.87333,2.12523,1.70881,1.5902,1.5902,1.40624,1.44069,1.40158,1.40938,1.34699,1.25287,1.7237},
    { 2.13422,2.13425,2.13427,2.1343,2.13433,2.13093,2.12547,2.12549,2.1233,2.12346,2.12153,2.12155,2.12772,2.12644,2.12427,2.12234,1.60884,1.60328,1.60063,2.18068,2.18076} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 21 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

