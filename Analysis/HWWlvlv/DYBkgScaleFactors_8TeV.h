static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.60365, 3.95191, 1.99512  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 24.9845,27.5802,28.1264,31.8166,34.2262,35.4147,35.4147,35.2563,38.0408,36.4729,31.9845,45.6018,12.5816,12.5816,9.32072,2.9575,1.76779,5.11761,6.46522,3.68979,3.9997},
    { 6.28581,8.88414,9.92761,8.50577,9.22781,7.28474,7.28474,9.38548,7.54937,8.61633,5.72657,30.8128,24.6981,24.6981,26.794,23.6775,23.3527,32.1944,29.9745,20.3736,14.3545},
    { 2.61503,2.85417,2.8524,2.85107,2.8502,3.08834,3.56793,3.56695,3.80403,3.79356,4.27367,4.26956,4.26157,4.49418,4.72838,4.96516,5.19363,5.65968,6.12308,5.24332,5.23514} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14042, 1.20401, 1.26796  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.2129,1.21099,1.21316,1.24113,1.29474,1.34826,1.34826,1.35707,1.40489,1.46235,1.46831,1.54259,1.65752,1.65752,2.39711,2.79064,4.66727,1.93668,1.63435,1.84455,1.5276},
    { 1.61568,1.59849,1.59511,1.62474,1.64091,1.66032,1.66032,1.65982,1.70721,1.69413,1.75582,1.47492,1.54322,1.54322,1.39871,1.43952,1.41674,1.45014,1.39692,1.28693,1.77082},
    { 2.10332,2.09977,2.09982,2.09985,2.09988,2.09688,2.09201,2.09203,2.0901,2.09026,2.08691,2.08696,2.09004,2.09108,2.08934,2.08776,1.8488,1.8455,1.84261,1.83764,1.83772} };
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

static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.60365, 3.95191, 1.99512  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 48.9953,49.9856,50.7261,52.7478,53.9665,58.6942,60.1205,60.1823,59.9366,60.6265,60.6259,82.1603,82.3483,82.8352,82.5085,81.9415,83.5118,86.388,88.1385,97.6885,95.2449},
    { 12.8259,13.9888,14.3669,14.7466,14.7847,16.6521,16.517,17.6275,17.6694,18.7818,18.8259,74.4455,79.8819,83.9252,86.201,91.3646,97.8663,108.154,114.174,119.643,119.967},
    { 2.61503,2.85417,2.8524,2.85107,2.8502,3.08834,3.56793,3.56695,3.80403,3.79356,4.27367,4.26956,4.26157,4.49418,4.72838,4.96516,5.19363,5.65968,6.12308,5.24332,5.23514} };
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

static Double_t DYBkgScaleFactorBDTKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14042, 1.20401, 1.26796  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.30233,1.3025,1.3024,1.30195,1.30182,1.28713,1.28681,1.28697,1.28709,1.28701,1.28701,1.139,1.13962,1.1401,1.14056,1.14174,1.14511,1.13744,1.15513,1.13488,1.1432},
    { 1.51137,1.50911,1.50866,1.50912,1.51221,1.52926,1.52982,1.52911,1.53023,1.52978,1.53259,1.29545,1.29513,1.29497,1.29489,1.29464,1.20224,1.2497,1.26113,1.20504,1.20469},
    { 2.10332,2.09977,2.09982,2.09985,2.09988,2.09688,2.09201,2.09203,2.0901,2.09026,2.08691,2.08696,2.09004,2.09108,2.08934,2.08776,1.8488,1.8455,1.84261,1.83764,1.83772} };
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

