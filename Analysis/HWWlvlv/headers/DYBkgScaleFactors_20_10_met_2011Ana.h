Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.42246, 5.22338, 2.25589  };
  Double_t DYBkgScaleFactorHiggsSelection[3][20] = { 
    { 1.31689,1.16132,1.00723,0.61172,0.433257,0.0376201,0.0376201,0.142152,1.21744,1.0328,1.55383,1.55383,0.367971,0.564857,0.126096,0.116427,0.295187,0.350533,0.335328,1.15977},
    { 2.15647,2.86144,3.17854,2.8198,2.66859,2.55787,2.55787,2.98856,2.60192,2.4859,2.32374,2.32374,2.00555,2.87193,2.23385,2.01814,2.92669,2.00687,1.22758,0.99461},
    { 2.42611,3.175,3.92501,4.3501,3.70744,4.25893,4.25893,4.146,4.26613,3.97405,3.43954,3.43954,3.0764,2.20484,1.80948,1.87926,2.91758,2.83742,2.44544,2.51545} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 20 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.68234, 1.40436, 1.23335  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][20] = { 
    { 1.80684,1.97282,2.18773,3.1403,4.19633,46.7921,46.7921,16.0209,2.85836,3.1472,2.08668,2.08668,3.51498,3.10052,2.87479,2.8558,3.54874,2.81413,1.98944,2.17471},
    { 1.64463,1.62788,1.6232,1.75224,1.89195,1.89555,1.89555,1.93335,2.01646,1.96485,1.93515,1.93515,2.03981,1.89707,1.91492,1.8388,1.66067,1.70828,1.51461,2.68245},
    { 1.42683,1.40387,1.38897,1.3621,1.40768,1.45081,1.45081,1.44774,1.41475,1.40887,1.42393,1.42393,1.30534,1.53245,1.52344,1.70171,1.76619,1.78632,1.43946,1.38292} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 20 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

