static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 10.428, 3.97728, 2.31165  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 32.9982,48.2528,56.5573,63.5464,69.6013,72.1489,72.1489,72.0178,76.6173,69.4068,63.0764,63.0764,22.3512,22.3512,18.1926,8.33596,9.30812,14.7848,13.8192,8.68079,9.33823},
    { 3.63342,5.01935,7.58654,6.92514,7.34249,7.47136,7.47136,8.29908,8.99327,8.1631,6.83231,6.83231,7.04443,7.04443,7.45027,6.65772,6.26723,14.3896,15.5788,10.3464,6.82593},
    { 1.74339,2.23746,2.98673,3.08366,3.03584,2.47652,2.47652,2.41014,2.7953,3.5573,3.62006,3.62006,3.05738,3.05738,2.39464,3.59308,4.51216,6.92569,5.33554,4.87756,1.65497} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.18856, 1.39418, 1.09509  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.34769,1.34134,1.33981,1.44328,1.3526,1.27976,1.27976,1.22058,1.21635,1.21312,1.17638,1.17638,1.49265,1.49265,1.36463,1.5907,2.11315,1.27823,1.29851,1.22423,3.01051},
    { 1.34206,1.29342,1.25353,1.26564,1.27089,1.28485,1.28485,1.29386,1.2916,1.29498,1.32378,1.32378,1.32446,1.32446,1.46206,1.45231,1.42665,1.27289,1.23228,1.17313,1.3209},
    { 1.56428,1.53525,1.50783,1.50143,1.53206,1.53543,1.53543,1.5411,1.56591,1.70953,1.71887,1.71887,1.77432,1.77432,1.98199,2.1616,2.12715,2.07797,2.07557,2.07714,2.28016} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 10.428, 3.97728, 2.31165  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,178.284,61.8694},
    { 57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,57.8674,34.3043},
    { 13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,13.1181,7.50134} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.18856, 1.39418, 1.09509  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.17853,1.18427},
    { 1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.40735,1.13083},
    { 1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.40093,1.9443} };
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

