Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 2.22112, 2.46277, 3.16458  };
  Double_t DYBkgScaleFactorHiggsSelection[3][18] = { 
    { 3.65304,4.26228,3.45614,3.1491,3.23441,2.62245,2.35818,2.72939,2.66236,2.40131,3.06974,1.03943,1.29796,0.739913,2.20456,0.819261,0.119594,1.17595},
    { 3.57456,4.07207,3.90905,2.67279,2.81023,2.8621,2.73723,2.49352,3.12179,3.29926,2.9183,3.7569,5.65156,5.83702,4.18396,4.22443,1.33643,2.49317},
    { 2.22503,2.81944,3.02105,3.29723,3.55867,3.89334,3.47649,3.96763,4.28986,4.57054,3.05996,2.66737,2.98253,3.76817,3.15128,3.56392,1.75131,2.70318} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.27618, 1.15686, 1.11473  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][18] = { 
    { 1.55041,1.54889,1.55036,1.4995,1.51898,1.7074,1.54051,1.52773,1.5221,1.52912,2.09665,2.06232,2.02073,2.80953,1.74448,2.9043,8.31972,9.87529},
    { 1.22811,1.2156,1.20924,1.24495,1.26774,1.29399,1.34166,1.34802,1.33718,1.33501,1.48966,1.4073,1.40074,1.38995,1.23222,1.19204,1.46612,1.47479},
    { 1.37588,1.37108,1.36773,1.38164,1.40561,1.42187,1.44855,1.47671,1.47512,1.47453,1.2856,1.32719,1.39477,1.26564,1.36928,1.25442,1.60005,1.237} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

