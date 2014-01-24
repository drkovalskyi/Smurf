static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.15149, 4.04835, 2.07817  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 18.3935,20.2845,19.971,23.1732,24.8238,26.0957,26.0957,23.0436,26.641,25.7657,22.4099,32.3725,8.9918,8.9918,5.14129,0.933906,1.13735,5.0475,5.63916,2.58201,2.97812},
    { 4.60894,6.01138,7.34225,6.88424,6.46078,6.87114,6.87114,9.55947,7.24745,7.70403,4.60553,20.8183,14.9642,14.9642,15.7298,14.146,14.4657,18.4781,16.4676,11.7541,8.5011},
    { 2.08035,2.31231,2.31104,2.31048,2.30983,2.54199,2.54198,2.54136,2.77187,2.76544,2.99804,2.99478,3.22639,3.69031,3.91964,3.91781,4.14265,4.60341,5.06742,5.53082,5.52369} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14843, 1.265, 1.25434  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.23135,1.22996,1.2359,1.26151,1.31612,1.36621,1.36621,1.39816,1.43985,1.5029,1.51232,1.48614,1.76285,1.76285,2.72896,2.17543,5.53423,1.96347,1.61576,1.9631,1.55686},
    { 1.69364,1.67464,1.66358,1.69517,1.72451,1.71797,1.71797,1.7036,1.76282,1.75971,1.86942,1.45498,1.58017,1.58017,1.4189,1.47571,1.4475,1.49126,1.43475,1.24508,1.90841},
    { 1.75011,1.74252,1.7426,1.74263,1.74267,1.73537,1.73537,1.7354,1.73018,1.73045,1.72588,1.726,1.72209,1.71569,1.7131,1.71314,1.71091,1.7069,1.70355,1.80167,1.66994} };
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

