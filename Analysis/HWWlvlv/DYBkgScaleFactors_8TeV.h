static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.15149, 4.04835, 1.78384  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 18.3935,20.2845,19.971,23.1732,24.8238,26.0957,26.0957,23.0436,26.641,25.7657,22.4099,32.3725,8.99181,8.99181,5.14129,0.933906,1.13735,5.0475,5.63916,2.58201,2.97812},
    { 4.60894,6.01138,7.34225,6.88424,6.46078,6.87114,6.87114,9.55947,7.24745,7.70403,4.60553,20.8183,14.9642,14.9642,15.7298,14.146,14.4657,18.4781,16.4676,11.7541,8.5011},
    { 4.12677,4.58692,4.5844,4.58328,4.58199,5.04253,5.04251,5.04128,5.49854,5.48578,5.94719,5.94072,6.40017,7.32045,7.77537,7.77174,8.21773,10.9627,12.0674,13.2507,13.2336} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14843, 1.265, 1.17705  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.23135,1.22996,1.2359,1.26151,1.31612,1.36621,1.36621,1.39816,1.43985,1.5029,1.51232,1.48614,1.76285,1.76285,2.72896,2.17543,5.53423,1.96347,1.61576,1.9631,1.55686},
    { 1.69364,1.67464,1.66358,1.69517,1.72451,1.71797,1.71797,1.7036,1.76282,1.75971,1.86942,1.45498,1.58017,1.58017,1.4189,1.47571,1.4475,1.49126,1.43475,1.24508,1.90841},
    { 2.13672,2.13172,2.13177,2.13179,2.13182,2.1277,2.1277,2.12772,2.12432,2.1245,2.12154,2.12161,2.11909,2.11497,2.1133,2.11333,2.1119,2.09527,2.09311,1.43546,1.43559} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.15149, 4.04835, 1.78384  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 36.5949,37.0593,37.1834,38.687,39.8833,43.4694,44.4702,44.3805,44.1897,44.1654,44.1649,57.8469,57.4806,58.1579,58.1671,56.9669,58.2259,59.9297,61.2625,69.0535,67.4314},
    { 9.00669,10.0544,10.5416,10.6822,11.1556,12.1497,12.0405,12.1723,12.2985,13.7324,14.3764,48.4077,51.8054,54.3776,56.1067,59.9348,64.6611,71.2794,74.9116,78.2718,77.9029},
    { 4.12677,4.58692,4.5844,4.58328,4.58199,5.04253,5.04251,5.04128,5.49854,5.48578,5.94719,5.94072,6.40017,7.32045,7.77537,7.77174,8.21773,10.9627,12.0674,13.2507,13.2336} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14843, 1.265, 1.17705  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.37545,1.37576,1.37586,1.37529,1.37504,1.35764,1.35729,1.35758,1.35771,1.35783,1.35783,1.14219,1.14359,1.14409,1.14459,1.14708,1.14924,1.14314,1.16,1.14193,1.15023},
    { 1.56824,1.56368,1.56201,1.56392,1.56583,1.59462,1.59549,1.59906,1.60077,1.59574,1.59637,1.36651,1.36617,1.36593,1.36582,1.36543,1.2672,1.31133,1.32233,1.26525,1.26535},
    { 2.13672,2.13172,2.13177,2.13179,2.13182,2.1277,2.1277,2.12772,2.12432,2.1245,2.12154,2.12161,2.11909,2.11497,2.1133,2.11333,2.1119,2.09527,2.09311,1.43546,1.43559} };
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

