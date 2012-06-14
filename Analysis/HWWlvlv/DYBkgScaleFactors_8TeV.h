static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.15149, 4.04835, 1.78384  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 18.3935,20.2845,19.971,23.1732,24.8238,26.0957,26.0957,23.0436,26.641,25.7657,22.4099,32.3725,8.98541,8.98541,5.13166,0.931934,1.13578,5.04296,5.63561,2.58176,2.97808},
    { 4.60894,6.01138,7.34225,6.88424,6.46078,6.87114,6.87114,9.55947,7.24745,7.70403,4.60553,20.8183,14.9642,14.9642,15.7298,14.146,14.4657,18.4781,16.4676,11.7533,8.50106},
    { 15.2272,20.146,24.5077,25.686,28.3518,27.6817,27.6817,25.8884,25.8055,26.3504,22.9476,22.9476,17.0993,17.0993,15.1585,13.663,14.4543,20.2691,19.5073,18.4452,17.4118} };
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
    { 1.23135,1.22996,1.2359,1.26151,1.31612,1.36621,1.36621,1.39816,1.43985,1.5029,1.51232,1.48614,1.76389,1.76389,2.73364,2.17535,5.54163,1.96494,1.61644,1.96326,1.55686},
    { 1.69364,1.67464,1.66358,1.69517,1.72451,1.71797,1.71797,1.7036,1.76282,1.75971,1.86942,1.45498,1.58017,1.58017,1.4189,1.47571,1.4475,1.49126,1.43475,1.24513,1.90841},
    { 1.37191,1.36724,1.36471,1.31493,1.25823,1.26935,1.26935,1.25877,1.26634,1.24357,1.24568,1.24568,1.59677,1.59677,1.87745,1.85105,2.00069,2.17279,2.07441,1.39849,1.41659} };
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

