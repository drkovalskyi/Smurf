static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.61887, 3.65061, 1.85388  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 11.4668,12.6231,12.9996,14.5267,15.2944,17.1661,17.1661,17.2433,18.8097,19.2494,19.1907,29.6489,12.7887,12.7887,9.15453,2.82527,1.70731,5.27327,6.52442,3.52161,4.0589},
    { 3.25848,4.34062,4.74342,4.17166,4.20841,2.9203,2.9203,3.98493,2.80638,3.16305,2.60621,21.9576,25.3049,25.3049,27.0314,23.8848,23.2025,32.3057,30.3995,20.8028,14.5749},
    { 2.02788,2.02712,2.02581,2.02523,2.02469,2.20723,2.57826,2.5775,2.76076,2.75424,3.12517,3.12416,2.93137,3.11237,3.29445,3.47793,3.65366,4.01949,4.20308,3.10804,3.10249} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31844, 1.15623, 1.24688  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.21215,1.2108,1.2124,1.21927,1.26154,1.32892,1.32892,1.34399,1.40083,1.39954,1.39994,1.5891,1.69134,1.69134,2.4668,2.94009,4.9057,1.97529,1.6662,1.96775,1.52467},
    { 2.07623,2.06738,2.0659,2.0781,2.08607,2.12305,2.12305,1.62132,1.74846,1.7213,1.82355,1.66495,1.54104,1.54104,1.39696,1.43759,1.41498,1.43867,1.39016,1.22478,1.77248},
    { 2.12369,2.12373,2.12378,2.1238,2.12382,2.1205,2.115,2.11502,2.11287,2.11302,2.10935,2.10937,2.11488,2.11606,2.11387,2.11192,1.67096,1.66608,1.66397,1.83575,1.83586} };
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

