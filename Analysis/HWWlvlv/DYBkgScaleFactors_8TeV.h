static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.05734, 4.03429, 1.93404  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 44.9973,62.7464,81.1692,84.216,84.477,91.9831,91.9831,92.4893,103.312,93.2308,80.4987,80.4987,36.9057,36.9057,18.7188,7.8172,4.23068,21.4992,19.9687,11.0874,10.9969},
    { 5.25788,8.86863,13.4163,12.7971,12.7168,14.65,14.65,14.2283,14.1597,15.2743,14.2771,14.2771,11.7669,11.7669,8.78178,6.46411,7.48024,15.6654,16.5559,12.0269,7.88368},
    { 2.89182,3.3636,4.08445,5.09761,4.823,4.17652,4.17652,4.02652,3.77719,3.96489,3.98785,3.98785,2.04681,2.04681,3.2972,5.26024,6.23631,6.15224,5.07565,6.95517,4.05555} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31616, 1.32001, 1.30172  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31851,1.31424,1.31148,1.31327,1.32068,1.33609,1.33609,1.32357,1.32436,1.32856,1.33011,1.33011,1.34378,1.34378,1.51559,1.90571,2.39661,1.39467,1.3856,1.36177,1.44605},
    { 1.43249,1.37288,1.34793,1.35744,1.36485,1.36007,1.36007,1.37061,1.3862,1.37411,1.3729,1.3729,1.36442,1.36442,1.45017,1.51757,1.47005,1.35721,1.34525,1.34748,1.39547},
    { 1.46832,1.45571,1.44135,1.43513,1.44597,1.46296,1.46296,1.47983,1.50247,1.49392,1.49647,1.49647,1.57182,1.57182,1.862,1.98324,1.59181,1.5104,1.4747,1.50299,1.60171} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.05734, 4.03429, 1.93404  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,213.234,71.91},
    { 82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,82.1928,50.9535},
    { 20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,20.0525,11.9807} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31616, 1.32001, 1.30172  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.31417,1.32215},
    { 1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.31939,1.32228},
    { 1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.40699,1.53399} };
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

