static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 7.20511, 3.9943, 1.94123  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 43.2603,60.0708,71.0979,75.6906,77.932,82.2983,82.2983,80.6712,88.8181,80.8726,68.5283,68.5283,31.3169,31.3169,19.3065,4.42766,3.45765,18.8293,18.7834,9.79293,9.43448},
    { 4.39633,6.00481,11.6174,10.728,10.7869,11.2787,11.2787,11.4989,12.855,12.4053,10.1596,10.1596,10.1044,10.1044,7.24066,6.63746,6.64194,15.413,15.0902,11.121,7.67461},
    { 2.52341,2.97478,3.8937,4.14038,4.48312,3.42972,3.42972,3.17666,2.76387,3.48411,3.50589,3.50589,3.45118,3.45118,4.05072,6.11393,6.50046,6.41141,5.47001,6.34334,3.12971} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31688, 1.32177, 1.30183  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31955,1.31504,1.31367,1.31636,1.32708,1.34917,1.34917,1.3279,1.32883,1.33291,1.33619,1.33619,1.36082,1.36082,1.54838,2.56248,2.63918,1.40447,1.38622,1.39002,1.45702},
    { 1.46301,1.41929,1.35352,1.3656,1.37262,1.37563,1.37563,1.3834,1.38293,1.3861,1.40201,1.40201,1.37314,1.37314,1.47193,1.48706,1.48279,1.35422,1.34735,1.3473,1.38924},
    { 1.48968,1.47533,1.45531,1.45007,1.45443,1.48199,1.48199,1.50464,1.5371,1.50691,1.50927,1.50927,1.57986,1.57986,1.82455,1.965,1.64279,1.54995,1.50268,1.45309,1.79193} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 7.20511, 3.9943, 1.94123  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,195.098,66.5719},
    { 74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,74.3488,45.2467},
    { 17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,17.651,10.1035} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31688, 1.32177, 1.30183  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.31485,1.32262},
    { 1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32056,1.32385},
    { 1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,1.60335,2.04004} };
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

