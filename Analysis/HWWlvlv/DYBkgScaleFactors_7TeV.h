static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 1.69739, 2.72072, 3.80435  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 6.42228,7.42272,8.0475,9.42917,9.44855,9.55061,9.55061,7.72613,7.55779,8.19898,8.44881,8.44881,8.37203,8.37203,2.65471,1.57684,0.873558,1.36461,1.78504,0.0116291,0.240766},
    { 5.49617,7.17249,9.33636,9.30408,9.26136,9.3302,9.3302,9.44101,9.35273,11.2673,11.9517,11.9517,10.199,10.199,8.47814,10.2831,10.858,14.7813,14.0457,6.66384,5.68448},
    { 0.150383,0.302265,0.758691,0.897296,0.987863,1.65388,1.65388,1.48783,1.48783,1.48292,1.48292,1.48292,0.844706,0.844706,1.96281e-05,1.96128e-05,3.96422e-05,1.28269,1.37093,0.000138051,9.82131e-05} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.44874, 1.32561, 1.32483  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.43061,1.42983,1.4309,1.43223,1.43837,1.46181,1.46181,1.51004,1.57367,1.56254,1.56498,1.56498,1.51741,1.51741,2.31998,2.45782,4.44743,2.43283,2.018,27.1431,9.88982},
    { 1.3595,1.35241,1.34737,1.37166,1.43207,1.43783,1.43783,1.41056,1.41201,1.40464,1.40269,1.40269,1.4284,1.4284,1.44515,1.42799,1.41852,1.37437,1.36395,1.3703,1.79729},
    { 2.61875,2.44917,2.33881,2.34391,2.34715,2.2855,2.2855,2.31056,2.31056,2.3112,2.3112,2.3112,2.52941,2.52941,1.78051,1.78103,1.58698,2.50246,2.4895,1.40431,1.4401} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 1.69739, 2.72072, 3.80435  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,4.90401},
    { 43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,24.9212},
    { 5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,3.26177} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.44874, 1.32561, 1.32483  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.58247},
    { 1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.33761},
    { 1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,2.19926} };
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

