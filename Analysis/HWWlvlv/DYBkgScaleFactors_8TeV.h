static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.02016, 3.17328, 2.02906  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 24.8852,35.2923,40.7324,45.1284,49.0589,51.7219,51.7219,50.9528,56.1193,48.9294,46.0007,46.0007,13.476,13.476,10.0915,2.38069,1.81084,6.56709,6.24789,3.98136,3.0386},
    { 3.02769,3.90479,6.73656,6.3827,6.96089,7.49644,7.49644,8.3552,9.06031,8.19321,6.86358,6.86358,5.86794,5.86794,4.80927,4.80855,4.02113,9.13528,9.48884,6.23621,3.77307},
    { 1.72083,2.20909,3.19652,3.34021,2.93548,2.39539,2.39539,2.78767,2.78448,3.53643,3.80284,3.80284,3.2231,3.2231,3.75676,5.26026,4.91395,6.08612,4.31848,5.19081,2.25656} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14955, 1.33684, 1.10891  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.30293,1.29366,1.29155,1.38544,1.29643,1.25382,1.25382,1.2689,1.25296,1.25101,1.23119,1.23119,1.67855,1.67855,1.56957,2.78752,3.3349,1.53264,1.49266,1.43186,3.1795},
    { 1.39413,1.35672,1.28574,1.29529,1.29719,1.30531,1.30531,1.3147,1.30359,1.30842,1.342,1.342,1.36699,1.36699,1.55584,1.53457,1.53826,1.32197,1.2847,1.24361,1.41547},
    { 1.64034,1.61483,1.58546,1.54179,1.57315,1.64366,1.64366,1.52593,1.55379,1.71077,1.73683,1.73683,1.79163,1.79163,2.02907,2.15851,2.12858,2.0852,2.08395,2.0723,2.20478} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.02016, 3.17328, 2.02906  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,113.486,37.0216},
    { 41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,41.56,23.9578},
    { 13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,8.65576} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14955, 1.33684, 1.10891  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.14642,1.23741},
    { 1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.34642,1.17061},
    { 1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.81016} };
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

