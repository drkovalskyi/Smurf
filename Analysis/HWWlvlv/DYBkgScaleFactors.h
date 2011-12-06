Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.05333, 3.77023, 8.26298  };
  Double_t DYBkgScaleFactorHiggsSelection[3][18] = { 
    { 3.62265,3.45818,3.68928,4.19747,4.70718,4.64437,3.68338,4.60246,4.96441,5.11402,3.53062,2.38553,2.78923,6.95108,22.814,21.144,0.164476,0.180441},
    { 1.9284,2.61853,3.3674,2.59525,2.58522,2.39913,2.37523,2.23619,2.86873,2.54084,4.1612,3.66422,6.43097,7.56567,7.91583,7.56728,2.70019,4.55992},
    { 4.61932,4.50239,5.48947,5.8775,5.99477,5.24129,4.92926,5.32894,4.9428,5.27937,4.94931,5.24715,6.06666,8.74202,8.77867,7.93401,6.95753,7.1839} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.62216, 1.195, 1.18877  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][18] = { 
    { 1.71408,1.71499,1.71602,1.5841,1.53852,1.46239,1.71253,1.67433,1.66573,1.66728,1.97881,2.94665,3.22312,3.46736,2.53005,2.06235,5.00285,9.83972},
    { 1.4439,1.43709,1.43326,1.50766,1.54834,1.51361,1.5246,1.51938,1.50944,1.50779,1.36303,1.41852,1.42939,1.42123,1.39771,1.35249,1.32435,1.32602},
    { 1.27235,1.26491,1.25952,1.29255,1.33194,1.33929,1.37152,1.40794,1.40588,1.40514,1.42581,1.56323,1.40958,1.28809,1.25688,1.26009,1.51671,1.46785} };
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

