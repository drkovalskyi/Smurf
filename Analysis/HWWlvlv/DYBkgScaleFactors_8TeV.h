static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.99679, 4.4907, 2.01994  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 45.3827,61.6894,79.367,80.253,80.2882,86.8727,86.8727,87.0773,97.4558,89.5298,77.9245,77.9245,35.245,35.245,20.9003,8.55292,5.6227,20.9976,19.9465,10.7175,10.4647},
    { 4.48433,7.78964,12.996,12.927,12.5658,14.2201,14.2201,14.2857,14.9363,15.985,14.5824,14.5824,11.6029,11.6029,8.253,6.39961,7.12643,15.2996,15.5534,11.7697,7.92836},
    { 2.81619,3.27687,4.21632,5.00098,4.51426,3.80142,3.80142,3.97154,3.71507,3.90031,3.92031,3.92031,2.19613,2.19613,4.78999,7.31615,7.86767,7.25684,5.95664,7.28447,4.14709} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31619, 1.32082, 1.30172  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31854,1.31468,1.31193,1.31441,1.32972,1.34611,1.34611,1.32614,1.32665,1.33041,1.33169,1.33169,1.34761,1.34761,1.49133,1.84207,2.06825,1.3932,1.38486,1.39441,1.45476},
    { 1.46315,1.38716,1.34975,1.35537,1.36513,1.36184,1.36184,1.36889,1.37872,1.36829,1.36931,1.36931,1.3635,1.3635,1.45127,1.50418,1.47213,1.35763,1.34791,1.3476,1.39261},
    { 1.46821,1.45552,1.43723,1.43487,1.45089,1.47031,1.47031,1.47985,1.50254,1.49393,1.50589,1.50589,1.79651,1.79651,1.80471,1.8907,1.68693,1.59529,1.61836,1.46397,1.59941} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.99679, 4.4907, 2.01994  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,210.929,71.6649},
    { 79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,79.8919,49.0057},
    { 19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,19.6145,12.25} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31619, 1.32082, 1.30172  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.31412,1.32195},
    { 1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.31992,1.32325},
    { 1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.37145,1.51736} };
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

