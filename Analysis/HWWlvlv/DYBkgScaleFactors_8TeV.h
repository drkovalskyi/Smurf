static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.29144, 3.36513, 1.88769  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 25.2273,35.7793,41.266,44.4779,47.893,50.5781,50.5781,49.8633,57.1387,51.0527,45.2996,45.2996,18.4594,18.4594,11.5058,2.84783,2.10473,8.10156,7.35707,5.34335,4.10627},
    { 3.56608,4.59719,7.94589,7.54577,8.26697,8.82283,8.82283,9.84729,10.8609,10.0042,8.37024,8.37024,7.23447,7.23447,5.51001,5.50789,4.52007,9.63808,9.62509,6.44241,4.20961},
    { 1.72124,2.20942,3.19694,3.34061,2.93585,2.39533,2.39533,2.78753,2.78438,3.53662,3.80298,3.80298,3.22265,3.22265,3.75574,5.25903,4.91325,6.08628,4.31891,5.19172,2.25686} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.63694, 1.44236, 1.20413  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.66591,1.66169,1.66077,1.68285,1.67118,1.69151,1.69151,1.83632,1.79842,1.85336,1.50636,1.50636,1.48043,1.48043,1.88071,2.87682,3.35548,1.64058,2.10529,2.04484,2.1721},
    { 1.43505,1.40101,1.33796,1.36734,1.78524,1.76087,1.76087,1.7587,1.61706,2.00699,2.03951,2.03951,2.02402,2.02402,2.07361,2.07276,2.08507,1.4104,1.29603,2.02723,2.05502},
    { 1.64023,1.61473,1.58537,1.54174,1.57309,1.64359,1.64359,1.52592,1.55393,1.71084,1.7369,1.7369,1.79164,1.79164,2.02908,2.15852,2.12858,2.08519,2.08393,2.07228,2.20474} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.29144, 3.36513, 1.88769  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,118.772,40.2733},
    { 45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,45.0404,26.0312},
    { 13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,13.0643,8.65336} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.63694, 1.44236, 1.20413  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,1.85715,2.01087},
    { 1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,1.61118,2.01305},
    { 1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.48454,1.81012} };
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

