static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.31188, 4.2367, 1.98627  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 45.0747,62.8625,81.3252,84.3904,84.6806,92.2172,92.2172,92.7455,103.606,93.5248,80.7221,80.7221,37.0137,37.0137,18.827,7.92401,4.34218,21.6306,20.1002,11.1747,11.1602},
    { 5.27753,8.89495,13.4501,12.8324,12.7545,14.6897,14.6897,14.2716,14.208,15.3243,14.3172,14.3172,11.7891,11.7891,8.80328,6.48572,7.50849,15.7061,16.6037,12.0837,7.94375},
    { 2.8919,3.3637,4.08452,5.09769,4.82299,4.17656,4.17656,4.02649,3.77713,3.96481,3.9877,3.9877,2.0469,2.0469,3.29733,5.26022,6.23648,6.15248,5.076,6.95594,4.05621} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31574, 1.31962, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31844,1.31419,1.31144,1.31321,1.3206,1.33601,1.33601,1.32344,1.32422,1.32838,1.32994,1.32994,1.34354,1.34354,1.5136,1.89467,2.3619,1.39356,1.38446,1.36068,1.44196},
    { 1.43162,1.37247,1.34768,1.35713,1.36447,1.35975,1.35975,1.3702,1.38565,1.37365,1.3725,1.3725,1.36419,1.36419,1.44953,1.51637,1.46895,1.35691,1.34496,1.34697,1.39401},
    { 1.46831,1.4557,1.44133,1.43512,1.44596,1.46295,1.46295,1.47982,1.50245,1.49391,1.49646,1.49646,1.57177,1.57177,1.86196,1.98322,1.5918,1.51036,1.47469,1.50296,1.60167} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.31188, 4.2367, 1.98627  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,214.757,72.5445},
    { 82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,82.7081,51.3039},
    { 20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,20.054,11.9819} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31574, 1.31962, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.31383,1.32157},
    { 1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.32185},
    { 1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.40707,1.5339} };
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

