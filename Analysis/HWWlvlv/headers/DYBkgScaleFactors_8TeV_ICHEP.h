static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.64965, 3.06983, 1.89593  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 12.1589,14.376,15.9424,17.9746,21.1091,20.7244,20.7244,18.7806,21.1629,19.3239,16.7333,16.7333,6.15002,6.15002,4.52802,2.91971,2.38975,2.39797,2.28673,1.39512,0.916324},
    { 1.43265,1.77909,2.43445,2.26498,2.74459,3.21998,3.21998,2.91896,3.11601,2.85103,2.55229,2.55229,1.32437,1.32437,2.73875,2.4371,1.58118,3.26534,4.30409,2.27453,1.35992},
    { 0.413953,0.620922,0.828922,0.743788,1.10894,0.877895,0.877895,0.998296,0.767418,1.15681,1.15629,1.15629,0.984367,0.984367,1.20111,1.19655,1.4848,1.47213,1.36815,0.688742,-0.0038552} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.34633, 1.37394, 1.30637  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.51167,1.51014,1.5084,1.50586,1.5821,1.72598,1.72598,1.53967,1.53429,1.58205,1.56481,1.56481,1.54546,1.54546,2.31416,2.51227,2.53668,2.05334,1.96705,1.99469,2.65562},
    { 1.62774,1.56512,1.51343,1.54998,1.51252,1.4881,1.4881,1.54249,1.55939,1.62226,1.65366,1.65366,1.95186,1.95186,1.58492,1.65132,1.88587,1.57607,1.47033,1.52402,1.7508},
    { 1.92119,1.82274,1.76839,1.82883,1.77704,1.83329,1.83329,1.83715,1.93964,1.84184,1.84199,1.84199,2.3872,2.3872,2.43164,2.433,2.36205,2.28728,2.26523,1.9035,-175.837} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.64965, 3.06983, 1.89593  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,47.2512,16.3752},
    { 16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,16.994,10.4995},
    { 3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,3.50782,1.94057} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.34633, 1.37394, 1.30637  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.34218,1.36607},
    { 1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.37711,1.39342},
    { 1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,1.8383,2.1852} };
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

