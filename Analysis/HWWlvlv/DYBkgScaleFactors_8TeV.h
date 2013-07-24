static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.33639, 4.23607, 1.98672  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 45.4406,63.2461,81.7321,84.8995,85.2648,92.8983,92.8983,93.5421,104.581,94.4633,81.6099,81.6099,37.4818,37.4818,18.7886,7.90582,4.31929,21.7962,20.2022,11.1644,11.2812},
    { 5.27896,8.89869,13.454,12.8352,12.7557,14.6898,14.6898,14.273,14.2101,15.3353,14.3398,14.3398,11.7974,11.7974,8.38979,6.09362,7.18342,15.4866,16.4279,11.9797,7.93912},
    { 2.64909,3.12137,3.84236,5.09832,4.82067,4.17455,4.17455,4.02478,3.77449,3.96189,3.98465,3.98465,2.04393,2.04393,3.29129,5.25043,6.23499,6.1524,5.07778,6.95649,4.05658} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31564, 1.31962, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31828,1.3141,1.31138,1.31314,1.31629,1.32869,1.32869,1.3232,1.32396,1.32804,1.32957,1.32957,1.3429,1.3429,1.51363,1.89495,2.35913,1.39152,1.38312,1.36112,1.44184},
    { 1.43163,1.37246,1.34769,1.35714,1.3645,1.35977,1.35977,1.37024,1.3857,1.37369,1.37254,1.37254,1.36426,1.36426,1.45932,1.53472,1.47945,1.3581,1.34565,1.34753,1.394},
    { 1.47655,1.4619,1.44561,1.43508,1.44596,1.46295,1.46295,1.4798,1.50244,1.49392,1.49647,1.49647,1.57345,1.57345,1.86163,1.98308,1.59167,1.51027,1.47462,1.50291,1.60163} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.33639, 4.23607, 1.98672  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,215.593,72.8273},
    { 82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,82.6936,51.2669},
    { 19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,19.6976,11.9118} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31564, 1.31962, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.31374,1.32151},
    { 1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.31904,1.32186},
    { 1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.41396,1.54036} };
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

