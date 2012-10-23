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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.326, 1.33469, 1.30277  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.33581,1.32737,1.3255,1.32829,1.33225,1.34461,1.34461,1.3487,1.3485,1.35695,1.35776,1.35776,1.40803,1.40803,1.63987,2.87682,3.19921,1.59514,1.55826,1.46715,1.68105},
    { 1.46298,1.43115,1.37323,1.38574,1.3858,1.38679,1.38679,1.38377,1.38697,1.39764,1.41302,1.41302,1.39854,1.39854,1.49258,1.49072,1.51708,1.38512,1.37099,1.38104,1.45063},
    { 1.5492,1.51925,1.48413,1.49262,1.52019,1.55196,1.55196,1.55387,1.57774,1.54154,1.69047,1.69047,1.7378,1.7378,2.00774,2.15852,1.89698,1.7426,1.70205,1.49381,1.81989} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.326, 1.33469, 1.30277  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.32301,1.33445},
    { 1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.33238,1.34099},
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

