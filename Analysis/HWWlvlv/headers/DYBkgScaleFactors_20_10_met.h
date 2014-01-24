static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.33382, 4.29494, 2.01531  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 44.7699,62.2535,80.7163,84.0196,84.2118,92.2428,92.2428,92.0498,102.776,92.7502,80.1407,80.1407,37.782,37.782,18.6794,7.86802,3.72146,21.9608,20.3559,11.2138,11.3676},
    { 5.26462,8.87387,13.4177,12.8065,12.7298,14.6684,14.6684,14.2541,14.1975,15.3008,14.2774,14.2774,11.7511,11.7511,8.78977,6.47567,7.50232,15.9373,16.7967,12.1844,8.05013},
    { 2.89458,3.36691,4.08792,5.10131,4.824,4.1781,4.1781,4.02485,3.77453,3.96193,3.98473,3.98473,2.04408,2.04408,3.29151,5.25066,6.23517,6.15282,5.07998,6.95872,4.05891} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31559, 1.31918, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31857,1.31434,1.31152,1.31799,1.32567,1.3377,1.3377,1.32363,1.32443,1.32864,1.33023,1.33023,1.34219,1.34219,1.51319,1.89575,2.54881,1.39049,1.38188,1.35957,1.4379},
    { 1.43173,1.37253,1.34773,1.35716,1.36449,1.35978,1.35978,1.37018,1.38556,1.37361,1.37247,1.37247,1.36426,1.36426,1.44968,1.51626,1.46861,1.35557,1.34404,1.34592,1.39087},
    { 1.4681,1.45551,1.44119,1.435,1.44585,1.4628,1.4628,1.47979,1.50243,1.4939,1.49645,1.49645,1.57341,1.57341,1.8616,1.98306,1.59166,1.51027,1.47458,1.5029,1.60158} };
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
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.33382, 4.29494, 2.01531  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,214.587,72.4564},
    { 83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,83.4252,51.7811},
    { 19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,19.9434,11.916} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31559, 1.31918, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.31369,1.32139},
    { 1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.31862,1.32134},
    { 1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.41372,1.54031} };
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

