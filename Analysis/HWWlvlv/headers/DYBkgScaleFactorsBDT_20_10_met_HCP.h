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

