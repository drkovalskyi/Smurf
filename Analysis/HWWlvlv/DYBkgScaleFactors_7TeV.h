static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 1.69739, 2.72072, 3.80435  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,16.685,4.90401},
    { 43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,43.9678,24.9212},
    { 5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,5.42167,3.26177} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.44874, 1.32561, 1.32483  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.43671,1.58247},
    { 1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.32434,1.33761},
    { 1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,1.98361,2.19926} };
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

