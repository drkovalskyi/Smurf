static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.41526, 4.28409, 2.00013  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,217.055,73.3652},
    { 83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,83.1877,51.579},
    { 19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,19.9389,11.9131} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.31535, 1.31933, 1.30171  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.31351,1.32111},
    { 1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.31877,1.32153},
    { 1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.41373,1.54032} };
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

