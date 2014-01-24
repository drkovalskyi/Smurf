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

