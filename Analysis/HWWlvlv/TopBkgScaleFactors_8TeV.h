Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.03291, 1.09622, 1.59376   };
  return TopBkgScaleFactor[jetBin];
}

Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.23788, 1.08813, 1.39466   };
  return TopBkgScaleFactorKappa[jetBin];
}

