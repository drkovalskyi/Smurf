Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.39523, 1.09418, 1.19448   };
  return TopBkgScaleFactor[jetBin];
}

Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.20999, 1.06078, 1.02772   };
  return TopBkgScaleFactorKappa[jetBin];
}

