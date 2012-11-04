static Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.15994, 1.06107, 1.19284   };
  return TopBkgScaleFactor[jetBin];
}

static Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.21538, 1.06222, 1.06535   };
  return TopBkgScaleFactorKappa[jetBin];
}

