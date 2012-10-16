static Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 0.980734, 1.08389, 1.08267   };
  return TopBkgScaleFactor[jetBin];
}

static Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.20604, 1.03266, 1.03313   };
  return TopBkgScaleFactorKappa[jetBin];
}

