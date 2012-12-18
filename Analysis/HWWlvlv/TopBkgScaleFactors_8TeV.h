static Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.0172, 1.04882, 1.04211   };
  return TopBkgScaleFactor[jetBin];
}

static Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.19517, 1.03188, 1.0259   };
  return TopBkgScaleFactorKappa[jetBin];
}

