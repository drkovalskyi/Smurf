static Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.05513, 1.07211, 1.07562   };
  return TopBkgScaleFactor[jetBin];
}

static Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.1926, 1.03155, 1.02575   };
  return TopBkgScaleFactorKappa[jetBin];
}

