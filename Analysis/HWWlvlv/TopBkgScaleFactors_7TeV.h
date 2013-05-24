static Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.24224, 1.057, 1.17214   };
  return TopBkgScaleFactor[jetBin];
}

static Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.21692, 1.06435, 1.06319   };
  return TopBkgScaleFactorKappa[jetBin];
}

