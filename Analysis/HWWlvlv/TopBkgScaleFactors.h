Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.39442, 1.09648, 1.19488   };
  return TopBkgScaleFactor[jetBin];
}

Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.21003, 1.06076, 1.25000   };
  return TopBkgScaleFactorKappa[jetBin];
}

