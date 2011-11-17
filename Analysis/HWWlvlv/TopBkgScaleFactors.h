Double_t TopBkgScaleFactor(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactor[3] = { 1.39426, 1.09693, 1.12014   };
  return TopBkgScaleFactor[jetBin];
}

Double_t TopBkgScaleFactorKappa(Int_t jetBin) {
  assert(jetBin >=0 && jetBin <= 2);
  Double_t TopBkgScaleFactorKappa[3] = { 1.21003, 1.06076, 1.34503   };
  return TopBkgScaleFactorKappa[jetBin];
}

