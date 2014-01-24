static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.31831, 3.48874, 1.61643  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 5.13531,5.31224,5.08014,5.52396,5.67646,5.80125,5.80125,3.75283,3.94043,4.53523,4.5078,9.75859,4.43294,4.43294,3.73741,0.639086,0.630674,2.42199,2.24642,1.35514,1.57414},
    { 1.35465,1.48052,1.3437,1.08587,1.04038,0.890331,0.890331,2.17386,2.04224,1.97811,1.94318,7.63844,9.73881,9.73881,8.80979,8.05391,9.26346,12.7806,11.2557,8.26549,5.21292},
    { 4.73648,6.67019,8.0844,8.5426,9.31266,8.60674,8.60674,7.91876,7.54063,8.3078,9.28001,9.28001,11.16,11.16,8.88319,7.62502,8.13152,12.3731,12.2241,11.3595,11.3376} };
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

static Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.27718, 1.2173, 1.20118  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.28971,1.30489,1.32098,1.34268,1.38344,1.47711,1.47711,1.74487,1.78984,1.72217,1.72575,1.6971,1.82299,1.82299,2.82071,2.32184,2.12589,2.03988,1.74619,1.70814,1.72615},
    { 2.14716,2.14327,2.16273,2.2688,2.32499,2.43698,2.43698,1.76233,1.80942,1.82295,1.8308,1.71885,1.61278,1.61278,1.46179,1.51926,1.48124,1.52579,1.46654,1.27457,1.94157},
    { 1.63669,1.63154,1.62934,1.55399,1.44583,1.45828,1.45828,1.44497,1.4334,1.43054,1.42754,1.42754,1.61197,1.61197,1.89574,1.86636,2.01161,2.1855,2.09765,1.41947,1.44177} };
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

