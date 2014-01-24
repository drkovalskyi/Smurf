static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.61075, 4.05844, 1.97568  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 48.7357,49.7302,50.4708,52.4925,53.7134,58.4354,59.8602,59.9224,59.6793,60.371,60.3704,82.6956,82.8985,83.3997,83.077,82.5213,84.1221,87.0376,88.451,98.0917,95.6686},
    { 12.3152,13.4888,13.8702,14.2537,14.3008,16.0986,15.9635,17.0819,17.1259,18.2565,18.317,75.3649,80.8538,84.9583,87.2601,92.2644,98.8382,109.228,115.31,120.853,121.189},
    { 2.6186,2.85763,2.85595,2.8551,2.85423,3.09092,3.57056,3.56958,3.80649,3.79746,4.27695,4.27281,4.26208,4.49496,4.72997,4.96688,5.1938,5.66673,6.14391,5.26751,5.25933} };
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

static Double_t DYBkgScaleFactorBDTKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.12416, 1.22258, 1.27555  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.29924,1.2994,1.29931,1.29884,1.29871,1.29202,1.2917,1.29187,1.29198,1.2919,1.2919,1.11592,1.11664,1.11719,1.11773,1.11912,1.1317,1.12792,1.13413,1.11887,1.12596},
    { 1.51539,1.51289,1.51239,1.51287,1.51613,1.54352,1.54411,1.54332,1.54447,1.5439,1.54675,1.31722,1.31692,1.31677,1.31669,1.31647,1.22302,1.26958,1.28065,1.22368,1.2233},
    { 2.10321,2.09969,2.09974,2.09976,2.09978,2.09686,2.09202,2.09203,2.09011,2.09025,2.08693,2.08698,2.09008,2.09112,2.08937,2.0878,1.84578,1.84263,1.83998,1.83431,1.83439} };
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

