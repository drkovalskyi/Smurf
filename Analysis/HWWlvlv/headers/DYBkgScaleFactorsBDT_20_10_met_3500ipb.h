static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.15149, 4.04835, 2.07817  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 36.595,37.0593,37.1835,38.6871,39.8833,43.4694,44.4702,44.3806,44.1897,44.1654,44.165,57.847,57.4806,58.1579,58.1672,56.9669,58.2259,59.9297,61.2625,69.0535,67.4314},
    { 9.00669,10.0544,10.5416,10.6822,11.1556,12.1497,12.0405,12.1723,12.2985,13.7324,14.3764,48.4077,51.8054,54.3776,56.1067,59.9348,64.6611,71.2794,74.9116,78.2718,77.9029},
    { 2.08035,2.31231,2.31104,2.31048,2.30983,2.54199,2.54198,2.54136,2.77187,2.76544,2.99804,2.99478,3.22639,3.69031,3.91964,3.91781,4.14265,4.60341,5.06742,5.53082,5.52369} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.14843, 1.265, 1.25434  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.37546,1.37576,1.37586,1.37529,1.37504,1.35764,1.35729,1.35758,1.35772,1.35783,1.35783,1.14219,1.14359,1.14409,1.14459,1.14708,1.14923,1.14314,1.16,1.14193,1.15023},
    { 1.56824,1.56368,1.56201,1.56392,1.56583,1.59462,1.59549,1.59906,1.60077,1.59574,1.59637,1.36651,1.36617,1.36593,1.36582,1.36543,1.2672,1.31133,1.32233,1.26525,1.26535},
    { 1.75011,1.74252,1.7426,1.74263,1.74267,1.73537,1.73537,1.7354,1.73018,1.73045,1.72588,1.726,1.72209,1.71569,1.7131,1.71314,1.71091,1.7069,1.70355,1.80167,1.66994} };
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

