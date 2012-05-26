Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.87405, 4.34857, 1.93195  };
  Double_t DYBkgScaleFactorHiggsSelection[3][20] = { 
    { 2.08701,1.90078,1.75245,1.06284,0.946468,1.42348,1.42348,0.944581,2.54336,2.25984,3.03866,4.11355,0.386191,0.57775,0.222977,0.156839,0.25182,0.955667,0.718017,1.45181},
    { 1.05686,1.00542,0.922407,1.10728,0.766804,0.887547,0.887547,1.85987,2.03567,1.92799,1.05525,3.38554,3.29051,3.68798,3.00994,2.44858,4.21775,3.01319,2.4733,1.42494},
    { 2.5341,3.07618,4.16976,4.37302,3.74306,4.25097,4.25097,3.85437,3.95711,3.73262,3.10781,3.10781,3.18907,2.18202,1.83206,2.27018,4.03328,3.73492,3.30268,3.40302} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 20 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.34384, 1.29294, 1.17997  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][20] = { 
    { 1.69288,1.80565,1.92072,2.78064,3.33587,3.05405,3.05405,4.69577,2.638,2.74563,2.03131,1.76299,3.20597,3.27949,3.09349,3.07967,3.96443,1.90403,1.71924,1.65491},
    { 1.93316,1.9582,2.00608,1.98544,2.16863,2.11934,2.11934,1.9253,1.9156,1.9351,2.09155,1.78377,1.74464,1.64561,1.64932,1.67705,1.62698,1.62801,1.48297,2.36516},
    { 1.46114,1.44675,1.42831,1.39075,1.38089,1.35728,1.35728,1.37746,1.389,1.40078,1.44212,1.44212,1.51803,1.94467,1.91022,2.07751,1.96989,1.96585,1.51077,1.42409} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 20 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

