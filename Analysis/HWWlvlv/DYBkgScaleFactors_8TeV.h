static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.61308, 3.50698, 1.99591  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 25.3867,36.0035,41.5532,44.8199,48.2903,51.051,51.051,50.4018,57.7277,51.6553,45.7915,45.7915,18.5242,18.5242,11.5986,2.92566,2.16665,8.12944,7.37409,5.36337,4.14031},
    { 3.59534,4.63688,7.99957,7.60355,8.32775,8.89052,8.89052,9.92173,10.9402,10.0783,8.43504,8.43504,7.24509,7.24509,5.52062,5.51539,4.52333,9.64039,9.63005,6.45008,4.21483},
    { 1.72083,2.20909,3.19652,3.34021,2.93548,2.39539,2.39539,2.78767,2.78448,3.53643,3.80284,3.80284,3.2231,3.2231,3.75676,5.26026,4.91395,6.08612,4.31848,5.19081,2.25656} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.62864, 1.43304, 1.20412  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.66,1.65579,1.65485,1.67689,1.66619,1.68644,1.68644,1.82755,1.79197,1.84575,1.50633,1.50633,1.47602,1.47602,1.87481,2.82915,3.30312,1.63335,2.10479,2.04447,2.16927},
    { 1.42969,1.39565,1.33307,1.36201,1.77371,1.75031,1.75031,1.74783,1.60865,1.99374,2.03891,2.03891,2.0212,2.0212,2.07335,2.07258,2.08497,1.41471,1.29596,2.02716,2.05488},
    { 1.64034,1.61483,1.58546,1.54179,1.57315,1.64366,1.64366,1.52593,1.55379,1.71077,1.73683,1.73683,1.79163,1.79163,2.02907,2.15851,2.12858,2.0852,2.08395,2.0723,2.20478} };
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

static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 6.61308, 3.50698, 1.99591  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,120.943,41.1068},
    { 45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,45.731,26.4772},
    { 13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,13.0666,8.65576} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.62864, 1.43304, 1.20412  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,1.84456,2.01029},
    { 1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,1.59956,2.00495},
    { 1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.48456,1.81016} };
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

