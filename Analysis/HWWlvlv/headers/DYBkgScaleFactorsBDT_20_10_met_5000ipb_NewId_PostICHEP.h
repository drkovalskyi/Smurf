static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 4.33919, 4.0774, 2.04371  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 48.4224,49.7045,50.2103,51.5031,52.7516,57.5051,58.9687,59.3379,59.0837,59.7927,59.7921,80.3989,81.1108,81.3952,81.324,80.7866,82.4461,85.4162,86.9596,95.6707,93.1892},
    { 11.1937,12.7531,13.1475,13.5383,13.5808,15.2655,15.1251,16.2901,16.3232,17.4962,17.5432,74.6896,80.2591,84.2077,86.5331,91.8345,99.11,109.5,115.959,121.614,122.477},
    { 2.35342,2.56863,2.5682,2.567,2.56622,2.78049,3.21208,3.2112,3.42712,3.41773,3.63364,3.62998,3.62349,4.04885,4.26116,4.47395,4.68297,5.10151,5.30124,4.42188,4.41471} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.11459, 1.22699, 1.2628  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.4192,1.41925,1.41924,1.41905,1.41894,1.44783,1.4476,1.44766,1.44774,1.44768,1.44768,1.12226,1.12277,1.12335,1.12381,1.12522,1.15,1.14448,1.14092,1.11248,1.1145},
    { 1.52485,1.51993,1.51922,1.51973,1.52365,1.53612,1.53689,1.53556,1.53702,1.53602,1.53951,1.30582,1.30551,1.30535,1.30529,1.30503,1.21151,1.26346,1.27639,1.22813,1.22771},
    { 2.10838,2.10484,2.10486,2.10489,2.10492,2.10194,2.09711,2.09712,2.09516,2.09532,2.09356,2.09362,2.09704,2.09614,2.0944,2.09283,1.81306,1.8096,1.80785,1.80517,1.80526} };
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

