static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 5.51023, 2.85191, 2.11109  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 13.0501,17.974,20.9027,23.0248,24.6685,26.7443,26.7443,27.4396,31.068,29.5663,29.4736,29.4736,18.4868,18.4868,11.5312,2.87135,2.12793,8.12077,7.37775,5.35196,4.11557},
    { 2.21755,2.81011,4.46016,4.37708,4.64597,4.93843,4.93843,5.52204,6.30857,6.09053,5.68431,5.68431,7.24986,7.24986,5.51577,5.51494,4.52715,9.6537,9.64104,6.45618,4.22035},
    { 1.27694,1.45582,2.18834,2.22465,2.20015,1.80191,1.80191,2.12214,2.09511,2.80005,2.79899,2.79899,3.22383,3.22383,3.75574,5.25903,4.91522,6.0891,4.32083,5.19304,2.25802} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.32653, 1.33608, 1.30357  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.33974,1.33207,1.32971,1.33247,1.33751,1.35318,1.35318,1.35422,1.35643,1.36431,1.36748,1.36748,1.40771,1.40771,1.63861,2.86101,3.1745,1.59403,1.55704,1.46681,1.67975},
    { 1.4688,1.44022,1.38911,1.39948,1.40439,1.40662,1.40662,1.40275,1.40263,1.41732,1.43858,1.43858,1.39804,1.39804,1.49209,1.49014,1.51633,1.38478,1.37071,1.38067,1.44988},
    { 1.69916,1.68648,1.65438,1.65142,1.68023,1.72416,1.72416,1.70357,1.76957,1.74105,1.74111,1.74111,1.7377,1.7377,2.00774,2.15852,1.89689,1.74253,1.70198,1.49378,1.81983} };
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

