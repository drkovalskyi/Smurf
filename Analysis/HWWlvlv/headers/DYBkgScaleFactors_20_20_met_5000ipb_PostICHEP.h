static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.63941, 3.53197, 1.87479  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 11.6403,12.7939,13.168,14.7557,15.5935,17.5842,17.5842,17.7756,18.7474,19.1871,19.1283,29.0371,12.5909,12.5909,9.32834,2.96232,1.77211,5.1244,6.47322,3.69342,4.00429},
    { 3.21704,4.28765,4.68528,4.11237,4.15785,2.87642,2.87642,3.93633,2.7632,3.1155,2.56236,21.1167,24.7052,24.7052,26.7979,23.6828,23.361,32.2053,29.9843,20.3801,14.3566},
    { 2.02624,2.02542,2.02406,2.0231,2.02256,2.20669,2.57778,2.57702,2.76046,2.75295,3.12446,3.12344,2.93243,3.1132,3.29468,3.47811,3.65509,4.0139,4.18398,3.08699,3.08147} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.30901, 1.13628, 1.23622  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.21329,1.21204,1.21363,1.21945,1.25091,1.3167,1.3167,1.33165,1.3877,1.38634,1.38677,1.68268,1.65745,1.65745,2.39703,2.79032,4.66687,1.93639,1.63412,1.84447,1.52724},
    { 2.07664,2.06766,2.06616,2.07861,2.08649,2.12425,2.12425,1.62771,1.75638,1.72898,1.8326,1.67799,1.54321,1.54321,1.3987,1.4395,1.41671,1.45013,1.39691,1.28692,1.77082},
    { 2.12361,2.12364,2.1237,2.12374,2.12376,2.12031,2.11479,2.11481,2.11265,2.11281,2.10912,2.10914,2.11464,2.11583,2.11365,2.11169,1.6741,1.66893,1.66632,1.84163,1.84174} };
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

