static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[21] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 3.93946, 3.57574, 2.00467  };
  Double_t DYBkgScaleFactorHiggsSelection[3][21] = { 
    { 8.49876,9.40819,9.28269,10.7997,11.5339,13.1146,13.1146,11.7585,12.8061,13.3285,13.2843,20.1848,8.99832,8.99832,5.14735,0.937346,1.14004,5.05158,5.64501,2.58441,2.98118},
    { 2.07267,2.60918,3.10622,2.90964,2.5991,2.66564,2.66564,4.04534,2.71638,2.61962,2.05102,13.4367,14.9688,14.9688,15.7318,14.1495,14.4711,18.486,16.4745,11.7597,8.50309},
    { 2.05776,2.05705,2.0558,2.05529,2.05477,2.28425,2.28424,2.28363,2.51147,2.50556,2.7355,2.73474,2.73386,3.19242,3.41937,3.41779,3.64019,4.09563,4.32431,4.5523,4.54534} };
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
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.28283, 1.19702, 1.2872  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][21] = { 
    { 1.23677,1.23632,1.24184,1.24832,1.28439,1.34581,1.34581,1.38545,1.4401,1.43553,1.43622,1.65378,1.76277,1.76277,2.72878,2.17079,5.53382,1.96333,1.61557,1.96302,1.55648},
    { 2.10906,2.09703,2.08998,2.10833,2.13236,2.14192,2.14192,1.64861,1.78546,1.80065,1.95906,1.70013,1.58015,1.58015,1.41888,1.47568,1.44746,1.49124,1.43473,1.24504,1.90841},
    { 1.76487,1.76492,1.76501,1.76505,1.76508,1.75761,1.75761,1.75765,1.75156,1.75185,1.74661,1.74664,1.74668,1.73842,1.73516,1.7352,1.73244,1.72757,1.7255,1.78429,1.64641} };
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

