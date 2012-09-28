static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.22779,1.22779,1.22779,1.22779,1.22779,1.22779,1.22779,1.22779,1.22783,1.23062,1.2277,1.2277,1.18554,1.18554,1.18561,1.18146,1.18259,1.17879,1.17869},
    { 0.950422,0.950422,0.950422,0.950422,0.950422,0.950422,0.950422,0.950422,0.950452,0.96038,0.955367,0.955367,0.953374,0.953374,0.951422,0.952295,0.954354,0.950802,0.945048} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

static Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047,1.33047},
    { 1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283,1.17283} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

static Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][19] = { 
    { 1.06304,1.06304,1.06304,1.06304,1.06304,1.06304,1.06304,1.06304,1.06304,1.06231,1.06205,1.06205,1.06525,1.06525,1.06528,1.06566,1.06586,1.06644,1.06687,},
    { 1.13346,1.13346,1.13346,1.13346,1.13346,1.13346,1.13346,1.13346,1.13346,1.13043,1.12986,1.12986,1.13244,1.13244,1.13268,1.13192,1.13165,1.13198,1.13321} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

static Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][19] = { 
    { 1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722,1.05722},
    { 1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882,1.11882} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 19 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

