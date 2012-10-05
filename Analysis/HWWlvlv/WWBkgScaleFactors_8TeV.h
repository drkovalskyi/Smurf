static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.13868,1.13868,1.13868,1.13868,1.13868,1.13868,1.13868,1.13868,1.13876,1.14172,1.13606,1.13606,1.10207,1.10207,1.10106,1.09681,1.09632,1.09352,1.09507},
    { 0.858159,0.858159,0.858159,0.858159,0.858159,0.858159,0.858159,0.858159,0.858171,0.868351,0.865085,0.865085,0.858239,0.858239,0.85691,0.857716,0.858292,0.856375,0.849614} };
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
    { 1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735,1.22735},
    { 1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848,1.06848} };
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
    { 1.06087,1.06087,1.06087,1.06087,1.06087,1.06087,1.06087,1.06087,1.06086,1.06014,1.05994,1.05994,1.06228,1.06228,1.06234,1.06269,1.06295,1.06342,1.06373,},
    { 1.14222,1.14222,1.14222,1.14222,1.14222,1.14222,1.14222,1.14222,1.14222,1.13864,1.13752,1.13752,1.13972,1.13972,1.13989,1.13893,1.13865,1.13887,1.14044} };
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
    { 1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674,1.05674},
    { 1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361,1.12361} };
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

