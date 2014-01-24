static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.2445,1.25335,1.25894,1.12633,1.10157,1.10157,1.0971,1.09753,1.08726,1.09006},
    { 0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.834473,0.826315,0.915547,0.843932,0.843932,0.83182,0.827822,0.825131,0.801824} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
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
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][18] = { 
    { 1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.24436,1.11864,1.11864,1.11864,1.11864,1.11864,1.15663,1.13094},
    { 0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.845645,0.943397,0.943397,0.943397,0.943397,0.943397,0.896942,0.857766} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
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
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][18] = { 
    { 1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10069,1.09958,1.10238,1.10694,1.10694,1.10731,1.10759,1.10861,1.1088,},
    { 1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.27126,1.27308,1.26243,1.29751,1.29751,1.30165,1.30339,1.30343,1.3122} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
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
  Int_t mHiggs[18] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][18] = { 
    { 1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10116,1.10307,1.10307,1.10307,1.10307,1.10307,1.10544,1.11355},
    { 1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.26745,1.25503,1.25503,1.25503,1.25503,1.25503,1.27274,1.31027} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

