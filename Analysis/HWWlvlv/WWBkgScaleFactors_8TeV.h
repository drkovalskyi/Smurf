static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.20693,1.20264,1.21626,1.15204,1.13442,1.13442,1.13175,1.1301,1.12523,1.12622,1.11102},
    { 1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01216,1.00792,1.03472,0.949683,0.949683,0.941703,0.932443,0.938111,0.914689,0.889615} };
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
    { 1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.2068,1.14871,1.14871,1.14871,1.14871,1.14871,1.1613,1.14682,1.18541},
    { 1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.01493,1.04803,1.04803,1.04803,1.04803,1.04803,1.08367,0.991001,0.973591} };
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
    { 1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09196,1.09014,1.08923,1.09192,1.09192,1.09212,1.09252,1.09305,1.09334,1.09504,},
    { 1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.19056,1.19075,1.19883,1.22629,1.22629,1.22817,1.23058,1.22807,1.23391,1.24068} };
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
    { 1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.09173,1.08984,1.08984,1.08984,1.08984,1.08984,1.09275,1.09897,1.10042},
    { 1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.18983,1.19631,1.19631,1.19631,1.19631,1.19631,1.1911,1.22677,1.24443} };
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

