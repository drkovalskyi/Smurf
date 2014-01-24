static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.09145,1.09145,1.09145,1.09145,1.09145,1.09145,1.09145,1.09145,1.09154,1.09432,1.09007,1.09007,1.07265,1.07265,1.07258,1.06957,1.06728,1.06288,1.06587},
    { 0.933518,0.933518,0.933518,0.933518,0.933518,0.933518,0.933518,0.933518,0.933528,0.933977,0.935819,0.935819,0.929446,0.929446,0.928979,0.926997,0.926689,0.921152,0.916185} };
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
    { 1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949,1.1949},
    { 1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759,1.08759} };
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
    { 1.05263,1.05263,1.05263,1.05263,1.05263,1.05263,1.05263,1.05263,1.05262,1.0522,1.05226,1.05226,1.0515,1.0515,1.05151,1.05174,1.05201,1.05244,1.05263,},
    { 1.10308,1.10308,1.10308,1.10308,1.10308,1.10308,1.10308,1.10308,1.10308,1.10199,1.10134,1.10134,1.10423,1.10423,1.10426,1.10426,1.10442,1.10504,1.10582} };
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
    { 1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482,1.04482},
    { 1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079,1.09079} };
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

