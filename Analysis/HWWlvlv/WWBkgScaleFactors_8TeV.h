static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.10788,1.10788,1.10788,1.10788,1.10788,1.10788,1.10788,1.10788,1.10797,1.11296,1.11236,1.11236,1.08295,1.08295,1.08274,1.07845,1.0764,1.07246,1.07548},
    { 0.873788,0.873788,0.873788,0.873788,0.873788,0.873788,0.873788,0.873788,0.873798,0.878975,0.883721,0.883721,0.87971,0.87971,0.879221,0.880722,0.88149,0.875204,0.870832} };
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
    { 1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364,1.17364},
    { 1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313,1.0313} };
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
    { 1.05821,1.05821,1.05821,1.05821,1.05821,1.05821,1.05821,1.05821,1.0582,1.05735,1.057,1.057,1.05934,1.05934,1.05937,1.05969,1.06002,1.06051,1.06074,},
    { 1.12612,1.12612,1.12612,1.12612,1.12612,1.12612,1.12612,1.12612,1.12612,1.12319,1.12121,1.12121,1.1234,1.1234,1.12344,1.12292,1.12298,1.12383,1.12474} };
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
    { 1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467,1.05467},
    { 1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216,1.11216} };
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

