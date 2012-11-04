static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.11765,1.11765,1.11765,1.11765,1.11765,1.11765,1.11765,1.11764,1.11764,1.11761,1.12126,1.12126,1.10662,1.10662,1.10885,1.10983,1.10582,1.11017,1.1063},
    { 1.06318,1.06318,1.06318,1.06318,1.06318,1.06318,1.06318,1.06321,1.06315,1.06887,1.0809,1.0809,1.05416,1.05416,1.05536,1.0636,1.05472,1.05393,1.03726} };
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
    { 1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175,1.06175},
    { 0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652,0.926652} };
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
    { 1.07412,1.07412,1.07412,1.07412,1.07412,1.07412,1.07412,1.07412,1.07412,1.07407,1.07373,1.07373,1.07604,1.07604,1.07593,1.07602,1.07633,1.07652,1.07711,},
    { 1.14595,1.14595,1.14595,1.14595,1.14595,1.14595,1.14595,1.14594,1.14596,1.14536,1.14341,1.14341,1.15358,1.15358,1.15352,1.15245,1.1537,1.15424,1.15663} };
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
    { 1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193,1.06193},
    { 1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275,1.15275} };
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

