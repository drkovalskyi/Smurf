static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.184,1.18517,1.19678,1.12148,1.07847,1.07847,1.07408,1.06922,1.06722,1.06919,1.06405},
    { 0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898501,0.901723,0.914984,0.949162,0.925648,0.925648,0.920216,0.916234,0.924806,0.911336,0.894628} };
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
    { 1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.11549,1.11549,1.11549,1.11549,1.11549,1.12952,1.12714,1.14179},
    { 0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.898488,0.934726,0.934726,0.934726,0.934726,0.934726,0.999845,0.973619,0.943467} };
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
    { 1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08349,1.08187,1.083,1.08658,1.08658,1.08688,1.08741,1.08777,1.08809,1.08905,},
    { 1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18063,1.17606,1.18689,1.19768,1.19768,1.19885,1.1999,1.19687,1.19982,1.20398} };
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
    { 1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08419,1.08419,1.08419,1.08419,1.08419,1.08578,1.08922,1.09216},
    { 1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.1917,1.1917,1.1917,1.1917,1.1917,1.17695,1.19565,1.21298} };
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

