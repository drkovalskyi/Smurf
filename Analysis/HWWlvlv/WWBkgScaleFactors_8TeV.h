static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18399,1.18517,1.19677,1.12147,1.07845,1.07845,1.07406,1.0692,1.0672,1.06918,1.06403},
    { 0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898475,0.901698,0.914958,0.949137,0.92562,0.92562,0.920189,0.916206,0.924778,0.911309,0.8946} };
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
    { 1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.18389,1.11547,1.11547,1.11547,1.11547,1.11547,1.12951,1.12713,1.14178},
    { 0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.898463,0.934702,0.934702,0.934702,0.934702,0.934702,0.999817,0.973597,0.943456} };
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
    { 1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08349,1.08187,1.083,1.08658,1.08658,1.08688,1.08741,1.08777,1.08809,1.08906,},
    { 1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18161,1.18064,1.17607,1.18689,1.19769,1.19769,1.19886,1.1999,1.19688,1.19983,1.20399} };
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
    { 1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.08392,1.0842,1.0842,1.0842,1.0842,1.0842,1.08578,1.08922,1.09216},
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

