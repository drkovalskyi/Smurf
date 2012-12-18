static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.08845,1.08845,1.08845,1.08845,1.08845,1.08845,1.08845,1.08845,1.08854,1.09366,1.09177,1.09177,1.0754,1.0754,1.07536,1.07122,1.06911,1.06495,1.0683},
    { 0.925531,0.925531,0.925531,0.925531,0.925531,0.925531,0.925531,0.925532,0.925543,0.930231,0.935417,0.935417,0.939054,0.939054,0.938622,0.935158,0.935024,0.929679,0.925004} };
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
    { 1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593,1.15593},
    { 1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296,1.04296} };
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
    { 1.05899,1.05899,1.05899,1.05899,1.05899,1.05899,1.05899,1.05899,1.05898,1.0582,1.05805,1.05805,1.0577,1.0577,1.05771,1.05806,1.05838,1.05886,1.05907,},
    { 1.10498,1.10498,1.10498,1.10498,1.10498,1.10498,1.10498,1.10498,1.10498,1.10265,1.10102,1.10102,1.10232,1.10232,1.10234,1.10247,1.1026,1.10317,1.10391} };
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
    { 1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248,1.05248},
    { 1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997,1.09997} };
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

