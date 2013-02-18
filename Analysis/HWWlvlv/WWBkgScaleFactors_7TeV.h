static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07634,1.07634,1.07634,1.07634,1.07634,1.07634,1.07634,1.07634,1.07634,1.0768,1.08309,1.08309,1.08086,1.08086,1.08396,1.08437,1.07995,1.08109,1.0763},
    { 1.01664,1.01664,1.01664,1.01664,1.01664,1.01664,1.01664,1.01667,1.01663,1.027,1.04242,1.04242,1.04531,1.04531,1.04675,1.05526,1.04674,1.04614,1.02996} };
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
    { 1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184,1.07184},
    { 0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248,0.950248} };
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
    { 1.07517,1.07517,1.07517,1.07517,1.07517,1.07517,1.07517,1.07517,1.07517,1.0751,1.07449,1.07449,1.07634,1.07634,1.07612,1.07627,1.07662,1.07697,1.0777,},
    { 1.14914,1.14914,1.14914,1.14914,1.14914,1.14914,1.14914,1.14913,1.14914,1.14764,1.14521,1.14521,1.15155,1.15155,1.15148,1.15044,1.15166,1.15221,1.15452} };
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
    { 1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575,1.0575},
    { 1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452,1.14452} };
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

