static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.10464,1.10464,1.10464,1.10464,1.10464,1.10464,1.10464,1.10464,1.10473,1.10962,1.10883,1.10883,1.07946,1.07946,1.07924,1.07507,1.0729,1.06873,1.0717},
    { 0.859665,0.859665,0.859665,0.859665,0.859665,0.859665,0.859665,0.859665,0.859675,0.864517,0.86868,0.86868,0.864577,0.864577,0.864106,0.86671,0.868389,0.862143,0.857489} };
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
    { 1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236,1.16236},
    { 1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668,1.00668} };
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
    { 1.05917,1.05917,1.05917,1.05917,1.05917,1.05917,1.05917,1.05917,1.05917,1.05828,1.05793,1.05793,1.06024,1.06024,1.06026,1.06057,1.06092,1.06145,1.06169,},
    { 1.13463,1.13463,1.13463,1.13463,1.13463,1.13463,1.13463,1.13463,1.13463,1.13145,1.12936,1.12936,1.13131,1.13131,1.13134,1.13033,1.13006,1.13095,1.13203} };
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
    { 1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799,1.05799},
    { 1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415,1.12415} };
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

