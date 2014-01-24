static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.07485,1.07485,1.07485,1.07485,1.07485,1.07485,1.07485,1.07484,1.07484,1.07531,1.08159,1.08159,1.0795,1.0795,1.08259,1.083,1.07857,1.07971,1.0749},
    { 1.00399,1.00399,1.00399,1.00399,1.00399,1.00399,1.00399,1.00402,1.00399,1.01437,1.02983,1.02983,1.03255,1.03255,1.03399,1.04248,1.03395,1.03329,1.01711} };
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
    { 1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034,1.07034},
    { 0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421,0.938421} };
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
    { 1.07539,1.07539,1.07539,1.07539,1.07539,1.07539,1.07539,1.07539,1.07539,1.07532,1.07471,1.07471,1.07654,1.07654,1.07632,1.07647,1.07682,1.07717,1.07791,},
    { 1.15132,1.15132,1.15132,1.15132,1.15132,1.15132,1.15132,1.15131,1.15132,1.14977,1.14728,1.14728,1.15373,1.15373,1.15366,1.15259,1.15385,1.15442,1.15679} };
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
    { 1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771,1.05771},
    { 1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675,1.14675} };
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

