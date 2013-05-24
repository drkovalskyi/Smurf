static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};
  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { 
    { 1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.0961,1.09619,1.09917,1.09507,1.09507,1.07438,1.07438,1.0743,1.06997,1.06776,1.06346,1.06656},
    { 0.923526,0.923526,0.923526,0.923526,0.923526,0.923526,0.923526,0.923527,0.923537,0.924018,0.92448,0.92448,0.920462,0.920462,0.920084,0.916269,0.915989,0.910471,0.905591} };
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
    { 1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038,1.19038},
    { 1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844,1.07844} };
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
    { 1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.0589,1.05889,1.05847,1.05858,1.05858,1.05862,1.05862,1.05864,1.059,1.05933,1.05983,1.06005,},
    { 1.10408,1.10408,1.10408,1.10408,1.10408,1.10408,1.10408,1.10408,1.10408,1.10299,1.10251,1.10251,1.10514,1.10514,1.10516,1.10539,1.10554,1.10617,1.10696} };
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
    { 1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864,1.04864},
    { 1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143,1.09143} };
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

