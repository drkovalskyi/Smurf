static Double_t TopVBFBkgScaleFactor(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgSF[31] = { 
1.33196,3.82802,7.06264,6.92777,4.93032,4.06937,4.39601,5.20212,5.20212,5.06923,4.69447,4.12564,4.1404,4.33762,4.07794,4.46293,4.37258,3.04052,2.13261,2.06215,1.47791,1.3903,1.37731,1.3608,1.36079,1.34454,1.34454,1.34454,1.34332,1.34332,1.34177};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFBkgSF[massIndex];
  } else { assert(0); }
  return 1.0;
}
static Double_t TopVBFBkgScaleFactorKappa(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgKappa[31] = { 
1.42976,1.70513,1.63978,1.63978,1.63978,1.63978,1.58708,1.53773,1.53773,1.53814,1.58063,1.58071,1.57864,1.51816,1.51816,1.47297,1.47297,1.50256,1.50269,1.52521,1.52486,1.44896,1.42407,1.42426,1.42426,1.42933,1.42933,1.42933,1.42933,1.42933,1.42971};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFBkgKappa[massIndex];
  } else { assert(0); }
  return 1.0;
}

static Double_t TopVBFOFYield(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgOF[31] = { 
1.10763,0.223707,11.2985,10.6657,4.87381,3.4851,3.30958,2.39922,2.39922,2.39922,2.39922,1.99511,1.99511,1.75473,1.6167,1.6167,1.58195,1.14094,0.770511,0.701396,0.655065,1.07468,1.15835,1.13174,1.13174,1.129,1.129,1.129,1.12697,1.12697,1.12388};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFBkgOF[massIndex];
  } else { assert(0); }
  return 1.0;
}
static Double_t TopVBFOFYieldKappa(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFOFBkgKappa[31] = { 
1.66468,1.30086,2.03238,2.03238,2.03238,2.03238,2.03238,2.03156,2.03156,2.03156,2.03156,2.03156,2.03156,2.022,2.022,2.022,2.022,2.022,2.02215,2.07606,2.07606,1.7812,1.66141,1.66141,1.66141,1.66301,1.66301,1.66301,1.66301,1.66301,1.66452};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFOFBkgKappa[massIndex];
  } else { assert(0); }
  return 1.0;
}

static Double_t TopVBFSFYield(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgSF[31] = { 
1.93462,4.99428,4.99428,4.99428,4.99428,4.99428,6.4564,11.9019,11.9019,10.9716,9.80073,9.80073,9.80073,12.1008,12.1008,12.5776,12.2668,9.32751,7.48313,7.10634,3.07894,2.05099,1.96576,1.96454,1.96452,1.93462,1.93462,1.93462,1.93462,1.93462,1.93462};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFBkgSF[massIndex];
  } else { assert(0); }
  return 1.0;
}
static Double_t TopVBFSFYieldKappa(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFSFBkgKappa[31] = { 
1.43862,1.71533,1.71533,1.71533,1.71533,1.71533,1.59825,1.5779,1.5779,1.5786,1.64677,1.64677,1.64677,1.55898,1.55898,1.49401,1.49803,1.54597,1.54598,1.48854,1.48854,1.43123,1.43192,1.43219,1.4322,1.43862,1.43862,1.43862,1.43862,1.43862,1.43862};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFSFBkgKappa[massIndex];
  } else { assert(0); }
  return 1.0;
}
