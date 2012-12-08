static Double_t TopVBFBkgScaleFactor(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgSF[31] = { 
0.811066,0.874913,0.874913,0.799659,1.03919,0.811796,0.783746,0.727536,0.727536,0.721623,0.692198,0.677879,0.543302,0.543302,0.913991,0.913991,1.02834,0.923313,1.13017,0.835928,0.698123,0.953945,1.10516,1.01651,1.09918,0.694384,0.849495,0.726843,0.653919,0.755674,0.949182};
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
1.19313,1.64392,1.64392,1.64392,1.54367,1.54367,1.54027,1.55354,1.55354,1.54793,1.54407,1.51029,1.55174,1.55174,1.50145,1.50145,1.45128,1.45101,1.36099,1.41037,1.4057,1.30951,1.30014,1.30085,1.3396,1.54803,1.53944,1.54277,1.67321,1.65126,1.65163};
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
0.938195,0.64058,0.64058,0.617545,1.15825,0.801719,0.757248,0.71628,0.71628,0.697326,0.687224,0.691344,0.401181,0.401181,0.37646,0.37646,0.367119,0.305115,0.846995,0.757541,0.724181,0.987158,1.36958,1.36463,1.78778,1.22779,1.65948,1.36625,1.1033,1.17134,1.56468};
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
1.24066,2.01881,2.01881,2.01881,1.69913,1.69913,1.69618,1.69618,1.69618,1.69618,1.6931,1.65538,1.8426,1.8426,1.97323,1.97323,1.97576,1.97314,1.50417,1.58943,1.50848,1.37576,1.34649,1.34682,1.36659,1.59996,1.54907,1.55342,1.70982,1.71058,1.71103};
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
0.618685,1.51147,1.51147,1.26426,1.01771,0.974423,0.962893,0.864683,0.864683,0.886304,0.801705,0.688662,0.817341,0.817341,1.96243,1.96243,2.19343,2.30333,1.84077,0.99471,0.656067,0.894179,0.575388,0.441805,0.265526,0.208866,0.120974,0.145857,0.191691,0.228532,0.265308};
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
1.26825,1.69554,1.69554,1.69554,1.72433,1.72433,1.733,1.78528,1.78528,1.76613,1.75795,1.75795,1.71903,1.71903,1.58196,1.58196,1.50644,1.50644,1.50727,1.51381,1.59445,1.4808,1.47754,1.47613,1.59223,1.7281,2.05861,2.00081,2.00081,2.00081,2.00081};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFSFBkgKappa[massIndex];
  } else { assert(0); }
  return 1.0;
}

