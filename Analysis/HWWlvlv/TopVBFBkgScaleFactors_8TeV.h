static Double_t TopVBFBkgScaleFactor(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgSF[31] = { 
0.921694,0.711883,1.58393,1.54286,1.50623,1.38614,1.46438,1.9186,1.88,1.83402,1.67754,1.5876,1.45405,1.27996,1.20154,1.2881,1.21417,1.03164,0.951844,0.931264,0.838354,0.981753,0.993205,0.982243,0.961054,0.946632,0.94426,0.944001,0.943915,0.944009,0.938259};
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
1.38046,1.69676,1.68349,1.68291,1.68229,1.68009,1.62318,1.55547,1.55514,1.55579,1.58868,1.58839,1.58307,1.54299,1.53977,1.49645,1.50574,1.54458,1.54296,1.50983,1.52254,1.39318,1.37715,1.37708,1.37698,1.38096,1.38096,1.38096,1.38089,1.3808,1.38123};
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
0.79391,0.0199048,1.51602,1.47344,1.43043,1.30369,1.26506,1.09424,1.07975,1.07221,1.06029,1.08438,1.02788,0.750231,0.699312,0.693273,0.664661,0.623325,0.576794,0.332953,0.315061,0.793953,0.843512,0.841723,0.819877,0.814635,0.811233,0.810856,0.810856,0.811293,0.809175};
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
1.57521,1.08501,2.02452,2.02324,2.02324,2.02088,2.0206,2.01029,2.00864,2.00792,2.00792,2.00476,1.99891,2.037,2.02965,2.0138,2.01213,2.00761,2.00502,2.49834,2.48921,1.63095,1.5755,1.57488,1.57488,1.57594,1.57594,1.57594,1.57594,1.57566,1.57702};
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
1.19097,2.09138,1.72679,1.69295,1.67238,1.56568,1.89692,3.37586,3.27337,3.10988,2.67025,2.45183,2.10855,2.22172,2.11409,2.34304,2.2451,1.78445,1.63371,1.90322,1.70895,1.34529,1.30333,1.27302,1.25336,1.22559,1.22547,1.22547,1.22465,1.22413,1.20997};
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
1.44155,1.70271,1.70075,1.69884,1.69563,1.69216,1.57827,1.61313,1.61312,1.6146,1.67178,1.67119,1.66795,1.58502,1.58202,1.52046,1.52482,1.59202,1.59148,1.48799,1.48658,1.43294,1.43382,1.43411,1.43388,1.44187,1.44187,1.44187,1.44187,1.44187,1.44187};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFSFBkgKappa[massIndex];
  } else { assert(0); }
  return 1.0;
}

