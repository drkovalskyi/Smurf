static Double_t TopVBFBkgScaleFactor(Int_t mH) {
  Int_t mHiggs[31] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800};
  Double_t TopVBFBkgSF[31] = { 
0.786126,0.0244211,0.0244211,0.0287184,0.0317483,0.0317217,0.905396,0.883289,0.883289,0.788989,0.793406,0.762665,0.706203,0.706203,1.6709,1.6709,1.81132,0.0159415,0.779477,0.403413,0.335134,0.246432,0.220865,0.28771,0.339341,0.478699,0.00332012,0.00342416,0.0036577,0.00362687,0.00761899};
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
1.34942,1.12665,1.12665,1.08001,1.09507,1.0928,1.98254,1.98302,1.98302,1.98413,1.98211,1.98065,1.98192,1.98192,2.04119,2.04119,2.04319,1.12279,1.98364,1.97798,2.06257,2.02502,2.02983,2.00936,2.01599,2.00696,1.09555,1.0952,1.09473,1.09404,1.07911};
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
0.832702,0.025263,0.025263,0.0223989,0.0333141,0.0291275,0.0418583,0.0393408,0.0393408,0.0325301,0.0367336,0.0366963,0.0325807,0.0325807,0.0306458,0.0306458,0.0269507,0.0219368,1.18689,0.63535,0.553061,0.325682,0.297513,0.371299,0.439538,0.671424,0.00216506,0.0030629,0.00394743,0.000188397,0.000502456};
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
1.38171,1.13058,1.13058,1.13013,1.14414,1.14414,1.11052,1.10954,1.10954,1.11384,1.1189,1.12028,1.11945,1.11945,1.13373,1.13373,1.12624,1.12681,1.98769,1.98424,1.97958,2.02212,2.02841,2.03644,2.03829,2.01129,1.14629,1.1527,1.1527,1.07982,1.12801};
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
0.705814,0.026884,0.026884,0.0401714,0.0320759,0.0378589,2.09071,2.09079,2.09079,1.97193,1.87698,1.78859,1.64614,1.64614,4.3345,4.3345,5.13518,0.016516,0.0106362,0.010083,0.00813074,0.00556603,0.00379035,0.035579,0.0388519,0.00930359,0.0135932,0.0129779,0.0147796,0.0190866,0.0232263};
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
1.7075,1.11678,1.11678,1.06305,1.06305,1.06703,2.00677,2.00674,2.00674,2.00674,2.00674,2.00614,2.00614,2.00614,2.0509,2.0509,2.05116,1.32404,1.32404,1.22581,1.30606,1.12801,1.09225,1.12217,1.12464,1.07812,1.07763,1.07754,1.0784,1.0784,1.0784};
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 31 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return TopVBFSFBkgKappa[massIndex];
  } else { assert(0); }
  return 1.0;
}

