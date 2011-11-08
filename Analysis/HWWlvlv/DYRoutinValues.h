Double_t RoutinValue(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t RoutinValuesWWPreselection[3] = { 0.183769, 0.153971, 0.165033  };
  Double_t RoutinValuesHiggsSelection[3][12] = { 
    { 0.206,0.206,0.699591,0.699591,0.261051,0.683625,0.336968,0.200775,0.155963,0.175754,0.0226749,0.0795888},
    { 0.108559,0.108559,0.247001,0.247001,0.180431,0.443684,0.379921,0.273375,0.178415,0.126694,0.0680514,0.0726171},
    { 0.11342,0.11342,0.219115,0.219115,0.159056,0.354693,0.358721,0.295036,0.201365,0.158846,0.0846966,0.094967} };
  if(mH == 0) return RoutinValuesHiggsSelection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return RoutinValuesHiggsSelection[jetBin][massIndex];
  } else {
    return RoutinValuesWWPreselection[jetBin];
  }
}

Double_t RoutinStatError(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t RoutinStatErrorWWPreselection[3] = { 0.0238961, 0.0117904, 0.0196004  };
  Double_t RoutinStatErrorHiggsSelection[3][12] = { 
    { 0.0414708,0.0414708,0.148771,0.148771,0.105408,0.365431,0.251062,0.160073,0.0766924,0.0615089,0.0100519,0.0406822},
    { 0.0132021,0.0132021,0.0287995,0.0287995,0.0273856,0.0814199,0.0750884,0.0498571,0.0264324,0.0185563,0.00938745,0.0132094},
    { 0.0205428,0.0205428,0.0390112,0.0390112,0.0356017,0.0884655,0.0896351,0.0672736,0.0391928,0.029574,0.0168193,0.0222032} };
  if(mH == 0) return RoutinStatErrorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return RoutinStatErrorHiggsSelection[jetBin][massIndex];
  } else {
    return RoutinStatErrorWWPreselection[jetBin];
  }
}

Double_t RoutinSystError(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t RoutinSystErrorWWPreselection[3] = { 0.0965501, 0.0236541, 0.0238791  };
  Double_t RoutinSystErrorHiggsSelection[3][12] = { 
    { 0.12568,0.12568,0.405953,0.405953,0.23945,1.0719,0.639724,0.349116,0.187602,0.0806421,0.0289556,0.0460485},
    { 0.0430623,0.0430623,0.116687,0.116687,0.0582758,0.12937,0.125818,0.0909856,0.0596441,0.0345044,0.0153562,0.0161262},
    { 0.0268188,0.0268188,0.0473847,0.0473847,0.0659056,0.140244,0.152258,0.110428,0.0582835,0.0457305,0.0449969,0.0206421} };
  if(mH == 0) return RoutinSystErrorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return RoutinSystErrorHiggsSelection[jetBin][massIndex];
  } else {
    return RoutinSystErrorWWPreselection[jetBin];
  }
}

