Double_t HiggsSignalQCDScaleKappa(TString nuisanceName, Int_t mH, Int_t jetBin) {
    if(mH == 110) mH = 115;
    if(mH == 118) mH = 120;
    if(mH == 122) mH = 120;
    if(mH == 124) mH = 120;
    if(mH == 125) mH = 120;
    if(mH == 126) mH = 130;
    if(mH == 128) mH = 130;
    if(mH == 135) mH = 140;
    if(mH == 145) mH = 140;
    if(mH == 155) mH = 150;
    if(mH == 700) mH = 600;
    if(mH == 800) mH = 600;
    if(mH == 900) mH = 600;
    if(mH ==1000) mH = 600;
    Int_t mHiggs[21] = {115,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
    Double_t kappa0_ggH[21] = {1.16,1.16,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.16,1.16,1.17,1.20,1.17,1.19,1.22,1.24,1.25};
    Double_t kappa0_ggH1in[21] = {0.92,0.92,0.91,0.91,0.90,0.90,0.89,0.89,0.88,0.88,0.88,0.87,0.87,0.86,0.84,0.83,0.82,0.81,0.80,0.78,0.78};
    Double_t kappa1_ggH1in[21] = {1.28,1.28,1.29,1.28,1.28,1.28,1.28,1.28,1.28,1.27,1.27,1.27,1.27,1.27,1.27,1.27,1.26,1.26,1.25,1.26,1.26};
    Double_t kappa1_ggH2in[21] = {0.97,0.97,0.98,0.97,0.97,0.96,0.96,0.96,0.96,0.96,0.96,0.96,0.95,0.96,0.95,0.95,0.95,0.95,0.95,0.95,0.94};
    Double_t kappa2_ggH2in[21] = {1.15,1.12,1.12,1.13,1.12,1.20,1.18,1.17,1.17,1.20,1.17,1.17,1.20,1.17,1.20,1.21,1.20,1.20,1.17,1.19,1.19};
    Int_t massIndex = -1;
    for (UInt_t m=0; m < 21 ; ++m) {
        if (mH == mHiggs[m]) massIndex = m;
    }
    assert(massIndex >= 0);
    if (nuisanceName.Data() == "QCDscale_ggH" && jetBin == 0) {
        return kappa0_ggH[massIndex];
    } else if (nuisanceName.Data() == "QCDscale_ggH1in" && jetBin == 0) {
        return kappa0_ggH1in[massIndex];
    } else if (nuisanceName.Data() == "QCDscale_ggH1in" && jetBin == 1) {
        return kappa1_ggH1in[massIndex];
    } else if (nuisanceName.Data() == "QCDscale_ggH2in" && jetBin == 1) {
        return kappa1_ggH2in[massIndex];
    } else if (nuisanceName.Data() == "QCDscale_ggH2in" && jetBin == 2) {
        return kappa2_ggH2in[massIndex];
    } else { 
        return 1.0;
    } 
    return 0; 
}
