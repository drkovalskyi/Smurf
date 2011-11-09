Double_t HiggsSignalPSUESystematics(Int_t mH, Int_t jetBin) {
    if(mH == 110) mH = 115;
    Int_t mHiggs[21] = {115,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
    Double_t kappa0[21] = {0.941,0.940,0.937,0.941,0.942,0.943,0.946,0.947,0.948,0.952,0.948,0.950,0.950,0.955,0.958,0.964,0.966,0.954,0.946,0.931,0.920};
    Double_t kappa1[21] = {1.128,1.110,1.113,1.104,1.093,1.084,1.075,1.067,1.068,1.055,1.061,1.061,1.061,1.058,1.061,1.068,1.078,1.092,1.102,1.117,1.121};
    Double_t kappa2[21] = {1.212,1.293,1.237,1.168,1.156,1.138,1.108,1.092,1.083,1.059,1.042,1.028,1.024,0.990,0.942,0.889,0.856,0.864,0.868,0.861,0.872};

    Int_t massIndex = -1;
    for (UInt_t m=0; m < 21 ; ++m) {
        if (mH == mHiggs[m]) massIndex = m;
    }
    assert(massIndex >= 0);
    if (jetBin == 0) {
        return kappa0[massIndex];
    } else if (jetBin == 1) {
        return kappa1[massIndex];
    } else if (jetBin == 2) {
        return kappa2[massIndex];
    }
    return 0; 
}
