Double_t PDFgHHSystematics(Int_t mH) {
    Int_t mHiggs[35] = {110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600,700,800,900,1000};
    Double_t kappa[35] = {1.072,1.072,1.071,1.072,1.072,1.072,1.072,1.072,1.072,1.072,1.072,1.071,1.071,1.072,1.073,1.073,1.075,1.075,1.075,1.075,1.076,1.075,1.075,1.075,1.078,1.081,1.082,1.085,1.088,1.091,1.091,1.097,1.102,1.107,1.117};

    Int_t massIndex = -1;
    for (UInt_t m=0; m < 35 ; ++m) {
        if (mH == mHiggs[m]) massIndex = m;
    }
    assert(massIndex >= 0);
    return kappa[massIndex];
}
