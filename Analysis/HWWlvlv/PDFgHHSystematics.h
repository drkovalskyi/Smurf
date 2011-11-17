Double_t PDFgHHSystematics(Int_t mH) {
    Int_t mHiggs[22] = {110,115,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
    Double_t kappa[22] = {1.078,1.078,1.078,1.079,1.080,1.081,1.081,1.082,1.082,1.082,1.082,1.082,1.083,1.083,1.084,1.087,1.091,1.092,1.095,1.098,1.101,1.105};

    Int_t massIndex = -1;
    for (UInt_t m=0; m < 22 ; ++m) {
        if (mH == mHiggs[m]) massIndex = m;
    }
    assert(massIndex >= 0);
    return kappa[massIndex];
}
