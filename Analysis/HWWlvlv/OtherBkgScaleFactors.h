
Double_t nVtxScaleFactorZtt(Int_t nvtx, Int_t period){
  double mynvtx = TMath::Min((double)nvtx,29.499);
  if    (period == 1){
     if    (mynvtx == 0) return    0.00000;
    else if(mynvtx == 1) return    0.05555;
    else if(mynvtx == 2) return    0.09257;
    else if(mynvtx == 3) return    0.14974;
    else if(mynvtx == 4) return    0.24219;
    else if(mynvtx == 5) return    0.38123;
    else if(mynvtx == 6) return    0.59133;
    else if(mynvtx == 7) return    0.93503;
    else if(mynvtx == 8) return    1.45941;
    else if(mynvtx == 9) return    2.26860;
    else if(mynvtx ==10) return    3.62841;
    else if(mynvtx ==11) return    5.31113;
    else if(mynvtx ==12) return    8.92745;
    else if(mynvtx ==13) return   11.87322;
    else if(mynvtx ==14) return   20.45361;
    else if(mynvtx >=15) return   26.80333;
  }
  else if(period == 2){
     if    (mynvtx == 0) return    0.00000;
    else if(mynvtx == 1) return    0.49871;
    else if(mynvtx == 2) return    0.51323;
    else if(mynvtx == 3) return    0.56341;
    else if(mynvtx == 4) return    0.61741;
    else if(mynvtx == 5) return    0.69943;
    else if(mynvtx == 6) return    0.81422;
    else if(mynvtx == 7) return    0.98780;
    else if(mynvtx == 8) return    1.24125;
    else if(mynvtx == 9) return    1.64361;
    else if(mynvtx ==10) return    2.34832;
    else if(mynvtx ==11) return    3.12128;
    else if(mynvtx ==12) return    4.87722;
    else if(mynvtx ==13) return    6.19816;
    else if(mynvtx ==14) return   10.39581;
    else if(mynvtx >=15) return   13.36622;
  }
  return 1.0;
}

Double_t ZttScaleFactor(Int_t nvtx, Int_t period) {
  return 0.019*nVtxScaleFactorZtt(nvtx, period);
}

Double_t ZttScaleFactorKappa() {
  return 1.10;
}

Double_t WGstarScaleFactor() {
  return 1.45;
}

Double_t WGstarScaleFactorSyst() {
  return 0.30;
}

Double_t WJetsMCScaleFactor() {
  return 2.00;
}
