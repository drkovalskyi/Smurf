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

Double_t ZttScaleFactor(Int_t nvtx, Int_t period, Double_t scale1fb) {
  if(period == 3) return 0.0165173*scale1fb;
  return 0.019*nVtxScaleFactorZtt(nvtx, period);
}

Double_t ZttScaleFactorKappa() {
  return 1.10;
}

Double_t WGstarScaleFactor() {
  return 1.60;
}

Double_t WGstarScaleFactorSyst() {
  return 0.30;
}

Double_t WJetsMCScaleFactor() {
  return 2.00;
}

Double_t reweightBToA(Int_t nvtx){
  double mynvtx = TMath::Min((double)nvtx,29.499);
   if	 (mynvtx == 0) return  0.00000;
  else if(mynvtx == 1) return 15.05738;
  else if(mynvtx == 2) return  9.65345;
  else if(mynvtx == 3) return  6.00292;
  else if(mynvtx == 4) return  3.71893;
  else if(mynvtx == 5) return  2.38827;
  else if(mynvtx == 6) return  1.52199;
  else if(mynvtx == 7) return  0.98081;
  else if(mynvtx == 8) return  0.63142;
  else if(mynvtx == 9) return  0.41253;
  else if(mynvtx ==10) return  0.27211;
  else if(mynvtx ==11) return  0.18064;
  else if(mynvtx ==12) return  0.11943;
  else if(mynvtx ==13) return  0.08065;
  else if(mynvtx ==14) return  0.05402;
  else if(mynvtx ==15) return  0.03806;
  else if(mynvtx ==16) return  0.02427;
  else if(mynvtx ==17) return  0.01958;
  else if(mynvtx ==18) return  0.01334;
  else if(mynvtx >=19) return  0.00667;
  return 1.0;
}
