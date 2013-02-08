Double_t ZttScaleFactor(Int_t period, Double_t scale1fb, Double_t trg, Double_t eff) {
  if(period != 2 && period != 3) assert(0);
  return 0.00573547*scale1fb*trg*eff;
}

Double_t ZttScaleFactorKappa() {
  return 1.10;
}

Double_t WGstarScaleFactor(Int_t type, Double_t met) {
  if(type == SmurfTree::ee && met > 130) return 0.0;
  return 1.50;
}

Double_t WGstarScaleFactorSyst() {
  return 0.40;
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

Double_t TopMCScaleFactor_VHqqll(int nJetsType) {
  if     (nJetsType == 2) return 4.670;
  else if(nJetsType == 1) return 1.703;
  else                    assert(0);
  return 1.0;
}

Double_t TopMCScaleFactor_VHqqll_Kappa(int nJetsType) {
  if     (nJetsType == 2) return 1.298;
  else if(nJetsType == 1) return 1.334;
  else                    assert(0);
  return 1.0;
}

Double_t WrongChargeScaleFactor_VHqqll() {
  return 0.9000;
}

Double_t WrongChargeScaleFactor_VHqqll_Kappa() {
  return 1.0346;
}

Double_t DYMCScaleFactor_VHqqll(int nJetsType) {
  if     (nJetsType == 2) return 1.700;
  else if(nJetsType == 1) return 1.400;
  else                    assert(0);
  return 1.0;
}

Double_t DYMCScaleFactor_VHqqll_Kappa(int nJetsType) {
  if     (nJetsType == 2) return 1.400;
  else if(nJetsType == 1) return 1.400;
  else                    assert(0);
  return 1.0;
}

Double_t WWVBFScaleFactor() {
  return 2.0;
}
