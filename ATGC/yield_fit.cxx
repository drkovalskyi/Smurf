{
  cpdf = new RooAddPdf("cpdf","combined pdf",RooArgList(*cSigPdf,*cBkgPdf));
  x_par->setConstant(0);
  y_par->setConstant(0);
  x_par->setVal(0);
  y_par->setVal(0);
  // cpdf->fitTo(*ds_data_pt1, RooFit::Minos());                                                                                                                                
  RooAbsReal* nll = cpdf->createNLL(*glb_data,RooFit::Extended(),RooFit::Constrain(RooArgSet(*n_top,*n_wjets,*n_wz,*n_zz,*n_dy)));
  // RooMinuit m(*nll);                                                                                                                                                         
  RooMinimizer m(*nll);
  m.setMinimizerType("Minuit2");
  m.setErrorLevel(0.5); //68% C.L.                                                                                                                                              
  m.migrad();
  m.hesse();
  m.minos();

  // Plot
  RooPlot* frame = var_pt1.frame(Title("Leading lepton Pt distribution in Data")) ;
  glb_data->plotOn(frame) ;
  cpdf->plotOn(frame) ;

}
