{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  ww_expected /= 5; // don't have enough events to do full stat
  setDefaults();
  setSigPdf_LZ_GZ();
  // setSigPdf_LZ_KG();
  // setBkgPdf();
  // sensitivity();
  // ww1DFits();
  wwATGC1DFits();
  // wwATGC1DLzKgFits();
  // fitData();
  // fitTop();
}
