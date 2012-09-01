// Pseudo-data experiments using aTGC samples
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  ww_expected /= 5; // don't have enough events to do full stat
  setDefaults();
  setSigPdf_LZ_GZ();
  wwATGC1DFits();
}
