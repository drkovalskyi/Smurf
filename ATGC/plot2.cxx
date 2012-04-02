// Pseudo-data fits on a sample with zero aTGC
// Signal only PDF
// LZ_KG space
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  setSigPdf_LZ_KG();
  // ww1DFits();
  ww1DFits("samples/qqww.root",true);
}
