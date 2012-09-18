// Pseudo-data fits on a sample with zero aTGC
// Signal only PDF
// LZ_GZ space
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  // setSigPdf_LZ_GZ();
  setSigPdf_LZ_KG();
  // setBkgPdf();
  // prepareForDataFits();
  // ww1DFits("smurf/qqww.root");
  // toy_wwATGC1DFits();
  toy_wwATGC1DLzKgFits();
}
