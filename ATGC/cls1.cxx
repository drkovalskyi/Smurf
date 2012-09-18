// Pseudo-data fits on a sample with zero aTGC
// Signal only PDF
// LZ_GZ space
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  setSigPdf_LZ_GZ();
  // setSigPdf_LZ_KG();
  setBkgPdf();
  prepareForDataFits();
  // toy_fullATGC1DFits();
  // toy_fullATGC1DLzKgFits();
//   x_par->setVal(0);
//   y_par->setVal(0);
//   makeCard("test");
  makeGrid();
}
