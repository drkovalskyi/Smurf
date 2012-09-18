// data fits
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  setSigPdf_LZ_GZ();
  setBkgPdf();
  fitData();
  yieldDataFit();
  c12->Print("lz_dgz_contourplot.pdf");
  c12->Print("lz_dgz_contourplot.root");
  cc1->Print("correlations.root");
  cc2->Print("lz_dgz_pdf_data.pdf");
  cc2->Print("lz_dgz_pdf_data.root");
}
