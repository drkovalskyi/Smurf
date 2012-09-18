// data fits
{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  setSigPdf_LZ_KG();
  setBkgPdf();
  fitData();
  yieldDataFit();
  c12->Print("lz_dkg_contourplot.pdf");
  c12->Print("lz_dkg_contourplot.root");
  cc2->Print("lz_dkg_pdf_data.pdf");
  cc2->Print("lz_dkg_pdf_data.root");
}
