{
  gSystem->CompileMacro("RooATGCPdf.C","k");
  gSystem->CompileMacro("wwfit.C","k");
  setDefaults();
  setSigPdf_LZ_GZ();
  // n_ww->setConstant(1);
  // simpleFitForATGC("smurf/qqww.root");
  // simpleFitForATGC("smurf/ww_mcnlo.root");
  // simpleFitForATGC("smurf/ww_mcnlo_up.root");
  simpleFitForATGC("smurf/ww_mcnlo_down.root");
  // simpleFitForATGC("samples/processed_data_WW_1j_kg100_lg0_gg100_kz175_lz0_gz175_fastsim386_v1.root",false);
}
