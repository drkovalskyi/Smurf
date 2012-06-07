{
  gSystem->CompileMacro("tauify.cc","k");
  gSystem->CompileMacro("dytt.C","k");
  dytt(10,1.6,
       "/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dyee.root",dy::MCee,
       "/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dytt.root");
}
