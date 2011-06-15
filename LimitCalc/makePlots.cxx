void makePlots(const char* pattern="output/limit-*", 
	       const char* title="H #rightarrow WW #rightarrow 2l2#nu + 0/1 jets")
{
  gSystem->Load("lands.so");
  gSystem->CompileMacro("PlotExpectedLimits.C","k");
  PlotExpectedLimits(pattern,title);
}
