{
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat("nemruosk");
   gStyle->SetPalette(1);
   if (gSystem->Getenv("CMSSW_VERSION")){
    cout << "loading..." <<endl;
    gSystem->Load("libCintex");
    Cintex::Enable();
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
    gSystem->Load("libDataFormatsFWLite.so");
    gSystem->SetIncludePath( "-I$ROOFITSYS/include" );
    gSystem->Load("libRooFit.so");
    using namespace RooFit;
    }
}
