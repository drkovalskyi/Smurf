#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"
#include "TStyle.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TF1.h"

#include "smurfproducer.C"

void zzskim()
{
  gROOT->ProcessLine(".L smurfproducer.C+");
  

  // data
  smurfproducer("/smurf/ksung/ntuples/", "data_2l.goodlumiRun2011A.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "data_2l.goodlumiRun2011A.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "data_2l.goodlumi.Run2011B.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "data_2l.goodlumi.Run2011B.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);


  // gg->H->ZZ
  smurfproducer("/smurf/ksung/ntuples/", "gfhzz250.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "gfhzz300.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "gfhzz350.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "gfhzz400.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "gfhzz500.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "gfhzz600.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);

  // qq->H->ZZ
  smurfproducer("/smurf/ksung/ntuples/", "vbfhzz250.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "vbfhzz300.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "vbfhzz350.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "vbfhzz400.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "vbfhzz500.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "vbfhzz600.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);

  // diboson backgrounds
  smurfproducer("/smurf/ksung/ntuples/", "zz2l-pythia.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "zz2l-madgraph.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "ww2l-pythia.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "ww-madgraph.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "wz3l-pythia.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "wz3l-madgraph.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);

  // top
  smurfproducer("/smurf/ksung/ntuples/", "ttbar2l-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "stop-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "ttop-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "wtop-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "stopb-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "ttopb-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "wtopb-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);

  // DY
  smurfproducer("/smurf/ksung/ntuples/", "zee-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "zmm-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "ztt-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples/", "zll-madgraph.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZ",3);

 
  // BTAG samples
  smurfproducer("/smurf/ksung/ntuples/", "ttbar2l-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "wtop-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "wtopb-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "ttop-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "ttopb-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "stop-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
  smurfproducer("/smurf/ksung/ntuples/", "stopb-powheg.root", "/smurf/yygao/data/HZZ2011A/ZZ/allj/", "ZZBTAG",3);
 
}
