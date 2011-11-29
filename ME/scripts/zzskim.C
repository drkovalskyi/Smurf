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
  smurfproducer("/smurf/ksung/ntuples_v2/", "data_2l.goodlumi.Run2011B.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "data_2l.goodlumiRun2011FinalB.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "data_2l.goodlumiFull2011.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "data_2l.goodlumiRun2011A.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "data_2l.goodlumi.Run2011B.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);

  // gg->H->ZZ
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz250.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz275.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz300.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz325.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz350.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz400.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz500.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "gfhzz600.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);

  // qq->H->ZZ
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz250.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz275.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz300.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz325.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz350.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz400.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz500.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "vbfhzz600.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);

  // diboson backgrounds
  smurfproducer("/smurf/ksung/ntuples_v2/", "zz2l-pythia.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "zz2l-madgraph.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "ww-madgraph.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "wz3l-pythia.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "wz3l-madgraph.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);

  // top
  smurfproducer("/smurf/ksung/ntuples_v2/", "ttbar2l-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "stop-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "ttop-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "wtop-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "stopb-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "ttopb-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "wtopb-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);

  // DY
  smurfproducer("/smurf/ksung/ntuples_v2/", "zee-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "zmm-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "ztt-powheg.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  smurfproducer("/smurf/ksung/ntuples_v2/", "zll-madgraph.root", "/smurf/yygao/data/HZZ2011Final/ZZ/allj/", "ZZ",3);
  
}
