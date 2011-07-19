// makeOverlay() function
// overlay LR histograms
// calculate the fom based on S/B
// create the histograms for the limits setting
 
#include "TFile.h"
#include "../TVar.hh"
#include "TString.h"

void getRefYields(int &mH, Float_t & nsig_CB, Float_t & nbkg_CB, Float_t & nsig_MVA, Float_t & nbkg_MVA);

void makePlots(int mH, TString outputDir, Float_t lumi, Float_t datalumi )
{
  gROOT->ProcessLine(".L makeOverlay_HZZ.C++"); 
  
  using namespace std;
  Float_t nsig_CB(0.0), nbkg_CB(0.0), nsig_MVA(0.0), nbkg_MVA(0.0);
  getRefYields(mH, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA);

  // this assumes if we want to draw the data points on top of the MC stacks
  // be aware, this is a temporarily solution,
  makeOverlay_HZZ(mH, outputDir, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA, lumi, datalumi);
}


// Get the reference final event counts based on MC
// Taken from V3 of the AN-2011/155
// http://cms.cern.ch/iCMS/jsp/db_notes/showNoteDetails.jsp?noteID=CMS%20AN-2011/155
void getRefYields(int &mH, Float_t & nsig_CB, Float_t & nbkg_CB, Float_t & nsig_MVA, Float_t & nbkg_MVA) {
  switch (mH) {

  case (115):
    nsig_CB = 5.7;
    nbkg_CB = 73.4;
    nsig_MVA = 5.7;
    nbkg_MVA = 65.8;
    break;

  case (120):
    nsig_CB = 10.6;
    nbkg_CB = 86.2;
    nsig_MVA = 10.2;
    nbkg_MVA = 79.1;
    //nsig_MVA = 2.18;
    //nbkg_MVA = 5.55;
    break;

  case (130):
    nsig_CB = 21.4;
    nbkg_CB = 92.6;
    nsig_MVA = 21.6;
    nbkg_MVA = 92.3;
    //nsig_MVA = 6.14;
    //nbkg_MVA = 12.15;
    break;

  case (140):
    nsig_CB = 26.9;
    nbkg_CB = 73.2;
    nsig_MVA = 24.6;
    nbkg_MVA = 56.6;
    // nsig_MVA = 10.82;
    // nbkg_MVA = 17.68;
    break;

  case (150):
    nsig_CB = 26.1;
    nbkg_CB = 44.0;
    nsig_MVA = 25.0;
    nbkg_MVA = 33.7;
    break;

  case (160):
    nsig_CB = 45.7;
    nbkg_CB = 40.5;
    nsig_MVA = 51.2;
    nbkg_MVA = 44.1;
    break;

  case (170):
    nsig_CB = 34.2;
    nbkg_CB = 24.9;
    nsig_MVA = 37.5;
    nbkg_MVA = 25.4;
    break;

  case (180):
    nsig_CB = 24.7;
    nbkg_CB = 28.0;
    nsig_MVA = 27.4;
    nbkg_MVA = 29.6;
    break;

  case (190):
    nsig_CB = 19.8;
    nbkg_CB = 41.7;
    nsig_MVA = 20.3;
    nbkg_MVA = 41.7;
    break;

  case (200):
    nsig_CB = 14.4;
    nbkg_CB = 43.0;
    nsig_MVA = 14.2;
    nbkg_MVA = 38.3;
    break;

  case (250):
    nsig_CB = 6.6;
    nbkg_CB = 46.2;
    nsig_MVA = 6.5;
    nbkg_MVA = 33.5;
    break;

  case (300):
    nsig_CB = 5.5;
    nbkg_CB = 39.1;
    nsig_MVA = 5.9;
    nbkg_MVA = 33.9;
    break;
    
  default:
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 0.0;
    nbkg_MVA = 0.0;
    break;
  }
  return;
}
