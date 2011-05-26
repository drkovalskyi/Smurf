// makeOverlay() function
// overlay LR histograms
// calculate the fom based on S/B
// create the histograms for the limits setting
 
#include "TFile.h"
#include "../TVar.hh"
#include "TString.h"

void getRefYields(int &mH, Float_t & nsig_CB, Float_t & nbkg_CB, Float_t & nsig_MVA, Float_t & nbkg_MVA);

void makePlots(int mH, TString outputDir)
{
  gROOT->ProcessLine(".L makeOverlay.C+"); 
  
  using namespace std;
  Float_t nsig_CB(0.0), nbkg_CB(0.0), nsig_MVA(0.0), nbkg_MVA(0.0);
  getRefYields(mH, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA);
  makeOverlay(mH, outputDir, nsig_CB, nbkg_CB, nsig_MVA, nbkg_MVA);

}


// Get the reference final event counts based on MC
// Taken from V3 of the AN-2011/155
// http://cms.cern.ch/iCMS/jsp/db_notes/showNoteDetails.jsp?noteID=CMS%20AN-2011/155
void getRefYields(int &mH, Float_t & nsig_CB, Float_t & nbkg_CB, Float_t & nsig_MVA, Float_t & nbkg_MVA) {
  switch (mH) {
  case (120):
    nsig_CB = 7.6;
    nbkg_CB = 67.9;
    nsig_MVA = 2.18;
    nbkg_MVA = 5.55;
    break;

  case (130):
    nsig_CB = 15.0;
    nbkg_CB = 68.9;
    nsig_MVA = 6.14;
    nbkg_MVA = 12.15;
    break;

  case (140):
    nsig_CB = 21.2;
    nbkg_CB = 62.1;
    nsig_MVA = 10.82;
    nbkg_MVA = 17.68;
    break;

  case (150):
    nsig_CB = 21.5;
    nbkg_CB = 35.2;
    nsig_MVA = 10.38;
    nbkg_MVA = 8.21;
    break;

  case (160):
    nsig_CB = 30.9;
    nbkg_CB = 21.9;
    nsig_MVA = 11.05;
    nbkg_MVA = 2.31;
    break;
  case (170):
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 17.72;
    nbkg_MVA = 6.23;
    break;
  case (180):
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 10.64;
    nbkg_MVA = 6.26;
    break;

  case (190):
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 8.11;
    nbkg_MVA = 8.08;
    break;

  case (200):
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 8.05;
    nbkg_MVA = 11.23;
    break;

  case (250):
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 2.35;
    nbkg_MVA = 7.34;
    break;

  case (300):
    nsig_CB = 0.0;
    nbkg_CB = 0.0;
    nsig_MVA = 2.20;
    nbkg_MVA = 8.05;
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
