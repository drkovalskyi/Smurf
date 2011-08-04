// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)


#ifndef WW_PHASE
#define WW_PHASE
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TFile.h>
#include <TF1.h>
#include "TRandom3.h"
#include <TVar.hh>
#include "TMCFM.hh"
#include "TUtil.hh"

#include "Generator.hh"
#include "KinematicSolver.hh"


using namespace std;


//-------------------------
void genDY    (double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList);

//-------------------------
void genMw1Mw2(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, double hmass);
void genMHiggsMw1(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMHiggs   (double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMHiggsYHiggs(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);

void genMw1Mw2(double* r,int SmearLevel,     event_type  cdf_event, event_type* sol, double* jacobian);
void genMw1Mw2(double* r,int SmearLevel, cdf_event_type* cdf_event, array_event_type* sol);

//--------------------------
void genMw1Mw2(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMw1Mz2(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList);

//---------------------------
void genMw_Wgamma(double* r,int SmearLevel, cdf_event_type cdf_event,mcfm_event_type* PSList);
void genMw_W1jet (double* r,int SmearLevel, cdf_event_type cdf_event,mcfm_event_type* PSList,  BoostHist boosthist);

//---------------------------
void genMzNu3(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMzNu3HZZ(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);

void genMwNu3WZ(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);



#endif
