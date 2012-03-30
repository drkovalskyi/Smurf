//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TEvtProb.hh"
#include "TVar.hh"
#include "TMatrixElement.hh"


ClassImp(TEvtProb)

    using namespace std;

    //-----------------------------------------------------------------------------
    // Constructors and Destructor
    //-----------------------------------------------------------------------------
    TEvtProb::TEvtProb() {
        _matrixElement = TVar::MCFM;
        TMatrixElement::inst()->SetMatrixElement(_matrixElement);
        _hwwPhaseSpace=TVar::MH;
        mcfm_init_();
    }

TEvtProb::~TEvtProb() {}

//-----------------------------------------------
//
// Importance Sampleing Integration
//
//-----------------------------------------------

void   TEvtProb::NeutrinoIntegrate(TVar::Process proc,
        const cdf_event_type &cdf_event,
        double *Xsec,
        double *XsecErr, TVar::VerbosityLevel verbosity){

    //Initialize Process
    SetProcess(proc);

    // This is bad.  It seems to be setting some kind of 
    // global BW function somewhere.
    // This function, My_choose, is defined in TUtil
    My_choose(proc);

    //delta(L1) delta(L2)
    int NDim = bveg1_mcfm_.ndim-6;
    //qx,qy
    if (_smearLevel >= 1) NDim+=2;
    //dE1,dE2
    if (_smearLevel >= 2) NDim+=2;


    if (verbosity >= TVar::DEBUG) {
        cout <<" [NeutrinoIntegrate]: Evaluate " << TVar::ProcessName(proc)
            <<" Ncalls " << _ncalls
            <<" npart._npart=" << npart_.npart
            <<" bveg1_mcfm_.ndim= " << bveg1_mcfm_.ndim
            <<" NDim= "<<NDim<<endl;
    }

    //Phase space generation
    TRandom3 myRandom;
    myRandom.SetSeed(_seed);

    double r  [22] = { 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
        0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
        0.5,0.5};


    int count_PS=0;
    double sumW=0,sumW2=0;

    double probAcceptanceEfficiency;

    if (proc == TVar::Wp_1jet || proc == TVar::Wm_1jet  )
        probAcceptanceEfficiency = getFakeRateProb(cdf_event, _effhist, _FRhist, proc, verbosity);
    else 
        probAcceptanceEfficiency = getProbAcceptanceEfficiency(cdf_event, _effhist);

    if(probAcceptanceEfficiency<=0) {
      cout <<"ERROR: probAcceptanceEfficiency =  " << probAcceptanceEfficiency <<endl;
        return;
    }

    for(int idx=0;idx<_ncalls;idx++){
        count_PS++;
        myRandom.RndmArray(NDim,r); // generate NDim random numbers and set the first NDim entries of r arrary
        double dXsec=0;
        dXsec=Integrand_NeutrinoIntegration(proc, cdf_event, r, NDim, 0, _boosthist, verbosity)*probAcceptanceEfficiency;
        
	if (dXsec<=0) continue;
        sumW  += dXsec;
        sumW2 += dXsec*dXsec;

    }

    *Xsec = sumW/count_PS;

    *XsecErr = fabs(sumW2/count_PS-sumW/count_PS*sumW/count_PS);
    if(*XsecErr>0.0) *XsecErr = sqrt(*XsecErr/count_PS);
    else             *XsecErr = -1;  


}

//=======================================
// Integrand
//=======================================
double TEvtProb::Integrand_NeutrinoIntegration(TVar::Process proc, const cdf_event_type &cdf_event, double * r, unsigned int NDim, void * param, BoostHist boosthist, TVar::VerbosityLevel verbosity){

    //constants
    double sqrts = 2.*EBEAM;
    double W=sqrts*sqrts;

    int NSol = 4;

    //Weight calculation
    double flux=1.;
    double dXsec=0.;
    double sumW =0.;

    // array of mcfm event type
    mcfm_event_type arr_mcfm_event[4];


    if (proc == TVar::WW)
    {
        NSol=4;         
        genMw1Mw2(r, _smearLevel, cdf_event , arr_mcfm_event, boosthist);
    }

    if (proc ==TVar::HWW)
    {
        if(_hwwPhaseSpace == TVar::MWMW)
        { 
            NSol=4;
            genMw1Mw2(r, _smearLevel, cdf_event, arr_mcfm_event, boosthist);
        }
        else if (_hwwPhaseSpace == TVar::MHMW)
        { 
            NSol=4;
            genMHiggsMw1(r, _smearLevel, cdf_event, arr_mcfm_event, boosthist);
        }
        else if (_hwwPhaseSpace == TVar::MHYH)
        { 
            NSol=2;
            genMHiggsYHiggs(r, _smearLevel, cdf_event, arr_mcfm_event, boosthist);
        }
        else if (_hwwPhaseSpace ==TVar::MH)
        { 
            NSol=2;
            genMHiggs(r,_smearLevel, cdf_event, arr_mcfm_event, boosthist);
        }
    }

    else if (proc == TVar::Wp_gamma || proc == TVar::Wm_gamma )
    { 
        NSol=2;  
        genMw_Wgamma(r,_smearLevel, cdf_event, arr_mcfm_event);
    }
    else if (proc == TVar::Wp_1jet || proc == TVar::Wm_1jet  )
    { 
        NSol=2; 
        genMw_W1jet(r,_smearLevel, cdf_event, arr_mcfm_event, boosthist);
    }
    else if (proc == TVar::ZZ)
    { 
        NSol=2;
        genMzNu3     (r,_smearLevel, cdf_event, arr_mcfm_event, boosthist);
    }
    else if (proc == TVar::WZ)
    { 
        NSol=2;
        genMwNu3WZ     (r,_smearLevel, cdf_event, arr_mcfm_event, boosthist);
    }
    else if (proc == TVar::HZZ)
    { 
        NSol=2;
        genMzNu3HZZ  (r,_smearLevel, cdf_event, arr_mcfm_event, boosthist);
    }
    else if (proc == TVar::Z_2l)
    { 
        NSol=1;    
        genDY(r,_smearLevel, cdf_event, arr_mcfm_event);
    } 

    //loop over solutions
    for(int jps=0;jps<NSol;jps++){

        mcfm_event_type& mcfm_event = arr_mcfm_event[jps];

        //      cout << "mcfm_event.pswt = " << mcfm_event.pswt << endl;
        if(mcfm_event.pswt<=0) continue;

        //Matrix Element evaluation in qX=qY=0 frame
        //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
        // 
        double qX=mcfm_event.p[0].Px()+mcfm_event.p[1].Px();
        double qY=mcfm_event.p[0].Py()+mcfm_event.p[1].Py();
        // cout<< "TEvtProb:: qX = " << qX<< "qY = "<<qY<<"\n";

        if((qX*qX+qY*qY)>0){
            double qE = mcfm_event.p[0].Energy()+mcfm_event.p[1].Energy();
            TVector3 boostV(qX/qE,qY/qE,0);
            for(int ipt=0;ipt<6;ipt++) mcfm_event.p[ipt].Boost(-boostV);
        }

        //event selections in Lab Frame
        if(My_eventcuts(proc, mcfm_event, cdf_event)){ mcfm_event.pswt=0; continue;}
        double flavor_msq[nmsq][nmsq];
        double msqjk = SumMatrixElementPDF(proc, &mcfm_event, flavor_msq, &flux);

        if(msqjk<=0){ mcfm_event.pswt=0; continue;}

	//	cout<<"msqjk = "<<msqjk<<"\n";
        //npart_ means number of final state particles. 
        //Don't forget the first two are 2 initial state particles
        //for(int ipt=2;ipt<npart_.npart+2;ipt++) msqjk = msqjk/mcfm_event.p[ipt].Energy();

        flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
	
	//	cout<<"pswt = "<<mcfm_event.pswt<<"\n";
        dXsec=msqjk*flux*mcfm_event.pswt;

        if (dXsec>0.0){
            sumW  += dXsec;
            mcfm_event.pswt=dXsec;
        }  
        else if (verbosity >= TVar::DEBUG)
	  {
	    cout <<" NeutrinoIntegrate Warning: dXsec " << dXsec
		 << " dXsec==dXsec " << (dXsec==dXsec)
		 << " dXsec>0.0 " << (dXsec>0.0)<<" "
		 <<" Msq="<<msqjk 
		 <<" flux="<<flux 
		 <<" wgt="<<mcfm_event.pswt
		 <<endl;
	  }
    }//loop solutions
    
    return sumW;

}

void TEvtProb::SetFRHist(TString elFRFile, TString elFRHist, TString muFRFile, TString muFRHist, TVar::VerbosityLevel verbosity) {
  TFile *elFRfile = TFile::Open(elFRFile, "READ");
  assert(elFRfile);
  gROOT->cd();
  if (verbosity >= TVar::DEBUG) 
    std::cout << "TEvtProb::SetFRHist set electron Fakerate from " << elFRFile << std::endl;
  _FRhist.els_fr = (TH2F*) elFRfile->Get(elFRHist)->Clone();
  elFRfile->Close();
  
  TFile *muFRfile = TFile::Open(muFRFile, "READ");
  assert(muFRfile);
  gROOT->cd();
  if (verbosity >= TVar::DEBUG) 
    std::cout << "TEvtProb::SetFRHist set muon Fakerate from " << muFRFile << std::endl;
  _FRhist.mus_fr = (TH2F*) muFRfile->Get(muFRHist)->Clone();
  muFRfile->Close();
}

void TEvtProb::SetMCHist(int proc, TString MCFileName, bool setFR, int njets, TVar::VerbosityLevel verbosity) {
  
  if (verbosity >= TVar::DEBUG) std::cout << "TEvtProb::SetMCHist for process " << TVar::ProcessName(proc) << " from " << MCFileName << std::endl;
    TFile *fUtil = TFile::Open(MCFileName, "READ");
    assert(fUtil);
    gROOT->cd();
    // lepton efficiencies are taken from the WW MC
    _effhist.els_eff_mc = (TH2F*) fUtil->Get("ww_heleEff")->Clone();
    _effhist.mus_eff_mc = (TH2F*) fUtil->Get("ww_hmuEff")->Clone();


    if (TVar::SmurfProcessName(proc) == "wjets") {
        _boosthist.kx = (TH1F*) fUtil->Get(TString("ww_kx"))->Clone();
        _boosthist.ky = (TH1F*) fUtil->Get(TString("ww_ky"))->Clone();
    }
    else 
      if (TVar::SmurfProcessName(proc) == "dyll") {}
    else
    {
      if ( njets == 0) {
        _boosthist.kx = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_kx"))->Clone();
        _boosthist.ky = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_ky"))->Clone();
      }
      if ( njets == 1) {
	_boosthist.kx = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_kx_1jet"))->Clone();
        _boosthist.ky = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_ky_1jet"))->Clone();
      }
      if ( njets >= 2) {
	_boosthist.kx = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_kx_2jet"))->Clone();
        _boosthist.ky = (TH1F*) fUtil->Get(TString(TVar::SmurfProcessName(proc)+"_ky_2jet"))->Clone();
      }
      
      
    }

    // fake-rate histograms
    if (setFR) {
      _FRhist.els_fr = (TH2F*) fUtil->Get("wjets_heleFR")->Clone();
      _FRhist.mus_fr = (TH2F*) fUtil->Get("wjets_hmuFR")->Clone();
    }
    
    if (  fUtil->Get("wjets_heleGenFR") != 0x0 ) 
      _FRhist.els_part_fo = (TH2F*) fUtil->Get("wjets_heleGenFR")->Clone();
      
    if ( fUtil->Get("wjets_hmuGenFR") != 0x0 )
      _FRhist.mus_part_fo = (TH2F*) fUtil->Get("wjets_hmuGenFR")->Clone();

    if ( fUtil->Get("wjets_heleFOResponse") != 0x0 ) 
      _FRhist.els_ptres = (TH2F*) fUtil->Get("wjets_heleFOResponse")->Clone();

    if ( fUtil->Get("wjets_hmuFOResponse") != 0x0 ) 
      _FRhist.mus_ptres = (TH2F*) fUtil->Get("wjets_hmuFOResponse")->Clone();

    fUtil->Close();
}

