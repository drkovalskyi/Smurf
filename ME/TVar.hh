#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

#define EBEAM 4000.00
#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define   eight_2Pi_to_5 7.83410393050320417e+04
#define    four_2Pi_to_2 39.478417604357432
class TVar{
  public:

    enum VerbosityLevel {
        ERROR = 0,
        INFO = 1,
        DEBUG = 2
    };

       enum MatrixElement{
          MCFM,
          MadGraph
       };


       enum Process{
          WW      = 0,
          Wp_1jet = 1,
          Wm_1jet = 2,
          Wp_gamma= 3,
          Wm_gamma= 4,
          HWW110  = 5,
          HWW115  = 6,
          HWW120  = 7,
          HWW125  = 8,
          HWW130  = 9,
          HWW135  = 10,
          HWW140  = 11,
          HWW145  = 12,
          HWW150  = 13,
          HWW155  = 14,
          HWW160  = 15,
          HWW170  = 16,
          HWW180  = 17,
          HWW190  = 18,
          HWW200  = 19,
          HWW250  = 20,
          HWW300  = 21,
          HWW350  = 22,
          HWW400  = 23,
          HWW450  = 24,
          HWW500  = 25,
          HWW550  = 26,
	  HWW600  = 27,
	  HWW700  = 28,
	  HWW800  = 29,
	  HWW900  = 30,
	  HWW1000 = 31,
	  HWW     = 32,
          HZZ200  = 33,
          HZZ250  = 34,
          HZZ300  = 35,
          HZZ350  = 36,
          HZZ400  = 37,
	  HZZ500  = 38,
	  HZZ600  = 39,
	  HZZ     = 40,
	  Z_2l    = 41,
          WpZ_lostW=42,
          WpZ_lostZ=43,
          ZZ_4l    =44,
          ttbar    =45,
          Wlnu     =46,
          WZ       =47,
          ZZ       =48,
	  WWj      =49,
	  HWWj     =50,
	  ggWW     =51,
	  Null
       };


       enum HWWPhaseSpace{
          MHMW=0,
          MH  =1,
          MWMW=2,
          MHYH=3
       };


  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){

       if(temp==TVar::WW      ) return TString("WW");
  else if(temp==TVar::Wp_gamma) return TString("Wp_gamma");
  else if(temp==TVar::Wm_gamma) return TString("Wm_gamma");
  else if(temp==TVar::Wp_1jet ) return TString("Wp_1jet");
  else if(temp==TVar::Wm_1jet ) return TString("Wm_1jet");
  else if(temp==TVar::ZZ      ) return TString("ZZ");
  else if(temp==TVar::WZ      ) return TString("WZ");
  else if(temp==TVar::HZZ200  ) return TString("HZZ200");
  else if(temp==TVar::HZZ250  ) return TString("HZZ250");
  else if(temp==TVar::HZZ300  ) return TString("HZZ300");
  else if(temp==TVar::HZZ350  ) return TString("HZZ350");
  else if(temp==TVar::HZZ400  ) return TString("HZZ400");
  else if(temp==TVar::HZZ500  ) return TString("HZZ500");
  else if(temp==TVar::HZZ600  ) return TString("HZZ600");
  else if(temp==TVar::ZZ_4l   ) return TString("ZZ_4l");
  else if(temp==TVar::Z_2l    ) return TString("Z_2l");
  else if(temp==TVar::WpZ_lostW) return TString("WpZ_lostW");
  else if(temp==TVar::WpZ_lostZ) return TString("WpZ_lostZ");
  else if(temp==TVar::Wlnu   ) return TString("Wlnu");
  else if(temp==TVar::ttbar   ) return TString("ttbar");
  else if(temp==TVar::HWW110  ) return TString("HWW110");
  else if(temp==TVar::HWW115  ) return TString("HWW115");
  else if(temp==TVar::HWW120  ) return TString("HWW120");
  else if(temp==TVar::HWW125  ) return TString("HWW125");
  else if(temp==TVar::HWW130  ) return TString("HWW130");
  else if(temp==TVar::HWW135  ) return TString("HWW135");
  else if(temp==TVar::HWW140  ) return TString("HWW140");
  else if(temp==TVar::HWW145  ) return TString("HWW145");
  else if(temp==TVar::HWW150  ) return TString("HWW150");
  else if(temp==TVar::HWW155  ) return TString("HWW155");
  else if(temp==TVar::HWW160  ) return TString("HWW160");
  else if(temp==TVar::HWW170  ) return TString("HWW170");
  else if(temp==TVar::HWW180  ) return TString("HWW180");
  else if(temp==TVar::HWW190  ) return TString("HWW190");
  else if(temp==TVar::HWW200  ) return TString("HWW200");
  else if(temp==TVar::HWW250  ) return TString("HWW250");
  else if(temp==TVar::HWW300  ) return TString("HWW300");
  else if(temp==TVar::HWW350  ) return TString("HWW350");
  else if(temp==TVar::HWW400  ) return TString("HWW400");
  else if(temp==TVar::HWW450  ) return TString("HWW450");
  else if(temp==TVar::HWW500  ) return TString("HWW500");
  else if(temp==TVar::HWW550  ) return TString("HWW550");
  else if(temp==TVar::HWW600  ) return TString("HWW600");
  else if(temp==TVar::HWW700  ) return TString("HWW700");
  else if(temp==TVar::HWW800  ) return TString("HWW800");
  else if(temp==TVar::HWW900  ) return TString("HWW900");
  else if(temp==TVar::HWW1000 ) return TString("HWW1000");
  else if(temp==TVar::WWj     ) return TString("WWj");
  else if(temp==TVar::HWWj    ) return TString("HWWj");
  else if(temp==TVar::HWW     ) return TString("HWW");

  return TString("UnKnown");
    
  };
 

 static TString SmurfProcessName(int temp){
   if(temp==TVar::WW      ) return TString("ww");
  else if(temp==TVar::Wp_gamma) return TString("wgamma");
  else if(temp==TVar::Wm_gamma) return TString("wgamma");
  else if(temp==TVar::Wp_1jet ) return TString("wjets");
  else if(temp==TVar::Wm_1jet ) return TString("wjets");
  else if(temp==TVar::ZZ      ) return TString("zz");
  else if(temp==TVar::WZ      ) return TString("wz");
  else if(temp==TVar::Z_2l    ) return TString("dyll");
  else if(temp==TVar::ttbar   ) return TString("ttbar");
  else if(temp==TVar::HWW110  ) return TString("hww110");
  else if(temp==TVar::HWW115  ) return TString("hww115");
  else if(temp==TVar::HWW120  ) return TString("hww120");
  else if(temp==TVar::HWW125  ) return TString("hww125");
  else if(temp==TVar::HWW130  ) return TString("hww130");
  else if(temp==TVar::HWW135  ) return TString("hww135");
  else if(temp==TVar::HWW140  ) return TString("hww140");
  else if(temp==TVar::HWW145  ) return TString("hww145");
  else if(temp==TVar::HWW150  ) return TString("hww150");
  else if(temp==TVar::HWW155  ) return TString("hww155");
  else if(temp==TVar::HWW160  ) return TString("hww160");
  else if(temp==TVar::HWW170  ) return TString("hww170");
  else if(temp==TVar::HWW180  ) return TString("hww180");
  else if(temp==TVar::HWW190  ) return TString("hww190");
  else if(temp==TVar::HWW200  ) return TString("hww200");
  else if(temp==TVar::HWW250  ) return TString("hww250");
  else if(temp==TVar::HWW300  ) return TString("hww300");
  else if(temp==TVar::HWW350  ) return TString("hww350");
  else if(temp==TVar::HWW400  ) return TString("hww400");
  else if(temp==TVar::HWW450  ) return TString("hww450");
  else if(temp==TVar::HWW500  ) return TString("hww500");
  else if(temp==TVar::HWW550  ) return TString("hww550");
  else if(temp==TVar::HWW600  ) return TString("hww600");
  else if(temp==TVar::HWW700  ) return TString("hww700");
  else if(temp==TVar::HWW800  ) return TString("hww800");
  else if(temp==TVar::HWW900  ) return TString("hww900");
  else if(temp==TVar::HWW1000 ) return TString("hww1000");
  else if(temp==TVar::HZZ200  ) return TString("hzz200");
  else if(temp==TVar::HZZ250  ) return TString("hzz250");
  else if(temp==TVar::HZZ300  ) return TString("hzz300");
  else if(temp==TVar::HZZ350  ) return TString("hzz350");
  else if(temp==TVar::HZZ400  ) return TString("hzz400");
  else if(temp==TVar::HZZ500  ) return TString("hzz500");
  else if(temp==TVar::HZZ600  ) return TString("hzz600");
  else
    return TString("UnKnown");
 };
 

  ClassDef(TVar,0)
};

struct branch_particle {
  int   PdgCode   ;
  int   Charge    ;
  double Px       ;
  double Py       ;
  double Pz       ;
  double E        ;
  double Eta      ;
  double Phi      ;

};
static const TString branch_format_particle =
 "PdgCode/I:"
 "Charge/I:"
 "Px/D:"
 "Py/D:"
 "Pz/D:"
 "E/D:"
 "Eta/D:"
 "Phi/D";

// in development
struct hzz4l_event_type{
  int PdgCode[4];
  TLorentzVector p[4];
  double Xsec   [10];
  double XsecErr[10];  
};
// in development
  
struct cdf_event_type{
  int PdgCode[2];
  TLorentzVector p[2];
  double MetX;
  double MetY;

  double Xsec   [60];
  double XsecErr[60];  
};
struct mcfm_event_type{
  int PdgCode[6];
  TLorentzVector p[6];
  double pswt;
};
struct event_type{
  TLorentzVector p1,p2,ep,em,nu,nb;
  double PSWeight;
};

struct array_event_type{
  int PdgCode[6];
  double p4[4][6];
  double PSWeight;
};

struct event_XSec{
  double dXSec, dXSecErr;
};

struct phasespace_Mw1Mw2{
	           double nuZ,nbZ,Mw1,Mw2,PSWeight;
};


struct phasespace_4D{
	   double x1,x2,costh_nu,phi_nu,PSWeight;
};


struct rand_type{
	   double r[10];
};

struct anomcoup{
	   double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z_,tevscale;
};

struct EffHist{
  TH2F* els_eff_mc;
  TH2F* mus_eff_mc;
};

struct FRHist{
  TH2F* els_fr;
  TH2F* mus_fr;
  TH2F* els_part_fo;
  TH2F* mus_part_fo;
  TH2F* els_ptres;
  TH2F* mus_ptres;
};

struct BoostHist{
  TH1F* kx;
  TH1F* ky;
};

#endif
