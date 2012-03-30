#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

#define EBEAM 3500.0
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
          HWW115  = 5,
          HWW120  = 6,
          HWW130  = 7,
          HWW140  = 8,
          HWW150  = 9,
          HWW160  =10,
          HWW170  =11,
          HWW180  =12,
          HWW190  =13,
          HWW200  =14,
          HWW210  =15,
          HWW220  =16,
          HWW230  =17,
          HWW250  =18,
          HWW300  =19,
          HWW350  =20,
          HWW400  =21,
          HWW450  =22,
          HWW500  =23,
          HWW550  =24,
	  HWW600  =25,
          HZZ200  =26,
          HZZ250  =27,
          HZZ300  =28,
          HZZ350  =29,
          HZZ400  =30,
	  HZZ500  =31,
	  HZZ600  =32,
	  HWW     =33,
          HZZ     =34,
	  Z_2l    =35,
          WpZ_lostW=36,
          WpZ_lostZ=37,
          ZZ_4l    =38,
          ttbar    =39,
          Wlnu     =40,
          WZ       =41,
          ZZ       =42,
	  WWj      =43,
	  HWWj     =44,
	  ggWW     =45,
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
  else if(temp==TVar::HWW     ) return TString("HWW");
  else if(temp==TVar::HZZ     ) return TString("HZZ");
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
  else if(temp==TVar::HWW115  ) return TString("HWW115");
  else if(temp==TVar::HWW120  ) return TString("HWW120");
  else if(temp==TVar::HWW130  ) return TString("HWW130");
  else if(temp==TVar::HWW140  ) return TString("HWW140");
  else if(temp==TVar::HWW150  ) return TString("HWW150");
  else if(temp==TVar::HWW160  ) return TString("HWW160");
  else if(temp==TVar::HWW170  ) return TString("HWW170");
  else if(temp==TVar::HWW180  ) return TString("HWW180");
  else if(temp==TVar::HWW190  ) return TString("HWW190");
  else if(temp==TVar::HWW200  ) return TString("HWW200");
  else if(temp==TVar::HWW210  ) return TString("HWW210");
  else if(temp==TVar::HWW220  ) return TString("HWW220");
  else if(temp==TVar::HWW230  ) return TString("HWW230");
  else if(temp==TVar::HWW250  ) return TString("HWW250");
  else if(temp==TVar::HWW300  ) return TString("HWW300");
  else if(temp==TVar::HWW350  ) return TString("HWW350");
  else if(temp==TVar::HWW400  ) return TString("HWW400");
  else if(temp==TVar::HWW450  ) return TString("HWW450");
  else if(temp==TVar::HWW500  ) return TString("HWW500");
  else if(temp==TVar::HWW550  ) return TString("HWW550");
  else if(temp==TVar::HWW600  ) return TString("HWW600");
  else if(temp==TVar::WWj      ) return TString("WWj");
  else if(temp==TVar::HWWj     ) return TString("HWWj");

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
  else if(temp==TVar::HWW115  ) return TString("hww115");
  else if(temp==TVar::HWW120  ) return TString("hww120");
  else if(temp==TVar::HWW130  ) return TString("hww130");
  else if(temp==TVar::HWW140  ) return TString("hww140");
  else if(temp==TVar::HWW150  ) return TString("hww150");
  else if(temp==TVar::HWW160  ) return TString("hww160");
  else if(temp==TVar::HWW170  ) return TString("hww170");
  else if(temp==TVar::HWW180  ) return TString("hww180");
  else if(temp==TVar::HWW190  ) return TString("hww190");
  else if(temp==TVar::HWW200  ) return TString("hww200");
  else if(temp==TVar::HWW210  ) return TString("hww210");
  else if(temp==TVar::HWW220  ) return TString("hww220");
  else if(temp==TVar::HWW230  ) return TString("hww230");
  else if(temp==TVar::HWW250  ) return TString("hww250");
  else if(temp==TVar::HWW300  ) return TString("hww300");
  else if(temp==TVar::HWW350  ) return TString("hww350");
  else if(temp==TVar::HWW400  ) return TString("hww400");
  else if(temp==TVar::HWW450  ) return TString("hww450");
  else if(temp==TVar::HWW500  ) return TString("hww500");
  else if(temp==TVar::HWW550  ) return TString("hww550");
  else if(temp==TVar::HWW600  ) return TString("hww600");
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
