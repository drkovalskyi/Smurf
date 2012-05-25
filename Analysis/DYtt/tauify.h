
///////////////////////////////////////////////////////////////////////////////
//
// --- Usage --- 
//
// *Electrons and Muons are treated as massless
//
// 1) Initialize the class once at the beginning of your looper:
// 
//    TauData* taudata = new TauData();
//
// 2) Call the "Tauify" method with the 4-momenta of the lepton
//    you want to use to model as a tau, and the polarization sign s:  
//
//    taudata->Tauify( p4_taulplus_, s );
//
// 3) Now you can access the 
//
//    myLep = taudata->GetLeptonFromTau();
//    myMET = taudata->GetMetFromTau();
//
// 4) Repeat 2 & 3 as necessary
//
///////////////////////////////////////////////////////////////////////////////
// Authors: D.Barge(UCSB), C.Compagnari(UCSB)

//////////////////
// C++ Includes //
//////////////////

#include "math.h"


///////////////////
// ROOT Includes //
///////////////////

#include "Math/Boost.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TRandom3.h"


//////////////
// Typedefs //
//////////////

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4d;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef ROOT::Math::Boost Boost; 


//////////////////////
// Class Definition //
//////////////////////

class TauData {

  public:
    
    //////////////////////////////
    // Constructor & Destructor //
    //////////////////////////////
  
    TauData(int seed = 4357){
    
      /////////////////////////////
      // Initialize Data Members //
      /////////////////////////////

      // Lepton Masses with Ridiculous Precision
      m_tau_ = 1.77682;
      m_mu_  = 0.10565;
      m_el_  = 0.00051;


      // Polarization of Tau
      polarization_sign_ = 0;


      // Randoms
      seed_              = seed;                // ROOT Default
      rand3_             = new TRandom3(seed_); // 
      randA_             = new TRandom3(seed_);      // 


      // Some Angles
      cosThetaCM_  = -999.0;
      cosThetaLAB_ = -999.0;


      // Parameters for ( d^2 / dx dcos(theta) ) Monte Carlo from which x, cos(theta) are drawn for Tau Decay
      xmin_  = 0;    // Minimum Value of x = E/Emax
      xmax_  = 1;    // Maximum Value of x = E/Emax
      ymin_  = -1;   // Minimum Value of y = cos(theta)
      ymax_  = 1;    // Maximum Value of y = cos(theta)
      fmax_  = 2.01; // Maximum Value of ( d^2 sigma / dx dy )


      // TH2 
      xbins_ = 100;  // Number of x bins for histogram used to record x used in tau decay
      ybins_ = 100;  // Number of y bins for histogram used to record y used in tau decay
      h2_xy_ = new TH2F("h2_xy", "h2_xy", xbins_, xmin_, xmax_, ybins_, ymin_, ymax_ );


      // TF2
      npx_         = 30;   // Number of x points used for TF2
      npy_         = 30;   // Number of y points used for TF2

      expr_left_   = Form( "pow(x,2) * ( 3 - 2*x + %d*( 2*x - 1 )*y )", -1 );
      expr_zero_   = Form( "pow(x,2) * ( 3 - 2*x + %d*( 2*x - 1 )*y )",  0 );
      expr_right_  = Form( "pow(x,2) * ( 3 - 2*x + %d*( 2*x - 1 )*y )",  1 );

      f2_xy_left_  = new TF2( "f2_xy", Form("%s", expr_left_  ), xmin_, xmax_, ymin_, ymax_ );
      f2_xy_zero_  = new TF2( "f2_xy", Form("%s", expr_zero_  ), xmin_, xmax_, ymin_, ymax_ );
      f2_xy_right_ = new TF2( "f2_xy", Form("%s", expr_right_ ), xmin_, xmax_, ymin_, ymax_ );

      //func_->SetTitleX("test");
      //func_->GetXaxis()->SetTitle("test");
      //func->SetNpx(npx);
      //func->SetNpy(npy);

      return;

    }

    ~TauData(void){ 
      delete rand3_;
      delete randA_;
      delete f2_xy_left_;
      delete f2_xy_zero_;
      delete f2_xy_right_;
      delete h2_xy_;
      return; 
    }

    /////////////////////////////
    // Public Member Functions //
    /////////////////////////////

    bool   Tauify ( LorentzVector, int );                                   // Decay the Given Lorentz Vector as Tau with the Polarization Given

    inline float CosThetaCM  (void) { return cosThetaCM_;  }                // Angle Between the CM Lepton and the Original Tau direction in the LAB
    inline float CosThetaLAB (void) { return cosThetaLAB_; }                // Angle Between the Lepton and Tau in the LAB

    inline LorentzVector GetLeptonFromTau (void) { return p4_lep_;      }   // 4 momenta of Lepton From Tau
    inline LorentzVector GetMetFromTau    (void) { return p4_met_;      }   // 4 momenta of MET From Tau
    inline TH2F*         GetTH2dxdy       (void) { return h2_xy_;       }   // Histogram of (x, y) values used to decay Tau
    inline TF2*          GetTF2dxdy_left  (void) { return f2_xy_left_;  }   // Analytic function of (x, y) - left  Polarization
    inline TF2*          GetTF2dxdy_zero  (void) { return f2_xy_zero_;  }   // Analytic function of (x, y) - no    Polarization
    inline TF2*          GetTF2dxdy_right (void) { return f2_xy_right_; }   // Analytic function of (x, y) - right Polarization


  private:

    //////////////////
    // Data Members //
    //////////////////

    int     xbins_;
    int     ybins_;
    int     npx_;
    int     npy_;
    int     polarization_sign_;

    float   m_tau_;
    float   m_mu_;
    float   m_el_;
    float   fmax_;
    float   xmin_;
    float   xmax_;
    float   ymin_;
    float   ymax_;

    double  cosThetaCM_;
    double  cosThetaLAB_;

    char*   expr_left_;
    char*   expr_zero_;
    char*   expr_right_;

    Boost         boost_CM_;
    Boost         boost_LAB_;

    LorentzVector p4_lep_;
    LorentzVector p4_met_;
    LorentzVector p4_lep_CM_;

    TRandom3*     rand3_;
    TRandom3*     randA_;

    TH2F* h2_xy_;

    TF2*  f2_xy_left_;
    TF2*  f2_xy_zero_;
    TF2*  f2_xy_right_;

    UInt_t        seed_;


    //////////////////////
    // Member Functions //
    //////////////////////

    void   DecayAsTau            ( LorentzVector );
    Vector GetRandomNormalVector ( Vector );


};




