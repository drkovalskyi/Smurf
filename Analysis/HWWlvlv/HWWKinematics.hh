#ifndef HWWKINEMATICS_H
#define HWWKINEMATICS_H

#include <TVector3.h>
#include <TLorentzVector.h>
#include <vector>

class HWWKinematics {
public:
  //! constructor at least wants leptons and MET
  HWWKinematics(TLorentzVector l1, TLorentzVector l2, TVector3 met);
  virtual ~HWWKinematics() {}
  //! set the jets
  void setJets(std::vector<TLorentzVector> jets) { JETS = jets; }
  
  //! kinematic fuctions
  double CalcMRNEW();
  double CalcDeltaPhiRFRAME();
  double CalcDoubleDphiRFRAME();
  double CalcMR();

private:
  TLorentzVector L1,L2;
  TVector3 MET;
  std::vector<TLorentzVector> JETS;
};

#endif // HWWKINEMATICS_H
