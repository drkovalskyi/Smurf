// -*- C++ -*-
// $Id$
// Description:
//   Smurf Tree producer with the following use cases:
//   * step-by-step synchronization
//   * Smurf ntuple producer
//   * a standard EDFilter to be used in CMSSW processing
//
#ifndef SMURF_PROCESSINGANDSKIMMING_SMURFTREEPRODUCER_H
#define SMURF_PROCESSINGANDSKIMMING_SMURFTREEPRODUCER_H

#include "Smurf/ProcessingAndSkimming/interface/Core.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TDatabasePDG.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/ProcessingAndSkimming/interface/Monitor.h"

namespace reco { class Vertex; }

class TFile;
class MuonMVAEstimator;
class PFPileUpAlgo;
class EGammaMvaEleEstimator;
class EcalClusterLazyTools;
class TransientTrackBuilder;

namespace Smurf{

  namespace muon {enum SelectionType {All,Iso,Id,Acceptance,d0,dz};}
  namespace ele {enum SelectionType {All,Iso,Id,Acceptance,d0,dz,Conversion,FO};}
  namespace jet {enum SelectionType {All};}

  class SmurfTreeProducer : public edm::EDFilter {
  public:
    
    explicit SmurfTreeProducer(const edm::ParameterSet&);
    ~SmurfTreeProducer();
  
  private:
    virtual void beginJob(){}
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    const reco::Vertex* primaryVertex( edm::Event& ) const;
    void fillCommonVariables(const edm::Event& iEvent);
    bool passedMuonSelection( edm::Event&, const reco::Muon&, muon::SelectionType type ) const;
    bool passedElectronSelection( edm::Event& iEvent, const reco::GsfElectron& ele, ele::SelectionType selection ) const;

    //   bool passedElectronSelection( const reco::Muon&, ele::SelectionType type ) const;
    //   bool softmuon( const edm::Event&, const SmurfTree& ww ) const;
    //   bool extralepton( const edm::Event&, const SmurfTree& ww ) const;
    //   bool metPreselection( const SmurfTree& ww ) const;
    //   bool metFull( const SmurfTree& ww ) const;
    //   bool metBad( const SmurfTree& ww ) const;
    //   template<typename T> bool jetVeto( const T& collection, const SmurfTree& ww,
    // 				     double maxPt, double maxEta ) const;
    //   template<typename T> T mostEnergeticJet( const std::vector<T>& collection, 
    // 					   const SmurfTree& ww, double maxEta ) const;
    
    //   double dxy(const reco::Candidate* cand, const math::XYZPoint& vertex) const;
    //   double dz( const reco::Candidate* cand, const math::XYZPoint& vertex) const;
    
    //   bool muAccept( const reco::Muon& muon ) const;
    //   bool elAccept( const reco::GsfElectron& el ) const;
    // ============================================ //
    edm::InputTag m_muons;
    edm::InputTag m_electrons;
    edm::InputTag m_jets;
    edm::InputTag m_primaryVertex;
    edm::InputTag m_pfCands;
    edm::InputTag m_rhoIsoAllInputTag;
    edm::InputTag m_rhoIsoNeutral2011InputTag;
    edm::InputTag m_conversionsInputTag;
    edm::InputTag m_beamSpotInputTag;

    //  TDatabasePDG m_pdg;
    bool m_matchCharge;
    double m_l1MinPt;
    double m_l2MinPt;
    double m_minMass;
    Monitor m_monitor;
    SmurfTree* m_tree;
    TFile* m_file;
    MuonMVAEstimator* m_reader_muonHZZ2012IsoRingsMVA;
    EGammaMvaEleEstimator* m_reader_egammaPOG2012MVA;
    PFPileUpAlgo* m_pfPileUpAlgo;
    EcalClusterLazyTools* m_clusterTools;
    const TransientTrackBuilder* m_ttBuilder;
  };
}
#endif
