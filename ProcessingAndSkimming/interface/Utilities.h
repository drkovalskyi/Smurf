#ifndef SMURF_PROCESSINGANDSKIMMING_UTILITIES_H
#define SMURF_PROCESSINGANDSKIMMING_UTILITIES_H

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <string>

namespace smurfutilities {

typedef math::XYZTLorentzVectorD LorentzVector;

// dump save tags for trigger name
void DumpSaveTags(const std::string triggerName,
        const HLTConfigProvider &hltConfig);

// does trigger object match offline object
// prescale is returned - 0 if no match.
unsigned int MatchTriggerObject(const edm::Event &iEvent, const edm::EventSetup &iSetup,
        const std::string triggerName, const std::string filterName,
        const std::string processName, const HLTConfigProvider &hltConfig, const edm::TriggerResults* triggerResults,
        const trigger::TriggerEvent *triggerEvent,
        const trigger::TriggerObjectCollection &allObjects,
        const LorentzVector &offlineObject);

// does trigger object match offline object
bool MatchTriggerObject(const trigger::TriggerObjectCollection &triggerObjects, const LorentzVector &offlineObject);

// match to gen particle
// of specified id and status
float MatchGenParticle(const reco::GenParticleCollection &genParticleCollection, const LorentzVector &object, 
        const int pdgId, const int status);

// calculate mt
float Mt(const float &pt1, const float &pt2, const float &dphi);

}

#endif

