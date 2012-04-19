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

#include <string>

namespace smurfutilities {

typedef math::XYZTLorentzVectorD LorentzVector;

// compare particle flow isolation definitions
void ValidatePFIsolation(const edm::Event& iEvent, const reco::GsfElectron &ele,  
        const reco::PFCandidateCollection &pfCands,
        const edm::Handle<reco::VertexCollection> &vertexHandle,
        const float &pfiso_ch, const float &pfiso_em, const float &pfiso_nh);

// dump save tags for trigger name
void DumpSaveTags(const std::string triggerName,
        const HLTConfigProvider &hltConfig);

// find trigger objects corresponding
// to trigger name and filter name
trigger::TriggerObjectCollection GetTriggerObjects(const std::string triggerName, const std::string filterName,
        const std::string procesName, const HLTConfigProvider &hltConfig, const edm::TriggerResults* triggerResults,
        const trigger::TriggerEvent *triggerEvent,
        const trigger::TriggerObjectCollection &allObjects);

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

}

#endif

