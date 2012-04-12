#ifndef SMURF_PROCESSINGANDSKIMMING_UTILITIES_H
#define SMURF_PROCESSINGANDSKIMMING_UTILITIES_H

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <string>

namespace smurfutilities {

// compare particle flow isolation definitions
void ValidatePFIsolation(const edm::Event& iEvent, const reco::GsfElectron &ele,  
        const reco::PFCandidateCollection &pfCands,
        const edm::Handle<reco::VertexCollection> &vertexHandle,
        const float &pfiso_ch, const float &pfiso_em, const float &pfiso_nh);

// dump save tags for trigger name
void DumpSaveTags(const std::string triggerName,
        const HLTConfigProvider &hltConfig);

}

#endif

