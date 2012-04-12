#include "Smurf/ProcessingAndSkimming/interface/Utilities.h"
#include "Smurf/ProcessingAndSkimming/interface/Selections.h"

#include "TString.h"
#include "TRegexp.h"

#include <vector>


void smurfutilities::ValidatePFIsolation(const edm::Event& iEvent, const reco::GsfElectron &ele,
        const reco::PFCandidateCollection &pfCands,
        const edm::Handle<reco::VertexCollection> &vertexHandle,
        const float &pfiso_ch, const float &pfiso_em, const float &pfiso_nh)
{
    float dle_ch = 0.0;
    float dle_em = 0.0;
    float dle_nh = 0.0;
    smurfselections::PFIsolation2012(ele, pfCands, vertexHandle, 0, 0.3, dle_ch, dle_em, dle_nh);
    if (fabs(dle_ch - pfiso_ch) > 0.01 || fabs(dle_em - pfiso_em) > 0.01 || fabs(dle_nh - pfiso_nh) > 0.01) {
        std::cout << "ERROR" << std::endl;
    }
    std::cout << iEvent.id().run() << ", " << iEvent.luminosityBlock() << ", " << iEvent.id().event() << std::endl;
    std::cout << "Electron pT, eta, phi : " << ele.pt() << ", " << ele.eta() << ", " << ele.phi() << std::endl;
    std::cout << "--- Florian (DLE) : " << pfiso_ch << "(" << dle_ch << ") , " 
              << pfiso_em << "(" << dle_em << ") , " << pfiso_nh << "(" << dle_nh << ")" << std::endl;

}

void smurfutilities::DumpSaveTags(const std::string triggerName, 
        const HLTConfigProvider &hltConfig)
{

    for (unsigned int i = 0; i < hltConfig.size(); i++) {

        // get name of ith trigger
        TString hltTrigName(hltConfig.triggerName(i));
        hltTrigName.ToLower();

        // pattern to match
        TString pattern(triggerName);
        pattern.ToLower();

        // match pattern
        TRegexp reg(Form("%s", pattern.Data()), true);

        // if trigger matches
        // then get the names of the save tags filters
        if (hltTrigName.Index(reg) >= 0) {

            std::cout << "[smurfutilities::DumpSaveTags] " << hltConfig.triggerName(i) << std::endl;
            const std::vector<std::string> &modules = hltConfig.saveTagsModules(i);
            for (size_t m = 0; m < modules.size(); ++m) {
                std::cout << "\t" << modules[m] << std::endl;
            }
        }
    }

}


