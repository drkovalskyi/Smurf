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

trigger::TriggerObjectCollection smurfutilities::GetTriggerObjects(const std::string triggerName, const std::string filterName,
        const std::string processName, const HLTConfigProvider &hltConfig, const edm::TriggerResults* triggerResults,
        const edm::Handle<trigger::TriggerEvent> &triggerEvent,
        const trigger::TriggerObjectCollection &allObjects)
{

    // objects matching trigger
    trigger::TriggerObjectCollection selectedObjects;

    // loop on triggers
    for (unsigned int i = 0; i < hltConfig.size(); i++) {

        // does this trigger pass
        if(!triggerResults->accept(i)) continue;

        // get name of ith trigger
        TString hltTrigName(hltConfig.triggerName(i));
        hltTrigName.ToLower();

        // pattern to match
        TString pattern(triggerName);
        pattern.ToLower();

        // match pattern
        TRegexp reg(Form("%s", pattern.Data()), true);

        // if trigger matches
        // then look for the objects corresponding to
        // the specified filter name
        if (hltTrigName.Index(reg) >= 0) {

            edm::InputTag filterNameTag(filterName, "", processName);
            size_t filterIndex = triggerEvent->filterIndex(filterNameTag);

            if (filterIndex < triggerEvent->sizeFilters()) {
                const trigger::Keys &keys = triggerEvent->filterKeys(filterIndex);
                for (size_t j = 0; j < keys.size(); j++) {
                    trigger::TriggerObject foundObject = allObjects[keys[j]];
                    selectedObjects.push_back(foundObject);
                }
            }
        }
    }

    return selectedObjects;

}





