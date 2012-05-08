#include "Smurf/ProcessingAndSkimming/interface/Utilities.h"
#include "Smurf/ProcessingAndSkimming/interface/Selections.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TString.h"
#include "TRegexp.h"

#include <vector>


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

unsigned int smurfutilities::MatchTriggerObject(const edm::Event &iEvent, const edm::EventSetup &iSetup,
        const std::string triggerName, const std::string filterName,
        const std::string processName, const HLTConfigProvider &hltConfig, const edm::TriggerResults* triggerResults,
        const trigger::TriggerEvent *triggerEvent,
        const trigger::TriggerObjectCollection &allObjects,
        const LorentzVector &offlineObject)
{

    unsigned int prescale = 0;

    // loop on triggers
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
        // then look for the objects corresponding to
        // the specified filter name
        if (hltTrigName.Index(reg) >= 0) {

            edm::InputTag filterNameTag(filterName, "", processName);
            if (filterName == "") {
                const std::vector<std::string> &modules = hltConfig.saveTagsModules(i);
                filterNameTag = edm::InputTag(modules.back(), "", processName);
            }
            else {
                filterNameTag = edm::InputTag(filterName, "", processName);
            }

            size_t filterIndex = triggerEvent->filterIndex(filterNameTag);

            if (filterIndex < triggerEvent->sizeFilters()) {
                const trigger::Keys &keys = triggerEvent->filterKeys(filterIndex);
                for (size_t j = 0; j < keys.size(); j++) {
                    trigger::TriggerObject foundObject = allObjects[keys[j]];
                    if (deltaR(foundObject.eta(), foundObject.phi(), offlineObject.eta(), offlineObject.phi()) < 0.2) {
                         prescale = hltConfig.prescaleValue(iEvent, iSetup, hltConfig.triggerName(i));
                         return prescale;
                    }
                }
            }
        }
    }

    assert(prescale == 0);
    return prescale;

}

bool smurfutilities::MatchTriggerObject(const trigger::TriggerObjectCollection &triggerObjects, const LorentzVector &offlineObject)
{
    for (unsigned int i = 0; i < triggerObjects.size(); ++i) {
        if (deltaR(triggerObjects[i].eta(), triggerObjects[i].phi(), offlineObject.eta(), offlineObject.phi()) < 0.2) 
            return true ;
    }
    return false;
}

float smurfutilities::MatchGenParticle(const reco::GenParticleCollection &genParticleCollection, const LorentzVector &object, 
        const int pdgId, const int status)
{

    float dRMin = 999.9;
    reco::GenParticleCollection::const_iterator genParticle;
    for (genParticle = genParticleCollection.begin(); genParticle != genParticleCollection.end(); ++genParticle) {
        if (abs(genParticle->pdgId()) != pdgId) continue;
        if (genParticle->status() != status)    continue;
        float dR = deltaR(genParticle->eta(), genParticle->phi(), object.eta(), object.phi());
        if (dR < dRMin) dRMin = dR;
    }

    return dRMin;

}


float smurfutilities::Mt(const float &pt1, const float &pt2, const float &dphi)
{
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

