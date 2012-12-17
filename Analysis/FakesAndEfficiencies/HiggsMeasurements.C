
#include "../../Core/LeptonTree.h"

// I feel very bad
const bool useTheBitsIfPossible = true;
//

bool passElectronFO2012(const LeptonTree &leptonTree);
bool passElectronID2012(const LeptonTree &leptonTree);
bool passElectronIso2012(const LeptonTree &leptonTree);

bool passMuonIso2012(const LeptonTree &leptonTree);
bool passMuonID2012(const LeptonTree &leptonTree);

//
// Mode
//

enum Mode {
    MuonTagAndProbe         =1,
    MuonFakeRate            =2,
    ElectronTagAndProbe     =3,
    ElectronFakeRate        =4,
    MuonTagAndProbeMC       =5,
    ElectronTagAndProbeMC   =6
};

//
// tag and probe / FR options
//

enum Option {

    //
    // Muon
    //

    MuonIsoDenominator2012=10,
    MuonIDDenominator2012=11,
    MuonFO0=120,
    MuonFO5=121,    
    MuonFO10=122, 
    MuonFO15=123, 
    MuonFO20=124, 
    MuonFO25=125, 
    MuonFO30=126, 
    MuonFO50=127,
    MuonIDIsoDenominator=13,
    MuonIDIsoRndTagDenominator=14,
    MuonIDIsoRndTagDblNoDzDenominator=15,

    MuonIsoNumerator2011=200,
    MuonIDNumerator2011=201,
    MuonIDIsoNumerator2011=202,
    MuonIsoNumerator2012=20,
    MuonIDNumerator2012=21,
    MuonIDIsoNumerator2012=22,
    MuonTrigSglNumerator2012=23,
    MuonTrigDblLeadNumerator2012=24,
    MuonTrigDblTrailNumerator2012=25,
    MuonTrigDblNumerator2012=26,


    //
    // Electron
    //

    ElectronRecoDenominator = 999,

    ElectronIsoDenominator2012=30,
    ElectronIDDenominator2012=31,
    ElectronFO15=320,
    ElectronFO20=321,
    ElectronFO25=322,
    ElectronFO30=323,
    ElectronFO35=324,
    ElectronFO40=325,
    ElectronFO45=326,
    ElectronFO50=327,
    ElectronIDIsoDenominator=33,
    ElectronIDIsoRndTagDenominator=34,
    ElectronIDIsoRndTagDblNoDzDenominator=35,
    ElectronIsoFODenominator2012=36,

    ElectronIsoNumerator2011=400,
    ElectronIDNumerator2011=401,
    ElectronIDIsoNumerator2011=402,
    ElectronIsoNumerator2012=40,
    ElectronIDNumerator2012=41,
    ElectronIDIsoNumerator2012=42,
    ElectronTrigSglNumerator2012=43,
    ElectronTrigDblLeadNumerator2012=44,
    ElectronTrigDblTrailNumerator2012=45,
    ElectronTrigDblNumerator2012=46,
    ElectronFONumerator=47,

};

//
// triggers
//

struct TriggerResults {
    UInt_t HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_;
    UInt_t HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_;
    UInt_t HLT_Mu17_Mu8_TrailingLeg_tag_;
    UInt_t HLT_Mu17_Mu8_TrailingLeg_probe_;
    UInt_t HLT_Mu17_Mu8_LeadingLeg_tag_;
    UInt_t HLT_Mu17_Mu8_LeadingLeg_probe_;
    UInt_t HLT_Mu17_Mu8_tag_;
    UInt_t HLT_Mu17_Mu8_probe_;
    UInt_t HLT_Mu17_TkMu8_TrailingLeg_tag_;
    UInt_t HLT_Mu17_TkMu8_TrailingLeg_probe_;
    UInt_t HLT_Mu17_TkMu8_LeadingLeg_tag_;
    UInt_t HLT_Mu17_TkMu8_LeadingLeg_probe_;
    UInt_t HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_;
    UInt_t HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_;
    UInt_t HLT_Mu17_TkMu8_tag_;
    UInt_t HLT_Mu17_TkMu8_probe_;
    UInt_t HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_;
    UInt_t HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_;
    UInt_t HLT_IsoMu24_eta2p1_tag_;
    UInt_t HLT_IsoMu24_eta2p1_probe_;
    UInt_t HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_;
    UInt_t HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_;
    UInt_t HLT_Ele17_Ele8_LeadingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_LeadingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_TrailingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_TrailingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_tag_;
    UInt_t HLT_Ele17_Ele8_probe_;
    UInt_t HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_;
    UInt_t HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_;
    UInt_t HLT_Ele27_WP80_tag_;
    UInt_t HLT_Ele27_WP80_probe_;
    UInt_t HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_;
    UInt_t HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_;
    UInt_t HLT_Ele20_SC4_Mass50_LeadingLeg_tag_;
    UInt_t HLT_Ele20_SC4_Mass50_LeadingLeg_probe_;
    UInt_t HLT_Ele20_SC4_Mass50_TrailingLeg_tag_;
    UInt_t HLT_Ele20_SC4_Mass50_TrailingLeg_probe_;
    UInt_t HLT_Ele32_SC17_Mass50_LeadingLeg_tag_;
    UInt_t HLT_Ele32_SC17_Mass50_LeadingLeg_probe_;
    UInt_t HLT_Ele32_SC17_Mass50_TrailingLeg_tag_;
    UInt_t HLT_Ele32_SC17_Mass50_TrailingLeg_probe_;
    UInt_t HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
    UInt_t HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_;
    UInt_t HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_;
    UInt_t HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_;
    UInt_t HLT_Ele8_CaloIdT_TrkIdVL_probe_;
    UInt_t HLT_Mu8_probe_;
    UInt_t HLT_Mu17_probe_;
};

//
// map branches
//

void mapExtraBranches(LeptonTree *leptonTree, TriggerResults &triggerResults) 
{
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag", &triggerResults.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe", &triggerResults.HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_TrailingLeg_tag", &triggerResults.HLT_Mu17_Mu8_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_TrailingLeg_probe", &triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_LeadingLeg_tag", &triggerResults.HLT_Mu17_Mu8_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_LeadingLeg_probe", &triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_tag", &triggerResults.HLT_Mu17_Mu8_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_probe", &triggerResults.HLT_Mu17_Mu8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLeg_tag", &triggerResults.HLT_Mu17_TkMu8_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLeg_probe", &triggerResults.HLT_Mu17_TkMu8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_LeadingLeg_tag", &triggerResults.HLT_Mu17_TkMu8_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_LeadingLeg_probe", &triggerResults.HLT_Mu17_TkMu8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag", &triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe", &triggerResults.HLT_Mu17_TkMu8_TrailingLegTrkFiltered_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_tag", &triggerResults.HLT_Mu17_TkMu8_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_TkMu8_probe", &triggerResults.HLT_Mu17_TkMu8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag", &triggerResults.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe", &triggerResults.HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_tag", &triggerResults.HLT_IsoMu24_eta2p1_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_probe", &triggerResults.HLT_IsoMu24_eta2p1_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_L1sL1DoubleEG137_tag", &triggerResults.HLT_Ele17_Ele8_L1sL1DoubleEG137_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_L1sL1DoubleEG137_probe", &triggerResults.HLT_Ele17_Ele8_L1sL1DoubleEG137_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_LeadingLeg_tag", &triggerResults.HLT_Ele17_Ele8_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_LeadingLeg_probe", &triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_TrailingLeg_tag", &triggerResults.HLT_Ele17_Ele8_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_TrailingLeg_probe", &triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_tag", &triggerResults.HLT_Ele17_Ele8_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_probe", &triggerResults.HLT_Ele17_Ele8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag", &triggerResults.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe", &triggerResults.HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_tag", &triggerResults.HLT_Ele27_WP80_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_probe", &triggerResults.HLT_Ele27_WP80_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_LeadingLeg_tag", &triggerResults.HLT_Ele17_Ele8_Mass50_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_LeadingLeg_probe", &triggerResults.HLT_Ele17_Ele8_Mass50_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_TrailingLeg_tag", &triggerResults.HLT_Ele17_Ele8_Mass50_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_Mass50_TrailingLeg_probe", &triggerResults.HLT_Ele17_Ele8_Mass50_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_LeadingLeg_tag", &triggerResults.HLT_Ele20_SC4_Mass50_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_LeadingLeg_probe", &triggerResults.HLT_Ele20_SC4_Mass50_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_TrailingLeg_tag", &triggerResults.HLT_Ele20_SC4_Mass50_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele20_SC4_Mass50_TrailingLeg_probe", &triggerResults.HLT_Ele20_SC4_Mass50_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_LeadingLeg_tag", &triggerResults.HLT_Ele32_SC17_Mass50_LeadingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_LeadingLeg_probe", &triggerResults.HLT_Ele32_SC17_Mass50_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_TrailingLeg_tag", &triggerResults.HLT_Ele32_SC17_Mass50_TrailingLeg_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele32_SC17_Mass50_TrailingLeg_probe", &triggerResults.HLT_Ele32_SC17_Mass50_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe", &triggerResults.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe", &triggerResults.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe", &triggerResults.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe", &triggerResults.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele8_CaloIdT_TrkIdVL_probe", &triggerResults.HLT_Ele8_CaloIdT_TrkIdVL_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu8_probe", &triggerResults.HLT_Mu8_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_probe", &triggerResults.HLT_Mu17_probe_);
}

//
// utilities
//

bool passElectronID2012(const LeptonTree &leptonTree)
{

    if (useTheBitsIfPossible) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleFOICHEP2012) 
            != LeptonTree::PassEleFOICHEP2012)       return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIDICHEP2012)
            != LeptonTree::PassEleIDICHEP2012)     return false;
    } else {
        if (!passElectronFO2012(leptonTree)) return false;
        float eta = fabs(leptonTree.sceta_);
        float mvaValue = leptonTree.egammaPOG2012MVA_;
        if (leptonTree.probe_.Pt() > 20.0) {
            if (eta < 0.8                 && mvaValue <= 0.94) return false;
            if (eta >= 0.8 && eta < 1.479 && mvaValue <= 0.85) return false;
            if (eta >= 1.479              && mvaValue <= 0.92) return false;
        } else {
            if (eta < 0.8                 && mvaValue <= 0.00) return false;
            if (eta >= 0.8 && eta < 1.479 && mvaValue <= 0.10) return false;
            if (eta >= 1.479              && mvaValue <= 0.62) return false;
        }
    }

    return true;
}

bool passElectronIso2012(const LeptonTree &leptonTree)
{

        float EffectiveArea = 0.0;
          if (fabs(leptonTree.sceta_) >= 0.0 && fabs(leptonTree.sceta_) < 1.0 ) EffectiveArea = 0.176;
          if (fabs(leptonTree.sceta_) >= 1.0 && fabs(leptonTree.sceta_) < 1.479 ) EffectiveArea = 0.206;
          if (fabs(leptonTree.sceta_) >= 1.479 && fabs(leptonTree.sceta_) < 2.0 ) EffectiveArea = 0.094;
          if (fabs(leptonTree.sceta_) >= 2.2 && fabs(leptonTree.sceta_) < 2.2 ) EffectiveArea = 0.172;
          if (fabs(leptonTree.sceta_) >= 2.3 && fabs(leptonTree.sceta_) < 2.3 ) EffectiveArea = 0.244;
          if (fabs(leptonTree.sceta_) >= 2.4 && fabs(leptonTree.sceta_) < 2.4 ) EffectiveArea = 0.333;
          if (fabs(leptonTree.sceta_) >= 2.4 ) EffectiveArea = 0.348;

    if (useTheBitsIfPossible) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIsoICHEP2012)
            != LeptonTree::PassEleIsoICHEP2012)     return false;
    } else {
        if ((leptonTree.pfchiso04_ + TMath::Max(double(0.0), leptonTree.pfemiso04_ + leptonTree.pfnhiso04_
            - EffectiveArea * TMath::Max(double(0.0), leptonTree.rhoIsoAll_)))/leptonTree.probe_.Pt() > 0.15) return false;
    }

    return true;
}

bool passElectronFO2012(const LeptonTree &leptonTree)
{

    if (useTheBitsIfPossible) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleFOICHEP2012)
                != LeptonTree::PassEleFOICHEP2012)          return false;

    } else {
        float pt = leptonTree.probe_.Pt();
        float d0 = leptonTree.d0vtx_;
        float dz = leptonTree.dzvtx_;
        if (leptonTree.trkiso_/pt               > 0.2)      return false; 
        if (leptonTree.hcaliso_/pt              > 0.2)      return false;
        if (fabs(d0) > 0.02)                                return false;
        if (fabs(dz) > 0.1)                                 return false;

        unsigned int mhits = leptonTree.mhit_;
        bool conv = leptonTree.vfitprob_;
        if (mhits > 0)      return false;
        if (conv)           return false;

        if (fabs(leptonTree.sceta_) < 1.479) {
            if (leptonTree.sieie_               > 0.01)  return false;
            if (fabs(leptonTree.detain_)        > 0.007) return false;
            if (fabs(leptonTree.dphiin_)        > 0.15)  return false;
            if (leptonTree.hoe_                 > 0.12)  return false;
            if ((leptonTree.ecaliso_ - 1.0)/pt  > 0.2)   return false;
        } else {
            if (leptonTree.sieie_               > 0.03)  return false;
            if (fabs(leptonTree.detain_)        > 0.009) return false;
            if (fabs(leptonTree.dphiin_)        > 0.10)  return false;
            if (leptonTree.hoe_                 > 0.10)  return false;
            if ((leptonTree.ecaliso_)/pt        > 0.2)   return false; 
        }
    }

    return true;

}

bool passMuonIso2012(const LeptonTree &leptonTree)
{

    if (useTheBitsIfPossible) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsoICHEP2012)
            != LeptonTree::PassMuIsoICHEP2012)                                 return false;
    } else {
        float mvaValue = leptonTree.muonHZZ2012IsoRingsMVA_;
        float eta = fabs(leptonTree.probe_.Eta());
        if (leptonTree.probe_.Pt() > 20.0) {
            if (eta < 1.479  && mvaValue <= 0.82)   return false;
            if (eta >= 1.479 && mvaValue <= 0.86)   return false;
        } else {
            if (eta < 1.479  && mvaValue <= 0.86)   return false;
            if (eta >= 1.479 && mvaValue <= 0.82)   return false;
        }
    }

    return true;
}

bool passMuonID2012(const LeptonTree &leptonTree)
{
    if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIDICHEP2012)
        != LeptonTree::PassMuIDICHEP2012)                                  return false;
    return true;
}

//
// probes
//

bool isProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults, const Mode mode)
{

    bool isMC = false;
    if (mode == ElectronTagAndProbeMC || mode == MuonTagAndProbeMC) isMC = true;

    //
    // electrons
    //

    if (option == ElectronRecoDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
    }

    if (option == ElectronFO15 || option == ElectronFO20 || option == ElectronFO25
            || option == ElectronFO30 || option == ElectronFO35 || option == ElectronFO40
            || option == ElectronFO45 || option == ElectronFO50 ) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleFOICHEP2012) 
                != LeptonTree::PassEleFOICHEP2012)          return false;
        if (triggerResults.HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_ == 0 
                && leptonTree.probe_.Pt() < 20.0)           return false;
        if (triggerResults.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_probe_ == 0 
                && leptonTree.probe_.Pt() >= 20.0)          return false;
        if (sqrt(pow(leptonTree.probe_.Eta() - leptonTree.jet1_.Eta(), 2) + pow(leptonTree.probe_.Phi() - leptonTree.jet1_.Phi(), 2)) < 1.0) return false;
    }
    if (option == ElectronFO15 && leptonTree.jet1_.Pt() < 15.0)                        return false;
    if (option == ElectronFO20 && leptonTree.jet1_.Pt() < 20.0)                        return false;
    if (option == ElectronFO25 && leptonTree.jet1_.Pt() < 25.0)                        return false;
    if (option == ElectronFO30 && leptonTree.jet1_.Pt() < 30.0)                        return false;
    if (option == ElectronFO35 && leptonTree.jet1_.Pt() < 35.0)                        return false;
    if (option == ElectronFO40 && leptonTree.jet1_.Pt() < 40.0)                        return false;
    if (option == ElectronFO45 && leptonTree.jet1_.Pt() < 45.0)                        return false;
    if (option == ElectronFO50 && leptonTree.jet1_.Pt() < 50.0)                        return false;


    if (option == ElectronIsoDenominator2012) {
        if (leptonTree.tag_.Pt() < 32.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                       return false;
        if (!passElectronIso2012(leptonTree))       return false;
    }

    if (option == ElectronIsoFODenominator2012) {
        if (leptonTree.tag_.Pt() < 32.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                                    return false;
        if (!passElectronFO2012(leptonTree))                                  return false;
        if (!passElectronIso2012(leptonTree))                                       return false;
    }

    if (option == ElectronIDDenominator2012) {
        if (leptonTree.tag_.Pt() < 32.0)                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                                    return false;
        if (!passElectronID2012(leptonTree))                                         return false;
    }

    if (option == ElectronIDIsoDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if (!passElectronIso2012(leptonTree))                                       return false;
        if (!passElectronID2012(leptonTree))                                         return false;
    }

    if (option == ElectronIDIsoRndTagDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if (!passElectronIso2012(leptonTree))                                       return false;
        if (!passElectronID2012(leptonTree))                                         return false;
        //if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
        //    (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if (leptonTree.tagAndProbeIsRandom_ != 1)   return false;
    }

    if (option == ElectronIDIsoRndTagDblNoDzDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if (!passElectronIso2012(leptonTree))                                       return false;
        if (!passElectronID2012(leptonTree))                                         return false;
        //if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
        //    (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if (leptonTree.tagAndProbeIsRandom_ != 1)   return false;
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0 && !isMC)        return false; 
    }

    //
    // muons
    //

    if (option == MuonFO5 || option == MuonFO10 || option == MuonFO15
            || option == MuonFO20 || option == MuonFO25 || option == MuonFO30 || option == MuonFO50) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuFOICHEP2012) != LeptonTree::PassMuFOICHEP2012)       return false;
        if (triggerResults.HLT_Mu8_probe_ == 0 && leptonTree.probe_.Pt() < 20.0)    return false;
        if (triggerResults.HLT_Mu17_probe_ == 0 && leptonTree.probe_.Pt() >= 20.0)  return false;
        if (sqrt(pow(leptonTree.probe_.Eta() - leptonTree.jet1_.Eta(), 2) + pow(leptonTree.probe_.Phi() - leptonTree.jet1_.Phi(), 2)) < 1.0) return false;
    }
    if (option == MuonFO5  && leptonTree.jet1_.Pt() < 5.0)                           return false;
    if (option == MuonFO10 && leptonTree.jet1_.Pt() < 10.0)                          return false;
    if (option == MuonFO15 && leptonTree.jet1_.Pt() < 15.0)                          return false;
    if (option == MuonFO20 && leptonTree.jet1_.Pt() < 20.0)                          return false;
    if (option == MuonFO25 && leptonTree.jet1_.Pt() < 25.0)                          return false;
    if (option == MuonFO30 && leptonTree.jet1_.Pt() < 30.0)                          return false;
    if (option == MuonFO50 && leptonTree.jet1_.Pt() < 50.0)                          return false;


    if (option == MuonIsoDenominator2012) {
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                   return false;
        if (!passMuonIso2012(leptonTree))                                           return false;
    }

    if (option == MuonIDDenominator2012) {
        if (leptonTree.tag_.Pt() < 30.0)                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                    return false;
        if (!passMuonID2012(leptonTree))                                           return false;
    }

    if (option == MuonIDIsoDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                           return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                  return false;
        if (!passMuonIso2012(leptonTree))                                           return false;
        if (!passMuonID2012(leptonTree))                                           return false;
    }

    if (option == MuonIDIsoRndTagDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                           return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                  return false;
        if (!passMuonIso2012(leptonTree))                                           return false;
        if (!passMuonID2012(leptonTree))                                           return false;
        //if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
        //    (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if (leptonTree.tagAndProbeIsRandom_ != 1)   return false;
    }

    if (option == MuonIDIsoRndTagDblNoDzDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                           return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                  return false;
        if (!passMuonIso2012(leptonTree))                                           return false;
        if (!passMuonID2012(leptonTree))                                           return false;
        //if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
        //    (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if (leptonTree.tagAndProbeIsRandom_ != 1)   return false;
        if ((triggerResults.HLT_Mu17_TkMu8_TrailingLeg_probe_ == 0
                && triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) && !isMC) return false;
    }

    return true;
}

bool isPassProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults, const Mode mode)
{

    bool isMC = false;
    if (mode == ElectronTagAndProbeMC || mode == MuonTagAndProbeMC) isMC = true;

    //
    // Electrons
    //

    if (option == ElectronIsoNumerator2012) {
        if (!passElectronIso2012(leptonTree))   return false;
    }

    if (option == ElectronIDNumerator2012) {
        if (!passElectronID2012(leptonTree))    return false;
    }

    if (option == ElectronIDIsoNumerator2012) {
        if (!passElectronID2012(leptonTree))    return false;
        if (!passElectronIso2012(leptonTree))   return false;
    }

    if (option == ElectronIsoNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleIso) != LeptonTree::PassEleIso)       return false;
    }

    if (option == ElectronIDNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassEleID) != LeptonTree::PassEleID)       return false;
    }

    if (option == ElectronIDIsoNumerator2011) {
        unsigned int IDIso = LeptonTree::PassEleIso | LeptonTree::PassEleID;
        if ((leptonTree.leptonSelection_ & IDIso) != IDIso)       return false;
    }

    if (option == ElectronTrigSglNumerator2012) {
        if (triggerResults.HLT_Ele27_WP80_probe_ == 0 && !isMC)   return false;
    }   
        
    if (option == ElectronTrigDblLeadNumerator2012) {
        if (triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_ == 0 && !isMC) return false;
    }   
        
    if (option == ElectronTrigDblTrailNumerator2012) {
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0 && !isMC) return false;
    }

    if (option == ElectronTrigDblNumerator2012) {
        if (triggerResults.HLT_Ele17_Ele8_probe_ == 0 && !isMC) return false;
    }

    if (option == ElectronFONumerator) {
        if (!passElectronFO2012(leptonTree))               return false;
    }

    //
    // Muons
    //

    if (option == MuonIsoNumerator2012) {
        if (!passMuonIso2012(leptonTree))   return false;
    }

    if (option == MuonIDNumerator2012) {
        if (!passMuonID2012(leptonTree))   return false;
    }

    if (option == MuonIDIsoNumerator2012) {
        if (!passMuonIso2012(leptonTree))   return false;
        if (!passMuonID2012(leptonTree))   return false;
    }

    if (option == MuonIsoNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIso) != LeptonTree::PassMuIso)       return false;
    }

    if (option == MuonIDNumerator2011) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuID) != LeptonTree::PassMuID)       return false;
    }

    if (option == MuonIDIsoNumerator2011) {
        unsigned int IDIso = LeptonTree::PassMuIso | LeptonTree::PassMuID;
        if ((leptonTree.leptonSelection_ & IDIso) != IDIso)       return false;
    }

    if (option == MuonTrigSglNumerator2012) {
        if (triggerResults.HLT_IsoMu24_eta2p1_probe_ == 0 && !isMC)   return false;
    }

    if (option == MuonTrigDblLeadNumerator2012) {
        if ((triggerResults.HLT_Mu17_TkMu8_LeadingLeg_probe_ == 0
                && triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_ == 0) && !isMC) return false;
    }

    if (option == MuonTrigDblTrailNumerator2012) {
        if ((triggerResults.HLT_Mu17_TkMu8_TrailingLeg_probe_ == 0
                && triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) && !isMC) return false;
    }

    if (option == MuonTrigDblNumerator2012) {
        if ((triggerResults.HLT_Mu17_TkMu8_probe_ == 0
                && triggerResults.HLT_Mu17_Mu8_probe_ == 0) && !isMC) return false;
    }

    return true;
}

