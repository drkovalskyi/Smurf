#! /usr/bin/env python

from ROOT import TFile
from array import array
import sys, os, math, argparse

from cutDefinitions import * 

# weights for MC sample  (no fake spillage correction)
def computeWeights(sample, mass=140): 

    # functions
    import scaleFactors

    # PU reweight factors
    puReweightFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/puReweighting.root') 
    puReweight = puReweightFile.Get('puWeights')
    
    # lepton scale factors
    # 41X
    leptonEfficiencyFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/efficiency_41X.root')
    # $2X
    #leptonEfficiencyFile = TFile('/data/smurf/mzanetti/data/42X/auxiliar/efficiency_42X.root')

    # trigger efficiencies
    triggerEfficiencyHistos = {'single_e':leptonEfficiencyFile.Get('h2_results_electron_single'),
                               'double_e':leptonEfficiencyFile.Get('h2_results_electron_double'),
                               'single_m':leptonEfficiencyFile.Get('h2_results_muon_single'),
                               'double_m':leptonEfficiencyFile.Get('h2_results_muon_double')}
    # clones for the moment
    triggerEfficiencyHistos['cross_e'] =  triggerEfficiencyHistos['double_e']  
    triggerEfficiencyHistos['cross_m'] =  triggerEfficiencyHistos['double_m']  
    # lepton efficiencies
    leptonEfficiencyHistos = {'muon':leptonEfficiencyFile.Get('h2_results_muon_selection'),
                              'electron':leptonEfficiencyFile.Get('h2_results_electron_selection')}

    # higgs pT factors
    higgsPtKFactorFile = TFile('/data/smurf/data/EPS/auxiliar/ggHWW_KFactors_PowhegToHQT.root')
    higgsPtKFactorHisto = higgsPtKFactorFile.Get('KFactor_PowhegToHQT_mH'+str(mass));
    
    # reference file
    refFile = TFile(sample)
    tree = refFile.Get('tree')
    
    # output File
    newFile = TFile(sample[:sample.find('.root')]+'_weighted'+'.root','RECREATE')
    newTree = tree.CloneTree(-1,'fast')

    # the variable to be added to the tree
    weight = array('f',[1.0]) 
    weightBranch = newTree.Branch('weight',weight,'weight/F') 

    for event in range(0, tree.GetEntries()):

        if event%10000==0: print event 

        tree.GetEntry(event)

        weight[0] = 1


        # PU (watch out the bins!)
        weight[0] = weight[0]*puReweight.GetBinContent(min(tree.nvtx,19))

        # leptons scale factor
        weight[0] = weight[0]*scaleFactors.leptonEfficiency(math.fabs(tree.lid1), math.fabs(tree.lep1.eta()), tree.lep1.pt(), leptonEfficiencyHistos)
        weight[0] = weight[0]*scaleFactors.leptonEfficiency(math.fabs(tree.lid2), math.fabs(tree.lep2.eta()), tree.lep2.pt(), leptonEfficiencyHistos)

        # trigger efficiency
        weight[0] = weight[0]*scaleFactors.triggerEfficiency([math.fabs(tree.lid1), math.fabs(tree.lid2)],
                                                             [math.fabs(tree.lep1.eta()), math.fabs(tree.lep2.eta())],
                                                             [tree.lep1.pt(), tree.lep2.pt()],
                                                             triggerEfficiencyHistos)

        # higgs pt 
        if 'hww' in sample:
            if event%10000==0: print 'correcting higgs Pt'            
            weight[0] = weight[0]*higgsPtKFactorHisto.GetBinContent( higgsPtKFactorHisto.GetXaxis().FindFixBin(tree.higgsPt));


        weightBranch.Fill()

    newFile.cd()
    newTree.Write()
    newFile.Close()



# compute fake weights
def computeFakeWeights(sample, isData, luminosity=1.092): 

    # functions
    import scaleFactors

    # PU reweight factors
    puReweightFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/puReweighting.root')
    puReweight = puReweightFile.Get('puWeights')
    
    
    # Lepton Scale Factors
    # 41X
    leptonEfficiencyFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/efficiency_41X.root')
    # $2X
    #leptonEfficiencyFile = TFile('/data/smurf/mzanetti/data/42X/auxiliar/efficiency_42X.root')

    # trigger efficiencies
    triggerEfficiencyHistos = {'single_e':leptonEfficiencyFile.Get('h2_results_electron_single'),
                               'double_e':leptonEfficiencyFile.Get('h2_results_electron_double'),
                               'single_m':leptonEfficiencyFile.Get('h2_results_muon_single'),
                               'double_m':leptonEfficiencyFile.Get('h2_results_muon_double')}
    # clones for the moment
    triggerEfficiencyHistos['cross_e'] =  triggerEfficiencyHistos['double_e']  
    triggerEfficiencyHistos['cross_m'] =  triggerEfficiencyHistos['double_m']  
    # lepton efficiencies
    leptonEfficiencyHistos = {'muon':leptonEfficiencyFile.Get('h2_results_muon_selection'),
                              'electron':leptonEfficiencyFile.Get('h2_results_electron_selection')}


    # Fake Rates
    # Si
    fakeRateFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/FakeRates_SmurfV6.root')
    fakeRateHistos = {'muon':fakeRateFile.Get('MuonFakeRate_M2_ptThreshold15_PtEta'),'electron': fakeRateFile.Get('ElectronFakeRate_V4_ptThreshold35_PtEta')}

    # Dima
    #fakeRateMuFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/ww_mu_fr.root')
    #fakeRateElFile = TFile('/data/smurf/mzanetti/data/EPS/auxiliar/ww_el_fr.root')
    #fakeRateHistos = {'muon':fakeRateMuFile.Get('mu_fr_m2_15'),'electron': fakeRateElFile.Get('el_fr_v4_35')}    


    # reference File
    refFile = TFile(sample)
    tree = refFile.Get('tree')
    
    # output File

    if isData: newFile = TFile(sample[:sample.rfind('/')+1] + 'fakes.root','RECREATE')
    else: newFile = TFile(sample[:sample.find('.root')]+'_spillage'+'.root','RECREATE')
    newTree = tree.CloneTree(-1,'fast')

    # the variable to be added to the tree
    weight = array('f',[1.0]) 
    weightBranch = newTree.Branch('weight',weight,'weight/F') 

    for event in range(0, tree.GetEntries()):

        if event%10000==0: print event 

        tree.GetEntry(event)

        weight[0] = 0

        # fake rate (both data and MC)
        nMuonFakes = [0,0]
        if (tree.cuts & cuts['Lep1LooseMuV2']) == cuts['Lep1LooseMuV2']  and (tree.cuts & cuts['Lep1FullSelection']) != cuts['Lep1FullSelection']: nMuonFakes[0]+=1
        if (tree.cuts & cuts['Lep2LooseMuV2']) == cuts['Lep2LooseMuV2']  and (tree.cuts & cuts['Lep2FullSelection']) != cuts['Lep2FullSelection']: nMuonFakes[1]+=1
        nElectronFakes = [0,0]
        if (tree.cuts & cuts['Lep1LooseEleV4']) == cuts['Lep1LooseEleV4'] and (tree.cuts & cuts['Lep1FullSelection']) != cuts['Lep1FullSelection']: nElectronFakes[0]+=1
        if (tree.cuts & cuts['Lep2LooseEleV4']) == cuts['Lep2LooseEleV4'] and (tree.cuts & cuts['Lep2FullSelection']) != cuts['Lep2FullSelection']: nElectronFakes[1]+=1

        # consider events with exactly one fake
        if (nMuonFakes[0]+nMuonFakes[1]+nElectronFakes[0]+nElectronFakes[1]) == 1:

            # negative to compute spillage
            weight[0] = 1 if isData else -1
            
            if nMuonFakes[0] == 1: weight[0] = weight[0]*scaleFactors.fakeRate(tree.lep1.eta(), tree.lep1.pt(), fakeRateHistos['muon'])
            if nMuonFakes[1] == 1: weight[0] = weight[0]*scaleFactors.fakeRate(tree.lep2.eta(), tree.lep2.pt(), fakeRateHistos['muon'])            
            if nElectronFakes[0] == 1: weight[0] = weight[0]*scaleFactors.fakeRate(tree.lep1.eta(), tree.lep1.pt(), fakeRateHistos['electron'])
            if nElectronFakes[1] == 1: weight[0] = weight[0]*scaleFactors.fakeRate(tree.lep2.eta(), tree.lep2.pt(), fakeRateHistos['electron'])            

            # this is only for MC
            if not isData: 
                # PU (watch out the bins!)
                weight[0] = weight[0]*puReweight.GetBinContent(min(tree.nvtx,19))

                # leptons scale factor
                weight[0] = weight[0]*scaleFactors.leptonEfficiency(math.fabs(tree.lid1), math.fabs(tree.lep1.eta()), tree.lep1.pt(), leptonEfficiencyHistos)
                weight[0] = weight[0]*scaleFactors.leptonEfficiency(math.fabs(tree.lid2), math.fabs(tree.lep2.eta()), tree.lep2.pt(), leptonEfficiencyHistos)

                # trigger efficiency
                weight[0] = weight[0]*scaleFactors.triggerEfficiency([math.fabs(tree.lid1), math.fabs(tree.lid2)],
                                                                 [math.fabs(tree.lep1.eta()), math.fabs(tree.lep2.eta())],
                                                                 [tree.lep1.pt(), tree.lep2.pt()],
                                                                 triggerEfficiencyHistos)

                # add luminosity in the case spillage is been computed
                weight[0] = weight[0]*luminosity
            
        weightBranch.Fill()

    newFile.cd()
    newTree.Write()
    newFile.Close()



