#! /usr/bin/env python

from ROOT import TCut
import sys, os, math, argparse

##########################
# ntuple bitmap definition
##########################
cuts = {}
cuts['BaseLine']= 1<<0  
cuts['ChargeMatch']= 1<<1  
cuts['Lep1FullSelection']= 1<<2  
cuts['Lep1LooseEleV1']= 1<<3  
cuts['Lep1LooseEleV2']= 1<<4  
cuts['Lep1LooseEleV3']= 1<<5  
cuts['Lep1LooseEleV4']= 1<<6  
cuts['Lep1LooseMuV1']= 1<<7  
cuts['Lep1LooseMuV2']= 1<<8  
cuts['Lep2FullSelection']= 1<<9  
cuts['Lep2LooseEleV1']= 1<<10 
cuts['Lep2LooseEleV2']= 1<<11 
cuts['Lep2LooseEleV3']= 1<<12 
cuts['Lep2LooseEleV4']= 1<<13 
cuts['Lep2LooseMuV1']= 1<<14 
cuts['Lep2LooseMuV2']= 1<<15 
cuts['FullMET']= 1<<16 
cuts['ZVeto']= 1<<17 
cuts['TopTag']= 1<<18 
cuts['TopVeto']= 1<<19 
cuts['OneBJet']= 1<<20 
cuts['TopTagNotInJets']= 1<<21 
cuts['ExtraLeptonVeto']= 1<<22 
cuts['Lep3FullSelection']= 1<<23 
cuts['Lep3LooseEleV1']= 1<<24 
cuts['Lep3LooseEleV2']= 1<<25 
cuts['Lep3LooseEleV3']= 1<<26 
cuts['Lep3LooseEleV4']= 1<<27 
cuts['Lep3LooseMuV1']= 1<<28 
cuts['Lep3LooseMuV2']= 1<<29  

selections = {}
# WW selections
selections['PreSelection'] = cuts['BaseLine']|cuts['ChargeMatch']|cuts['Lep1FullSelection']|cuts['Lep2FullSelection']|cuts['ExtraLeptonVeto']
selections['WWSelection'] = selections['PreSelection']|cuts['ZVeto']|cuts['FullMET']|cuts['TopVeto'] 
selections['WWSelectionNoMet'] = selections['PreSelection']|cuts['ZVeto']|cuts['TopVeto'] 

# WW fake selection
selections['PreSelectionFake'] = cuts['BaseLine']|cuts['ChargeMatch']|cuts['ExtraLeptonVeto']
selections['WWSelectionFake'] = selections['PreSelectionFake']|cuts['ZVeto']|cuts['FullMET']|cuts['TopVeto'] 
selections['WWSelectionNoMetFake'] = selections['PreSelectionFake']|cuts['ZVeto']|cuts['TopVeto']

# DY selections
selections['ZSelection'] = selections['PreSelection']|cuts['TopVeto']
selections['AntiZSelection'] = selections['PreSelection']|cuts['ZVeto']|cuts['TopVeto']

# DY selections Fake
selections['ZSelectionFake'] = selections['PreSelectionFake']|cuts['TopVeto']
selections['AntiZSelectionFake'] = selections['PreSelectionFake']|cuts['ZVeto']|cuts['TopVeto']

# top selections
selections['topSelection'] = selections['PreSelection']|cuts['ZVeto']|cuts['FullMET']
selections['topDenumerator'] = selections['topSelection']|cuts['OneBJet']
selections['topNumerator'] = selections['topDenumerator']|cuts['TopTagNotInJets']
selections['topTagging'] = selections['topSelection']|cuts['TopTagNotInJets']#cuts['TopTag']

# top selections fakes
selections['topSelectionFake'] = selections['PreSelectionFake']|cuts['ZVeto']|cuts['FullMET']
selections['topDenumeratorFake'] = selections['topSelectionFake']|cuts['OneBJet']
selections['topNumeratorFake'] = selections['topDenumeratorFake']|cuts['TopTagNotInJets']
selections['topTaggingFake'] = selections['topSelectionFake']|cuts['TopTagNotInJets']#cuts['TopTag']


higgsCuts = {}
higgsCuts[120] = {'ptMax':20, 'ptMin':10,'mll':40, 'dPhi':115, 'mt':[70,120]} 
higgsCuts[130] = {'ptMax':25, 'ptMin':10,'mll':45, 'dPhi':90, 'mt':[75,125]} 
higgsCuts[140] = {'ptMax':25, 'ptMin':15,'mll':45, 'dPhi':90, 'mt':[80,130]} 
higgsCuts[150] = {'ptMax':27, 'ptMin':25,'mll':50, 'dPhi':90, 'mt':[80,150]} 
higgsCuts[160] = {'ptMax':30, 'ptMin':25,'mll':50, 'dPhi':60, 'mt':[90,160]} 
higgsCuts[170] = {'ptMax':34, 'ptMin':25,'mll':50, 'dPhi':60, 'mt':[20,170]} 
higgsCuts[180] = {'ptMax':36, 'ptMin':25,'mll':60, 'dPhi':70, 'mt':[20,180]} 
higgsCuts[190] = {'ptMax':38, 'ptMin':25,'mll':80, 'dPhi':90, 'mt':[20,190]} 
higgsCuts[200] = {'ptMax':40, 'ptMin':25,'mll':90, 'dPhi':100, 'mt':[20,200]} 
higgsCuts[250] = {'ptMax':55, 'ptMin':25,'mll':150, 'dPhi':140, 'mt':[20,250]} 
higgsCuts[300] = {'ptMax':70, 'ptMin':25,'mll':200, 'dPhi':175, 'mt':[20,300]} 
higgsCuts[350] = {'ptMax':80, 'ptMin':25,'mll':250, 'dPhi':175, 'mt':[20,350]} 
higgsCuts[400] = {'ptMax':90, 'ptMin':25,'mll':300, 'dPhi':175, 'mt':[20,400]} 
higgsCuts[450] = {'ptMax':110,'ptMin':25,'mll':350, 'dPhi':175, 'mt':[20,450]} 
higgsCuts[500] = {'ptMax':120,'ptMin':25,'mll':400, 'dPhi':175, 'mt':[20,500]} 
higgsCuts[550] = {'ptMax':130,'ptMin':25,'mll':450, 'dPhi':175, 'mt':[20,550]} 
higgsCuts[600] = {'ptMax':140,'ptMin':25,'mll':500, 'dPhi':175, 'mt':[20,600]} 


