#!/bin/bash

#mergeJSON.py RemoteGlidein_DoubleMu_Run2012A-13Jul2012-v1_AOD_190456_193621/res/lumiSummary.json          \
#    RemoteGlidein_DoubleMu_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/res/lumiSummary.json  \
#    RemoteGlidein_DoubleMu_Run2012B-13Jul2012-v4_AOD_193834_196531/res/lumiSummary.json          \
#    RemoteGlidein_DoubleMu_Run2012C-24Aug2012-v1_AOD_198022_198523/res/lumiSummary.json          \
#    RemoteGlidein_DoubleMu_Run2012C-PromptReco-v2_AOD_198934_203755/res/lumiSummary.json         \
#    RemoteGlidein_DoubleMu_Run2012D-PromptReco-v1_AOD_203773_208299/res/lumiSummary.json > DoubleMu.json

#compareJSON.py --and DoubleMu.json ../runlists/HCP.json > DoubleMu_good.json

mergeJSON.py RemoteGlidein_DoubleElectron_Run2012A-13Jul2012-v1_AOD_190456_193621/res/lumiSummary.json          \
    RemoteGlidein_DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/res/lumiSummary.json  \
    RemoteGlidein_DoubleElectron_Run2012B-13Jul2012-v1_AOD_193834_196531/res/lumiSummary.json          \
    RemoteGlidein_DoubleElectron_Run2012C-24Aug2012-v1_AOD_198022_198523/res/lumiSummary.json          \
    RemoteGlidein_DoubleElectron_Run2012C-PromptReco-v2_AOD_198934_203755/res/lumiSummary.json         \
    RemoteGlidein_DoubleElectron_Run2012D-PromptReco-v1_AOD_203773_208299/res/lumiSummary.json > DoubleElectron.json

compareJSON.py --and DoubleElectron.json ../runlists/HCP.json > DoubleElectron_good.json

lumiCalc2.py -i DoubleElectron_good.json recorded --hltpath HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*
lumiCalc2.py -i DoubleElectron_good.json recorded --hltpath HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v*

