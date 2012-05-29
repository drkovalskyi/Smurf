#!/bin/bash

#crab -create -submit -cfg crab_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.cfg

#crab -create -submit -cfg crab_DoubleElectronRun2012APromptV1.cfg
crab -create -submit -cfg crab_DoubleMu2012APromptV1.cfg
crab -create -submit -cfg crab_SingleElectronRun2012APromptV1.cfg
crab -create -submit -cfg crab_SingleMuRun2012APromptV1.cfg
crab -create -submit -cfg crab_PhotonsRun2012APromptV1.cfg

crab -create -submit -cfg crab_DoubleElectronRun2012BPromptV1.cfg
crab -create -submit -cfg crab_DoubleMuRun2012BPromptV1.cfg
crab -create -submit -cfg crab_SingleElectronRun2012BPromptV1.cfg
crab -create -submit -cfg crab_SingleMuRun2012BPromptV1.cfg
crab -create -submit -cfg crab_PhotonsRun2012BPromptV1.cfg

