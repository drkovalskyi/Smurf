#!/bin/bash

RUNLIST=Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON.txt
JOB=$1

compareJSON.py --and ${JOB}/res/lumiSummary.json ${RUNLIST} > processed_${JOB1}.json
lumiCalc2.py -i processed_${JOB1}.json overview

