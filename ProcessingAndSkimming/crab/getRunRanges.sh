#!/bin/bash

# input parameter validation
if [ ! $# -eq 3 ]; then
    echo "./getRunRanges.sh FIRSTRUN  DSET  GTAG
    FIRSTRUN - run number of the first run. Enter zero to use first run in dataset
    DSET - e.g. Run2012A-13Jul2012-v1
    GTAG - e.g. FT_53_V6_AN2::All"
    exit 1
fi

# input parameters
FIRSTRUN=$1
DSET=$2
GTAG=$3

#DSETS="/SingleElectron/${DSET}/AOD /SingleMu/${DSET}/AOD /DoubleElectron/${DSET}/AOD /DoubleMu/${DSET}/AOD"
DSETS="/DoubleElectron/${DSET}/AOD /DoubleMu/${DSET}/AOD"
for DATASET in $DSETS
do

    # annoying
    if [ "${DATASET}" == "/DoubleMu/Run2012B-13Jul2012-v1/AOD" ]; then
        DATASET="/DoubleMu/Run2012B-13Jul2012-v4/AOD"
    fi

    # get available run range
    DATASET_NOSLASH=`echo $DATASET | sed 's/\///' | sed 's/\//_/g'`
    echo "doing... python das_client.py --limit=0 --query=\"run dataset=${DATASET}\""
    python das_client.py --limit=0 --query="run dataset=${DATASET}" | sort > runs_${DATASET_NOSLASH}.txt
    if [ ${FIRSTRUN} == 0 ]; then
        FIRSTRUN=`cat runs_${DATASET_NOSLASH}.txt | head -n 1`
    fi
    LASTRUN=`cat runs_${DATASET_NOSLASH}.txt | tail -n 1`

    # determine pset to use for CMSSW
    PSET=muontreemaker_Data2012_cfg.py
    if [[ "$DATASET" == *Single* ]]; then
        if [[ "$DATASET" == *Muon* ]]; then
            PSET=muontreemaker_Data2012_cfg.py
        elif [[ "$DATASET" == *Electron* ]]; then
            PSET=muontreemaker_Data2012_cfg.py
        fi
    elif [[ "$DATASET" == *Double* ]]; then
        if [[ "$DATASET" == *Muon* ]]; then
            PSET=muontreemaker_Data2012FR_cfg.py
        elif [[ "$DATASET" == *Electron* ]]; then
            PSET=electrontreemaker_Data2012FR_cfg.py
        fi
    fi

    # print crab config
    echo "
[CRAB]
jobtype   = cmssw
scheduler = remoteGlidein
use_server = 0
#scheduler = glidein
#use_server = 1

[CMSSW]
datasetpath             = $DATASET
pset                    = configs/pset_${DATASET_NOSLASH}_${FIRSTRUN}_${LASTRUN}_cfg.py
total_number_of_lumis   = -1
lumis_per_job           = 200
runselection            = ${FIRSTRUN}-${LASTRUN}
output_file             = leptonTree.root

[USER]
return_data             = 1
copy_data               = 0
#ui_working_dir          = ${DATASET_NOSLASH}_${FIRSTRUN}_${LASTRUN}
ui_working_dir          = RemoteGlidein_${DATASET_NOSLASH}_${FIRSTRUN}_${LASTRUN}

[GRID]
maxtarballsize = 50
se_black_list = T2_US_Caltech,T2_US_UCSD" \
    > crab_RemoteGlidein_${DATASET_NOSLASH}_${FIRSTRUN}_${LASTRUN}.cfg

    # change gtag in CMSSW config
    if [ ! -d "configs" ]; then
        mkdir configs
    fi
    sed s/GR_R_52_V7::All/${GTAG}/ ../test/${PSET} > configs/pset_${DATASET_NOSLASH}_${FIRSTRUN}_${LASTRUN}_cfg.py

    # submit
    crab -create -submit -cfg crab_RemoteGlidein_${DATASET_NOSLASH}_${FIRSTRUN}_${LASTRUN}.cfg

done

