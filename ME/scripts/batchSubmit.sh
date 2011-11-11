#!/bin/bash

#
# usage
#

if [ ! $# -eq 6 ]; then
    echo "USAGE: ./batchSubmit.sh   ME_CODE_LOCATION PROCESS   ME_NTUPLE_LOCATION   NSECTIONS   NEVT_PER_SECTION
	ME_CODE_LOCATION   - location of ME code (e.g. /uscms/home/dlevans/CMSSW_3_11_3_ME_Smurf/src)
	PROCESS            - name of input smurf ntuple (e.g. ww)
	ME_NTUPLE_LOCATION - location of input smurf ntuple (e.g. /uscms/home/dlevans/smurf/data/Run2011_Spring11_SmurfV3/tas-zerojet/)
	NSECTIONS          - number of jobs to submit
	NEVT_PER_SECTION   - number of events per job
	Mode               - WW or ZZ"
    exit 1
fi

ME_CODE_LOCATION=$1
PROCESS=$2
ME_NTUPLE_LOCATION=$3
NSECTIONS=$4
NEVT_PER_SECTION=$5
MODE=$6

#
# make a tar of the ME code
# to run on the worker node
#

echo making tar file of ME code to send to worker nodes
if [ -f ME_tarball.tgz ]; then
	rm ME_tarball.tgz
fi
tar --exclude "scripts" -czf ME_tarball.tgz ../*

#
# make the commands file
# that describes what to send to the worker node
# and what to do with it
#

echo writting commands_${PROCESS}.cmd
echo "
# -*- sh -*- # for font lock mode
# variable definitions
- env =  cd $ME_CODE_LOCATION; . /uscmst1/prod/sw/cms/bashrc prod; eval \`export SCRAM_ARCH=slc5_amd64_gcc434; scramv1 runtime -sh\`; cd -
- tag =
- output = outputFile=
- tagmode = none
- tarfile = ME_tarball.tgz
- untardir = tardir
- copycommand = cp

# Sections listed
# arguments of cafRun.sh: section number, process, input directory, number of events to run mode" > commands_${PROCESS}.cmd

for (( SECTION=1; SECTION<=$NSECTIONS; SECTION++)) 
do
  echo "output_\$(JID)         tardir/cafRun.sh $SECTION $PROCESS $ME_NTUPLE_LOCATION $NEVT_PER_SECTION $MODE" >> commands_${PROCESS}.cmd
done

#
# now submit the jobs
# described in the output commands.cmd file
# from above
#

echo python runManySections.py --submitCondor commands_${PROCESS}.cmd
python runManySections.py --submitCondor commands_${PROCESS}.cmd

echo HINT: monitor your jobs with condor_q -submitter username

