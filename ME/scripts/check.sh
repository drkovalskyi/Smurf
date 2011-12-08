#!/bin/bash

# input is process name to check
# This the script to compare the number of jobs submitted and the output to make sure no job is missing

PROC=$1

if [ ! $# -eq 1 ]; then
    echo "USAGE: ./check.sh proc
        proc - the process name e.g. hww120"
    exit 1
fi

if [ -f  commands_${PROC}_resubmit.cmd ]; then
    rm -f  commands_${PROC}_resubmit.cmd
fi

while read LINE; do 
    TEST=`echo $LINE | grep "scramv1"`
    if [ ! -z "$TEST" ]; then
	export ENV=$TEST
    fi
done < commands_${PROC}.cmd

echo writting commands_${PROC}_resubmit.cmd
echo "
# -*- sh -*- # for font lock mode
# variable definitions
$ENV
- tag =
- output = outputFile=
- tagmode = none
- tarfile = ME_tarball.tgz
- untardir = tardir
- copycommand = cp

# Sections listed
# arguments of cafRun.sh: section number, process, input directory, number of events to run" > commands_${PROC}_resubmit.cmd


export RESUBMIT=0
while read LINE; do 
	TEST=`echo $LINE | grep "tardir/cafRun"`
	if [ ! -z "$TEST" ]; then
        	SECTION=`echo $LINE | awk '{print $3}'`
        	FILE=`ls | grep ^${PROC}_ME_${SECTION}.root`
		echo "$FILE"
        	if [ -z "$FILE" ]; then
		    echo $LINE >>  commands_${PROC}_resubmit.cmd
		    echo Found missing files in ${PROC} in Job ${SECTION}..
		    RESUBMIT=1
	        fi
	fi

done < commands_${PROC}.cmd

if [ ${RESUBMIT} == "1" ]; then
    echo python runManySections.py --submitCondor commands_${PROC}_resubmit.cmd
#    python runManySections.py --submitCondor commands_${PROC}_resubmit.cmd
fi

