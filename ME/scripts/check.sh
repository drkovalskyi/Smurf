#!/bin/bash

# input is process name to check
# This the script to compare the number of jobs submitted and the output to make sure no job is missing

PROC=$1

if [ ! $# -eq 1 ]; then
    echo "USAGE: ./check.sh proc
        proc - the process name e.g. hww120"
    exit 1
fi


while read LINE; do 
	TEST=`echo $LINE | grep "tardir/cafRun"`
	if [ ! -z "$TEST" ]; then
        	SECTION=`echo $LINE | awk '{print $3}'`
        	FILE=`ls | grep ${PROC}_ME_${SECTION}.root`
        	if [ -z "$FILE" ]; then
        	        echo ${PROC}_ME_${SECTION}.root is missing
	        fi
	fi

done < commands_${PROC}.cmd

