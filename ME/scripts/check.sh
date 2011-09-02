#!/bin/bash

# input is process name to check

PROC=$1

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

