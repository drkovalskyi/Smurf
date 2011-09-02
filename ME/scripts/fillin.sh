#!/bin/bash

for LINE in `grep root bad.txt | awk '{print $9}'`; do

	mv $LINE bad/
	SECTION=`echo $LINE | sed 's/\(.*\)_ME_\(.*\).root/\2/'`
	PROC=`echo $LINE | sed 's/\(.*\)_ME_\(.*\).root/\1/'`

	HYP=0
	if [ ! -z `echo $PROC | grep hww` ]; then
                HYP=`echo $PROC | sed 's/hww\(.*\)/\1/g'`
	fi

	echo "output_\$(JID)         tardir/cafRun.sh $SECTION $PROC /uscms/home/dlevans/smurf/data/Run2011_Spring11_SmurfV5/mitf-zerojet/ 250 $HYP"
done

