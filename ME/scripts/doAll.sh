#!/bin/bash

# first run the mva

for MASS in 115 120 130 140 150 160 170 180 190 200 250 300; do
    echo doing $MASS
    nohup ./run_mva.sh  0 $MASS >& output/mva_$MASS.log
done

# now we have everything for all mass points
# and we will make the text file from which to draw the bands

export LIMITDIR=output/limits/
export LUMI=1000
export NJET=0

for MVA in BDT ME; do
   echo $MVA method
    if [ -f ${LIMITDIR}/results_${NJET}j_${LUMI}pb_$MVA.txt ]; then rm -f ${LIMITDIR}/results_${NJET}j_${LUMI}pb_$MVA.txt; fi
    touch ${LIMITDIR}/results_${NJET}j_${LUMI}pb_$MVA.txt
    for LOG in `ls ${LIMITDIR} | grep $MVA.res`; do
        # get the data from this log file and print it to the results file
	    OBSERVED=`cat ${LIMITDIR}/$LOG | grep "Observed Upper Limit" | awk '{print $12}'`
	    BANDS=`cat ${LIMITDIR}/$LOG | grep "BANDS" | awk '{print  $2"  "$3"  "$4"  "$5" "$6}'`
	    echo $OBSERVED $BANDS >> ${LIMITDIR}/results_${NJET}j_${LUMI}pb_$MVA.txt
    done

   # now call the root script to draw the bands
   root -l -q plotBand2.C\(${NJET},\"${LIMITDIR}\",\"${LUMI}pb_${MVA}\"\); 
done
