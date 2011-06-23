#!/bin/bash

# usage ./doAll.sh NJET LUMI DATALUMI

export NJETS=$1;
export LUMI=$2;
export DATALUMI=$3;

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./doAll.sh njets lumi
        njet - njet bin e.g. 0, 1 or 2
	projlumi  - the luminosity to which the limit corresponds to
	datalumi  - the luminosity of data "
    exit 1
fi


for MASS in 115 120 130 140 150 160 170 180 190 200 250 300; do
    echo doing $MASS
    nohup ./run_mva.sh  0 $MASS $LUMI $DATALUMI >& output/mva_$MASS.log 
done

# Make the text file from the lands outpouts and draw the limits

export LIMITDIR=output/limits/1000pb_smallwjets/

for MVA in BDT ME; do
   echo $MVA method
    if [ -f ${LIMITDIR}/results_${NJETS}j_${LUMI}pb_$MVA.txt ]; then rm -f ${LIMITDIR}/results_${NJETS}j_${LUMI}pb_$MVA.txt; fi
    touch ${LIMITDIR}/results_${NJETS}j_${LUMI}pb_$MVA.txt
    for LOG in `ls ${LIMITDIR} | grep $MVA.res`; do
        # get the data from this log file and print it to the results file
            OBSERVED=`cat ${LIMITDIR}/$LOG | grep "Observed Upper Limit" | awk '{print $12}'`
            BANDS=`cat ${LIMITDIR}/$LOG | grep "BANDS" | awk '{print  $2"  "$3"  "$4"  "$5" "$6" "$7}'`
            echo $OBSERVED $BANDS >> ${LIMITDIR}/results_${NJETS}j_${LUMI}pb_$MVA.txt
    done

   # now call the root script to draw the bands
#    root -l -q plotBand2.C\(${NJETS},\"${LIMITDIR}\",\"${LUMI}pb_${MVA}\"\); 
done
