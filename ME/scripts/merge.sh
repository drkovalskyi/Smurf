#!/bin/bash

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./merge.sh   FILEPATH PROCES
	FILEPATH - path to ME files to merge
	PROCESS - name of process, e.g. ww, hww120 etc.
	OUTPUT - directory to put merged files"
    exit 1
fi

FILEPATH=$1
PROCESS=$2
OUTPUT=$3

. /uscmst1/prod/sw/cms/bashrc prod
eval `export SCRAM_ARCH=slc5_ia32_gcc434; scramv1 runtime -sh`

if [ -f merge.C ]; then
	rm -f merge.C
fi

echo -e "{\tTChain s(\"tree\");" > merge.C
for FILE in $FILEPATH/${PROCESS}_ME_*.root; do

        if [ .$FILE = ."" ]; then
                echo "ERROR : File ${PROCESS}_ME_$i.root not found, skip"  
        else
                echo -e "\ts.Add(\"$FILE\");" >> merge.C
        fi
done
echo -e "\ts.SetMaxTreeSize(1e9);" >> merge.C
echo -e "\ts.Merge(\"${OUTPUT}/${PROCESS}_ME.root\");" >> merge.C
echo "}" >> merge.C

echo "Merging $PROCESS"
root -l -q merge.C

