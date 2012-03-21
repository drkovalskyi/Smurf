#!/bin/sh 

# arguments
# 1 - grid node 
# 2 - filename
# 3 - input path
# 4 - number of events to run
# 5 - higgs mass argument

SECTION=$1
ROOTFILE=$2.root
INPUTPATH=$3
NEV=$4
EVSTART=$((NEV*($1-1)))
MODE=$5

cd tardir

echo "Event Start:" $EVSTART
echo "Copying $INPUTPATH/$ROOTFILE to the current folder"
cp $INPUTPATH/$ROOTFILE ./

echo "{" > temp.C

if [ "${MODE}" == 'WW' ]; then
    echo "gSystem->CompileMacro(\"runME_test.C\");" >> temp.C
    echo "runME_test(\"./\", \"$ROOTFILE\", \"./\", 10, 2, 100000, 1.0, $NEV, $EVSTART);" >> temp.C
fi

if [ "${MODE}" == 'ZZ' ]; then
    echo "gSystem->CompileMacro(\"runME_HZZ.C\");" >> temp.C
    echo "runME_HZZ(\"./\", \"$ROOTFILE\", \"./\", 10, 1, 100000, 1.0, $NEV, $EVSTART);" >> temp.C
fi
echo "}" >> temp.C
cat temp.C

echo `date`
root -l -b -q temp.C
mv $2_ME.root ../$2_ME_$1.root
find . -type f -name "*.root" ! -name '*ME*.root' -execdir rm {} +
echo `date`

