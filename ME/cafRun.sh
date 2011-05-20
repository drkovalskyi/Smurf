#!/bin/sh 

# arguments
# 1 - grid node 
# 2 - filename
# 3 - input path
# 4 - number of events to run

SECTION=$1
ROOTFILE=$2.root
INPUTPATH=$3
NEV=$4
EVSTART=$((NEV*($1-1)))

cd tardir

echo "Event Start:" $EVSTART
echo "Copying $INPUTPATH/$ROOTFILE to the current folder"
cp $INPUTPATH/$ROOTFILE ./

echo "{" > temp.C
echo "gSystem->CompileMacro(\"runME_test.C\");" >> temp.C
echo "runME_test(\"./\", \"$ROOTFILE\", \"./\", 10, 1, 100000, 1.0, $NEV, $EVSTART);" >> temp.C
echo "}" >> temp.C
cat temp.C

root -l -b -q temp.C
mv $2_ME.root ../$2_ME_$1.root
find . -type f -name "*.root" ! -name '*ME*.root' -execdir rm {} +

