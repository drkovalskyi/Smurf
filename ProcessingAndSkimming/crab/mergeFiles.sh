#!/bin/bash

DIR=$1

echo "
void merge() {
  TChain *chain = new TChain(\"leptons\");
  chain->SetMaxTreeSize(5000000000LL); //default is 100000000000LL
" > merge.C

for FILE in `ls $DIR | grep leptonTree`; do
    echo $FILE
    echo "chain->Add(\"$DIR/${FILE}\");" >> merge.C
done

echo "
chain->Merge(\"$DIR/merged.root\", \"fast\");
}" >> merge.C

#root -b -q merge.C

