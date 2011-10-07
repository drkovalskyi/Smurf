#!/bin/bash

export INPUTSMURFFDIR=$1
export MEFDIR=$2
export LUMI=1143

if [ ! $# -eq 2 ]; then
    echo "USAGE: ./addLRHZZ.sh inputSmurfFDir meFDir
       inputSmurfFDir - the directory containing the input smurf ntuples for dXsec calculation. E.g. /smurf/yygao/data/EPSHZZV0/ZZ/allj/
        meFDir        - the directory containing the dXsec variables. Eg. /smurf/yygao/data/EPSHZZV0/ZZ/allj/ME/"
    exit 1
fi

rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data_2l.goodlumi1092ipb
ww2l
ggww
ttbar
wtop
ttop
wz
zz
zee
zmm
ztt
EOF

# ===========================================
# Fill the smurfntuples with LR
# ===========================================

for MH in 250 300 350 400 500 600; do
    rm -f log/add_lr_$MH.log;
    echo doing $MH
    # first add the LR for the higgs sample
    root -q -l LR_HZZ.C+\(${MH},\"gfhzz${MH}\",\"${INPUTSMURFFDIR}\",\"${MEFDIR}\",-1,${LUMI}\) >> log/add_lr_$MH.log
    mv ${MEFDIR}/gfhzz${MH}_ME.root ${MEFDIR}/gfhzz${MH}.root 
    # add the LR for the backgrounds
    for i in `cat list_samples.txt` ; do
	sample=${i%%,*};
	echo filling ME LR in sample:  $sample
	root -q -l LR_HZZ.C+\(${MH},\"$sample\",\"${INPUTSMURFFDIR}\",\"${MEFDIR}\",-1,${LUMI}\) >> log/add_lr_$MH.log 
	mv ${MEFDIR}/${sample}_ME.root ${MEFDIR}/${sample}.root 
    done
done

rm -f list_samples.txt;
