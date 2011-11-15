#!/bin/bash

export INPUTSMURFFDIR=$1
export MEFDIR=$2
export LUMI=1545;

if [ ! $# -eq 2 ]; then
    echo "USAGE: ./addLR.sh inputSmurfFDir meFDir
       inputSmurfFDir - the directory containing the input smurf ntuples for dXsec calculation. E.g. /smurf/yygao/data/LP2011/WW/0j/
        meFDir        - the directory containing the dXsec variables. Eg. /smurf/yygao/data/LP2011/WW/0j/ME"
    exit 1
fi

rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data
data-emb-tau123
wgamma
wg3l
qqww
ggww
ww_mcnlo
ww_mcnlo_up
ww_mcnlo_down
wjets
wjets_data
wjets_PassFail
ttbar
ttbar_mg
tw
tw_ds
wz
zz_py
dyee
dymm
dytt
dyee_LooseMET
dymm_LooseMET
EOF

# ===========================================
# Fill the smurfntuples with LR
# ===========================================

#for MH in 115 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600; do 
#for MH in 115; do
for MH in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600; do 
    rm -f log/add_lr_$MH.log;
    echo doing $MH
    # first add the LR for the higgs sample
    root -q -l LR.C+\(${MH},\"hww${MH}\",\"${INPUTSMURFFDIR}\",\"${MEFDIR}\",-1,${LUMI}\) >> log/add_lr_$MH.log
    mv ${MEFDIR}/hww${MH}_ME.root ${MEFDIR}/hww${MH}.root 
    # add the LR for the backgrounds
    for i in `cat list_samples.txt` ; do
	sample=${i%%,*};
	echo filling ME LR in sample:  $sample
	root -q -l LR.C+\(${MH},\"$sample\",\"${INPUTSMURFFDIR}\",\"${MEFDIR}\",-1,${LUMI}\) >> log/add_lr_$MH.log 
	mv ${MEFDIR}/${sample}_ME.root ${MEFDIR}/${sample}.root 
    done
done

rm -f list_samples.txt;
