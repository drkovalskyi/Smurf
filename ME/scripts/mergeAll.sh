#!/bin/bash

export OUTPUTDIR=/uscms_data/d2/ygao/Run2011_Summer11_SmurfV7_42X/mitf-alljets/WW/0j/ME/
for PROC in data data-emb-tau123 hww115 hww120 hww130 hww140 hww150 hww160 hww170 hww180 hww190 hww200 hww250 hww300 hww350 hww400 hww450 hww500 hww500 hww550 hww600 qqww ggww wjets wjets_data wjets_PassFail ttbar tw dyee dymm dytt wgamma wz zz_py ttbar_mg tw_ds wg3l dyee_LooseMET dymm_LooseMET ww_mcnlo ww_mcnlo_up ww_mcnlo_down; do
    echo "Merging " $PROC
    ./merge.sh ./ $PROC $OUTPUTDIR
done