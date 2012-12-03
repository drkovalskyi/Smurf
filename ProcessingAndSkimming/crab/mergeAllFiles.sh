#!/bin/bash

OUTDIR=/smurf/dlevans/LeptonTree
TAG=V00-02-07

if [ ! -d "${OUTDIR}/${TAG}" ]; then
    mkdir ${OUTDIR}/${TAG}
fi


for DIR in `ls | grep ^RemoteGlidein`;
do
    DATASET=`echo ${DIR} | sed 's/RemoteGlidein_//'`
    if [ ! -d "${OUTDIR}/${TAG}/${DATASET}" ]; then
        mkdir ${OUTDIR}/${TAG}/${DATASET}
    fi

    echo "doing $DIR"
    ./mergeFiles.sh $DIR/res
    cp $DIR/res/merged.root ${OUTDIR}/${TAG}/${DATASET}/merged.root

done

