#!/bin/bash
TASK=$1
for STDOUT in `ls ${TASK}/res/*stdout`; do
    FAIL=`grep "Executable failed  8028" ${STDOUT}`
    if [ "${FAIL}" != "" ]; then
        HOST=`grep "Job submitted on" ${STDOUT}`
        ERROR=`grep -A 1 "Exception Message" ${STDOUT} | tail -n 1`
        echo ${HOST} ---- ${ERROR}
    fi
done

