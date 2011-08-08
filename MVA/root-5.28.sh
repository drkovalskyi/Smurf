#!/bin/sh                                                                                                                                                                         
export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00/slc4_ia32_gcc34/root/

if [ x$LD_LIBRARY_PATH != x ]; then
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
else
    LD_LIBRARY_PATH=${ROOTSYS}/lib
fi
export LD_LIBRARY_PATH

if [ x$PATH != x ]; then
    PATH=${PATH}:${ROOTSYS}/bin
else
    PATH=${ROOTSYS}/bin
fi
export PATH

exec ${ROOTSYS}/bin/root $*


