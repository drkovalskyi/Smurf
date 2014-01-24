#!/bin/tcsh

foreach n(46 47 48 49 50 51 52 53 53 54 55 100)
./runAllJobsOneJetBin.sh $1 $n;
end
