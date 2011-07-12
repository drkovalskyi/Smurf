#!/bin/sh
# ./GetExpectLimits.pl eps_final/limits_0j_cut.txt
# ./GetExpectLimits.pl eps_final/limits_1j_cut.txt
# ./GetExpectLimits.pl eps_final/limits_2j_cut.txt
# ./GetExpectLimits.pl eps_final/limits_nj_cut.txt
# echo "wait five minutes before sending a new set of jobs"
# sleep 300
# ./GetExpectLimits.pl eps_final/limits_0j_shape.txt
# ./GetExpectLimits.pl eps_final/limits_1j_shape.txt
# ./GetExpectLimits.pl eps_final/limits_nj_shape.txt
# echo "wait five minutes before sending a new set of jobs"
# sleep 300
# ./GetExpectLimits.pl eps_final/limits_of_0j_shape.txt
# ./GetExpectLimits.pl eps_final/limits_of_1j_shape.txt
# ./GetExpectLimits.pl eps_final/limits_of_nj_shape.txt
# echo "wait five minutes before sending a new set of jobs"
# sleep 300
# ./GetExpectLimits.pl eps_final/limits_sf_0j_shape.txt
# ./GetExpectLimits.pl eps_final/limits_sf_1j_shape.txt
# ./GetExpectLimits.pl eps_final/limits_sf_nj_shape.txt
# echo "all jobs are submitted"

./GetExpectLimits.pl eps_final/limits_of_0j_cut.txt
./GetExpectLimits.pl eps_final/limits_of_1j_cut.txt
./GetExpectLimits.pl eps_final/limits_of_nj_cut.txt
echo "wait five minutes before sending a new set of jobs"
sleep 300
./GetExpectLimits.pl eps_final/limits_sf_0j_cut.txt
./GetExpectLimits.pl eps_final/limits_sf_1j_cut.txt
./GetExpectLimits.pl eps_final/limits_sf_nj_cut.txt
# echo "all jobs are submitted"