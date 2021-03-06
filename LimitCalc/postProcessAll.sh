#!/bin/sh
# ./postprocessing.pl eps_final2/limits_0j_cut.txt "0-jets (cut)"
# ./postprocessing.pl eps_final2/limits_1j_cut.txt "1-jet (cut)"
# ./postprocessing.pl eps_final2/limits_2j_cut.txt "2-jets (cut)"
# ./postprocessing.pl eps_final2/limits_nj_cut.txt "0/1/2-jets (cut)"

./postprocessing.pl eps_final2/limits_0j_shape.txt "0-jets (shape)"
./postprocessing.pl eps_final2/limits_1j_shape.txt "1-jet (shape)"
./postprocessing.pl eps_final2/limits_nj_shape.txt "0/1/2-jets (shape)"

# ./postprocessing.pl eps_final2/limits_nj_shape_drop_1jll.txt "0/1-jets (without SF 1-jet)"

./postprocessing.pl eps_final2/limits_of_0j_shape.txt "0-jets OF"
./postprocessing.pl eps_final2/limits_of_1j_shape.txt "1-jet OF"
./postprocessing.pl eps_final2/limits_of_nj_shape.txt "0/1-jets OF"
./postprocessing.pl eps_final2/limits_sf_0j_shape.txt "0-jets SF"
./postprocessing.pl eps_final2/limits_sf_1j_shape.txt "1-jet SF"
./postprocessing.pl eps_final2/limits_sf_nj_shape.txt "0/1-jets SF"

# ./postprocessing.pl eps_final2/limits_of_0j_cut.txt "0-jets OF"
# ./postprocessing.pl eps_final2/limits_of_1j_cut.txt "1-jet OF"
# ./postprocessing.pl eps_final2/limits_of_nj_cut.txt "0/1-jets OF"
# ./postprocessing.pl eps_final2/limits_sf_0j_cut.txt "0-jets SF"
# ./postprocessing.pl eps_final2/limits_sf_1j_cut.txt "1-jet SF"
# ./postprocessing.pl eps_final2/limits_sf_nj_cut.txt "0/1-jets SF"
