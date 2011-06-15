Here are basic instruction how to get LandS run in parallel.

0) Install LandS

cvs co -d LandS UserCode/mschen/LandS
cvs co -d Smurf/LimitCalc UserCode/Smurf/LimitCalc
cd Smurf/LimitCalc
make lands


1) create input text file with cards to be used. Example:

# Format:
# mass <card1> [<card2> ...]

115 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH115.txt
120 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH120.txt
130 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH130.txt
140 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH140.txt
150 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH150.txt
160 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH160.txt
170 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH170.txt
180 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH180.txt
190 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH190.txt
200 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH200.txt
250 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH250.txt
300 inputLimits/inputs_1j_1000pb_shape/hww-SM-mH300.txt

2) Modify GetExpectLimits.pl if you need to. Basically you can tell it
   if you want to start multiple jobs on a local machine or submit to
   LSF. You may also want to change number of jobs per mass point as
   well as adjust paramets of LandS jobs

3) Collect outputs in output/. Just let it finish. Nothing fancy.

4) do post processing to get a final plot:

   postprocess.pl <card_name> "n jets"

./postprocessing.pl hww-1000pb_0_shape.txt "0 jet"

What you get is hww-1000pb_0_shape.pdf - expected upper limits bands.


