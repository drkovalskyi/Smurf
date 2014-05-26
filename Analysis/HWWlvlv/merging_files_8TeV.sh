source ~/EVAL_SH65 5_3_14;
export myDir=$PWD;
export OUTPUT=/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_mva;
cd ${OUTPUT};

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_SMH.root      ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_hww125.root
hadd -f ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_SMH.root      ntuples2012_MultiClass_125train_1jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_1jets_hww125.root

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_0M.root       ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_xww0m125.root
rm -f   ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_0M.root
ln -s   ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_0M.root       ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_0M.root

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_1M.root       ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_xww1m125.root
rm -f   ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_1M.root
ln -s   ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_1M.root       ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_1M.root

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_1P.root       ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_xww1p125.root
rm -f   ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_1P.root
ln -s   ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_1P.root       ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_1P.root

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_1MIX.root     ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_xww1mix125.root
rm -f   ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_1MIX.root
ln -s   ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_1MIX.root     ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_1MIX.root

hadd -f ntuples2012_MultiClass_125train_0jets_xww2pqq0p5125.root              ntuples2012_MultiClass_125train_0jets_xww2pqq125.root ntuples2012_MultiClass_125train_0jets_xww2p125.root
rm -f   ntuples2012_MultiClass_125train_1jets_xww2pqq0p5125.root
ln -s   ntuples2012_MultiClass_125train_0jets_xww2pqq0p5125.root              ntuples2012_MultiClass_125train_1jets_xww2pqq0p5125.root

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_2P.root       ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_xww2p125.root
rm -f   ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_2P.root
ln -s   ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_2P.root       ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_2P.root

hadd -f ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_2PQQ.root     ntuples2012_MultiClass_125train_0jets_backgroundA_skim6.root ntuples2012_MultiClass_125train_0jets_xww2pqq125.root
rm -f   ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_2PQQ.root
ln -s   ntuples2012_MultiClass_125train_0jets_backgroundA_skim6_2PQQ.root     ntuples2012_MultiClass_125train_1jets_backgroundA_skim6_2PQQ.root

cd ~/releases/CMSSW_5_3_14/src/;
root -l -q -b Smurf/Analysis/Zll/WGammaStarScaleFactor.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww2pqq0p5125.root","${OUTPUT}/output12.root",1.035018)'; mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww2pqq0p5125.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2bp-v19.root","${OUTPUT}/output12.root")';   mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2bp-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2hm-v19.root","${OUTPUT}/output12.root")';   mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2hm-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2hp-v19.root","${OUTPUT}/output12.root")';   mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2hp-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph2-v19.root","${OUTPUT}/output12.root")';  mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph2-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph3-v19.root","${OUTPUT}/output12.root")';  mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph3-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph6-v19.root","${OUTPUT}/output12.root")';  mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph6-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph7-v19.root","${OUTPUT}/output12.root")';  mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph7-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2mh9-v19.root","${OUTPUT}/output12.root")';  mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2mh9-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2mh10-v19.root","${OUTPUT}/output12.root")'; mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2mh10-v19.root;

root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0p-vbf-v19.root","${OUTPUT}/output12.root",74,10001)';       mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0p-vbf-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0m-vbf-v19.root","${OUTPUT}/output12.root",74,10001)';       mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0m-vbf-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-vbf-v19.root","${OUTPUT}/output12.root",74,10001)';      mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-vbf-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-vbf-v19.root","${OUTPUT}/output12.root",74,10001)'; mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-vbf-v19.root;

root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0p-wh-v19.root","${OUTPUT}/output12.root",74,26)';       mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0p-wh-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0m-wh-v19.root","${OUTPUT}/output12.root",74,26)';       mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0m-wh-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0ph-wh-v19.root","${OUTPUT}/output12.root",74,26)';      mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0ph-wh-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-wh-v19.root","${OUTPUT}/output12.root",74,26)'; mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-wh-v19.root;

root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0p-zh-v19.root","${OUTPUT}/output12.root",74,24)';       mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0p-zh-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0m-zh-v19.root","${OUTPUT}/output12.root",74,24)';       mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0m-zh-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-zh-v19.root","${OUTPUT}/output12.root",74,24)';      mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-zh-v19.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-zh-v19.root","${OUTPUT}/output12.root",74,24)'; mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-zh-v19.root;

cd ${OUTPUT};
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2bp.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2bp-v19.root    ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2bp-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2hm.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2hm-v19.root    ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2hm-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2hp.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2hp-v19.root    ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2hp-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2ph2.root  ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2ph2-v19.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph2-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2ph3.root  ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2ph3-v19.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph3-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2ph6.root  ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2ph6-v19.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph6-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2ph7.root  ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2ph7-v19.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2ph7-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2mh9.root  ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2mh9-v19.root   ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2mh9-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_2mh10.root ntuples2012_MultiClass_125train_0jets_xww125p6_s12-g125ww4l-2mh10-v19.root  ntuples2012_MultiClass_125train_0jets_xww125p6_s12-gqq125ww4l-2mh10-v19.root

hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_s12-0p.root       ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0pm-v19.root      ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0p-vbf-v19.root       ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0p-wh-v19.root       ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0p-zh-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_s12-0m.root       ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mt-v19.root      ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0m-vbf-v19.root       ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0m-wh-v19.root       ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0m-zh-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_s12-0ph.root      ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-v19.root      ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-vbf-v19.root      ntuples2012_MultiClass_125train_0jets_xww125p6_s12-h125ww4l-0ph-wh-v19.root      ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0ph-zh-v19.root
hadd -f ntuples2012_MultiClass_125train_0jets_xww125p6_s12-0mf05ph0.root ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-v19.root ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-vbf-v19.root ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-wh-v19.root ntuples2012_MultiClass_125train_0jets_xww125p6_s12-x125ww4l-0mf05ph0-zh-v19.root

rm -f ntuples2012_MultiClass_125train_1jets_xww125p6_s12-0p.root;
ln -s ntuples2012_MultiClass_125train_0jets_xww125p6_s12-0p.root ntuples2012_MultiClass_125train_1jets_xww125p6_s12-0p.root;

cd $myDir;
