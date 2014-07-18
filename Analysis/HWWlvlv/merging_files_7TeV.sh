source ~/EVAL_SH65 5_3_14;
export myDir=$PWD;
export OUTPUT=/data/smurf/data/Run2011_Fall11_SmurfV9_42X/mitf-alljets_mva;
cd ${OUTPUT};

rm -f ntuples_126train_1jets_backgroundA_skim6.root
ln -s ntuples_126train_0jets_backgroundA_skim6.root            ntuples_126train_1jets_backgroundA_skim6.root

rm -f ntuples_126train_1jets_hww110.root \
      ntuples_126train_1jets_hww115.root \
      ntuples_126train_1jets_hww118.root \
      ntuples_126train_1jets_hww120.root \
      ntuples_126train_1jets_hww122.root \
      ntuples_126train_1jets_hww124.root \
      ntuples_126train_1jets_hww126.root \
      ntuples_126train_1jets_hww128.root \
      ntuples_126train_1jets_hww130.root \
      ntuples_126train_1jets_hww135.root \
      ntuples_126train_1jets_hww140.root \
      ntuples_126train_1jets_hww150.root \
      ntuples_126train_1jets_hww160.root \
      ntuples_126train_1jets_hww170.root \
      ntuples_126train_1jets_hww180.root \
      ntuples_126train_1jets_hww190.root \
      ntuples_126train_1jets_hww200.root \
      ntuples_126train_1jets_hww250.root \
      ntuples_126train_1jets_hww300.root \
      ntuples_126train_1jets_hww350.root \
      ntuples_126train_1jets_hww400.root \
      ntuples_126train_1jets_hww450.root \
      ntuples_126train_1jets_hww500.root \
      ntuples_126train_1jets_hww550.root \
      ntuples_126train_1jets_hww600.root

ln -s ntuples_126train_0jets_hww110.root ntuples_126train_1jets_hww110.root
ln -s ntuples_126train_0jets_hww115.root ntuples_126train_1jets_hww115.root
ln -s ntuples_126train_0jets_hww118.root ntuples_126train_1jets_hww118.root
ln -s ntuples_126train_0jets_hww120.root ntuples_126train_1jets_hww120.root
ln -s ntuples_126train_0jets_hww122.root ntuples_126train_1jets_hww122.root
ln -s ntuples_126train_0jets_hww124.root ntuples_126train_1jets_hww124.root
ln -s ntuples_126train_0jets_hww126.root ntuples_126train_1jets_hww126.root
ln -s ntuples_126train_0jets_hww128.root ntuples_126train_1jets_hww128.root
ln -s ntuples_126train_0jets_hww130.root ntuples_126train_1jets_hww130.root
ln -s ntuples_126train_0jets_hww135.root ntuples_126train_1jets_hww135.root
ln -s ntuples_126train_0jets_hww140.root ntuples_126train_1jets_hww140.root
ln -s ntuples_126train_0jets_hww150.root ntuples_126train_1jets_hww150.root
ln -s ntuples_126train_0jets_hww160.root ntuples_126train_1jets_hww160.root
ln -s ntuples_126train_0jets_hww170.root ntuples_126train_1jets_hww170.root
ln -s ntuples_126train_0jets_hww180.root ntuples_126train_1jets_hww180.root
ln -s ntuples_126train_0jets_hww190.root ntuples_126train_1jets_hww190.root
ln -s ntuples_126train_0jets_hww200.root ntuples_126train_1jets_hww200.root
ln -s ntuples_126train_0jets_hww250.root ntuples_126train_1jets_hww250.root
ln -s ntuples_126train_0jets_hww300.root ntuples_126train_1jets_hww300.root
ln -s ntuples_126train_0jets_hww350.root ntuples_126train_1jets_hww350.root
ln -s ntuples_126train_0jets_hww400.root ntuples_126train_1jets_hww400.root
ln -s ntuples_126train_0jets_hww450.root ntuples_126train_1jets_hww450.root
ln -s ntuples_126train_0jets_hww500.root ntuples_126train_1jets_hww500.root
ln -s ntuples_126train_0jets_hww550.root ntuples_126train_1jets_hww550.root
ln -s ntuples_126train_0jets_hww600.root ntuples_126train_1jets_hww600.root

rm -f ntuples_126train_1jets_data_skim6.root \
      ntuples_126train_1jets_hww_syst_skim6.root \
      ntuples_126train_1jets_xww0p125.root \
      ntuples_126train_1jets_xww0m125.root \
      ntuples_126train_1jets_xww2pqq125.root \
      ntuples_126train_1jets_xww2p125.root

ln -s ntuples_126train_0jets_data_skim6.root     ntuples_126train_1jets_data_skim6.root    
ln -s ntuples_126train_0jets_hww_syst_skim6.root ntuples_126train_1jets_hww_syst_skim6.root
ln -s ntuples_126train_0jets_xww0p125.root       ntuples_126train_1jets_xww0p125.root	   
ln -s ntuples_126train_0jets_xww0m125.root   	 ntuples_126train_1jets_xww0m125.root
ln -s ntuples_126train_0jets_xww2pqq125.root 	 ntuples_126train_1jets_xww2pqq125.root 
ln -s ntuples_126train_0jets_xww2p125.root	 ntuples_126train_1jets_xww2p125.root

hadd -f ntuples_126train_0jets_backgroundA_skim6_SMH.root      ntuples_126train_0jets_backgroundA_skim6.root ntuples_126train_0jets_hww126.root
rm -f   ntuples_126train_1jets_backgroundA_skim6_SMH.root
lns -s  ntuples_126train_0jets_backgroundA_skim6_SMH.root      ntuples_126train_1jets_backgroundA_skim6_SMH.root

hadd -f ntuples_126train_0jets_backgroundA_skim6_0M.root       ntuples_126train_0jets_backgroundA_skim6.root ntuples_126train_0jets_xww0m125.root
rm -f   ntuples_126train_1jets_backgroundA_skim6_0M.root
ln -s   ntuples_126train_0jets_backgroundA_skim6_0M.root       ntuples_126train_1jets_backgroundA_skim6_0M.root

hadd -f ntuples_126train_0jets_backgroundA_skim6_2P.root       ntuples_126train_0jets_backgroundA_skim6.root ntuples_126train_0jets_xww2p125.root
rm -f   ntuples_126train_1jets_backgroundA_skim6_2P.root
ln -s   ntuples_126train_0jets_backgroundA_skim6_2P.root       ntuples_126train_1jets_backgroundA_skim6_2P.root

hadd -f ntuples_126train_0jets_backgroundA_skim6_2PQQ.root     ntuples_126train_0jets_backgroundA_skim6.root ntuples_126train_0jets_xww2pqq125.root
rm -f   ntuples_126train_1jets_backgroundA_skim6_2PQQ.root
ln -s   ntuples_126train_0jets_backgroundA_skim6_2PQQ.root     ntuples_126train_1jets_backgroundA_skim6_2PQQ.root

hadd -f ntuples_126train_0jets_backgroundA_skim6_2PQQ0p5.root  ntuples_126train_0jets_backgroundA_skim6.root ntuples_126train_0jets_xww2pqq125.root ntuples_126train_0jets_xww2p125.root
rm -f   ntuples_126train_1jets_backgroundA_skim6_2PQQ0p5.root
ln -s   ntuples_126train_0jets_backgroundA_skim6_2PQQ0p5.root  ntuples_126train_1jets_backgroundA_skim6_2PQQ0p5.root

hadd -f ntuples_126train_0jets_xww2pqq0p5125.root              ntuples_126train_0jets_xww2pqq125.root ntuples_126train_0jets_xww2p125.root
rm -f   ntuples_126train_1jets_xww2pqq0p5125.root
ln -s   ntuples_126train_0jets_xww2pqq0p5125.root              ntuples_126train_1jets_xww2pqq0p5125.root

cd ~/releases/CMSSW_5_3_14/src/;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples_126train_0jets_xww125p6_gqq125ww4l-2hp.root","${OUTPUT}/output12.root")';    mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples_126train_0jets_xww125p6_gqq125ww4l-2hp.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples_126train_0jets_xww125p6_gqq125ww4l-2mh9.root","${OUTPUT}/output12.root")';   mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples_126train_0jets_xww125p6_gqq125ww4l-2mh9.root;
root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples_126train_0jets_xww125p6_gqq125ww4l-2ph7.root","${OUTPUT}/output12.root")';   mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples_126train_0jets_xww125p6_gqq125ww4l-2ph7.root;

root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples_126train_0jets_xww125p6_x125ww4l-0mf05ph0-vbf.root","${OUTPUT}/output12.root",74,10001)';    mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples_126train_0jets_xww125p6_x125ww4l-0mf05ph0-vbf.root;

root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples_126train_0jets_xww125p6_h125ww4l-0m-wh.root","${OUTPUT}/output12.root",74,26)';              mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples_126train_0jets_xww125p6_h125ww4l-0m-wh.root;

root -l -q -b Smurf/Analysis/Zll/changeDataType.C+'("${OUTPUT}/ntuples_126train_0jets_xww125p6_h125ww4l-0m-zh.root","${OUTPUT}/output12.root",74,24)';              mv ${OUTPUT}/output12.root ${OUTPUT}/ntuples_126train_0jets_xww125p6_h125ww4l-0m-zh.root;

cd ${OUTPUT};
hadd -f ntuples_126train_0jets_xww125p6_2bp.root                                                              ntuples_126train_0jets_xww125p6_g125ww4l-2bp.root;
hadd -f ntuples_126train_0jets_xww125p6_2mh10.root                                                            ntuples_126train_0jets_xww125p6_g125ww4l-2mh10.root;
hadd -f ntuples_126train_0jets_xww125p6_2ph3.root                                                             ntuples_126train_0jets_xww125p6_g125ww4l-2ph3.root;
hadd -f ntuples_126train_0jets_xww125p6_2ph6.root                                                             ntuples_126train_0jets_xww125p6_g125ww4l-2ph6.root;
hadd -f ntuples_126train_0jets_xww125p6_2hp.root      ntuples_126train_0jets_xww125p6_gqq125ww4l-2hp.root;
hadd -f ntuples_126train_0jets_xww125p6_2mh9.root     ntuples_126train_0jets_xww125p6_gqq125ww4l-2mh9.root    ntuples_126train_0jets_xww125p6_g125ww4l-2mh9.root;
hadd -f ntuples_126train_0jets_xww125p6_2ph7.root     ntuples_126train_0jets_xww125p6_gqq125ww4l-2ph7.root;

hadd -f ntuples_126train_0jets_xww125p6-0p.root       ntuples_126train_0jets_xww125p6_x125ww4l-0pm.root;
hadd -f ntuples_126train_0jets_xww125p6-0m.root                                                               ntuples_126train_0jets_xww125p6_h125ww4l-0m-wh.root;
hadd -f ntuples_126train_0jets_xww125p6-0ph.root      ntuples_126train_0jets_xww125p6_x125ww4l-0ph.root;
hadd -f ntuples_126train_0jets_xww125p6-0mf05ph0.root ntuples_126train_0jets_xww125p6_x125ww4l-0mf05ph0.root  ntuples_126train_0jets_xww125p6_x125ww4l-0mf05ph0-vbf.root;

rm -f ntuples_126train_1jets_xww125p6-0p.root;
ln -s ntuples_126train_0jets_xww125p6-0p.root ntuples_126train_1jets_xww125p6-0p.root;

cd $myDir;
