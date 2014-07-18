#!/bin/tcsh -f

foreach file (`cat samples_to_run_7TeV.txt | cut -d' ' -f1 `)

cp runAll-default_7TeV.sh runAll-default_7TeV_spin.sh;
chmod a+x runAll-default_7TeV_spin.sh;

hadd -f data2011/ntuples_126train_0jets_backgroundA_skim6_$file.root \
	data2011/ntuples_126train_0jets_backgroundA_skim6.root \
	data2011/ntuples_126train_0jets_$file.root;

rm -f data2011/ntuples_126train_1jets_backgroundA_skim6_$file.root;
cd data2011;
ln -s ntuples_126train_0jets_backgroundA_skim6_$file.root ntuples_126train_1jets_backgroundA_skim6_$file.root;
cd -;

sed -i 's/backgroundA_skim6.root/backgroundA_skim6_'"${file}"'.root/' runAll-default_7TeV_spin.sh;

./runAll-default_7TeV_spin.sh 0 126 6;
./runAll-default_7TeV_spin.sh 1 126 6;

#rm -rf /data/smurf/ceballos/inputLimits/ana_spin/${file};
mkdir -p /data/smurf/ceballos/inputLimits/ana_spin/${file};

mv output/histo_limits_ntuples_126train_0jets_0j_chan6_mh126_spin_7TeV.txt /data/smurf/ceballos/inputLimits/ana_spin/${file}/hwwof_0j_7TeV.txt;
mv output/histo_limits_ntuples_126train_1jets_1j_chan6_mh126_spin_7TeV.txt /data/smurf/ceballos/inputLimits/ana_spin/${file}/hwwof_1j_7TeV.txt;
mv hwwof_0j.input_7TeV.root /data/smurf/ceballos/inputLimits/ana_spin/${file}/hwwof_0j.input_7TeV.root;
mv hwwof_1j.input_7TeV.root /data/smurf/ceballos/inputLimits/ana_spin/${file}/hwwof_1j.input_7TeV.root;

rm -f data2011/ntuples_126train_0jets_backgroundA_skim6_$file.root
rm -f data2011/ntuples_126train_1jets_backgroundA_skim6_$file.root

end

rm -f runAll-default_7TeV_spin.sh;
