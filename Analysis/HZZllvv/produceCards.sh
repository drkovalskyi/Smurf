#!/bin/sh

for i in 250 300 400
do
  root -l -b -q Smurf/Analysis/HZZllvv/MakeHZZRes.C+\(0,${i},\"\",\"/data/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/hzz${i}.root\",\"/data/smurf/data/EPS/mitf/background42x.root\",\"/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data_photons.root\",\"/data/smurf/data/EPS/mitf/data_2l.root\",0,true,0\)
  root -l -b -q Smurf/Analysis/HZZllvv/MakeHZZRes.C+\(1,${i},\"\",\"/data/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/hzz${i}.root\",\"/data/smurf/data/EPS/mitf/background42x.root\",\"/data/smurf/sixie/data/Run2011_Spring11_SmurfV6/mitf-alljets/data_photons.root\",\"/data/smurf/data/EPS/mitf/data_2l.root\",0,true,0\)

done
