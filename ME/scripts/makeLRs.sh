#!/bin/sh

mass=$1

root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"zz\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"wz\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"qqww\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"ggww\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"ttbar\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"tw\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"dyee\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"dymm\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"wgamma\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"data\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q LR_HZZ.C++\(${mass},\"ntuple\",\"hzz${mass}\",\"nt_0J/\",\"ME_${mass}/\",-1\)
root -l -b -q makePlots.C\(${mass},\"ME_${mass}/\",1092,109\)
