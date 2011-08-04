#! /usr/bin/env python

import sys, os, math, argparse
from computeWeights import computeWeights, computeFakeWeights

datasets = ['qqww.root',
            'ggww.root',
            'wgamma.root',
            'ttbar.root',
            'tw.root',
            'zz.root',
            'wz.root',
            'dymm.root',
            'dyee.root',
            'dytt.root',
            ]



parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),description="analysis")
parser.add_argument('-s',dest='sample',action='store',required=False,help='sample')
parser.add_argument('-m',dest='mass',action='store',required=False,help='higgs mass')
parser.add_argument('-l',dest='luminosity',action='store',required=False,help='equivalent luminosity for MC')
parser.add_argument('-f',dest='computeFakes',action='store',required=False,help='compute fakes?')
args=parser.parse_args()

# signal mass
mass = 140
if args.mass: mass = int(args.mass)

# computeFakes
computeFakes = False
if args.computeFakes: computeFakes = bool(args.computeFakes)

# luminosity
luminosity = 1.092#*1.06
if args.luminosity: luminosity = float(args.luminosity)


# single sample
if args.sample:
    if computeFakes == True:
        print 'Computing Fakes'
        if 'data' in args.sample:
            print 'Processing as DATA'
            computeFakeWeights(args.sample,True)
        else:
            print 'Processing as MC'
            computeFakeWeights(args.sample,False,luminosity)
    else:
        print 'Computing ScaleFactor and Trigger Weights'
        computeWeights(args.sample, mass)

# loop over all samples
else:
    for sample in datasets:
        print 'processing: ', sample
        if computeFakes == True:
            print 'Computing Fakes'
            if 'data' in sample:
                print 'Processing as DATA'
                computeFakeWeights(sample,True)
            else:
                print 'Processing as MC'
                computeFakeWeights(sample,False,luminosity)
        else:
            print 'Computing ScaleFactor and Trigger Weights'
            computeWeights(sample, mass)



