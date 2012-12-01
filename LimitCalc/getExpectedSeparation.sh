#!/bin/bash

# ------------------------------------------------
# a script to generate toy MC according to card M0
# and perform fit based on M0 and M1
# and generate toy MC according to card M1
# and perform fit based on M0 and M1
# ------------------------------------------------
# You need to provide the actual card names 
# By default the M0 is for the SM Higgs
#                M1 is for spin 2
# ------------------------------------------------

if [ ! $# -eq 3 ]; then
    echo "
USAGE: ./getExpectedSeparation.sh TASK NJOBS NTOYS
    TASK  - Unique name for this task
    NJOBS - The number of jobs
    NTOYS - The number of pseudoexperiments per job"
    exit 1
fi

TASK=$1

if [ -d ${TASK} ]; then
	echo "Task ${TASK} already exists"
	exit 1
else
	mkdir -p ${TASK}/log
	mkdir ${TASK}/output
fi

NJOBS=$2
NTOYS=$3

#
# queue configuration
#

QUEUE="1nd"

#
# job configuration
#

WORKDIR=`pwd`
OUTPUTHOSTNAME="lxplus.cern.ch"
LANDS="../../../../LandS/test/lands.exe"
WRAPPER="process_job.sh"

#
# Define the cards and names of the two hypothesis
# make sure use the SMHiggs for the M0
# this is essential in calculating the hypothesis separations 
#

M0NAME=SMHiggs
M0="'$WORKDIR/hwwof_0j_shape_8TeV.txt $WORKDIR/hwwof_1j_shape_8TeV.txt'"
#M0="'$WORKDIR/hwwof_0j_shape_8TeV.txt'"

M1NAME=Graviton
M1="'$WORKDIR/xwwof_0j_shape_8TeV.txt $WORKDIR/xwwof_1j_shape_8TeV.txt'"
#M1="'$WORKDIR/xwwof_0j_shape_8TeV.txt'"

#
# make the wrapper
#

cat > ${WRAPPER} << EOF
#!/bin/sh

#
# script to generate toys from M0
# and fit toys for M0 and M1

# set up parameters
REMOTEDIR=\$1
SEED=\$2
NTOYS=\$3
M0=\$4
M1=\$5
M0NAME=\$6
M1NAME=\$7
echo LINE 84 \${M0NAME}
echo LINE 85 \${M1NAME}
WORKDIR=${WORKDIR}
REMOTEHOST=${OUTPUTHOSTNAME}
LANDS=\${WORKDIR}/${LANDS}
TMPDIR=\`pwd\`
SCRATCH=\`mktemp -d\`
LOGFILE=log_\${M0NAME}_\${M1NAME}_\${SEED}
echo LINE 92 \$LOGFILE
# get the environment
cd \$WORKDIR
eval \`scram runtime -sh\`
cd \$SCRATCH

# generate the toys
LIB="-L \$CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so"
TOYCOMMAND="-M Hybrid -m 125 --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF --bWriteToys 1 -rMin 0 -rMax 10  --freq -singlePoint 1"
\${LANDS} -d \${M0} \${TOYCOMMAND} \${LIB} -n \${M0NAME} --nToysForCLsb \${NTOYS} --nToysForCLb 1 --seed \${SEED} > \${LOGFILE} 2>&1

# fit the toys
FITCOMMAND="-M MaxLikelihoodFit -rMin 0 -rMax 10 -m 125 --NoErrorEstimate --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF --doExpectation 1"
\${LANDS} -d \${M0} \${FITCOMMAND} \${LIB} --loadToysFromFile \${M0NAME}_PseudoData_sb_seed\${SEED}.root -n LL_toy\${M0NAME}_fit\${M0NAME}_seed\${SEED} >> \${LOGFILE} 2>&1
\${LANDS} -d \${M1} \${FITCOMMAND} \${LIB} --loadToysFromFile \${M0NAME}_PseudoData_sb_seed\${SEED}.root -n LL_toy\${M0NAME}_fit\${M1NAME}_seed\${SEED} >> \${LOGFILE} 2>&1

# return output and tidy up
gzip \${LOGFILE}
scp * \${REMOTEHOST}:\${REMOTEDIR} && rm *
EOF
chmod +x ${WRAPPER}

#
# submit to batch
#

echo "Will submit two sets of $NJOBS jobs,
  - one set for ${M0} signal + background toys fitted for ${M0} and ${M1},
  - one set for ${M1} signal + background toys fitted for ${M0} and ${M1}."

JOB=0
while [ $JOB -lt $NJOBS ]; do
	SEED=$RANDOM
        # generate toys from M0 and fit for M0 and M1
	bsub -q ${QUEUE} -o ${TASK}/log/joblog_${M0NAME}_${M1NAME}_${SEED}.lsf "${WRAPPER} $WORKDIR/${TASK}/output ${SEED} ${NTOYS} \"${M0}\" \"${M1}\" ${M0NAME} ${M1NAME}"
        # generate toys from M1 and fit for M0 and M1
	bsub -q ${QUEUE} -o ${TASK}/log/joblog_${M1NAME}_${M0NAME}_${SEED}.lsf "${WRAPPER} $WORKDIR/${TASK}/output ${SEED} ${NTOYS} \"${M1}\" \"${M0}\" ${M1NAME} ${M0NAME}"
	let JOB=$JOB+1
done

