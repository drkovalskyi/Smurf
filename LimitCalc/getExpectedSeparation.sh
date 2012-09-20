#!/bin/bash

#
# a script to generate toy MC according to card M0
# and perform fit based on M0 and M1
# and generate toy MC according to card M1
# and perform fit based on M0 and M1
#

if [ ! $# -eq 5 ]; then
    echo "
USAGE: ./getExpectedSeparation.sh TASK M0 M1 NJOBS NTOYS
    TASK  - Unique name for this task
    M0    - The card for model0
    M1    - The card for model1
    NJOBS - The number of jobs
    NTOYS - The number of pseudoexperiments per job
"
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

M0=$2
M1=$3
NJOBS=$4
NTOYS=$5

#
# queue configuration
#

QUEUE="1nd"

#
# job configuration
#

WORKDIR=`pwd`
OUTPUTHOSTNAME="lxplus.cern.ch"
LANDS="../LandS/test/lands.exe"
WRAPPER="process_job.sh"

#
# make the wrapper
#

cat > ${WRAPPER} << EOF
#!/bin/sh

#
# script to generate toys from M0
# and fit toys for M0 and M1
#

# set up parameters
REMOTEDIR=\$1
SEED=\$2
NTOYS=\$3
M0=\$4
M1=\$5
M0NAME=\`echo \${M0} | sed 's/\.txt//'\`
M1NAME=\`echo \${M1} | sed 's/\.txt//'\`
WORKDIR=${WORKDIR}
REMOTEHOST=${OUTPUTHOSTNAME}
LANDS=\${WORKDIR}/${LANDS}
TMPDIR=\`pwd\`
SCRATCH=\`mktemp -d\`
LOGFILE=log_\${M0NAME}_\${M1NAME}_\${SEED}

# get the environment
cd \$WORKDIR
eval \`scram runtime -sh\`
cd \$SCRATCH

# generate the toys
LIB="-L \$CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so"
TOYCOMMAND="-M Hybrid -m 125 --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF --bWriteToys 1 -rMin 0 -rMax 10  --freq -singlePoint 1"
\${LANDS} -d \${WORKDIR}/\${M0} \${TOYCOMMAND} \${LIB} -n \${M0NAME} --nToysForCLsb \${NTOYS} --nToysForCLb 1 --seed \${SEED} > \${LOGFILE} 2>&1

# fit the toys
FITCOMMAND="-M MaxLikelihoodFit -rMin 0 -rMax 10 -m 125 --NoErrorEstimate --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF --doExpectation 1"
\${LANDS} -d \${WORKDIR}/\${M0} \${FITCOMMAND} \${LIB} --loadToysFromFile \${M0NAME}_PseudoData_sb_seed\${SEED}.root -n LL_toy\${M0NAME}_fit\${M0NAME}_seed\${SEED} >> \${LOGFILE} 2>&1
\${LANDS} -d \${WORKDIR}/\${M1} \${FITCOMMAND} \${LIB} --loadToysFromFile \${M0NAME}_PseudoData_sb_seed\${SEED}.root -n LL_toy\${M0NAME}_fit\${M1NAME}_seed\${SEED} >> \${LOGFILE} 2>&1

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
    bsub -q ${QUEUE} -o ${TASK}/log/joblog_${M0}_${M1}_${SEED}.lsf "${WRAPPER} $WORKDIR/${TASK}/output ${SEED} ${NTOYS} ${M0} ${M1}"
	# generate toys from M1 and fit for M0 and M1
    bsub -q ${QUEUE} -o ${TASK}/log/joblog_${M1}_${M0}_${SEED}.lsf "${WRAPPER} $WORKDIR/${TASK}/output ${SEED} ${NTOYS} ${M1} ${M0}"
	let JOB=$JOB+1
done

