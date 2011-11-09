#!/bin/bash

#
# script to skim smurf ntuples
# it outputs 3 different selections
# WW pre-selections
# ZZ pre-selections
# PassFail samples for the WW analysis
#

INPUTDIR=$1
OUTPUTDIR=$2
SELECTION=$3

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./skim.sh   INPUTDIR OUTPUTDIR
        INPUTDIR - location of smurf ntuples to skim (e.g. /smurf/data/Run2011_Spring11_SmurfV3/mitf-alljets/)
        OUTPUTDIR - location to output skimmed ntuples
        SELECTION - selection, choose from ZZ, WW, PassFail, LooseMET, Wgamma"
    exit 1
fi


# check that directories exist

if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir doesnt exist!
        exit 2
fi

if [ ! -d $OUTPUTDIR ]; then
    echo Error: Output dir doesnt exist!
    exit 3
fi

# loop over root files in input dir
# and do the skim root script

if [ "$SELECTION" == 'ZZ' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
gfhzz250.root
gfhzz300.root
gfhzz350.root
gfhzz400.root
gfhzz500.root
gfhzz600.root
zz.root
wz.root
ttbar.root
stop.root
ttop.root
wtop.root
ww2l.root
ggww.root
EOF
fi

zee.root
zmm.root
ztt.root
data_2l.goodlumi1092ipb.root






if [ "$SELECTION" == 'WW' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data.root
zz.root
wz.root
ttbar.root
tw.root
qqww.root
ggww.root
wjets.root
wgamma.root    
dyee.root
dymm.root
dytt.root
hww115.root
hww120.root
hww130.root
hww140.root
hww150.root
hww160.root
hww170.root
hww180.root
hww190.root
hww200.root
hww250.root
hww300.root
hww350.root
hww400.root
hww450.root
hww500.root
hww550.root
hww600.root
EOF
fi

if [ "$SELECTION" == 'PassFail' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data.root
zz.root
wz.root
ttbar.root
tw.root
qqww.root
ggww.root
wgamma.root
EOF
fi

if [ "$SELECTION" == 'LooseMET' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
dyee.root
dymm.root
EOF
fi


# Do the skimming...
for FILE in `cat list_samples.txt` ; do
    for JETBIN in 0 1 2 ; do 
	outputdir=$OUTPUTDIR/$SELECTION/${JETBIN}j/
	if [ "$SELECTION" == 'PassFail' ]; then
	    outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	fi
	if [ "$SELECTION" == 'LooseMET' ]; then
	    outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	fi
	
	mkdir -p $outputdir
	echo doing "root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",$JETBIN\);"
	root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",$JETBIN\);
    done
    # for ZZ do a special alljet bin category
    if [ "$SELECTION" == 'ZZ' ]; then
	outputdir=$OUTPUTDIR/$SELECTION/allj/
	mkdir -p $outputdir
	echo doing "root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",3\);"
	root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",3\);
    fi
done

# if the selection is the PassFail merge all the files
if [ "$SELECTION" == 'PassFail' ]; then
    for JETBIN in 0 1 2 ; do 
	outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	rm -f ${outputdir}/wjets_PassFail.root
	rm -f ${outputdir}/wgamma_lgamma_PassFail.root
	rm  merge.C
	touch merge.C
	echo -e "{\tTChain s(\"tree\");" >> merge.C
	for fn in $outputdir/*_PassFail.root; do
	        if [ .$fn = ."" ]; then 
		    echo "ERROR : File _PassFail.root not found, skip"  
		        else 
		    echo -e "\ts.Add(\"$fn\");" >> merge.C
		        fi
		    i=$((i+1))
		    done
	
	echo -e "\ts.SetMaxTreeSize(1e9);" >> merge.C
	echo -e "\ts.SetBranchStatus(\"quadlep\", 0);" >> merge.C
	echo -e "\ts.Merge(\"$outputdir/wjets_data.root\");" >> merge.C
	echo "}" >> merge.C
	echo "Merging $PROCESS"
	root -l -q merge.C
	rm -f ${outputdir}/*_PassFail.root
    done
fi


# now add the special pass and fail cases for wjets MC 
if [ "$SELECTION" == 'PassFail' ]; then  	    
    for JETBIN in 0 1 2 ; do 
	outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"wjets.root\",\"$outputdir\",\"$SELECTION\",$JETBIN\);
	root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"wgamma_lgamma.root\",\"$outputdir\",\"$SELECTION\",$JETBIN\);
    done
fi
 
