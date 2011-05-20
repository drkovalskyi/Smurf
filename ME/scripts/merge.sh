
#!/bin/sh

# Merges files in XXX directories 

if [ $# -lt 1 ]; then
	echo "USAGE : $0 <path> <process>"
	exit
fi

nfiles=`ls $2/*$1*.root | wc -l`

i=1

rm  merge.C
touch merge.C
echo -e "{\tTChain s(\"tree\");" >> merge.C

for fn in $2/$1_ME_*.root; do

	if [ .$fn = ."" ]; then 
		echo "ERROR : File $1_ME_$i.root not found, skip"  
	else 
		echo -e "\ts.Add(\"$fn\");" >> merge.C
	fi
	i=$((i+1))
done

echo -e "\ts.SetMaxTreeSize(1e9);" >> merge.C
echo -e "\ts.Merge(\"$1_ME_merged.root\");" >> merge.C
echo "}" >> merge.C

echo "Merging $1"
root -l -q merge.C

