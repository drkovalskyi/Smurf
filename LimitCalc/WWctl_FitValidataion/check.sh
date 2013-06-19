#!/bin/sh 

for x in {0..99}
do 

    #export newdefault=fitoutput/log_newdefault_${x}.txt
    #export CR1=fitoutput/log_CR1_${x}.txt
    #export CR2=fitoutput/log_CR2_${x}.txt

    export mu_newdefault=`grep POI fitoutput/log_newdefault_${x}.txt | awk '{print $4 "+" $6 "-" $8}'`
    export mu_CR1=`grep POI fitoutput/log_CR1_${x}.txt | awk '{print $4 "+" $6 "-" $8}'`
    export mu_CR2=`grep POI fitoutput/log_CR2_${x}.txt | awk '{print $4 "+" $6 "-" $8}'`

    echo -e " $x :: \t $mu_newdefault \t $mu_CR1 \t $mu_CR2"
done 
