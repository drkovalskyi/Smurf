#!/bin/bash


echo "\documentclass{cmspaper}
\usepackage{graphicx}
\begin{document}
\title{Efficiency Measurement Summary}
\tableofcontents
\clearpage" > summary.tex


DIR="/smurf/dlevans/Efficiencies/V00-02-09/"

#
# ID and Iso efficiency
#

echo "\section{ID and Iso efficiencies}
The ID and Iso efficiencies are measured using the N-1 method.  
The ID is masured with respect to probes that pass Iso and 
Iso is measured with respect to probes that pass ID.
The signal is modelled by a MC template convoluted with
a gaussian to model resolution differences between data and MC.
The background in the failing probe sample is modelled 
with an exponential times an error function.  The efficiency 
is extracted by simultaneously fitting the samples of passing
and failing probes in the mass window $ 60<M_{TP}<120 $ GeV.
Scale factors are derived by comparing the efficiency measured in data
with the efficiency measured in MC by counting passing and failing
probes in the same mass window." >> summary.tex

for MEASUREMENT in `ls ${DIR} | grep NM1Eff_I | egrep -v MC`; 
do

    NAME=`echo ${MEASUREMENT} | sed 's/_/ /g'`
    N=4
    HEADER="\begin{tabular}{c|c|c|c}
    \hline & $ 0 < |\eta| < 0.8 $ & $ 0.8 < |\eta| < 1.2 $ & $ 1.2 < |\eta| < 2.4 $  \\\\"
    if [[ "$MEASUREMENT" == *Electron* ]]; then
        N=5
        HEADER="\begin{tabular}{c|c|c|c|c}
        \hline & $ 0 < |\eta| < 0.8$ & $ 0.8 < |\eta| < 1.479$ & $ 1.479 < |\eta| < 2 $ & $ 2 < |\eta| < 2.5 $  \\\\"
    fi

    echo "\subsection{${NAME}}
    \vspace{50pt}
    \begin{table}[!ht]
    \begin{center}
    ${HEADER}
    \hline
    \multicolumn{${N}}{c} {N-1 Efficiencies in data} \\\\
    \hline" >> summary.tex
    head -n 11 ${DIR}/${MEASUREMENT}/extra/dat_eff_table.tex | tail -n 6 >> summary.tex
    echo "  \hline
    \multicolumn{${N}}{c} {N-1 Efficiencies in simulation} \\\\
    \hline" >> summary.tex
    head -n 11 ${DIR}/${MEASUREMENT}/extra/mc_eff_table.tex | tail -n 6 >> summary.tex
    echo "  \hline
    \multicolumn{${N}}{c} {Simulation-to-data scale factors} \\\\
    \hline" >> summary.tex
    head -n 11 ${DIR}/${MEASUREMENT}/extra/sf_table.tex | tail -n 6 >> summary.tex
    echo "  \end{tabular}
    \caption{${NAME}}
    \end{center}
    \end{table}

    \clearpage

    " >> summary.tex

done

#
# Trigger efficiency
#

echo "\section{Trigger efficiencies}
Trigger efficiencies are measured with respect to the full offline selection
using couting in the mass window $ 81<M_{TP}<101 $ GeV.
In the case of the double triggers, the efficiency of each leg is measured 
separately, and the dZ cut is measured with respect to both legs passing
the entire trigger except for the dZ cut." >> summary.tex

for MEASUREMENT in `ls ${DIR} | grep NM1Eff_Trig | egrep -v MC`; 
do

    NAME=`echo ${MEASUREMENT} | sed 's/_/ /g'`
    echo "\subsection{${NAME}}
    \vspace{50pt}" >> summary.tex

    cat ${DIR}/${MEASUREMENT}/extra/dat_eff_table.tex >> summary.tex
    echo "\clearpage" >> summary.tex


done

echo "\end{document}" >> summary.tex

pdflatex summary.tex
pdflatex summary.tex

