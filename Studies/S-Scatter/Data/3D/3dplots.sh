#!/bin/bash

gnuplot -e "IN='p0-a7-r0--H0T0RN.dat'; OUT='h0t0t1--FH.tex'; B=0.5" h0t0t1.gnu
gnuplot -e "IN='p0-a7-r0--H0T0RN.dat'; OUT='h0t0h1--FH.tex'; B=0.5" h0t0h1.gnu

gnuplot -e "IN='p0-a7-r1--H0T0RN.dat'; OUT='h0t0t1--IH.tex'; B=1.0" h0t0t1.gnu
gnuplot -e "IN='p0-a7-r1--H0T0RN.dat'; OUT='h0t0h1--IH.tex'; B=1.0" h0t0h1.gnu

#rm 3DPlot.pdf
rm *eps* 
#qpdf --empty --pages "h0t0t1--FH.pdf h0t0h1--FH.pdf h0t0t1--IH.pdf h0t0h1--IH.pdf" -- 3DPlot.pdf
qpdf --empty --pages *pdf -- 3DPlot.pdf
zathura 3DPlot.pdf 
#gnuplot -e "'IN=p0-a7-r0--H0T0RN.dat'; OUT='h0t0h1--FH.tex'" h0t0h1.gnu
#gnuplot -e "'IN=p0-a7-r1--H0T0RN.dat'; OUT='h0t0t1--IH.tex'" h0t0t1.gnu
#gnuplot -e "'IN=p0-a7-r1--H0T0RN.dat'; OUT='h0t0h1--IH.tex'" h0t0h1.gnu

