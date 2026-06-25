set terminal epslatex 12 color colortext standalone size 8,4

set output 'cunt.tex'

#set size 1, 1 
set size ratio 1
set multiplot layout 1,2

# Left panel
set xrange [0:0.5]
set yrange [-90:90]
set xlabel "$h_0$"
set ylabel "$\\theta_0$"
set ytics -90,30,90
set logscale cb

plot 'p0-a7-r0--ST.dat' using 1:2:3 with image t ''

# Right panel
reset
set size ratio 1
set xrange [0:1]
set yrange [-90:90]
set xlabel "$h_0$"
set ylabel "$\\theta_0$"
set ytics -90,30,90
plot 'p0-a7-r1--ST.dat' using 1:2:3 with image t ''


unset multiplot

set output # finish the current output file
system('latex cunt.tex && dvips cunt.dvi && ps2pdf cunt.ps')
