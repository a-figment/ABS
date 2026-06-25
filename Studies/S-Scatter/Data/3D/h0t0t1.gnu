set terminal epslatex 12 color colortext standalone size 10,10

set output OUT

set xrange [0:B]
set yrange [-(pi/2):(3*pi/2)]
set zrange [-(pi/2):(3*pi/2)]

set xlabel "$h_0$"
set ylabel "$\\theta_0$"
set zlabel "$\\theta_1$"
set xyplane 0
set view 60,105,0.9
splot IN using 1:2:3 with dots notitle

set output 
system(sprintf("pdflatex %s", OUT))

