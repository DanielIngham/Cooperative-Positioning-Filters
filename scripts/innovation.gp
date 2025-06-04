set grid

set parametric
r = 3
plot [0:2*pi] r*sin(t), r*cos(t), \
"./output/normalised_innovation.dat" using 1:2 with points pointsize 0.1
  
pause -1

unset object 1 
set parametric
set urange [0:2*pi]
set vrange [0:pi]
set view equal xyz
set isosamples 50,50
r = 3
plot [0:2*pi] r*sin(t), r*cos(t), \
   "./output/normalised_estimation.dat" using 1:3 with points pointsize 0.1
# splot  3*sin(v)*cos(u), 3*sin(v)*sin(u), 3*cos(v) with lines notitle,\
#   "./output/normalised_estimation.dat" using 1:2:3 with points pointsize 0.1
pause -1
