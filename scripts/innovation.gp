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
# set xrange [-5:5]
# set yrange [-5:5]
plot [0:2*pi] r*sin(t), r*cos(t), \
   "./output/normalised_estimation.dat" using 4:5 with points pointsize 0.1


unset xrange 
unset yrange 
plot [0:2*pi] r*sin(t), r*cos(t), \
   "./output/normalised_estimation.dat" using 3:5 with points pointsize 0.1

pause -1

plot   "./output/normalised_estimation.dat" using 0:4 with linespoints pointsize 0.1

pause -1

plot   "./output/normalised_estimation.dat" using 0:5 with linespoints pointsize 0.1

pause -1
# set xrange [-5:5]
# set yrange [-5:5]


# splot  3*sin(v)*cos(u), 3*sin(v)*sin(u), 3*cos(v) with lines notitle,\
#   "./output/normalised_estimation.dat" using 1:2:3 with points pointsize 0.1
# pause -1
