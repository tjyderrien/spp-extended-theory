reset
set output 'plot3dPeriod.eps'
set terminal postscript eps enhanced color font "Helvetica, 24"

unset log xz
set log y

set format "%g"
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'

# set pm3d interpolate 8,8
# set isosamples 100,100

set grid
splot "PeriodWithThickness.csv" u 1:2:3 w p notitle
