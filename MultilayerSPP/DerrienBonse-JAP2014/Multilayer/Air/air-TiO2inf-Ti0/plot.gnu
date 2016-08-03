#!gnuplot

reset

set output 'TiO2-Ti-Period.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set log x

xscale=1E9*2. #2 orignates from the Maple sheet, where a = 0.5 * thickness of the layer
yscale=1E9
set xlabel 'Thickness (nm)'
set ylabel 'SPP period (nm)'

plot "Result.tmp" u ($2*xscale):($3*yscale) w l lw 3 notitle


