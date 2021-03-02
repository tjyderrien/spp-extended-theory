#!gnuplot

reset

set output 'Al203-Al-Period.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set log x

xscale=1E9*2. #2 orignates from the Maple sheet, where a = 0.5 * thickness of the layer
yscale=1E9

set xlabel 'Thickness (nm)'
set ylabel 'SPP period (nm)'

plot "515nm/Result.tmp" u ($2*xscale):($3*yscale) w l lw 3 lc 1 t '515 nm', \
"800nm/Result.tmp" u ($2*xscale):($3*yscale) w l lw 3 lc 3 t '800 nm', \
"1030nm/Result.tmp" u ($2*xscale):($3*yscale) w l lw 3 lc 7 t '1030 nm'


