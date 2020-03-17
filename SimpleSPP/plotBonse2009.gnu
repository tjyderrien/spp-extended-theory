#!gnuplot

filename="Si-1030nm-SipeDrude.csv"

set output '1030nm-Bonse2009.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel 'CB electron density (cm^{-3})'
set ylabel 'LIPSS period (nm)'
set y2label 'Efficacy factor (u.a.)'

set key bottom right

set ytics nomirror
set y2tics nomirror

plot filename u ($1*1E-6):($3) w lp lw 3 t 'Period (nm)', \
     filename u ($1*1E-6):($4) w lp lw 3 axis x1y2 t 'Efficacy factor'
