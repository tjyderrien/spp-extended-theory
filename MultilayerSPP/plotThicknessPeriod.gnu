#!gnuplot
reset

set output 'plotThicknessPeriod.eps'
set terminal postscript eps enhanced color font "Helvetica, 24"

set grid

datafile="PeriodWithThickness.csv"

#=================== PLOTTING WAVE NUMBERS AS FUNCTION OF THICKNESS (Berini style)

set output 'plotThicknessWaveNumbers.eps'
set terminal postscript eps enhanced color font "Helvetica, 24"
# set log x
# set format "%g"
set ytics format "%g"
set y2tics format "10^{%L}"
set xtics format "%g"
set xlabel 'Thickness (nm)'

set ytics nomirror
set y2tics nomirror

set xrange [:22]
set yrange [:1.7]
# set log x
# set log y
set log y2

set ylabel  'Re({/Symbol b}_{SPP})'
set y2label 'Im({/Symbol b}_{SPP})'

wavelength = 30e-9

c=3E8
omega = 2.*pi*c/wavelength
beta0 = omega / c

xscale=1E9
plot datafile u ($2*xscale):(2.*pi/$3 / beta0 ) w p ps 2 lw 5 lc 1           t 'Re {/Symbol b}_{SPP}', \
     datafile u ($2*xscale):(0.5/$4 / beta0 )   w p ps 2 lw 5 lc 2 axis x1y2 t 'Im {/Symbol b}_{SPP}'
     
#=================== PLOTTING PERIOD AND LSPP AS FUNCTION OF THICKNESS (Derrien style)
     
set output 'plotThicknessPeriod.eps'
set ylabel 'Period (m)'
set y2label 'Mean-free path L_{SPP} (m)'

set ytics format "10^{%L}"

set log y
set log y2

set yrange [:1e-7]

plot datafile u ($2*xscale):3 w p ps 1 lw 3 lc 1           t '{/Symbol L}_{SPP}', \
     datafile u ($2*xscale):4 w p ps 1 lw 3 lc 2 axis x1y2 t 'L_{SPP}'
     


