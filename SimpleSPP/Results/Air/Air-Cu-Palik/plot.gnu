#!gnuplot

reset
set output 'SummaryOfAnalysis.eps' 
set terminal postscript eps enhanced color font 'Helvetica, 24'

fileSPP="SimpleSPPmodel.csv"
fileSipe="ResultsZhenyaCu.csv"

set xlabel 'Laser wavelength (nm)'
set ylabel 'SPP period (nm)'

#set key out above center
set key bottom right

set xtics 200
set ytics 200

plot x w l lc 7 lw 2 lt 2 t '{/Symbol L} = {/Symbol l}', \
fileSPP u 1:2 w l lt 1 lc 1 lw 2 t 'Ideally-flat SPP model', \
fileSipe u 1:4 w p lc 3 lw 5 t  'Roughness-coupled SPP model'
# x w l lc 7 lw 2 t '{/Symbol L} = {/Symbol l}'
