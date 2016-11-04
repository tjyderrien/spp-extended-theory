#!gnuplot

reset
set output 'SummaryOfAnalysis.eps' 
set terminal postscript eps enhanced color font 'Helvetica, 24'

fileSPP="SimpleSPPmodel.csv"
fileSipe="ResultsZhenyaCu.csv"
fileExp="MultiwavelenthExperiment.csv"


set xlabel 'Laser wavelength {/Symbol l} (nm)'
set ylabel 'SPP period {/Symbol L} (nm)'
set xrange [0:]
#set key out above center
set key bottom right

set xtics 200
set ytics 200

plot fileSPP u 1:2 w l lt 1 lc 1 lw 2 t 'Ideally-flat SPP model', \
fileSipe u 1:4 w p lc 3 pt 4 ps 1.5 lw 3 t  'Roughness-coupled SPP model', \
fileExp u 1:2:3 w yerrorbars lc 7 lt 1 pt 7 ps 1.5 lw 3 t 'Experiment', \
fileSipe u 1:5 w p lc 7 pt 5 ps 1.5 lw 3 t 'air/liquid Cu SPP model'
# x w l lc 7 lw 2 lt 2 t '{/Symbol L} = {/Symbol l}', \
