#!gnuplot

reset
load '/media/thibault/PHD/Travail/Berlin/Calculs/SingleSPPfile.gnu'

lambda=800e-9
xscale=1E6

set output '20141119-SPPvalidityRange.eps'
set terminal postscript eps enhanced color font 'Times, 26'

set key bottom left

set xtics format "10^{%L}"
set ytics nomirror format "%g"

set xrange [1e18:1e22]
set y2range [1e-3:1e2]

set log x
unset log y2

set xlabel 'Carrier density (cm^{-3})'
set ylabel 'Dielectric permittivity {/Symbol e}' 

plot real(epsilon(lambda,xscale*x,epsilonSi0,meffSi,nuSi)) w l lw 5 lc 1 t 'Re({/Symbol e})', \
imag(epsilon(lambda,xscale*x,epsilonSi0,meffSi,nuSi)) w l lw 5 lc 3 t 'Im({/Symbol e})'

###

set output '20141119-SPPapproximationValidity.eps'
set xrange [1e18:1e22]

set ytics nomirror format "%g"

set ylabel 'Approximation validity {/Symbol z} (%)'
plot 100e0*abs(1e0-imag(epsilon(lambda,x*xscale,epsilonSi0,meffSi,nuSi)) / real(epsilon(lambda,x,epsilonSi0,meffSi,nuSi))) w l lw 5 lc 7 t 'Perfect metal approximation' axis x1y2
