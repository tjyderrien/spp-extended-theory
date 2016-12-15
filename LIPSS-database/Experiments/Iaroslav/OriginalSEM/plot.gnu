#!gnuplot

reset

fileThibs='orientation.csv'
fileYoann='' #'10um.txt'

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Repeat.eps'

set xlabel 'Angle'
set ylabel 'Amplitude (un.ar.)'

set samples 10000

deltaAngle=1E0
angleZero=45
amplitude=3500E0
off=1E0

# sigmaAngle(angle) = angle / 2.*sqrt(2.*log(2.))
f(x) = amplitude*exp(-0.5e0*((x-angleZero)/(2e0*deltaAngle / (2e0*sqrt(2e0*log(2e0)))) )**2e0 ) + off
# try(x) = amplitude*exp(-0.5e0*((x-angleZero)/deltaAngle)**2.)+offset

fit f(x) fileThibs u 1:2 via angleZero, amplitude, deltaAngle, off
# fit try(x) fileThibs u 1:2 via angleZero, amplitude, deltaAngle, offset 

plot fileYoann u 1:2 w l t '10 um', \
fileThibs u 1:2 w l t '35 um', \
f(x) w l lw 3 t 'Fit'

