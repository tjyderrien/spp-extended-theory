#!/bin/bash

# Plotting complex acos with GNUplot

set output 'ComplexAcos.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )
set cbrange [-pi : pi]
set cbtics ("0" -pi, "2{/Symbol p}}" pi)
set cblabel "Phase Angle" rotate offset -2,0

Hue(x,y) = (pi + atan2(-y,-x)) / (2*pi)
phase(x,y) = hsv2rgb( Hue(x,y), sqrt(x**2+y**2), 1. )

set xrange [-pi/2. : pi/2.]
set yrange [-pi/2. : pi/2.]
set urange [-pi/2. : pi/2.]
set vrange [-pi/2. : pi/2.]
set xtics ("-{/Symbol p}/2" -pi/2., "-{/Symbol p}/4" -pi/4., "0" 0, "{/Symbol p}/4" pi/4., "{/Symbol p}/2" pi/2.)
set ytics ("-{/Symbol p}/2" -pi/2., "-{/Symbol p}/4" -pi/4., "0" 0, "{/Symbol p}/4" pi/4., "{/Symbol p}/2" pi/2.)

set view map; set size square; unset key
set isosamples 100,100

rp(x,y) = real(f(x,y))
ip(x,y) = imag(f(x,y))
color(x,y) = hsv2rgb( Hue( rp(x,y), ip(x,y) ), abs(f(x,y)), 1. )

f(x,y) = acos(x + y*{0,1})
splot '++' using 1:2:(color($1,$2)) with pm3d lc rgb variable
