#!gnuplot
author="Burke"
year="1986"
material="Ag"

set output 'Burke1986-SPP-result.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set size 1.3, 1.0

set format "%g"
wavelength=633e-9 #Ag
k0 = 2*pi/wavelength
ImBeta(Lspp) = 1./(2. * Lspp)
ReBeta(period) = 2. * pi / period

set colors podo

filename =  author."-SPPmodes.csv" # summary of everything
filename0 = author."-SPPmodes-branch0.csv" #branch (k2,k3) --
filename1 = author."-SPPmodes-branch1.csv" # branch -+
filename2 = author."-SPPmodes-branch2.csv" # branch +-
filename3 = author."-SPPmodes-branch3.csv" # branch ++

set xlabel 'Thickness (nm)'
# set ylabel 'SPP period (nm)'
set ylabel 'Re({/Symbol b})/k_0'

#unset key 
set key outside right Left reverse

xscale = 1E9
yscale = 1E0 #1E9
size_of_points = 1.3

set yrange [1.9:3.5] #like in the paper
unset log y
set ytics format "%3.2g" 

print "SPP period..."
plot filename0 u (xscale*$1):(ReBeta($7)/k0) w p lc rgb "red"    ps 1.3 pt 12 t 'Branch -, -', \
     filename1 u (xscale*$1):(ReBeta($7)/k0) w p lc rgb "blue"   ps 1.0 pt 9 t 'Branch -, +', \
     filename2 u (xscale*$1):(ReBeta($7)/k0) w p lc rgb "web-green"  ps 0.9 pt 5 t 'Branch +, -', \
     filename3 u (xscale*$1):(ReBeta($7)/k0) w p lc rgb "orange" ps size_of_points-0.5 pt 7 t 'Branch +, +'
#      wavelength*yscale w l lc rgb "black" lw 2 dt 3 notitle #t 'Laser wavelength {/Symbol l} (nm)'

# set output author.year.'-SPP-subwavelength.eps'
# set yrange [0:1200]
# replot

# ### PRINTING EPSILON(oxide)
# 
set output author.year.'-Lspp.eps'
# set ylabel 'L_{SPP} (m)'
set ylabel 'Im({/Symbol b})/k_0'
set ytics format "10^{%L}" 10
print "Lspp..."
set log y
set yrange [1e-5:1]   
plot filename0 u (xscale*$1):(ImBeta($8)/k0) w p lc rgb "red"    ps 1.3 pt 12 t 'Branch -, -', \
     filename1 u (xscale*$1):(ImBeta($8)/k0) w p lc rgb "blue"   ps 1.0 pt 9 t 'Branch -, +', \
     filename2 u (xscale*$1):(ImBeta($8)/k0) w p lc rgb "web-green"  ps 0.9 pt 5 t 'Branch +, -', \
     filename3 u (xscale*$1):(ImBeta($8)/k0) w p lc rgb "orange" ps size_of_points-0.5 pt 7 t 'Branch +, +'

