#!gnuplot

homepath="/home/thibault/Documents/Codes/SPP-extended-theory/FDTD-validation/GsVit"

reset
#
### analytical case
## set output '20150731-AnalyticalCase.eps'
## set terminal postscript eps enhanced color font 'Helvetica, 24'
#
#c=3E8
#lambda=633e-9
#
#meshSize=3e-8
#sampleSize=60
#
#r=3e-2
#k=2e0*pi*c/lambda
#p=sampleSize*meshSize
#q=sampleSize*meshSize
#A=1e0
## x(r,theta,phi)=1e0 #r*sin(phi)*cos(theta)
## y(r,theta,phi)=1e0 #r*sin(phi)*sin(theta)
## z(r,theta,phi)=1e0 #r*cos(phi)
## intensity(x,y)=(A*p*q * (sin(k*p*x)/(2e0*r))/(k*p*x/(2e0*r)) * ( (sin(k*q*y/2e0/r)) / (k*q*y/2e0/r) ) )**2
#
## intensity(r,theta,phi)=(A*p*q*(sin(k*p*x(r,theta,phi)/(2e0*r))/(k*p*x(r,theta,phi)/(2e0*r)) )*( (sin(k*q*y(r,theta,phi)/2e0/r)) / (k*q*y(r,theta,phi)/2e0/r) ) )**2
#
#set format "%g"
#set xlabel 'X'
#set ylabel 'Y'
## set xrange [0:lambda]
## set yrange [0:lambda]
#set xrange [0e0:360e0]
#set yrange [0e0:360e0]
#
#set log cb
#set cbrange [1e-30:1e0]
#set sample 200,200
#set isosample 200
# 
#set view map
#splot intensity(x,y) w pm3d 


# PRE-TREATMENT
# system "bash ../../calcDiffractionPatternSphere.sh"	



# reset
fileR=system('ls xnff_set001_* -1 | tail -n1')
fileT=system('ls xnff_set001_* -1 | head -n1')

# plot pm3d with Diffraction.dat
# set terminal x11 2

set output '20150731-ReproduceTransmissionPattern.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel '{/Symbol q} index'
set ylabel '{/Symbol f} index'
set format "%g"
set zrange [:]
set palette rgb 34,35,36;
set view map 
set xrange [:]
set size square
# set pm3d interpolate 4,4
set log z
unset log cb
splot "< awk -f ".homepath."/pm3d.awk Diffraction.dat" u 1:2:4 w pm3d t 'Ey^2'
# splot "Diffraction.dat" u 1:2:5 w p t 'Ez'
# splot "Diffraction.dat" u 1:2:3 w p t 'Ex'

# reset
# calculating total scattered energy (summed on all calculated points 
system "echo 'Total intensity out of the simulation box'"
system "echo '/!\\ must contain enough points in resolution to be accurate'"
# system "awk -f '../sumLines.awk' xnff_set001_*"

## 
## plot r(t) for various angles
# plotting brut data

set output '20150815-RawData.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

#set yrange [1e-5:1e-1]
print fileR; print fileT
set log y
set xlabel 'Time step'
set ylabel 'Far field amplitude'
set y2label 'Near field amplitude'
set ytics nomirror
set y2tics nomirror 

plot fileT u 1:($2**2+$3**2+$4**2) w l t 'T', \
fileR u 1:($2**2+$3**2+$4**2) w l t 'R', \
"Field.dat" u 1:2 w l t 'Field amplitude' axis x1y2

