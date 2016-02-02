#!gnuplot
reset

load '/media/PHD/Travail/Berlin/Codes/GNUplot/try2.plt'
set style increment userstyles

n1=1e0; k1=0e0; 

Snell(inc,N,K) = 180e0/pi * ( asin( sin(inc) / ({1,0}*N+{0,1}*K) ) )
Snell2(inc,N,K) = 180e0/pi * ( asin( real(sin(inc) / ({1,0}*N+{0,1}*K) ) ) )
angle(x) = 180e0/pi * atan(x) 
SnellKovalenko(inc,N,K)=180e0/pi * asin(sqrt(2e0)*sin(inc)/sqrt( N**2 - K**2 + sin(inc)**2 + sqrt( (N**2-K**2-sin(inc)**2)**2 + 4e0*N**2 * K**2 ) ))
Euler(inc, N, K) = 180e0/pi * sqrt( asin ( 0.5e0 * ( (1e0+n1**2*sin(inc)**2/(N**2+K**2))-sqrt( (1e0+n1**2*sin(inc)**2/(N**2+K**2))**2 - (4e0*N**2*n1**2*sin(inc)**2/(N**2+K**2)**2)) ) ) )
Chang(inc,N,K) = 180e0/pi*asin( sin(inc)/N * (1e0 - (N**2 * sin(inc)**2/(2e0*(N**2-sin(inc)**2)**2))*(K/N)**2) )
angleInc=45e0*pi/180e0

ReEps(x,y)=real(({1,0}*x+{0,1}*y)**2)
ImEps(x,y)=imag(({1,0}*x+{0,1}*y)**2)

#

######## refraction angle with complex indexes
##
##set output '20130610-RefractionOfnk.eps'
##set terminal postscript eps enhanced color font 'Helvetica, 18'
#

set view 67,224
set xlabel 'N'
set ylabel 'K'
splot "TreatedDatas.dat" u 1:2:(angle(($3/$5))) w p t 'Effective refraction angle'


"TreatedDatas.dat" u 1:2:(SnellKovalenko(angleInc,$1,$2)) w p lc 2 t 'Kovalenko', \
"TreatedDatas.dat" u 1:2:(Chang(angleInc,$1,$2)) w p lc 3 t 'Chang'

#"TreatedDatas.dat" u 1:2:(real(Snell(45e0*pi/180e0,$1,$2))) w p lc 3 t 'Real part of Snell law'

######### refraction angle with dielectric indexes
#set output '20130610-RefractionOfReImEps.eps'
#set terminal postscript eps enhanced color font 'Helvetica, 18'
#
#set view 67,224

set xlabel 'Re({/Symbol e})'
set ylabel 'Im({/Symbol e})'
set yrange [:1]
splot "TreatedDatas.dat" u (ReEps($1,$2)):(ImEps($1,$2)):(angle(abs($3/$5))) w p lw 1 t 'Refraction angle', \
"TreatedDatas.dat" u (ReEps($1,$2)):(ImEps($1,$2)):(SnellKovalenko(45e0*pi/180e0,$1,$2)) w p lw 1 t 'Kovalenko', \
"TreatedDatas.dat" u (ReEps($1,$2)):(ImEps($1,$2)):(Chang(45e0*pi/180e0,$1,$2)) w p lw 1 t 'Chang'

#"TreatedDatas.dat" u (ReEps($1,$2)):(ImEps($1,$2)):(real(Snell(45e0*pi/180e0,$1,$2))) w p t 'Real part of Snell law'



############ refraction angle in a more clear plot

set output '20130611-RefractionOfNk0.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'

ymin=0e0
ymax=45e0
set yrange [ymin:ymax]
set y2range [ymin:ymax]

set title 'Refraction in absorbing medium'
set xlabel 'Re(n)'
set ylabel 'Refraction angle (deg)'
# set y2label 'Snell law refractiona ngle (deg)'
# set ytics nomirror textcolor rgbcolor "red"
# set y2tics nomirror textcolor rgbcolor "blue"

plot "< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(angle(($3/$5))) w p lc 1 t 'FDTD, k=0', \
"< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(Snell(angleInc,$1,$2)) w l lc 1 lw 3 t 'Snell, k=0', \
"< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(Snell2(angleInc,$1,$2)) w l lc 1 lw 3 t 'Snell (real), k=0', \
"< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(Euler(angleInc,$1,$2)) w l lc 1 t 'Euler, k=0', \
"< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(SnellKovalenko(angleInc,$1,$2)) w l lc 1 lw 3 t 'Kovalenko, k=0', \
"< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(Chang(angleInc,$1,$2)) w l lc 1 t 'Chang, k=0.0'


set output '20130611-RefractionOfNk0.5.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'

plot "< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(angle(($3/$5))) w p lc 2 t 'FDTD, k=0.5', \
"< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(Snell(angleInc,$1,$2)) w l lc 1 lw 3 t 'Snell, k=0.5', \
"< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(Snell2(angleInc,$1,$2)) w l lc 1 lw 3 t 'Snell (real), k=0.5', \
"< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(Euler(angleInc,$1,$2)) w l lc 2 lw 3 t 'Euler, k=0.5', \
"< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(SnellKovalenko(angleInc,$1,$2)) w l lc 2 t 'Kovalenko, k=0.5', \
"< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(Chang(angleInc,$1,$2)) w l lc 3 t 'Chang, k=0.5'

set output '20130611-RefractionOfNk1.0.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'

plot "< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(angle(($3/$5))) w p lc 3 t 'FDTD, k=1', \
"< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(Snell(angleInc,$1,$2)) w l lc 1 lw 3 t 'Snell, k=1', \
"< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(Snell2(angleInc,$1,$2)) w l lc 1 lw 3 t 'Snell (real), k=1', \
"< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(Euler(angleInc,$1,$2)) w l lc 3 lw 3 t 'Euler, k=1.0', \
"< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(SnellKovalenko(angleInc,$1,$2)) w l lc 3 t 'Kovalenko, k=1.0', \
"< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(Chang(angleInc,$1,$2)) w l lc 3 t 'Chang, k=1.0'




# "< awk '{ if($2==0) print }' TreatedDatas.dat" u 1:(real(Snell(angleInc,$1,$2))) w l lc 1 t 'Re(Snell), k=0', \
# "< awk '{ if($2==0.5) print }' TreatedDatas.dat" u 1:(real(Snell(angleInc,$1,$2))) w l lc 3 t 'Re(Snell), k=0.5', \
#"< awk '{ if($2==1.0) print }' TreatedDatas.dat" u 1:(real(Snell(angleInc,$1,$2))) w l lc 4 t 'Re(Snell), k=1.0'
