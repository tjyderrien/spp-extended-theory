 #!/usr/bin/gnuplot -persist

 
SPPpath=system("echo $SPPextPath")
localPath=SPPpath."/MultilayerSPP/DerrienBonse-JAP2014"

############# plot a contour map

reset
# 
# set output '20130701-Contour.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

set style increment userstyles
set contour surface
set cntrparam levels discrete 50,100,200,400,600,800 
set cntrparam cubicspline
set key out right
# set view map
# unset surface

set xlabel "Free-carrier density N_e (cm^{-3})" 
set xtics format "10^{%L}"
# set xrange [ 6 : 100 ] noreverse nowriteback  # (currently [:1.00000e+23] )
set xrange [1.3e20:1.07e23]

# set hidden3d

set ylabel "Excited layer thickness t (nm)" 
# set yrange [ * : 15 ] noreverse nowriteback  # (currently [12.0000:0.00000] )
set yrange [0.1:2274]
set ytics 5

set title "SPP period {/Symbol L} (nm)" 
set ztics 10
# set zrange [ : 800 ]

set pm3d interpolate 16,16
# set palette rgbformulae 30, 31, 32

set log xy

scale=1e6

# plot "cont.dat" w l notitle
splot "Result.awked.tmp" u ($1*1e-6):($2*1e9*2e0):($3*1e9) w l lw 1 lc rgb "black" notitle
# "Result.awked.tmp" u ($1/scale):(1e3):(period( beta( epsilonWater, epsilon( $1,epsilonSi0,meffSi,nuSi )) )*1e9) w l lc rgb "black" t 'water/Si*'
# "Result.awked.tmp" u ($1/scale):(100e0):((real(epsilon($1)) < -1e0) ? 0.5e0*period( beta( epsilonWater, epsilon( $1 )) )*1e9 : 1/0) w l lc rgb "black" t 'Single SPP case w/ standing wave'

set output
set terminal x11
replot

################# plot the case of water/a-Si*/Si and compare with Miyaji

reset
set output '20130730-ValidationMyaji1.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

# load localPath.'/try2.plt'

# set style increment userstyles
# set log cb
# set pm3d interpolate 8,8
# set style fill solid

load localPath.'/SingleSPPfile.gnu'

#set xrange [4:100]

set xlabel "Free-carrier density N_e (cm^{-3})" 
set xtics format "10^{%L}"
# set ylabel "Excited layer thickness t (nm)" 
set ylabel 'Periodicity {/Symbol L} (nm)'

# set parametric
# n=10 #number of curves
set key top right font 'Helvetica, 16' spacing 0.7

set log x

scale=1e6

plot "Result.awked.tmp" u ($1/scale):(period( beta( epsilonWater, epsilon( $1, epsilonaSi0, meffaSi, nuaSi)) )*1e9) w p t 'SPP Water/a-Si*', \
"Result.awked.tmp" u ($1/scale):(period( beta( epsilon( $1, epsilonaSi0, meffaSi, nuaSi), epsilonSi0) )*1e9) w p t 'SPP a-Si*/Si', \
"Result.awked.tmp" u ($1/scale):(period( beta( epsilon( $1, epsilonaSi0, meffaSi, nuaSi), epsilon( $1, epsilonSi0, meffSi, nuSi)) )*1e9) w p t 'SPP a-Si*/Si*', \
"< awk '{ if($2==.1000000003e-5) print }' Result.awked.tmp" u ($1/scale):($3*1e9) w l t 'Multilayer SPP, t=1 um'

##################### SINGLE SPP periodicity considering SPP conditions ######################
reset

load localPath.'/try2.plt'
load localPath.'/SingleSPPfile.gnu'

set style increment userstyles

set output '20130731-SPPperiodAmorphousSi.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

set xlabel 'Free-carrier density N_e (m^{-3})'
set ylabel 'SPP periodicity {/Symbol L} (nm)'

set samples 1000
set key font 'Helvetica, 16' spacing 0.7

# water0 / a-Si interface

set xrange [1e26:1E29]
set xtics format "10^{%L}" 10e0
set log x

plot (ConditionPerfect(epsilonWater, epsilon(x,epsilonaSi0, meffaSi, nuaSi))) ? (period( beta( epsilonWater, epsilon(x,epsilonaSi0,meffaSi,nuaSi )) )*1e9) : (1/0) w l lw 2 t 'water/a-Si*, perfect', \
(ConditionAbs(epsilonWater, epsilon(x,epsilonaSi0, meffaSi, nuaSi))) ? (period( beta( epsilonWater, epsilon(x,epsilonaSi0,meffaSi,nuaSi )) )*1e9) : (1/0) w l lw 1 t 'water/a-Si*, absorbing', \
(ConditionPerfect(epsilon(x,epsilonaSi0, meffaSi, nuaSi),epsilonSi0)) ? (period( beta( epsilon(x,epsilonaSi0,meffaSi,nuaSi ), epsilonSi0) )*1e9) : (1/0) w l lw 2 t 'Sa-Si*/Si, perfect', \
(ConditionAbs(epsilon(x,epsilonaSi0, meffaSi, nuaSi),epsilonSi0)) ? (period( beta( epsilon(x,epsilonaSi0,meffaSi,nuaSi ), epsilonSi0) )*1e9) : (1/0) w l lw 1 ls 4 t 'a-Si*/Si, absorbing'


set output '20130731-SPPperiodExcitedSi.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'
set key bottom right font 'Helvetica, 16' spacing 0.7

scale=1e6

ncrLocal=Ncr(2e0*pi*c/lambda,meffSi,nuSi,epsilonSi0)
set arrow 1 from ncrLocal, 200 to ncrLocal, 900 nohead

plot "< awk '{ if($2==.9999999942e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 5 t 'air/Si*/Si, t=2 um', \
"< awk '{ if($2==.9999999956e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 5 t 'air/Si*/Si, t=200 nm', \
"< awk '{ if($2==.6309573419e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 5 t 'air/Si*/Si, t=120 nm', \
"< awk '{ if($2==.3981071690e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 5 t 'air/Si*/Si, t=80 nm', \
"< awk '{ if($2==.2511886423e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 5 t 'air/Si*/Si, t=50 nm', \
(period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonAir) )*1e9) w l lw 2 t 'SPP air/Si*', \
(ConditionAbs(epsilon(x,epsilonSi0, meffSi, nuSi),epsilonAir)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonAir) )*1e9) : (1/0) w l lw 10 t 'SPP air/Si*, restricted', \
1e9*lambda/real(sqrt(epsilon(x,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/n'

# "< awk '{ if($2==.2511886426e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 5 t 'air/Si*/Si, t=5 nm'

set output '20130731-SPPperiodExcitedSiWater.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'
# set key out top center horizontal font 'Helvetica, 16' spacing 0.7
set key inside bottom right vertical font 'Helvetica, 16' spacing 0.7

# set xrange [1E27:]
set yrange [:800]
set ytics 200
scale=1e6

ncrLocal=Ncr(2e0*pi*c/lambda,meffSi,nuSi,epsilonSi0)
set arrow 1 from ncrLocal, 200 to ncrLocal, 800 nohead

plot (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonWater) )*1e9) w l lw 2 t 'water/Si*', \
(ConditionAbs(epsilon(x,epsilonSi0, meffSi, nuSi),epsilonWater)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonWater) )*1e9) : (1/0) w l lw 10 t 'water/Si*, restricted', \
"< awk '{ if($2==1e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Si*/Si, t=1 um', \
"< awk '{ if($2==500e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Si*/Si, t=500 nm', \
"< awk '{ if($2==200e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Si*/Si, t=200 nm', \
"< awk '{ if($2==100e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Si*/Si, t=100 nm', \
"< awk '{ if($2==50e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Si*/Si, t=50 nm', \
1e9*lambda/real(sqrt(epsilon(x,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/n'


# "< awk '{ if($2==10e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Si*/Si, t=10 nm', \

set output '20130731-SPPperiodSiExcitedWaterSimple.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'
set key inside top right font 'Helvetica, 16' spacing 0.7
# set xrange [1E27:]
unset log y
set yrange [:2000]
scale=1e6

plot (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonWater) )*1e9) w l lw 2 t 'Water/Si*', \
(ConditionAbs(epsilon(x,epsilonSi0,meffSi,nuSi),epsilonWater)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonWater) )*1e9) : (1/0) w l lw 10 t 'Water/Si*, restricted', \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSi0) )*1e9) w l lw 2 t 'Water*/Si', \
(ConditionAbs(epsilon(x,epsilonWater,meffWater,nuWater),epsilonSi0)) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSi0) )*1e9) : (1/0) w l lw 10 t 'Water*/Si, restricted', \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) w l lw 2 t 'Water*/Si*', \
(ConditionAbs(epsilon(x,epsilonWater,meffWater,nuWater),epsilon(x,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 t 'Water*/Si*, restricted'

reset
set output '20130731-SPPperiodSiExcitedWater.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

set style increment userstyles
set key inside bottom right font 'Helvetica, 16' spacing 0.7
set xrange [1E26:1e29]
set log x
set xtics format "10^{%L}" 10
set log y
set yrange [10:3e3]
scale=1e6


ncrLocal=Ncr(2e0*pi*c/lambda,meffWater,nuWater,epsilonWater)
set arrow 1 from ncrLocal, 10 to ncrLocal, 3e3 nohead

plot (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonWater) )*1e9) w l lw 2 t 'Water/Water*', \
ConditionAbs(epsilon(x,epsilonWater,meffWater,nuWater),epsilonWater) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonWater) )*1e9) : (1/0) w l lw 10 t 'Water/Water*, restricted', \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSi0) )*1e9) w l lw 2 t 'Water*/Si', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilonSi0)) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSi0) )*1e9) : (1/0) w l lw 10 t 'Water*/Si, restricted', \
"< awk '{ if($2==100e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=100 um', \
"< awk '{ if($2==10e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=10 um', \
"< awk '{ if($2==1e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=1 um', \
"< awk '{ if($2==500e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=500 nm', \
"< awk '{ if($2==200e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=200 nm'



"< awk '{ if($2==200e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=200 nm', \
"< awk '{ if($2==100e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=100 nm', \
"< awk '{ if($2==50e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=50 nm', \
"< awk '{ if($2==10e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water/Water*/Si, t=10 nm'

################ water/Si*/Si

reset
lambda=790e-9
set output '20130731-SPPperiodSiExcitedSi.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

set style increment userstyles
# set key left Left font 'Helvetica, 10' spacing 0.8
set key inside bottom right
set xlabel 'Free-carrier density N_e (cm^{-3})'

set log x
set xtics format "%1.0l x 10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
set log y
set yrange [10:1.0e3]
set ytics nomirror 10

set y2label 'Normalized SPP period {/Symbol L}/{/Symbol l}'
set y2range [10/(lambda*1E9):1E3/(lambda*1E9)]
set log y2
set y2tics nomirror 10

xscale=1e6

ncrLocal=Ncr(lambda,meffSi,nuSi,epsilonSi0)
set arrow 1 from ncrLocal/xscale, 10 to ncrLocal/xscale, 1.0e3 nohead

set xrange [ncrLocal/xscale:5e28/xscale]

plot (period( beta(lambda, epsilonWater, epsilon(lambda,xscale*x,epsilonSi0,meffSi,nuSi)) )*1e9) w l lc 7 lw 5 t 'Water/Si*', \
(period( beta(lambda, epsilonSi0, epsilon(lambda,x*xscale,epsilonSi0,meffSi,nuSi)) )*1e9) w l lc 3 lw 5 t 'Si*/Si', \
"< awk '{ if($2==.100e-6 && $3>400e-9 && $1>2e27) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w p smooth bezier lw 5 t 'Water/Si*/Si, t=200 nm', \
"< awk '{ if($2==.50e-7 && $3>300e-9) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w p smooth bezier lw 5 t 'Water/Si*/Si, t=100 nm', \
"< awk -v ncr=ncrLocal '{ if($2==0.5e-8 && (($1<5e27 && $3>500e-9) || ($1>5e27) ) ) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lc 1 lw 5 t 'Water/Si*/Si, t=10 nm', \
lambda*1E9 w l t 'Laser wavelength {/Symbol l}_0' lc 7 lw 2

set table 'Values.dat'
set xtics format "%g"
set ytics format "%g"
replot
unset table


1E9*lambda/real(sqrt(epsilon(x,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/sqrt({/Symbol e}) Si*'


"< awk '{ if($2==0.4e-7) print }' LargerTickness/Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=80 nm', \
(ConditionAbs( epsilonWater, epsilon(x,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilonWater, epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 t 'Water/Si*, restricted', \
(ConditionAbs( epsilonSi0, epsilon(x,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilonSi0, epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 t 'Si*/Si, restricted', \
"< awk '{ if($2==.300e-6) print }' LargerTickness/Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=600 nm', \
"< awk '{ if($2==.200e-6) print }' LargerTickness/Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=400 nm', \
"< awk '{ if($2==.150e-6) print }' LargerTickness/Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=300 nm', \
"< awk '{ if($2==0.2e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=40 nm', \
"< awk '{ if($2==0.1e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=20 nm', \
"< awk '{ if($2==0.25e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w p  lw 3 t 'Water/Si*/Si, t=5 nm', \

################ air/Si*/Si

reset
set output '20130731-SPPperiodSiExcitedSi.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

lambda=790e-9

set style increment userstyles
set key inside bottom right
# set key out
set xlabel 'Free-carrier density N_e (cm^{-3})'

set log x
set xtics format "%1.0l x 10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
set log y
set yrange [10:1.2e3]
set ytics nomirror

set log y2
set y2label 'Normalized SPP period {/Symbol L}/{/Symbol l}'
set y2tics nomirror 10
set y2range [10e0/(lambda*1E9):1.2e3/(lambda*1E9)]

xscale=1e6

ncrLocal=Ncr(lambda,meffSi,nuSi,epsilonSi0)
set arrow 1 from ncrLocal/xscale, 10 to ncrLocal/xscale, 1.2e3 nohead

set xrange [ncrLocal/xscale:5E28/xscale]

plot (period( beta(lambda, epsilonAir, epsilon(lambda, xscale*x,epsilonSi0,meffSi,nuSi)) )*1e9) w l lw 4 lc 7 t 'Air/Si*', \
(period( beta( lambda, epsilonSi0, epsilon(lambda, xscale*x,epsilonSi0,meffSi,nuSi)) )*1e9) w l lw 5 lc 3 t 'Si*/Si', \
"< awk '{ if($2==.100e-6) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 7 t 'Air/Si*/Si, t=200 nm', \
"< awk '{ if($2==.50e-7) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 7 t 'Air/Si*/Si, t=100 nm', \
"< awk '{ if($2==0.5e-8) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 5 lc 1 t 'Air/Si*/Si, t=10 nm', \
lambda*1E9 w l t 'Laser wavelength {/Symbol l}_0' lc 7 lw 2

set table 'Values.dat'
set xtics format "%g"
set ytics format "%g"
replot
unset table


"< awk '{ if($2==0.4e-7 && $3>700e-9) print }' LargerThickness/Result.awked.tmp" u ($1):($3*1e9) w p smooth bezier lw 3 t 'Air/Si*/Si, t=80 nm', \


(ConditionAbs( epsilonSi0, epsilon(x,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilonSi0, epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 t 'Si*/Si, restricted', \
(ConditionAbs( epsilonAir, epsilon(x,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilonAir, epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 t 'air/Si*, restricted', \
"< awk '{ if($2==.300e-6) print }' LargerThickness/Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Air/Si*/Si, t=600 nm', \
"< awk '{ if($2==.200e-6) print }' LargerThickness/Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Air/Si*/Si, t=400 nm', \
"< awk '{ if($2==.150e-6) print }' LargerThickness/Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Air/Si*/Si, t=300 nm', \


"< awk '{ if($2==0.2e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Air/Si*/Si, t=40 nm', \
"< awk '{ if($2==0.1e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Air/Si*/Si, t=20 nm', \
"< awk '{ if($2==0.25e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Air/Si*/Si, t=5 nm'


# 1E9*lambda/real(sqrt(epsilon(x,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/sqrt({/Symbol e}) Si*'

# "Experiment.dat" u 1:2:3 w yerrorbar t 'Experiment: this work', \
# "Costache.dat" u 1:2:3 w yerrorbar t 'Experiment: Costache 2004'


################ air/SiO2*/SiO2

reset
set output '20130731-SPPperiodSiExcitedSi.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

lambda=800e-9

set style increment userstyles
set key inside bottom right font 'Helvetica, 20'
# set key out
set xlabel 'Free-carrier density N_e (cm^{-3})'

set log x
set xtics format "%1.0l x 10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
set log y
set yrange [10:1.2e3]
set ytics nomirror

set log y2
set y2label 'Normalized SPP period {/Symbol L}/{/Symbol l}'
set y2tics nomirror 10
set y2range [10/(lambda*1E9):1.2e3/(lambda*1E9)]

xscale=1e6

ncrLocal=Ncr(2e0*pi*c/lambda,meffSi,nuSi,epsilonSi0)
set arrow 1 from ncrLocal/xscale, 10 to ncrLocal/xscale, 1.2e3 nohead

set xrange [ncrLocal/xscale:4E28/xscale]

plot (period( beta( epsilonAir, epsilon(xscale*x,epsilonSiO2,meffSiO2,nuSiO2)) )*1e9) w l lw 4 lc 7 t 'Air/SiO_2*', \
(period( beta( epsilonSi0, epsilon(xscale*x,epsilonSiO2,meffSiO2,nuSiO2)) )*1e9) w l lw 5 lc 3 t 'SiO_2*/SiO_2', \
"< awk '{ if($2==100e-9) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 7 t 'Air/SiO_2*/SiO_2, t=200 nm', \
"< awk '{ if($2==50e-9) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 7 t 'Air/SiO_2*/SiO_2, t=100 nm', \
"< awk '{ if($2==20e-9) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 7 t 'Air/SiO_2*/SiO_2, t=40 nm', \
"< awk '{ if($2==5e-9) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l smooth bezier lw 5 lc 1 t 'Air/SiO_2*/SiO_2, t=10 nm', \
lambda*1E9 w l t 'Laser wavelength {/Symbol l}_0' lc 7 lw 2

set table 'Values.dat'
set xtics format "%g"
set ytics format "%g"
replot
unset table



################ air/SiO2*/Si

reset
set output '20130731-SPPperiodSiExcitedSiO2.eps'
set terminal postscript eps enhanced monochrome font 'Times, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

lambda=790e-9

set style increment userstyles
set key left Left font 'Times, 16' spacing 0.6
set xlabel 'Free-carrier density N_e (m^{-3})'
set xrange [1E26:1e29]
set log x
set xtics format "10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
unset log y
set yrange [10:1.2e3]

scale=1e6

ncrLocal=Ncr(2e0*pi*c/lambda,meffSiO2,nuSiO2,epsilonSiO2)
set arrow 1 from ncrLocal, 10 to ncrLocal, 1.2e3 nohead

plot (period( beta( epsilonAir, epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) w l lw 2 t 'air/Si*', \
(ConditionAbs( epsilonAir, epsilon(x,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilonAir, epsilon(x,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 t 'air/Si*, restricted', \
(period( beta( epsilonSi0, epsilon(x,epsilonSiO2,meffSiO2,nuSiO2)) )*1e9) w l lw 2 t 'SiO_2*/Si', \
(ConditionAbs( epsilonSi0, epsilon(x,epsilonSiO2,meffSiO2,nuSiO2))) ? (period( beta( epsilonSi0, epsilon(x,epsilonSiO2,meffSiO2,nuSiO2)) )*1e9) : (1/0) w l lw 10 t 'SiO_2*/Si, restricted', \
(period( beta( epsilonAir, epsilon(x,epsilonSiO2,meffSiO2,nuSiO2)) )*1e9) w l lw 2 t 'air/SiO_2*', \
(ConditionAbs( epsilonAir, epsilon(x,epsilonSiO2,meffSiO2,nuSiO2))) ? (period( beta( epsilon(x,epsilonSiO2,meffSiO2,nuSiO2), epsilonAir) )*1e9) : (1/0) w l lw 10 t 'air/SiO_2*, restricted', \
"< awk '{ if($2==0.1e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'air/SiO_2*/Si, t=100 nm', \
"< awk '{ if($2==10e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'air/SiO_2*/Si, t=20 nm', \
"< awk '{ if($2==7e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'air/SiO_2*/Si, t=14 nm', \
"< awk '{ if($2==5e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'air/SiO_2*/Si, t=10 nm', \
1E9*lambda/real(sqrt(epsilon(x,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/sqrt({/Symbol e}) Si*', \
"Experiment.dat" u 1:2:3 w yerrorbar t 'Experiment: this work', \
"Costache.dat" u 1:2:3 w yerrorbar t 'Experiment: Costache 2004'


################ Water* / Silicon oxide ###########

reset
set output '20130731-SPPperiodSiExcitedWater.eps'
set terminal postscript eps enhanced monochrome font 'Times, 26' size 20cm, 10cm

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

set style increment userstyles
# set key bottom left Left font 'Times, 16' spacing 0.6
set key out
set xlabel 'Free-carrier density N_e (m^{-3})'
set xrange [1E26:1e29]
set log x
set xtics format "10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
unset log y
set yrange [10:2.0e3]

scale=1e6

nSi=4.6e27

ncrLocal=Ncr(2e0*pi*c/lambda,meffWater,nuWater,epsilonWater)
nsppLocal=Nspp(2e0*pi*c/lambda,meffWater,nuWater,epsilonWater)
set arrow 1 from ncrLocal, 10 to ncrLocal, 2.0e3 nohead
# set arrow 2 from nsppLocal, 10 to nsppLocal, 2.0e3 nohead


plot (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSiO2) )*1e9) w l lw 2 t 'Water*/SiO_2', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilonSiO2)) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSiO2) )*1e9) : (1/0) w l lw 10 notitle, \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilon(nSi,epsilonSi0,meffSi,nuSi)) )*1e9) w l lw 2 t 'Water*/Si*', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilon(nSi,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilon(nSi,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l lw 10 notitle, \
"< awk '{ if($2==50e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Water*/SiO_2/Si, t=100 nm', \
"< awk '{ if($2==5e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Water*/SiO_2/Si, t=10 nm'

# 1E9*lambda/real(sqrt(epsilonSiO2)) w l t '{/Symbol l}/n, SiO_2'
# 1E9*lambda/real(sqrt(epsilon(x,epsilonWater,meffWater,nuWater))) w l t '{/Symbol l}/n*, Water*'
#t 'Water*/Si*, restricted', \
#t 'Water*/SiO_2, restricted'

# "< awk '{ if($2==20e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Water*/SiO_2/Si, t=40 nm', \
# "< awk '{ if($2==10e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Water*/SiO_2/Si, t=20 nm', \
# "< awk '{ if($2==1e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Water*/SiO_2/Si, t=2000 nm', \
# "../Experiment.dat" u 1:2:3 w yerrorbar t 'Experiment Bonse et al'

# "< awk '{ if($2==2.5e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w p lw 3 t 'Water*/SiO_2/Si, t=5 nm', \

set table 'Values.dat'
set xtics format "%g"
set ytics format "%g"
replot
unset table

"< awk '{ if($2==1e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=2 nm', \
"< awk '{ if($2==7e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=14 nm', \
"< awk '{ if($2==.10e-4) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=10 um', \
"< awk '{ if($2==.9999999942e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=2 um', \
"< awk '{ if($2==.9999999956e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=200 nm', \
"< awk '{ if($2==.9999999971e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=20 nm', \
"< awk '{ if($2==.9999999985e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=2 nm', \


################ Water* / Silicon oxide / Silicon* ###########

reset
set output '20131127-SPPperiodExcitedSiExcitedWater.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

lambda=790e-9

set style increment userstyles

set xlabel 'Free-carrier density N_e (cm^{-3})'

set log x
set xtics format "%1.0l x 10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
set log y
set yrange [10:3e3]
set ytics nomirror

set y2label 'Normalized SPP period {/Symbol L}/{/Symbol l}'
set log y2
set y2range [10/(lambda*1E9):3E3/(1E9*lambda)]
set y2tics nomirror 10

set key bottom center font 'Helvetica, 16' lmargin 0

xscale=1e6
neFixed=2.08E27 #4.6e27

ncrWater=Ncr(lambda,meffWater,nuWater,epsilonWater)
nsppWater=Nspp(lambda,meffWater,nuWater,epsilonWater)
ncrSi=Ncr(lambda,meffSi,nuSi,epsilonSi0)
nsppSi=Nspp(lambda,meffSi,nuSi,epsilonSi0)
ncrSiO2=Ncr(lambda,meffSiO2,nuSiO2,epsilonSiO2)
nsppSiO2=Nspp(lambda,meffSiO2,nuSiO2,epsilonSiO2)
set arrow 1 from ncrWater, 10 to ncrWater, 3e3 nohead
# set arrow 2 from nsppWater, 10 to nsppWater, 2.0e3 nohead
# set arrow 3 from ncrSiO2, 10 to ncrSiO2, 2.0e3 nohead
# set arrow 4 from nsppSiO2, 10 to nsppSiO2, 2.0e3 nohead

set xrange [ncrWater/xscale:5e28/xscale]

plot (ConditionAbs( epsilon(lambda,x*xscale,epsilonWater,meffWater,nuWater),epsilonSiO2)) ? (period( beta( lambda, epsilon(lambda,x*xscale,epsilonWater,meffWater,nuWater), epsilonSiO2) )*1e9) : (1/0) w l lc 7 lw 4 smooth bezier t 'Water*/SiO_2', \
"< awk '{ if($2==10e-9 && ($1>1.5e27)) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l lw 4 t 'Water*/SiO_2/Si*, t=20 nm', \
"< awk '{ if($2==5e-9 && ($1>1.5e27)) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l lw 5 t 'Water*/SiO_2/Si*, t=10 nm', \
"< awk '{ if($2==2.5e-9 && ($1>1.5e27)) print }' Result.awked.tmp" u ($1/xscale):($3*1e9) w l lw 4 lc 1 t 'Water*/SiO_2/Si*, t=5 nm', \
lambda*1E9 w l lc 7 lw 3 t 'Laser wavelength {/Symbol l}_0'

plot "< awk '{ if($2==5e-9 && ($1>1.5e27)) print }' Damping.awked.tmp" u ($1):($3*1e0) w l lw 4 t 'Damping, t=10 nm'

# "< awk '{ if($2==25e-9 && ($1<1.3e27 || $1>1.85e27)) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 lc 2 t 'Water*/SiO_2/Si*, t=50 nm', \
# (ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilon(neFixed,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilon(neFixed,epsilonSi0,meffSi,nuSi)) )*1e9) : (1/0) w l smooth bezier lw 4 lc 7 t 'Water*/Si*', \
# (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilon(neFixed,epsilonSi0,meffSi,nuSi)) )*1e9) w l smooth bezier lw 2 lc 7 t 'Water*/Si*', \
# (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSiO2) )*1e9) w l lc 7 lw 2 smooth bezier t 'Water*/SiO_2', \

set table 'Values.dat'
set xtics format "%g"
set ytics format "%g"
replot
unset table

1E9*lambda/real(sqrt(epsilon(x,epsilonWater,meffWater,nuWater))) w l t '{/Symbol l}/n*, Water*', \
1E9*lambda/real(sqrt(epsilon(neFixed,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/n*, Si*', \

"< awk '{ if($2==100e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w lp smooth bezier lw 3 t 'Water*/SiO_2/Si*, t=200 nm', \


# 1E9*lambda/real(sqrt(epsilon(x,epsilonWater,meffWater,nuWater))) w l t '{/Symbol l}/n*, Water*'

# "< awk '{ if($2==2.5e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w lp lw 3 smooth bezier t 'Water*/SiO_2/Si*, t=5 nm'

# "< awk '{ if($2==20e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w lp smooth bezier lw 3 t 'Water*/SiO_2/Si*, t=40 nm', \
# "< awk '{ if($2==10e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w lp smooth bezier lw 3 t 'Water*/SiO_2/Si*, t=20 nm', \
# 1E9*lambda/real(sqrt(epsilon(neFixed,epsilonSi0,meffSi,nuSi))) w l t '{/Symbol l}/n, Si*', \

# "Experiment.dat" u 1:2:3 w yerrorbar t 'Experiment Bonse et al'

# "< awk '{ if($2==2.5e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w lp smooth bezier lw 3 t 'Water*/SiO_2/Si*, t=5 nm', \

1E9*lambda/real(sqrt(epsilonSiO2)) w l t '{/Symbol l}/n, SiO_2', \
(ConditionAbs( epsilon(x,epsilonSi0,meffSi,nuSi),epsilonSiO2)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonSiO2) )*1e9) : (1/0) w l lw 4 lc 7 t 'SiO_2/Si*', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilonSiO2)) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater),epsilonSiO2)*1e9)) : (1/0) w l lw 4 lc 7 t 'Water*/SiO_2', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilon(neFixed,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater),epsilon(neFixed,epsilonSi0,meffSi,nuSi) )*1e9)) : (1/0) w l lw 4 lc 7 t 'Water*/Si*', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilon(neFixed,epsilonSi0,meffSi,nuSi))) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater),epsilon(neFixed,epsilonSi0,meffSi,nuSi) )*1e9)) : (1/0) w l lw 4 lc 7 t 'Water*/Si*', \
(ConditionAbs( epsilon(x,epsilonSi0,meffSi,nuSi),epsilonSiO2)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonSiO2) )*1e9) : (1/0) w l lw 10 t 'SiO_2/Si*, restricted', \
(ConditionAbs( epsilon(x,epsilonSi0,meffSi,nuSi),epsilonSiO2)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonSiO2) )*1e9) : (1/0) w l lw 10 t 'SiO_2/Si*, restricted', \
"< awk '{ if($2==1e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=2 nm', \
"< awk '{ if($2==7e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=14 nm', \
"< awk '{ if($2==.10e-4) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=10 um', \
"< awk '{ if($2==.9999999942e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=2 um', \
"< awk '{ if($2==.9999999956e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=200 nm', \
"< awk '{ if($2==.9999999971e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=20 nm', \
"< awk '{ if($2==.9999999985e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/SiO_2/Si, t=2 nm', \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonSiO2) )*1e9) w l lw 2 t 'Water*/SiO_2', \

############ AIR / TiO2* / Ti ##################

reset
set output '20130731-SPPperiodTiExcitedTiO2.eps'
set terminal postscript eps enhanced monochrome font 'Times, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

set style increment userstyles
set key bottom right font 'Times, 16' spacing 0.6
set xlabel 'Free-carrier density N_e (m^{-3})'
set xrange [2E27:1e30]
set log x
set xtics format "10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
set log y
set yrange [10:2.5e3]

scale=1e6

ncrLocal=Ncr(2e0*pi*c/lambda,meffTiO2,nuTiO2,epsilonTiO2)
set arrow 1 from ncrLocal, 10 to ncrLocal, 2.5e3 nohead

plot (period( beta( epsilonAir, epsilon(x,epsilonTiO2,meffTiO2,nuTiO2)) )*1e9) w l lw 2 t 'air/TiO_2*', \
(ConditionAbs( epsilonAir,epsilon(x,epsilonTiO2,meffTiO2,nuTiO2))) ? (period( beta( epsilonAir, epsilon(x,epsilonTiO2,meffTiO2,nuTiO2)) )*1e9) : (1/0) w l lw 10 t 'air/TiO_2*, restricted', \
(period( beta( epsilon(x,epsilonTiO2,meffTiO2,nuTiO2), epsilonTi) )*1e9) w l lw 2 t 'TiO_2*/Ti', \
(ConditionAbs( epsilon(x,epsilonTiO2,meffTiO2,nuTiO2),epsilonTi)) ? (period( beta( epsilon(x,epsilonTiO2,meffTiO2,nuTiO2), epsilonTi) )*1e9) : (1/0) w l lw 10 t 'TiO_2*/Ti, restricted', \
(period( beta( epsilonAir, epsilonTi) )*1e9) w l lw 2 t 'Air/Ti', \
(ConditionAbs( epsilonAir,epsilonTi)) ? (period( beta( epsilonAir, epsilonTi) )*1e9) : (1/0) w l lw 10 t 'Air/Ti, restricted', \
"< awk '{ if($2==0.1e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Air/TiO_2*/Ti, t=100 nm', \
"< awk '{ if($2==20e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Air/TiO_2*/Ti, t=40 nm', \
"< awk '{ if($2==10e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Air/TiO_2*/Ti, t=20 nm', \
"< awk '{ if($2==7e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Air/TiO_2*/Ti, t=14 nm', \
"< awk '{ if($2==5e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Air/TiO_2*/Ti, t=10 nm', \
"< awk '{ if($2==1e-9) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Air/TiO_2*/Ti, t=2 nm', \
1e9*lambda/real(sqrt(epsilon(x,epsilonTiO2,meffTiO2,nuTiO2))) w l t '{/Symbol l}/n* TiO_2', \
"Experiment.dat" u 1:2:3 w yerrorbar lw 5 t 'Experiment Bonse et al'


############ Water* / TITANIUM oxide ##################

reset
set output '20130731-SPPperiodTiExcitedWater.eps'
set terminal postscript eps enhanced monochrome font 'Times, 26'

load localPath.'/SingleSPPfile.gnu'
load localPath.'/try2.plt'

set style increment userstyles
set key bottom right font 'Times, 16' spacing 0.6
set xlabel 'Free-carrier density N_e (m^{-3})'
set xrange [2E27:1e30]
set log x
set xtics format "10^{%L}" 10

set ylabel 'SPP period {/Symbol L} (nm)'
set log y
set yrange [10:2.5e3]

scale=1e6

ncrLocal=Ncr(2e0*pi*c/lambda,meffWater,nuWater,epsilonWater)
set arrow 1 from ncrLocal, 10 to ncrLocal, 2.5e3 nohead

plot (period( beta( epsilonTiO2, epsilonTi) )*1e9) w l lw 2 t 'TiO_2/Ti', \
(ConditionAbs( epsilonTiO2,epsilonTi)) ? (period( beta( epsilonTiO2, epsilonTi) )*1e9) : (1/0) w l lw 10 t 'TiO_2/Ti, restricted', \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonTiO2) )*1e9) w l lw 2 t 'Water*/TiO_2', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilonTiO2)) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonTiO2) )*1e9) : (1/0) w l lw 10 t 'Water*/TiO_2, restricted', \
(period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonTiO22) )*1e9) w l lw 2 t 'Water*/TiO_2 extra', \
(ConditionAbs( epsilon(x,epsilonWater,meffWater,nuWater),epsilonTiO22)) ? (period( beta( epsilon(x,epsilonWater,meffWater,nuWater), epsilonTiO22) )*1e9) : (1/0) w l lw 10 t 'Water*/TiO_2 extra, restricted', \
"< awk '{ if($2==.9999999942e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/TiO_2/Ti, t=2 um', \
"< awk '{ if($2==.9999999971e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/TiO_2/Ti, t=20 nm', \
"< awk '{ if($2==.6309573428e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/TiO_2/Ti, t=13 nm', \
"< awk '{ if($2==.3981071696e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/TiO_2/Ti, t=8 nm', \
"< awk '{ if($2==.2511886426e-8) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 3 t 'Water*/TiO_2/Ti, t=5 nm', \
1e9*lambda/real(sqrt(epsilon(x,epsilonWater,meffWater,nuWater))) w l t '{/Symbol l}/n*', \
"Experiment.dat" u 1:2:3 w yerrorbar lw 5 t 'Experiment Bonse et al'

# "< awk '{ if($2==.9999999956e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water*/TiO_2/Ti, t=200 nm', \
# "< awk '{ if($2==0.1e-5) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water*/TiO_2/Ti, t=2 um', \
# "< awk '{ if($2==0.1e-6) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water*/TiO_2/Ti, t=200 nm', \
# "< awk '{ if($2==0.50e-7) print }' Result.awked.tmp" u ($1):($3*1e9) w l lw 5 t 'Water*/TiO_2/Ti, t=20 nm', \


#
#(ConditionPerfect(epsilon(x,epsilonSi0, meffSi, nuSi),epsilonWater)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonWater) )*1e9) : (1/0) w l lw 2 t 'water/Si*, perfect', \
#(ConditionAbs(epsilon(x,epsilonSi0, meffSi, nuSi),epsilonWater)) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilonWater) )*1e9) : (1/0) w l lw 1 t 'water/Si*, absorbing', \
#(ConditionPerfect(epsilon(x,epsilonSi0, meffSi, nuSi),epsilon(x,epsilonWater,meffWater,nuWater))) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilon(x,epsilonWater,meffWater,nuWater)) )*1e9) : (1/0) w l lw 2 t 'water*/Si*, perfect', \
#(ConditionAbs(epsilon(x,epsilonSi0,meffSi,nuSi),epsilon(x,epsilonWater,meffWater,nuWater))) ? (period( beta( epsilon(x,epsilonSi0,meffSi,nuSi), epsilon(x,epsilonWater,meffWater,nuWater)) )*1e9) : (1/0) w l lw 1 t 'water*/Si*, absorbing'
#
#
