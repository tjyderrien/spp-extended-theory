#!/bin/bash

# Gather datas from files
# 

rm *.tmp

cat NeList.dat | tr "\t" "\n" > NeList.tmp
cat aList.dat | tr "\t" "\n"  > aList.tmp
cat Period.dat | tr "\t" "\n" > Period.tmp
cat DielectricTiO2.dat | tr "\t" "\n" > Epsilon1.tmp
cat SPPdepth.dat | tr "\t" "\n" > Damping.tmp
cat beta.dat | tr "\t" "\n" > beta.tmp
# cat SPPdecayDepth.dat | tr "\t" "\n" > Depth.tmp

TotalM=`wc -l NeList.tmp | awk '{ print  $1}' ` #30
TotalN=`wc -l aList.tmp | awk '{ print  $1}' ` #30

indexI=0; 
# pour chaque ligne de NeList.dat, 
for i in `cat NeList.tmp` 
do
	indexJ=0
#	epsilon=`head -n$indexI Epsilon1.tmp | tail -n1`
	let "indexI=indexI+1"	
	epsilon=`head -n$indexI Epsilon1.tmp | tail -n1`
	# pour chaque ligne de aList.dat 
	for j in `cat aList.tmp`
	do
		let "indexJ=indexJ+1"
		let "k=(indexI-1)*TotalN+indexJ"
#		let "tailSize=TotalM*TotalN-k"
		result=`head -n$k Period.tmp | tail -n1 `
		epsilon=`head -n$indexI Epsilon1.tmp | tail -n1 `
		dampingdepth=`head -n$indexI Damping.tmp | tail -n1 `
		beta=`head -n$k beta.tmp | tail -n1 `
		#  concatener NeList(i), aList(j), Results(i,j)
#		echo $indexI, $indexJ, $k
		echo $i $j $result $epsilon >> Result.tmp
		echo $i $j $dampingdepth >> DampingLength.tmp
		echo $i $j $beta >> betaShaped.tmp
	done
done

cat betaShaped.tmp | tr '+' ' ' > betaShaped2.tmp
cat betaShaped2.tmp | tr '*I' ' ' > betaShaped3.tmp

# awk '{ if($3<10e-6) print; else {print $1,$2,"?", "?"}}' Result.tmp | awk -f '/home/thibault/Documents/Codes/Scripts/pm3d.awk' > Result.awked.tmp
awk '{ if($3<1e10) print; else {print $1,$2,"?", "?"}}' Result.tmp | awk -f '/home/thibault/Documents/LaAPT/Scripts/pm3d.awk' > Result.awked.tmp
awk '{ if($3<1e10) print; else {print $1,$2,"?", "?"}}' DampingLength.tmp | awk -f '/home/thibault/Documents/LaAPT/Scripts/pm3d.awk' > Damping.awked.tmp
awk '{ if($3<1e10) print; else {print $1,$2,"?", "?"}}' betaShaped3.tmp > Beta.awked.tmp
# | awk -f '/home/thibault/Documents/LaAPT/Scripts/pm3d.awk' > Beta.awked.tmp
