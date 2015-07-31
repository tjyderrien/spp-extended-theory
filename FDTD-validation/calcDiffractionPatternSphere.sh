#!/bin/bash
rm Diffraction.dat
echo "# theta, phi, Ix, Iy, Iz" >> Diffraction.dat
for i in `ls xnff* -1` 
do
	# collect the theta, phi coordinates from filename
	theta=`echo $i | tr "_" " " | awk '{ print $3}' `
	phi=`echo $i | tr "_" " " | awk '{ print $4}' `
#	echo theta=$theta, phi=$phi
	
	# calculate the total field for each coordinate
	value=`awk -f '../sumLines.awk' $i`
#	echo $value
	
	#export to a file
	echo "$theta $phi $value i" >> Diffraction.dat
done
