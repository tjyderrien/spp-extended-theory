#!/bin/bash

set LC_ALL=C
prefix="Datas/"
rm *.tmp
# extraite la ligne 20 000 de chaque fichier et les reunir dans un seul
# 1: indice N, 2: indice K, 3: ligne 20000, 4: arctan(Z/X), 5: deg
cd $prefix
for i in `ls *.vtk` 
do
	echo $i | tr '_' ' ' | awk '{ print $2*1E0}' >> RealIndexes.tmp
	echo $i | tr '_' ' ' | awk '{ print $3 }' | awk 'BEGIN {FS="."} { print $1+0.1*$2 }' >> ImagIndexes.tmp ## | tr '_' ' ' | tr -d 'vtk' | awk '{ print $3}'
	cat $i | tail -n +12800 | head -n1 | awk '{ print $1,$2,$3}' >> RawDatas.tmp
done
cd -
cp $prefix/*.tmp .
paste RealIndexes.tmp ImagIndexes.tmp RawDatas.tmp > TreatedDatas.tmp
awk '{ printf("%10.5lf \t %lf \t %lf \t %lf \t %lf \n", $1, $2, $3, $4, $5) }' TreatedDatas.tmp | sort -n > TreatedDatas.dat
