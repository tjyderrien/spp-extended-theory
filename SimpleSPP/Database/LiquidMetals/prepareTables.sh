#!/bin/bash

echo "Generate code-compatible data for this database."

file="Cu-liq-Miller"

if [ ! -e $file ]; then
	echo "** Error: file $file does not exist. "
	exit 1
fi

awk '{ print $1, $2 }' $file > $file-ReEps.csv
awk '{ print $3, $4 }' $file > $file-ImEps.csv
