#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: fixCsvFile.sh <Name of the CSV file to fix>"
	echo "Description: replaces comas <,> with tabs <tab>."
	exit -1
fi

mv $1 $1.coma
cat $1.coma | tr "," "\t" > $1.tab
mv $1.tab $1

echo "Initial file was saved in $1.coma."
