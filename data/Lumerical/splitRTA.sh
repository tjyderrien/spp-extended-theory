#!/bin/bash
# Splits file by blocks. 
# Assumption: we have want to split in 3 files.
# Usage: ./split.sh <filename>

if [ -z $1 ]; then 
  echo -e "Splits file by equal blocks into 3 new files. \nUsage: ./split.sh <filename>"
else
  filename=$1
  NumLines=$( wc -l $filename | awk '{ print $1 }' )
  echo $NumLines
  NumLinesPerBlock=$( python2 -c "print ( ${NumLines} / 3 + 1 )"  ) 
  echo "NumLinesPerBlock=${NumLinesPerBlock}"
  split --suffix-length=1 --lines=${NumLinesPerBlock} --elide-empty-files --additional-suffix=.csv --numeric-suffixes=1 $filename $filename
fi
