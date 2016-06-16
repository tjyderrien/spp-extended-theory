#!/bin/bash

set -e

material=$1

if [ -z ${material} ]; then
	echo "Usage: ./importFromPalik.sh <Material symbol to include (Si, Be, ...) >"
	exit
else 
	echo "Ok let's add material $1 to our databasis."
fi

# exit

searchString () {
#   if [ -n $1 ]; then
  find | grep "$1" -i
#   else
#     exit
#   fi
}

echo "Searching for ${material} in available files..."
search=$( searchString "${material}")
searchnum=$( searchString "${material}" | wc -l )
echo "${searchnum} result(s) found."
echo "${search}" 

if [ ${searchnum} -ge 2  ]; then
  echo "** Search: Can you be more specific ?"
  echo "** Results"
  echo "${search}"
  exit
else
  searchMod=$( echo ${search} | tr "_" " " | tr "/" " " | tr '.' ' ' )
  material=$( echo ${searchMod} | awk '{ print $2 }' )
  searchPath=$( echo ${search} | tr "/" " " )
  datafile=$( echo ${searchPath} | awk '{ print $3 }' )
  source="Palik"
fi

filename="${material}-${source}"

echo "Build file ${filename}..."
echo "[Building] Headers..."
cat << EOF > ${filename}.1
# Optical constants for $1
# 'Handbook of Optical Constants of Solids', Ed. by Edward D. Palik,
# Academic Press, Inc., 1985.
# Lambda (A) n k
EOF

echo "[Building] Formatting data file..."
awk '{ print $1" "$2" "$3 }' Palik/${datafile} > ${filename}.2

echo "[Building] Stacking header and n,k file..."
cat ${filename}.1 ${filename}.2 > ${filename}

echo "[Cleaning]"
rm ${filename}.1 ${filename}.2

ls -lhtr ${filename}

echo "Done."
