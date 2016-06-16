#!/bin/bash

if [ -e $1 ]; then
	echo "Usage: ./importFromPalik.sh <Material symbol to include (Si, Be, ...) >"
	exit
else 
	echo "Ok let's add material $1 to our databasis."
fi

# exit

echo "Searching for $1 in available files..."

search () {
	if [ -e $1 ]; then
		find | grep "${research}"
	else
		echo "[Search] No argument was given."
		exit
	fi
}

searchnum=$( search "$1" | wc -l )
echo "${searchnum} result(s) found."

if [ ${searchnum} -ge 2  ]; then
	echo "** Search: Can you be more specific ?"
	echo "** Results"
	search "$1"
	exit
fi

echo "Creating headers..."
cat << EOF > ${1}_Palik.1
# Optical constants for $1
# 'Handbook of Optical Constants of Solids', Ed. by Edward D. Palik,
# Academic Press, Inc., 1985.
# Lambda (A) n k
EOF

echo "Stacking header and n,k file..."
cat "${1}_Palik.1" Palik/$1_palik.nk > ${1}_Palik

echo "Cleaning..."
rm ${1}_Palik.1

echo "Done."
