#!/bin/bash

awk -f eVtoMicrons.awk l-Zr-Krishnan-eV-n.csv > l-Zr-Krishnan-um-n.csv
awk -f eVtoMicrons.awk l-Zr-Krishnan-eV-k.csv > l-Zr-Krishnan-um-k.csv

sort -n l-Zr-Krishnan-um-n.csv > l-Zr-Krishnan-n.csv
sort -n l-Zr-Krishnan-um-k.csv > l-Zr-Krishnan-k.csv

#rm Zr-Krishnan-um-n.csv Zr-Krishnan-um-k.csv