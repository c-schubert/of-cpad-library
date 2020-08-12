#!/bin/bash
filename="0_cpa.txt"

fromline=$(grep -in "^1.81228$" $filename | awk -F: '{print $1}')
fromline=$(($fromline-1))
linesoffile=$(wc -l $filename | awk '{print $1}')

toline=$linesoffile
newfilename="./testextract/${filename}"

echo "Extracting $fromline-$toline of file $filename to file $newfilename"

# echo "$fromline,${toline}p"
#sed -n "$fromline,${toline}p" $filename > $newfilename