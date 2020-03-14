#!/bin/bash
echo "converting images to gif"
cd WWURay
ls -v video >> $1.txt
cd video
convert -delay xx @../$1.txt ../results/$1.gif
rm s*
cd ..
rm $1.txt
echo "conversion done"

