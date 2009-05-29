#!/bin/bash
# This generates shared list for OpenMP
# This shared list must at least include all global variables accessed (those that weren't passed to be local).
# Shared should also include any local function variables common to loop, which should be done by user as appended with a comma after OPENMPSHAREDLIST

echo "0"
cpp $1 > cpp.openmp.temptemptemp.1
beginline=`grep -n BEGINOPENMPSHAREDLIST cpp.openmp.temptemptemp.1 | sed 's/:.*//g'`
truebeginline=$(($beginline+1))
tail -n +$truebeginline cpp.openmp.temptemptemp.1 > cpp.openmp.temptemptemp.2
endline=`grep -n ENDOPENMPSHAREDLIST cpp.openmp.temptemptemp.2 | sed 's/:.*//g'`
trueendline=$(($endline-1))
head -$trueendline cpp.openmp.temptemptemp.2 > cpp.openmp.temptemptemp.3
echo "1"
sed -e 's/(FTYPE.*);//g' -e 's/(double.*);//g'  -e 's/(int.*);//g'  cpp.openmp.temptemptemp.3 >  cpp.openmp.temptemptemp.3.5
echo "2"
sed -e 's/\[[\ _a-ZA-Z0-9\+\-\/\*?><:]*\]//g' -e 's/\*//g'  -e 's/{//g' -e 's/(//g' -e 's/)//g' -e 's/;//g' cpp.openmp.temptemptemp.3.5 > cpp.openmp.temptemptemp.4
echo "3"
sed -e 's/double //g' -e 's/signed char //g'  -e 's/char //g'  -e 's/int //g'  -e 's/float //g' -e 's/FILE //g' -e 's/long //g'  -e 's/void //g' -e 's/struct Ccoordparams //g' -e 's/struct blink //g' -e 's/SFTYPE //g' -e 's/FTYPE //g' -e 's/struct of_geom //g' -e 's/struct of_state //g' -e 's/=&geomdontuse//g'  cpp.openmp.temptemptemp.4 >  cpp.openmp.temptemptemp.5
echo "4"
sed -e '/^$/d' -e '/^#/d' cpp.openmp.temptemptemp.5 > cpp.openmp.temptemptemp.6
echo "5"
cat cpp.openmp.temptemptemp.6 | tr '\n' ',' > cpp.openmp.temptemptemp.7
echo "6"
sed -e 's/,,/,/g'  cpp.openmp.temptemptemp.7 >  cpp.openmp.temptemptemp.8
echo "7"
sed -e 's/\[[\ _a-ZA-Z0-9\+\-\/\*?><:]*\]//g' cpp.openmp.temptemptemp.8 >  $1.joncppoutput

rm -rf cpp.openmp.temptemptemp.1 cpp.openmp.temptemptemp.2 cpp.openmp.temptemptemp.3 cpp.openmp.temptemptemp.3.5 cpp.openmp.temptemptemp.4 cpp.openmp.temptemptemp.5 cpp.openmp.temptemptemp.6 cpp.openmp.temptemptemp.7 cpp.openmp.temptemptemp.8

echo "Now for main.c, create global.openmpsharedlist.h and put #define OPENMPSHAREDLIST in front of the result of this script and remove the comma at the end"
