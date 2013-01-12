#!/bin/bash

# run using "bash" not "sh"

# then do (and look for non-standard warnings)
# grep -i "warning" make.log | grep -v "unrecognized" | grep -v "subscript out of range" | grep -v "may lose significant bits" |less -S

####################
# get list of extern declarations to add to global.h to check declaration consistency

# 1) EITHER (for all c files -- even unused): cat *.c *.h > temptempallc.txt
# 1) OR: cat $totallist > temptempallc.txt

# auto-extraction( requires #END at end of thing to be put into export)
sed 's/\\//g' maketail.harm.inc | sed -e "s/\=/\=\'/g" | sed -e "s/\#END/\'/g" | sed -e 's/OBJS/export OBJS/g' | sed -e 's/SRCS/export SRCS/g' >  temptempmyexport.sh
source temptempmyexport.sh

# get all code in one file
export totallist=`echo $SRCSMANYPOINTS $SRCSONEPOINT`
cat $totallist > temptempallc.txt

sed 's/\/\/.*//g' temptempallc.txt > temptempallc2.txt # remove C 1-line comments
sed ':a;N;$!ba;s/\n/ /g' temptempallc2.txt > temptempallc3.txt # replace new lines with space
# to force simple interpretation of ;'s as indicating lines (so functions not broken on multiple lines)
sed 's/;/;\n/g' temptempallc3.txt > temptempallc35.txt  # replace ';' with ';<CTRL-ALT-J>'
sed 's/#define/\n#define/g' temptempallc35.txt > temptempallc4.txt  # replace '#define' with '<CTRL-ALT-J>#define'
# Replace all 'extern' with '<CTRL-ALT-J>extern ' to isolate lines that start with extern
sed 's/extern/\n extern /g' temptempallc4.txt > temptempallc5.txt
# get extern declarations
grep -e '.*\?extern .*\?(.*\?) *\?;' temptempallc5.txt > temptempfinalallc.txt


# 2) uncomment bottom #include from global.h
# 3) make and see if any conflicts:

# superclean so redo dependencies that depend upon filename: temptempfinalallc.txt
make superclean ; make prep ;  make grmhd &> make.log

#
# clean up
#
# CAREFUL DON'T CHANGE LINE or add star
rm -rf temptempall.txt temptempallc5.txt temptempallc4.txt temptempallc3.txt temptempallc2.txt temptempallc.txt
# didn't remove temptempfinalallc.txt  so user can look at what was included
#
#
#
#
#
#




