#!/bin/bash
source ~/.bashrc
for fil in `ls *.c`
do
    ./extractnonfunc < $fil | grep -v "^void " | grep -v "^//"| grep -v "^/\*" | grep -v "^#" | grep -v "FUNCBEGIN"
done
