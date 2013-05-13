#!/bin/bash
alias mydiff='diff -bBdp'
#source ~/.bashrc
#for fil in `ls -d *.c *.h`; do echo $fil; diff -bBdpy -W 600 --suppress-common-lines $fil $1/$fil; done
for fil in `ls initboundcode/*.c initboundcode/*.h *.c *.h make*`; do echo $fil; mydiff $fil $1/$fil; done
#for fil in `ls *.c *.h` makefile; do echo $fil; diff -q $fil $1/$fil; done
