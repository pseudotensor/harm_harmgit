#!/bin/bash
#source ~/.bashrc
#for fil in ls *.c *.h makefile; do echo $fil; mydiff $fil $1/$fil; done
for fil in ls *.c *.h makefile; do echo $fil; diff -bBdp $fil $1/$fil; done
