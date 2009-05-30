#!/bin/bash
for fil in `ls . | grep grmhd-256x128-hor-.05-a-0`
do
  rsh metric "mkdir $fil" ; rcp $fil/*.fli metric:$fil
done
