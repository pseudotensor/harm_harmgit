#!/bin/bash

xset r on
xset r on
xset r on
xset r on

END=255
for i in $(seq 0 $END)
do
    echo $i
    xset r $i
done




