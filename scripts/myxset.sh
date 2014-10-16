#!/bin/bash

xset r on
xset r on
xset r on
xset r on

END=255
for i in $(seq 0 $END)
do
    echo $i
    xset -r $i
done

# xev tells you which keycodes for each key for a given keyboard
xset r 111
xset r 113
xset r 114
xset r 116



