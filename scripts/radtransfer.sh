#!/bin/bash

alias lssdir='ls -apt | grep / | sed "s/\///"'

lssdir | grep rad > list.txt

for fil in `cat list.txt`
do
    echo $fil
    ssh pseudotensor@cli.globusonline.org scp -D -r -s 2 xsede#kraken:/lustre/scratch/jmckinne/$fil/ xsede#ranch:/home/01014/tg802609/
    sleep 360
done

