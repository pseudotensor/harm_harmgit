#!/bin/bash
for i in 1 2 3 4 7 8 9
do
    echo "Doing test #${i}..."
    cd ~/Research/code/harm
    dirname=tests/komi${i}
    mkdir -p tests
    mkdir -p $dirname
    cat initboundcode/init.koral.h | sed "s/WHICHKOMI [0-9]*/WHICHKOMI $i/g" > $dirname/init.koral.h
    cp $dirname/init.koral.h initboundcode/init.koral.h
    echo "Compiling test #${i}..."
    make superclean &> $dirname/compile_log.txt
    make prep &> $dirname/compile_log.txt
    make -j 4 &> $dirname/compile_log.txt
    cp grmhd $dirname
    cd $dirname
    echo "Running test #${i}..."
    ./grmhd 1 1 1  &> run_log.txt
done