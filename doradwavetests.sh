#!/bin/bash
name="radwave"
for i in 1 10 11 101 102 103 104 1001 1002 1003 1004 1101 1102 1103 1104
do
    echo "Doing test #${i}..."
    for n in 8 16 32 64 128 256 512
    do
        cd ~/code/harm
        dirname=tests/$name${i}_${n}
        mkdir -p tests
        mkdir -p $dirname
	#choose test
	cat initboundcode/init.koral.c | sed "s/^\([ \t]*\)RADWAVE_NUMERO=[0-9]*\;/\1RADWAVE_NUMERO=${i}\;/g" > $dirname/init.koral.c
        cp $dirname/init.koral.c initboundcode/init.koral.c
	#choose resolution
        cat initboundcode/init.koral.h | sed "s/N1 [0-9 ]*\/\/RADWAVE/N1 $n   \/\/RADWAVE/g" > $dirname/init.koral.h
        cp $dirname/init.koral.h initboundcode/init.koral.h
        echo "Compiling test #${i} using ${n} cells: $name${i}_${n} ..."
        make superclean &> $dirname/compile_log.txt
        make prep &> $dirname/compile_log.txt
        make -j 16 &> $dirname/compile_log.txt
        cp grmhd $dirname
        cd $dirname
        echo "Running (in the background) test #${i} using ${n} cells: $name${i}_${n} ..."
        ./grmhd 1 1 1  &> run_log.txt &
    done
done
echo "Waiting for all background processes to finish..."
#wait for all tests at different resolutions to complete
wait 
echo "Done!"