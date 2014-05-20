for fil in `cat scripts/doclist.txt`
do
#fil=./docs//pnmhd/2.txt
    echo $fil
    mv $fil $fil.temp
    echo "/*! \file $fil" > $fil.temp2
    echo "~~~~" > $fil.temp4
    echo "*/" > $fil.temp3
    cat $fil.temp2 $fil.temp4 $fil.temp $fil.temp4 $fil.temp3 > $fil
    rm -rf $fil.temp2 $fil.temp3 $fil.temp $fil.temp4
done
