head -1 gdump.bin-col0000 > gdump.head

list=`ls gdump.bin-col????`

for fil in $list
do
    bin2txt 1 2 0 1 3 288 128 128 1 $fil $fil.txt d 1
done

list=`ls gdump.bin-col????.txt`

#paste  $list | column -s $' ' -t > gdump.txt.data

paste -d " " $list > gdump.txt.data

#pr -t -m firstfile secondfile

cat gdump.head gdump.txt.data > gdump.txt

bin2txt 2 1 0 -1 3 288 128 128 1 gdump.txt gdump.bin d 126





