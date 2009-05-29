NUMTESTS=$3
a=$2
while [ $a -le "$NUMTESTS" ]
  do 
  (sh dotestone.sh $1 $a)
  a=$(($a+1))
done
