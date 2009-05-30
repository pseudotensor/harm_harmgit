#
nx=128
ny=128
nz=1

igrid=0
inx=128
iny=128
inz=128
ixmax=100
iymax=100
izmax=100

# not necessary, will get from header, just fill up
iRin=4.84
iRout=100
iR0=4.37
ihslope=0.6552
idefcoord=0


# from fieldline files
bin2txt 1 2 0 -1  2 128 128 1 1 fieldline0100.bin fieldline0100.txt f 11
head -1 fieldline0100.txt > head.txt
#
tail +2 fieldline0100.txt > nohead.txt
awk '{print " "$9" "$10" "$11" "}'  nohead.txt > B0100.txt
rm -f nohead.txt
#
tail +2 dump0000 > nohead.txt
# 6+NPR+1+4*4+4+.  .=35
awk '{print " "$1" "$2" "}'  nohead.txt > grid1.txt
awk '{print " "$36" "}'  nohead.txt > grid2.txt

cp -f head.txt blob0100.txt
paste grid1.txt B0100.txt grid2.txt >> blob0100.txt

# need to print out header and then columns with
# ti tj B1 B2 B3 gdet
~/sm/smcalc 1 1 0 $nx $ny blob0100.txt aphi.txt

# now interpolate
~/sm/iinterp 1 2 1 1 $nx $ny $nz  1.0 0.0 0.5 $igrid  $inx $iny $inz  -$ixmax $ixmax -$iymax $iymax -$izmax $izmax  $iRin $iRout $iR0 $ihslope $idefcoord < aphi.txt > iaphi.txt
head -1 iaphi.txt > myhead.txt
tail +2 iaphi.txt > myiaphi.txt
