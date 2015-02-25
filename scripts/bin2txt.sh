cd dumps

dnum=0000
#dnum=0023
#dnum=0113
#dnum=0237

bin2txt 1 2 0 -1 2 128 64 1 1 gdump.bin gdump d 126
bin2txt 1 2 0 -1 2 128 64 1 1 dump${dnum}.bin dump${dnum} d 61

bin2txt 1 2 0 -1 2 128 64 1 1 debug${dnum}.bin debug${dnum} d 288
ln -s debug${dnum} debugdump${dnum}
bin2txt 1 2 0 -1 2 128 64 1 1 eosdump${dnum}.bin eosdump${dnum} d 54
bin2txt 1 2 0 -1 2 128 64 1 1 failfloordudump${dnum}.bin failfloordudump${dnum} d 13
bin2txt 1 2 0 -1 2 128 64 1 1 raddump${dnum}.bin raddump${dnum} d 44
bin2txt 1 2 0 -1 2 128 64 1 1 vpotdump${dnum}.bin vpotdump${dnum} d 4

cd ..
