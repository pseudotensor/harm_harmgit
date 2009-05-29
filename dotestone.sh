# setup code for a test
testnumber=$2
codename=$1
# assumes sauron.qsub inside unicode
cp -a $codename $codename.$testnumber
cd $codename.$testnumber
sed -e 's/#define TESTNUMBER 0/#define TESTNUMBER '$testnumber'/g'  init.h > init.h.new
mv init.h.new init.h
make cleanall
rm -rf grmhd
nice -n 10 nohup make >make.out 2>make.err
mkdir run
cp grmhd run
cd run
nice -n 10 nohup ./grmhd >grmhd.out 2>grmhd.err
export DUMPFILE=`ls dumps/dump* | tail -n 1` #get last dump filename

mkdir -p /mnt/data2/atchekho/autoanalysis
mkdir -p /mnt/data2/atchekho/autoanalysis/dumps
cp $DUMPFILE /mnt/data2/atchekho/autoanalysis/dumps/$codename.$testnumber
cp dumps/gdump /mnt/data2/atchekho/autoanalysis/dumps/$codename.$testnumber.gdump
#cp coordparms.dat /raid5/atchekho/autoanalysis/dumps/$codename.$testnumber.params

#sed -e 's/code.TESTNUMBER/'$codename.$testnumber'/g'  sauron.qsub > sauron.qsub.temp
#sed -e 's/TESTNUMBER/'$testnumber'/g'  sauron.qsub.temp > sauron.qsub.new
#qsub sauron.qsub.new >pbs.out 2>pbs.err
