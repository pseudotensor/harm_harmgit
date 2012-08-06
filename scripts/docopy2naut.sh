
# from kraken directly (thickdisk except thickdisk7) to ki-jmck
alias ls='ls'
cd /lustre/scratch/rblandfo
fildirs=dirsnewfullrest.txt
bdir=`pwd`/listsfull/
mkdir $bdir
#ls | egrep 'thickdisk*|rtf*' > $fildirs
alias lsdir='ls -lrt | egrep "^d"'
alias lssdir='ls -prt | grep / | sed "s/\///"'
alias lssdir2='ls -prt | grep / | tail -n +3 | sed "s/\///"'
#diit1=`lssdir | grep thickdisk3`
#diit2=`lssdir | grep thickdiskhr3`
#echo $diit1 > $fildirs
#echo $diit2 >> $fildirs
diit=`lssdir | egrep 'thickdiskr7'`
echo $diit > $fildirs


for fil in `cat $fildirs`
    do  
    echo $fil
    export dirorig=`pwd`
    cd $fil/dumps/
    find . \( -name "fieldline*.bin" -o -name "dump0000.bin" -o -name "rdump-0.bin" -o -name "gdump.bin" \) -print > $bdir/listcopy$fil.txt
    #
    rm -rf $bdir/listcopy$fil.globuslist.txt
    numfiles=0
    for realfile in `cat $bdir/listcopy$fil.txt`
    do
	veryrealfile=`basename $realfile`
	echo "xsede#kraken/lustre/scratch/rblandfo/$fil/dumps/$veryrealfile xsede#nautilus/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/$fil/dumps/$veryrealfile" >> $bdir/listcopy$fil.globuslist.txt
        # give permission for jmckinne to access files
	chmod a+x /lustre/scratch/rblandfo
	chmod a+rx /lustre/scratch/rblandfo/$fil
	chmod a+rx /lustre/scratch/rblandfo/$fil/dumps
	chmod a+rx /lustre/scratch/rblandfo/$fil/dumps/$veryrealfile
	numfiles=$(($numfiles+1))
    done
    chmod a+rx /lustre/scratch/rblandfo/$fil/nprlistinfo.dat
    chmod a+rx /lustre/scratch/rblandfo/$fil/coordparms.dat
    #
    #
    ssh pseudotensor@cli.globusonline.org mkdir xsede#nautilus:/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/$fil/
    ssh pseudotensor@cli.globusonline.org mkdir xsede#nautilus:/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/$fil/dumps/
    #
    #
    if [ -e /lustre/scratch/rblandfo/$fil/nprlistinfo.dat ]
    then
	ssh pseudotensor@cli.globusonline.org scp xsede#kraken:/lustre/scratch/rblandfo/$fil/nprlistinfo.dat xsede#nautilus:/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/$fil/
    fi
    if [ -e /lustre/scratch/rblandfo/$fil/coordparms.dat ]
	then
	ssh pseudotensor@cli.globusonline.org scp xsede#kraken:/lustre/scratch/rblandfo/$fil/coordparms.dat xsede#nautilus:/lustre/medusa/jmckinne/data1/jmckinne/jmckinne/$fil/
    fi
    #
    # perform transfer
    #
    if [ $numfiles -gt 0 ]
    then
	ssh pseudotensor@cli.globusonline.org transfer -s 1 < $bdir/listcopy$fil.globuslist.txt
    fi
    #
    #
    # now that done, take away permission for jmckinne to access files
    # no, above exits immediately, so can't remove access.  Also, other copies for all files to ranch may need access
    if [ 1 -eq 0 ]
    then
	for realfile in `cat $bdir/listcopy$fil.txt`
	do
	    veryrealfile=`basename $realfile`
	    chmod og-x /lustre/scratch/rblandfo
	    chmod og-rx /lustre/scratch/rblandfo/$fil
	    chmod og-rx /lustre/scratch/rblandfo/$fil/dumps
	    chmod og-rx /lustre/scratch/rblandfo/$fil/dumps/$veryrealfile
	done
	chmod og-rx /lustre/scratch/rblandfo/$fil/nprlistinfo.dat
	chmod og-rx /lustre/scratch/rblandfo/$fil/coordparms.dat
    fi
    #
    #
    #
    cd $dirorig
done


########################################################################

