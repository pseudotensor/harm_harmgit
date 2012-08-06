
# from kraken directly (thickdisk except thickdisk7) to ki-jmck
alias ls='ls'
cd /lustre/scratch/jmckinne
fildirs=dirsnewfull.txt
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
diit=`lssdir | egrep 'thickdisk*'`
echo $diit >> $fildirs


for fil in `cat $fildirs`
    do  
    echo $fil
    export dirorig=`pwd`
    cd $fil/dumps/
    find . \( -name "fieldline*.bin" -o -name "dump0000.bin" -o -name "rdump-0.bin" -o -name "gdump.bin" \) -print > $bdir/listcopy$fil.txt
    ssh jmckinne@ki-jmck.slac.stanford.edu "mkdir -p /data2/jmckinne/$fil/dumps/"
    scp ../npr* ../coord* jmckinne@ki-jmck.slac.stanford.edu:/data2/jmckinne/$fil/
    # -I $bdir/listcopy$fil.txt
    filestocopy=`ls dump0000.bin* gdump.bin* rdump-0.bin* fieldline*.bin`
    #filestocopy="rdump-0.bin fieldline*.bin"
    ~/bin/bbcp -a -k -f -p -r -z -s 20 -b +500 -l $dirorig/stderr_direct.$fil   -P 5 -V -T 'ssh -x -a -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' -S 'ssh -x -a -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' $filestocopy jmckinne@ki-jmck.slac.stanford.edu:/data2/jmckinne/$fil/dumps/
    cd $dirorig
done


########################################################################

