##################
# shows how to use interpolation routine to take fieldline data and get back interpolation for each quantity desired.

#
# 1) make program

cd /lustre/ki/pfs/jmckinne/harmgit_jon2interp/
# setup for reduced code set
rm -rf init.c init.h
touch init.h


# 2) ensure PRINTHEADER and SCANHEADER in global.jon_interp.h are correct for older/newer simulations (i.e. THETAROT in new only)
# Do this by setting OLDERHEADER 1 if non-tilted runs.  Else set to 0.

# 3) make program itself (need Intel MKL -- modify makefile if path needs to be changed -- currently setup for ki-rh39)

make superclean ; make prepiinterp ; make iinterp &> make.log

# also make bin2txt program:

make superclean ; make prepbin2txt ; make bin2txt
# check makefile and setup for ki-rh39/orange/etc.

# ensure no errors during compile or link (need lapack!)

##############
4) copy programs to your path

cp iinterp ~/bin/iinterp.orange.thickdisk7
cp bin2txt ~/bin/bin2txt.orange

###############
# 5) do interpolation (directly read-in binary fieldline file and output full single file that contains interpolated data)

# 0=OLDER header with 30 entries (thickdisk/sasha sims)  1=NEWER header with 32 entries (tilted sims)
newheader=0

# get 3 times so can compute temporal derivative for (e.g.) current density at same spatial/temporal location as dump
dumpnum=5437
dumpnumm1=$(($dumpnum-1))
if [ -e dumps/fieldline$dumpnumm1.bin ]
then
    dumpnumm1=$(($dumpnum-1))
else
    dumpnumm1=$(($dumpnum))
fi
dumpnump1=$(($dumpnum+1))
if [ -e dumps/fieldline$dumpnump1.bin ]
then
    dumpnump1=$(($dumpnum+1))
else
    dumpnump1=$(($dumpnum))
fi	
#
# get times of dumps
time0=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $1}'`
timem1=`head -1 dumps/fieldline$dumpnumm1.bin |awk '{print $1}'`
timep1=`head -1 dumps/fieldline$dumpnump1.bin |awk '{print $1}'`
#
# get original resolution
nt=1
nx=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $2}'`
ny=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $3}'`
nz=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $4}'`
if [ $newheader -eq 0 ]
then
    numcolumns=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $30}'`
else
    numcolumns=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $32}'`
fi
R0=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $14}'`
Rin=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $15}'`
Rout=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $16}'`
hslope=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $17}'`
defcoord=`head -1 dumps/fieldline$dumpnum.bin |awk '{print $19}'`
#
# Note that iinterp has x->xc y->zc z->yc since originally was doing 2D in x-z
# That is:
# So box(n)x refers to true x
# So box(n)y refers to true z
# So box(n)z refers to true y
#
# True x is V5D's x and iinterp's x
# True z is V5D's y and iinterp's y
# True y is V5D's z and iinterp's z
#
#
# Vectors are still in order as columns as vx, vy, vz for true x,y,z, respectively
#
boxnt=1
boxnx=100
boxny=100
boxnz=100
boxxl=-10
boxyl=-10
boxzl=-10
boxxh=10
boxyh=10
boxzh=10

# set docurrent=0 if want quick result with no current density
# this will change number of output columns
docurrent=1

cd /lustre/ki/pfs/jmckinne/thickdisk7/
IDUMPDIR=/lustre/ki/pfs/jmckinne/thickdisk7/idumps/
# ensure coordparms.dat exists here -- required to read in harm internal grid parameters
mkdir $IDUMPDIR
#
#
#
whichoutput=14
outfilename=$IDUMPDIR/fieldline$dumpnum.cart.bin.$boxnx.$boxny.$boxnz
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN $nt $nx $ny $nz -refine 1.0 -filter 0 -grids 1 0 -nN $boxnt $boxnx $boxny $boxnz -ibox $time0 $time0 $boxxl $boxxh $boxyl $boxyh $boxzl $boxzh -coord $Rin $Rout $R0 $hslope -defcoord $defcoord -dofull2pi 1 -docurrent $docurrent -tdata $timem1 $timep1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double -infile dumps/fieldline$dumpnum.bin -infilem1 dumps/fieldline$dumpnumm1.bin -infilep1 dumps/fieldline$dumpnump1.bin -outfile $outfilename


# as a test, one can do just 1 variable (the density)
if [ 1 -eq 0 ]
then
    outfilename=$IDUMPDIR/fieldline$dumpnum.cart.bin.densityonly.$boxnx.$boxny.$boxnz
    ~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype 1 -itype 1 -head 1 1 -oN $nt $nx $ny $nz -refine 1.0 -filter 0 -grids 1 0 -nN $boxnt $boxnx $boxny $boxnz -ibox $time0 $time0 $boxxl $boxxh $boxyl $boxyh $boxzl $boxzh -coord $Rin $Rout $R0 $hslope -defcoord $defcoord -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -infile dumps/fieldline$dumpnum.bin -outfile $outfilename
fi


# The whichoutput==14 case results in a file with 1 line text header, line break, then data
# The binary data block of *floats* (4 bytes) is ordered as:
# fastest index: column or quantity list
# next fastest index: i (associated with true x)
# next fastest index: j (associated with true y)
# next fastest index: k (associated with true z)
# slowest index: h (time, but only single time here)

# Note that while array access of quantity is out of order in iinterp result, the 3D spatial grid is still written as a right-handed coordinate system.
# So, if you face the screen showing positive z-axis pointing up (increasing j) and positive x-axis pointing to the right (increasing i), then the positive y-axis points into the screen (increasing k).

# The columns or quantities in the list are ordered as:
# rho0,ug,vx,vy,vz,Bx,By,Bz,FEMrad,Bphi,Jt,Jx,Jy,Jz  (J's only exist if -docurrent 1 was set)
# FEMrad is the radial energy flux
# Bphi is the poloidal enclosed current density
# J\\mu=J^\\mu is the current density.
#
#
# Derived quantities:
#
# 3-velocity magnitude: v^2=vx^2+vy^2+vz^2 (i.e. just spatials are squared)
# Lorentz factor: ut = 1/sqrt(1-v^2)
# 4-velocity: u={ut,ut*vx,ut*vy,ut*vz)
# Lab current squared: J^2=J.J = -Jt*Jt + Jx*Jx + Jy*Jy + Jz*Jz (full space-time square)
# Comoving current density: j^\\nu = J^\\mu h_\\mu^\\nu = J^\\mu (\\delta_\\mu^\\nu + u_\\mu u^\\nu) = J^\\nu + (J^\\mu u_\\mu)u^\\nu
# Comoving square current density: j^2 = j^\\nu j_\\nu = J^2 + (J.u)^2 = J^2 + ut*(Jt + Jx*vx + Jy*vy + Jz*vz)
# field along flow: u.B = ux*Bx + uy*By + uz*Bz (i.e. Bt=0)
# comoving magnetic field: b^\mu = (B^\mu + (u.B)u^\mu)/ut
# comoving mag energy: b^2/2 = 0.5*(B^2 + (u.B)^2)/ut^2



# can check how looks in text by doing:

file=$outfilename
if [ $newheader -eq 0 ]
then
    numoutputcols=`head -1 $file  |awk '{print $30}'`
else
    numoutputcols=`head -1 $fil |awk '{print $32}'`
fi
newnx=`head -1 $file |awk '{print $2}'`
newny=`head -1 $file |awk '{print $3}'`
newnz=`head -1 $file |awk '{print $4}'`
bin2txt.orange 1 2 0 -1 3 $newnx $newny $newnz 1 $file $file.txt f $numoutputcols
less -S $file.txt

# should look reasonable.







#
################################### DONE!


















############
# bad way to create 3-time file (too slow for read-in to seek all over the place)


# another beginning of slow way:
# setup new compact 3-time file read-in as 1-time file
bin2txt.orange 1 2 0 -1 3 $nx $ny $nz 1 dumps/fieldline$dumpnum.bin $IDUMPDIR/fieldline$dumpnum.txt f $numcolumns
bin2txt.orange 1 2 0 -1 3 $nx $ny $nz 1 dumps/fieldline$dumpnumm1.bin $IDUMPDIR/fieldline$dumpnumm1.txt f $numcolumns
bin2txt.orange 1 2 0 -1 3 $nx $ny $nz 1 dumps/fieldline$dumpnump1.bin $IDUMPDIR/fieldline$dumpnump1.txt f $numcolumns
#... using "join" or "pr" or similar -- but too slow to create text files and pointless.




# setup new 3-time file:
if [ 1 -eq 0 ]
then
    head -1 ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.bin.head
    tail -n +2 ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.bin.data
    tail -n +2 ./dumps/fieldline$dumpnumm1.bin > $IDUMPDIR/fieldline$dumpnumm1.bin.data
    tail -n +2 ./dumps/fieldline$dumpnump1.bin > $IDUMPDIR/fieldline$dumpnump1.bin.data
    cat $IDUMPDIR/fieldline$dumpnum.bin.head $IDUMPDIR/fieldline$dumpnumm1.bin.data $IDUMPDIR/fieldline$dumpnum.bin.data $IDUMPDIR/fieldline$dumpnump1.bin.data > $IDUMPDIR/fieldline$dumpnum.bin.3time
fi



###########
# NOTES:

# note that files needed by iinterp and bin2txt are only:

CFILES1="jon_interp.c jon_interp_computepreprocess.c jon_interp_filter.c jon_interp_interpolationitself.c jon_interp_mnewt.c jon_interp_newt.c nrutil2.c coord.c nrutil.c ludcmp.c lubksb.c lnsrch.c fdjac.c fmin.c broydn.c qrdcmp.c rsolv.c qrupdt.c rotate.c ranc.c bcucof.c bcuint.c metric.tools.c mpi_fprintfs.c gaussj.c tetrad.c"
CFILES2="tensor.c smcalc.c generatenprs.c bin2txt.c initbase.defaultnprlists.c"

CFILES="$CFILES1 $CFILES2"

# create grep -v command:
GREPCOMMAND=""
for fil in $CFILES ; do GREPCOMMAND=$GREPCOMMAND"|grep -v $fil " ; done


ls -art *.c <PASTE echo $GREPCOMMAND here>




#######################################
# OLDER THINGS -- don't use -- just use above if reading in field lines

# interpolate to get rest-mass density
whichoutput=1000
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.rho.bin
#
# interpolate to get thermal gas internal energy density
whichoutput=1001
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.ug.bin
#
#
# interpolate to get orthonormal v^x
whichoutput=1003
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.vx.bin
#
# interpolate to get orthonormal v^y
whichoutput=1004
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.vy.bin
#
# interpolate to get orthonormal v^z
whichoutput=1005
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.vz.bin
#
# interpolate to get orthonormal B^x
whichoutput=1007
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.Bx.bin
#
# interpolate to get orthonormal B^y
whichoutput=1008
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.By.bin
#
# interpolate to get orthonormal B^z
whichoutput=1009
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.Bz.bin
#
# interpolate to get radial energy flux per unit sin(\theta)
whichoutput=1011
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.FE.bin
#
# interpolate to get poloidal current density (B_\phi)
whichoutput=1012
~/bin/iinterp.orange.thickdisk7 -binaryinput 1 -binaryoutput 1 -inFTYPE float -outFTYPE float -dtype $whichoutput -itype 1 -head 1 1 -oN 1 272 128 256 -refine 1.0 -filter 0 -grids 1 0 -nN 1 $boxnx $boxny $boxnz -ibox 0 0 -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvaluetype 0 -gdump ./dumps/gdump.bin -gdumphead 1 1 -binaryinputgdump 1 -inFTYPEgdump double < ./dumps/fieldline$dumpnum.bin > $IDUMPDIR/fieldline$dumpnum.cart.Bphi.bin

# merge interpolation results




#####################################
# NOTES, OLDER THINGS -- DONT USE!




# MTB12 simulations have the following headers in fieldline and gdump files:
tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2], dx[3], localrealnstep,gam,a,R0,Rin,Rout,hslope,localdt,defcoord,MBH,QBH,EP3,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns


# new tilted simulation header has extra THETAROT
tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2], dx[3], localrealnstep,gam,a,R0,Rin,Rout,hslope,localdt,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns);



############ SKIP THIS SECTION AND GO TO IINTERP SECTION
# convert fieldline????.bin into text format.
cd /lustre/ki/pfs/jmckinne/harmgit_jon2interp/
make superclean ; make prepbin2txt ; make bin2txt

~/bin/bin2txt.orange  1 2 0 -1 3 272 128 256 1 fieldline5437.bin fieldline5437.txt f 11

# now break into separate columns
head -1 fieldline5437.txt > fieldline5437.txt.head
tail -n +2 fieldline5437.txt > fieldline5437.txt.data
awk '{print $1}' fieldline5437.txt.data > fieldline5437.txt.data.rho0
awk '{print $2}' fieldline5437.txt.data > fieldline5437.txt.data.ug

# get u^i from u^t and v^i = u^i/u^t
awk '{print $4}' fieldline5437.txt.data > fieldline5437.txt.data.uux0
awk '{print $5*$4}' fieldline5437.txt.data > fieldline5437.txt.data.uux1
awk '{print $6*$4}' fieldline5437.txt.data > fieldline5437.txt.data.uux2
awk '{print $7*$4}' fieldline5437.txt.data > fieldline5437.txt.data.uux3

################

# -coord <Rin> <Rout> <R0> <hslope>
#~/bin/iinterp.orangenew -binaryinput 0 -binaryoutput 0 -inFTYPE float -outFTYPE float -dtype 1 -itype 1 -head 1 1 -oN 272 128 256 -refine 1.0 -filter 0 0 -grids 1 0 -nN $boxnx $boxny $boxnz -ibox -$boxx $boxx -$boxy $boxy -$boxz $boxz -coord 1.15256306940633 26000 0 1.04 -defcoord 1401 -dofull2pi 1 -extrap 1 -defaultvalue 0 -gdump ./gdump.txt < fieldline5437.rho.txt > fieldline5437.cart.rho.txt
