###################
#
# see makefile.notes for detailed comments/notes
#
###################


#####################
#
# PICK SYSTEM TO COMPILE ON
#
#####################
# WATCH OUT FOR SPACES, etc. AFTER ASSIGNMENTS!!!!
USEBG=0
#QUEENBEE (avoids making directories or forking in general)
USEQB=0
USEUB=0
USEOSX=0
# USEABE -> USEICCGENERIC is used
USEICCGENERIC=0
USEICCINTEL=0
USEGCC=0
USECCC=0
USEORANGE=0
USENERSC=0
USEKRAKEN=0
USEKRAKENICC=0
USETACCLONESTAR4=0
USETACCLONESTAR=0
USETACCRANGER=0
USEPFE=0
USEPFESGIMPT=0
USEORB=1
USENAU=0

#####################
#
# GSL CHOICE (needed for limotorus)
#
#####################
USEGSL=0


########################
#
# PERFORMANCE CHECKS WITH GPROF
#
########################
# default should be 0 since has large performance hit.
USEGPROF=0


########################
#
# Debug information in binary
#
########################
# Generally ok to have, no performance hit
# Required for using pfmon with intelligent output
USEDEBUG=1



########################
#
# BEGIN setup what compiling
#
########################


########################
#
# PREPARE COMMAND NAMES
#
########################
PREP=prep
PREPLD=prepgrmhdldouble
PREPLIAISON=prepliaison
PREPJONINTERP=prepiinterp
PREPBIN2TXT=prepbin2txt
PREPSMCALC=prepsmcalc

########################
#
# BINARY NAMES
#
########################
DEFTARGET=
CMD=grmhd
CMD2=grmhdldouble
LIAISONCMD=liaison
JONINTERPCMD=iinterp
BIN2TXTCMD=bin2txt
SMCALCCMD=smcalc


# different target sets
MPITARGETS=$(DEFTARGET) $(CMD) $(CMD2) $(LIAISONCMD) $(PREP) $(PREPLD) $(PREPLAISON)
LAPACKTARGETS=$(DEFTARGET) $(CMD) $(CMD2) $(LIAISONCMD) $(JONINTERPCMD) $(PREP) $(PREPLD) $(PREPLAISON) $(PREPJONINTERP)
BINARYTARGETS=$(CMD) $(CMD2) $(LIAISONCMD)

########################
#
# Set default goal
#
########################
.DEFAULT_GOAL := $(CMD)

$(warning default goal is $(.DEFAULT_GOAL))
$(warning goal is $(MAKECMDGOALS))
$(warning MPITARGETS is $(MPITARGETS))

# check for default make call
ifeq ($(MAKECMDGOALS),)
MYMAKECMDGOALS=$(.DEFAULT_GOAL)
else
MYMAKECMDGOALS=$(MAKECMDGOALS)
endif



########################
#
# DONE setup what compiling
#
########################



########################
#
# MPI CHOICE
#
########################

#ifeq ($(MYMAKECMDGOALS),$(filter $(MYMAKECMDGOALS), $(MPITARGETS)))
ifneq ($(filter $(MPITARGETS), $(MYMAKECMDGOALS)),)
# choice for below
# i.e. MPI is only for grmhd, grmhdldouble, and liaison codes
USEMPI=0

else

# no choice for below
# i.e. no MPI for iinterp, bin2txt, or smcalc
USEMPI=0

endif




########################
#
# OpenMP CHOICE
#
########################
USEOPENMP=0



# turns on -mp and -pc64 for USEICC?????=1
# Note that if using -pc64 -mp that error in inversion seems to be limited for doubles to 1E-11 instead of 1E-15
# causes code to be less accurate: e.g. inversion doesn't reach to 1E-15, but rather 1E-11
# for some simple linear wave tests with small amplitude can cause significant problems if any  noise exists in intial solution
# With OpenMP and icc 10.1 20080602 this feature with:
#      -O2 -finline -finline-functions -ip -fno-alias -unroll -Wall -Wcheck -Wshadow -w2 -wd=1419,869,177,310,593,810,981,1418 $(DFLAGS) -g
# caused the fluxcalc_fluxctstag_emf_1d() function fail to work properly.  It would cause NaN's to random appear in different places. I tried to find any associated bug, but could not.  Reverting to -O0 leads to no problems.
ENFORCEPRECISION=0


###################
# LAPACK
# can easily see if MKL already setup by doing:
# echo ${MKLROOT}
# in bash and see if defined
#
# otherwise must do "source" command below before compiling with lapack:
# on 64-bit emulated machines (most intel/amd machines):
# source /opt/intel/mkl/10.0.3.020/tools/environment/mklvarsem64t.sh
# For true 64-bit machines and 32-bit machines different source
# true 64-bit:
# source /opt/intel/mkl/10.0.3.020/tools/environment/mklvars64.sh
# true 32-bit:
# source /opt/intel/mkl/10.0.3.020/tools/environment/mklvars32.sh
# and recall the directory will be different on each machine
#
# on ki-rh39 and orange:
# source /u/ki/jmckinne/intel/mkl/10.0.5.025/tools/environment/mklvarsem64t.sh
# or:
# source /u/ki/jmckinne/intel/mkl/10.0.3.020/tools/environment/mklvarsem64t.sh
#
# also see:
# http://www.netlib.org/eispack/
# http://www.netlib.org/lapack/double/dsyev.f
#
ifneq ($(filter $(LAPACKTARGETS), $(MYMAKECMDGOALS)),)
# choice for below
USELAPACK=1
# below is for new form of library -l line
USELAPACKNEW=1
else
# NO choice for below
USELAPACK=0
USELAPACKNEW=0
endif
#
####################










##########################
#
# report
#
##########################

$(warning USEMPI is $(USEMPI))
$(warning USEOPENMP is $(USEOPENMP))
$(warning USELAPACK is $(USELAPACK))
$(warning USELAPACKNEW is $(USELAPACKNEW))
$(warning ENFORCEPRECISION is $(ENFORCEPRECISION))



#END


