###################
#
# see makefile.notes for detailed comments/notes
#
###################

#####################
#
# AUTO PICK SYSTEM TO COMPILE ON
#
#####################

# for this to work, should put in .bashrc or at command-line before running make:
#export HOSTNAME


$(warning HOSTNAME is $(HOSTNAME))
$(warning CFLAGS is $(CFLAGS))


#############################
#
# Default is all compilers are not used, so can reference (e.g.) OSX being set to 0 or not.
# DO NOT CHANGE THESE ZEROS!  Just use these as reference for what to add or change in the next section when you pick your system or compiler.
#
#############################

# WATCH OUT FOR SPACES, etc. AFTER ASSIGNMENTS!!!!
USEGCC=0
USEICCGENERIC=0
USEICCINTEL=0
USEBG=0
#QUEENBEE (avoids making directories or forking in general)
USEQB=0
USEUB=0
USEOSX=0
# USEABE -> USEICCGENERIC is used
USECCC=0
USEORANGE=0
USENERSC=0
USEKRAKEN=0
USEKRAKENICC=0   # best to use for Kraken
USETACCLONESTAR4=0
USETACCLONESTAR=0
USETACCRANGER=0
USEPFE=0
USENAU=0
USESTAMPEDE=0


#####################
#
# PICK SYSTEM TO COMPILE ON.  Only should change "HOSTNAME->COMPILER" section below, not default!
#
#####################

hoststriped=$(shell echo $$HOSTNAME | sed 's/[0-9]*//g')
$(warning hoststriped is $(hoststriped))


ifeq ($(HOSTNAME),)
#############
# default compiler is GCC that will always be present on linux
#############
USEGCC=1
$(warning USEcompiler is set arbitrarily in makehead.inc.  Put your system into makehead.inc to get specific compiler.)


else
#####################
# HOSTNAME->COMPILER: system name (HOSTNAME) chooses compiler to use.  So don't have to have git repo keep changing for each user.
####################

ifeq ($(HOSTNAME),physics-179.umd.edu)
USEICCINTEL=1
#USEGCC=1
$(warning USEICCINTEL is $(USEICCINTEL))
else ifeq ($(hoststriped),krakenpf)
USEKRAKENICC=1
$(warning USEKRAKENICC is $(USEKRAKENICC))
else ifeq ($(HOSTNAME),ki-rh42)
USEICCINTEL=1
$(warning USEICCINTEL is $(USEICCINTEL))
else ifeq ($(HOSTNAME),ki-rh39)
USEICCINTEL=1
$(warning USEICCINTEL is $(USEICCINTEL))
else ifeq ($(HOSTNAME),ki-jmck)
USEICCINTEL=1
$(warning USEICCINTEL is $(USEICCINTEL))
else ifeq ($(hoststriped),login.stampede.tacc.utexas.edu)
USESTAMPEDE=1
$(warning USESTAMPEDE is $(USESTAMPEDE))
else ifeq ($(hoststriped),pfe)
USEPFE=1
$(warning USEPFE is $(USEPFE))
else ifeq ($(HOSTNAME),ubuntu)
USEICCINTEL=1
$(warning USEICCINTEL is $(USEICCINTEL))
else ifeq ($(HOSTNAME),bh00)
USEICCINTEL=1
$(warning USEICCINTEL is $(USEICCINTEL))
else ifeq ($(HOSTNAME),bhjet)
USEOSX=1
$(warning USEOSX is $(USEOSX))
else
#$(error USEcompiler never set!.  Need to choose for your system: $(HOSTNAME) in makehead.inc .)
USEICCINTEL=1
$(warning USEICCINTEL is $(USEICCINTEL))
endif


############
# end else if having HOSTNAME choose compiler.
############
endif





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
PREPF=prepf
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
USEMPI=1

else

# no choice for below
# i.e. no MPI for iinterp, bin2txt, or smcalc
USEMPI=1

endif




########################
#
# OpenMP CHOICE
#
# With icc and koralinsert and bounds.koral.c, issue with compile with OpenMP.  Can turn off OpenMP for just that one file by trying compile and copying compile like for bounds.c and remove USINGOPENMP and -fopenmp
#
########################
USEOPENMP=0


########################
#
# Pthread CHOICE
#
########################
USEPTHREAD=1



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


# choice for below
# For Kraken, choose USEGSL=0 and comment-out #include "gslincludes.h" in global.general.h : Kraken doesn't have newer GSL that's required (missing gsl/gsl_odeiv2.h) and anyways has compile time issue for some reason.
USEGSL=0
#
########################









##########################
#
# report
#
##########################

$(warning USEMPI is $(USEMPI))
$(warning USEOPENMP is $(USEOPENMP))
$(warning USEPTHREAD is $(USEPTHREAD))
$(warning USELAPACK is $(USELAPACK))
$(warning USELAPACKNEW is $(USELAPACKNEW))
$(warning ENFORCEPRECISION is $(ENFORCEPRECISION))



#END


