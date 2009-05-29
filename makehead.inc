#
#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"





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
USEMPI=0

endif


########################
#
# OpenMP CHOICE
#
########################
USEOPENMP=1


$(warning USEMPI is $(USEMPI))
$(warning USEOPENMP is $(USEOPENMP))


#END


