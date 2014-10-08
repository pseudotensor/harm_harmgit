########################
#
# START of maketail.inc
#
########################

# for long doubles -m96bit-long-double 
##CFLAGSGCCOTHER= -O2 -m128bit-long-double -Wall -Wuninitialized -fno-builtin -DUSINGICC=0
##CFLAGSOTHER= -O3 -long-double -unroll -Wall -Wuninitialized -fno-builtin -DUSINGICC=0

# set some constants
ONE=1
ZERO=0


########################
#
# Setup CC and CC2 and MMCOMP
#
########################

ifeq ($(USEMPI),1)

ifeq ($(USEMCCSWITCH),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
endif

ifeq ($(USEMCCSWITCH),0)
CC=$(MCC) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
endif

ifeq ($(USEMCCSWITCHFORGCC),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) -cc=gcc $(GCCCFLAGS)
endif

ifeq ($(USEMCCSWITCHFORGCC),0)
CC2=gcc $(GCCCFLAGS)
endif

MMCOMP=$(CC)
#MMCOMP=$(MCC) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
#MMCOMP=$(COMP)

endif
############## endif usempich==1



ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=gcc  $(GCCCFLAGS)

MMCOMP=$(COMP)

endif
########### endif usempich==0




########################
#
# How to compile generatenprs
#
########################

ifeq ($(USESPECIAL4GENERATE),1)
CCGENERATEFINAL=$(CCGENERATE)
endif

ifeq ($(USESPECIAL4GENERATE),0)
CCGENERATEFINAL=$(CC)
endif


########################
#
# REPORTS
#
########################

$(warning COMP is $(COMP))
$(warning MMCOMP is $(MMCOMP))
$(warning CC is $(CC))
$(warning CC2 is $(CC2))
$(warning CCGENERATEFINAL is $(CCGENERATEFINAL))


########################
#
# INCLUDES
#
########################

# when adding C files, only need to modify the below
include maketail.harm.inc

# below is for long double code (not HARM)
include maketail.ldouble.inc

# below is for liaison code
include maketail.liaison.inc

# below is for joninterp code
include maketail.joninterp.inc

# below is for bin2txt code
include maketail.bin2txt.inc

# below is for smcalc code
include maketail.smcalc.inc



# "make prep" must be called first so these are set correctly
#MYCOMMANDS1 := $(shell cat objsonepoint.cmd.list ) 
#MYCOMMANDS2 := $(shell cat objsmanypoints.cmd.list ) 
#MYCOMMANDS3 := $(shell cat objspuredouble.cmd.list ) 


########################
#
# PRIMARY MAKE CALL:
#
########################

#all:	$(PREP) $(CMD) $(FINISH)
all:	$(CMD) $(FINISH)


########################
#
# CLEANUP:
#
########################
# CAREFUL with spaces below -- don't put space in front of "*":

cleanbackup:
	rm *~

clean:
	rm -rf *.o

cleani: clean
	rm -rf *.il

cleanall: cleani
	rm -rf *~

superclean: cleanall
	rm -rf $(CMD) $(CMD2) $(LIAISONCMD) $(JONINTERPCMD) $(BIN2TXTCMD) $(SMCALCCMD) $(PREP) $(PREPLD) $(PREPLIAISON) $(PREPJONINTERP) deps.onepoint.list deps.manypoints.list deps.puredouble.list objsmanypoints.cmd.list objsonepoint.cmd.list objslong.cmd.list objspuredouble.cmd.list global.defnprs.h generatenprs libml.a
	rm -rf objsonepoint.cmd.list objsmanypoints.cmd.list objspuredouble.cmd.list objslong.cmd.list objsdepsliaison.cmd.list objsdepsjoninterp.cmd.list
	rm -rf deps.manypoints.list deps.onepoint.list deps.puredouble.list deps.long.list deps.depsliaison.list deps.depsjoninterp.list
	rm -rf objsonepoint objsmanypoints objspuredouble onepoint manypoints puredouble
	rm -rf decs.h mpidecs.h superdecs.h superdecs.pointers.h decs.user.h decs.general.h supermpidecs.h mpidecs.mpi_grmhd_grray_liaison.h rancdecs.h
	rm -rf freespace.o tau_neededbyharm.o testfpp.o
	rm -rf depssmcalc depsliaison depsjoninterp deps.depssmcalc.list deps.depsbin2txt.list depsbin2txt
# below are just links
	rm -rf global.h defs.h
	rm -rf extractnonfunc gprof-helper.so
	rm -rf objslong objsdepssmcalc.cmd.list objsdepssmcalc objsdepsliaison objsdepsjoninterp objsdepsbin2txt.cmd.list objsdepsbin2txt
	rm -rf prepbin2txt prepsmcalc
	rm -rf initbase.defaultnprlists.c
	rm -rf *.deps.puredouble.list *.deps.manypoints.list *.deps.onepoint.list

superfclean: cleanall superclean
	rm -rf $(PREPF)
	
# removes in addition the links to user code
# useful when wanting to compare with SVN that has no such files
ultraclean: superclean
	    rm -rf init.h init.c bounds.c
	

########################
#
# Setup how to create global.defnprs.h using generatenprs
#
########################

	# make defnprs.h
global.defnprs.h : global.realdef.h global.nondepmnemonics.h definit.h init.h generatenprs.c
	$(CCGENERATEFINAL) -o generatenprs generatenprs.c ; ./generatenprs


########################
#
# Prepare to make GRMHD code
#
########################

PREPDEP=makefile makedecs.h.sh makehead.inc maketail.harm.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc defs.grmhd.h defs.general.h mpidefs.h rancdefs.h superdefs.h superdefs.pointers.h supermpidefs.h global.defnprs.h kazfulleos.defsglobalprivate.h  kazfulleos.superdefs.h mpidefs.mpi_grmhd_grray_liaison.h


$(PREP): $(PREPDEP)
 	############################################
	#below is for GRMHD code
 	############################################
	#
	#
	# can't make below 2 background
	rm -rf defs.h ; ln -s defs.grmhd.h defs.h
	rm -rf global.h ; ln -s global.grmhd.h global.h
	#
	sh ./makealldecs.sh
	#
 	############################################
	#below is for code working per point
 	############################################
	echo $(OBJSONEPOINT) > objsonepoint
	# remove and touch (so exists) files to create
	rm -rf objsonepoint.cmd.list ; rm -rf onepoint ; touch objsonepoint.cmd.list ; touch onepoint
	# create list of prefixes
	
	for fil in `cat objsonepoint` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> onepoint ; done
	# create list of commands
	for fil in `cat onepoint` ; do echo ${ECHOSWITCH} "$$fil.o :" >> objsonepoint.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS)" >> objsonepoint.cmd.list ; done
	# remove old depenedencies file and create
	rm -rf deps.onepoint.list ; touch deps.onepoint.list
	# create new dependencies file
	for fil in `cat onepoint` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c > $$fil.deps.onepoint.list & done; wait; 
	for fil in `cat onepoint` ; do  (cat $$fil.deps.onepoint.list >> deps.onepoint.list) ; done ; wait
	#
	#
 	############################################
	#below is for code working on multiple point at once
 	############################################
	# below is only different in CFLAGS part
	echo $(OBJSMANYPOINTS) > objsmanypoints
	# remove and touch (so exists) files to create
	rm -rf objsmanypoints.cmd.list ; rm -rf manypoints ; touch objsmanypoints.cmd.list ; touch manypoints
	# create list of prefixes
	for fil in `cat objsmanypoints` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> manypoints ; done
	# create list of commands
	for fil in `cat manypoints` ; do echo ${ECHOSWITCH} "$$fil.o :" >> objsmanypoints.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGS) -c $$fil.c $(LDFLAGS)" >> objsmanypoints.cmd.list ; done
	# remove old depenedencies file and create
	rm -rf deps.manypoints.list ; touch deps.manypoints.list
	# create new dependencies file
	for fil in `cat manypoints` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c > $$fil.deps.manypoints.list & done; wait; 
	for fil in `cat manypoints` ; do  (cat $$fil.deps.manypoints.list >> deps.manypoints.list) ; done
	echo	$(OBJSMANYPOINTS) > $(PREP)
	#
	# wait for background processes
	wait
	#
 	############################################
	#below is for double-only code (that can't be long double converted or compiled with -long_double
 	############################################
	echo $(OBJSPUREDOUBLE) > objspuredouble
	# remove and touch (so exists) files to create
	rm -rf objspuredouble.cmd.list ; rm -rf puredouble ; touch objspuredouble.cmd.list ; touch puredouble
	# create list of prefixes
	
	for fil in `cat objspuredouble` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> puredouble ; done
	# create list of commands
	for fil in `cat puredouble` ; do echo ${ECHOSWITCH} "$$fil.o :">> objspuredouble.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGSNONPRECISENOLONG) -c $$fil.c $(LDFLAGS)" >> objspuredouble.cmd.list ; done
	# remove old depenedencies file and create
	rm -rf deps.puredouble.list ; touch deps.puredouble.list
	# create new dependencies file
	for fil in `cat puredouble` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c > $$fil.deps.puredouble.list & done; wait; 
	for fil in `cat puredouble` ; do  (cat $$fil.deps.puredouble.list >> deps.puredouble.list) ; done
	#
$(PREPF): $(PREPDEP)
 	############################################
	#below is for GRMHD code
 	############################################
 	############################################
	#below is for fortran bits
 	############################################
	#Sasha added the two seds that cleaned up after fpp (remove stray characters "&" and "f" at the beginning of lines)
ifeq ($(USEKRAKEN),$(ZERO))
	fpp -P test.f > testfpp.temp1.f ; sed 's/^f.*$$//g' testfpp.temp1.f | sed 's/^\&.*$$//g' > testfpp.temp.f; f2c -Wall -f -P testfpp.temp.f ; sed 's/static//g' testfpp.temp.c > testfpp.c ;sed 's/static//g' testfpp.temp.P > testfpp.P
endif
    # otherwise assume created testfpp.P and testfpp.c on another computer and copied to Kraken
	echo "DONE" > $(PREPF)
    #
    

########################
#
# Prepare to make long double version of GRMHD code
#
########################

 	#
$(PREPLD):  $(PREPDEP) maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc
 	############################################
	# below is for long double code
 	############################################
	#
	( rm -rf defs.h ; ln -s defs.grmhd.h defs.h )
	( rm -rf global.h ; ln -s global.grmhd.h global.h )
	( sh ./makedecs.h.sh defs.h )
	( sh ./makedecs.h.sh mpidefs.h )
	( sh ./makedecs.h.sh rancdefs.h )
	( sh ./makedecs.h.sh kazfulleos.defsglobalprivate.h )
	( sh ./makedecs.h.sh kazfulleos.superdefs.h )
	( sh ./makedecs.h.sh superdefs.h )
	( sh ./makedecs.h.sh superdefs.pointers.h )
	( sh ./makedecs.h.sh superdefs.rad.h )
	( sh ./makedecs.h.sh superdefs.pointers.rad.h )
	( sh ./makedecs.h.sh defs.user.h )
	( sh ./makedecs.h.sh defs.general.h )
	( sh ./makedecs.h.sh supermpidefs.h )
	( sh ./makedecs.h.sh mpidefs.mpi_grmhd_grray_liaison.h )
 	#
	( echo $(OBJSOTHER) > objslong )
	# remove and touch (so exists) files to create
	( rm -rf objslong.cmd.list ; rm -rf long ; touch objslong.cmd.list ; touch long)
	# create list of prefixes
	( for fil in `cat objslong` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> long ; done )
	# create list of commands
	( for fil in `cat long` ; do echo ${ECHOSWITCH} "$$fil.o :" >> objslong.cmd.list ; echo ${ECHOSWITCH} "	$(CC2) $(GCCCFLAGS) -c $$fil.c $(LDFLAGS)" >> objslong.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.long.list ; touch deps.long.list )
	# create new dependencies file
	( for fil in `cat long` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c >> deps.long.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSMANYPOINTS) > $(PREPLD) )



########################
#
# Prepare to make liaison code
#
########################

$(PREPLIAISON): $(PREPDEP) maketail.liaison.inc defs.liaison.h defs.liaison.h liaison.defs.h
 	############################################
	#below is for liaison code
 	############################################
	#
	#
	#
	( rm -rf defs.h ; ln -s defs.liaison.h defs.h )
	( rm -rf global.h ; ln -s global.liaison.h global.h )
	( sh ./makedecs.h.sh defs.h )
	( sh ./makedecs.h.sh mpidefs.h )
	( sh ./makedecs.h.sh rancdefs.h )
	( sh ./makedecs.h.sh kazfulleos.defsglobalprivate.h )
	( sh ./makedecs.h.sh kazfulleos.superdefs.h )
	( sh ./makedecs.h.sh superdefs.h )
	( sh ./makedecs.h.sh superdefs.pointers.h )
	( sh ./makedecs.h.sh superdefs.rad.h )
	( sh ./makedecs.h.sh superdefs.pointers.rad.h )
	( sh ./makedecs.h.sh supermpidefs.h )
	( sh ./makedecs.h.sh liaison.defs.h )
	( sh ./makedecs.h.sh mpidefs.mpi_grmhd_grray_liaison.h )
	#
	( echo $(OBJSLIAISON) > objsdepsliaison )
	# remove and touch (so exists) files to create
	( rm -rf objsdepsliaison.cmd.list ; rm -rf depsliaison ; touch objsdepsliaison.cmd.list ; touch depsliaison)
	# create list of prefixes
	( for fil in `cat objsdepsliaison` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> depsliaison ; done )
	# create list of commands
	( for fil in `cat depsliaison` ; do echo ${ECHOSWITCH} "$$fil.o :" >> objsdepsliaison.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS)" >> objsdepsliaison.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.depsliaison.list ; touch deps.depsliaison.list )
	# create new dependencies file
	( for fil in `cat depsliaison` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c >> deps.depsliaison.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSLIAISON) > $(PREPLIAISON) )
	#


########################
#
# Prepare to make iinterp code
#
########################


$(PREPJONINTERP): $(PREPDEP) defs.jon_interp.h global.jon_interp.h
 	############################################
	#below is for jon_interp code
 	############################################
	#
	#
	#
	( rm -rf defs.h ; ln -s defs.jon_interp.h defs.h )
	( rm -rf global.h ; ln -s global.jon_interp.h global.h )
	( sh ./makedecs.h.sh rancdefs.h )
	( sh ./makedecs.h.sh defs.h )
	#
	( echo $(OBJSJONINTERP) > objsdepsjoninterp )
	# remove and touch (so exists) files to create
	( rm -rf objsdepsjoninterp.cmd.list ; rm -rf depsjoninterp ; touch objsdepsjoninterp.cmd.list ; touch depsjoninterp)
	# create list of prefixes
	( for fil in `cat objsdepsjoninterp` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> depsjoninterp ; done )
	# create list of commands
	( for fil in `cat depsjoninterp` ; do echo ${ECHOSWITCH} "$$fil.o :">> objsdepsjoninterp.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS)" >> objsdepsjoninterp.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.depsjoninterp.list ; touch deps.depsjoninterp.list )
	# create new dependencies file
	( for fil in `cat depsjoninterp` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c >> deps.depsjoninterp.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSJONINTERP) > $(PREPJONINTERP) )
	#

########################
#
# Prepare to make bin2txt code
#
########################

$(PREPBIN2TXT): makefile makehead.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc
 	############################################
	#below is for bin2txt code
 	############################################
	#
	( echo $(OBJSBIN2TXT) > objsdepsbin2txt )
	# remove and touch (so exists) files to create
	( rm -rf objsdepsbin2txt.cmd.list ; rm -rf depsbin2txt ; touch objsdepsbin2txt.cmd.list ; touch depsbin2txt)
	# create list of prefixes
	( for fil in `cat objsdepsbin2txt` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> depsbin2txt ; done )
	# create list of commands
	( for fil in `cat depsbin2txt` ; do echo ${ECHOSWITCH} "$$fil.o :" >> objsdepsbin2txt.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(BIN2TXTLIBS1) $(LDFLAGS) $(BIN2TXTLIBS2)" >> objsdepsbin2txt.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.depsbin2txt.list ; touch deps.depsbin2txt.list )
	# create new dependencies file
	( for fil in `cat depsbin2txt` ; do  echo $$fil.c ; $(MMCOMP) $(BIN2TXTLIBS1) -c -MM $$fil.c >> deps.depsbin2txt.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSBIN2TXT) > $(PREPBIN2TXT) )
	#

########################
#
# Prepare to make smcalc code
#
########################

$(PREPSMCALC): makefile makehead.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc
 	############################################
	#below is for smcalc code
 	############################################
	#
	( echo $(OBJSSMCALC) > objsdepssmcalc )
	# remove and touch (so exists) files to create
	( rm -rf objsdepssmcalc.cmd.list ; rm -rf depssmcalc ; touch objsdepssmcalc.cmd.list ; touch depssmcalc)
	# create list of prefixes
	( for fil in `cat objsdepssmcalc` ; do prefa=`basename $$fil .o` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo ${ECHOSWITCH} $$pref >> depssmcalc ; done )
	# create list of commands
	( for fil in `cat depssmcalc` ; do echo ${ECHOSWITCH} "$$fil.o :" >> objsdepssmcalc.cmd.list ; echo ${ECHOSWITCH} "	$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS)" >> objsdepssmcalc.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.depssmcalc.list ; touch deps.depssmcalc.list )
	# create new dependencies file
	( for fil in `cat depssmcalc` ; do  echo $$fil.c ; $(MMCOMP) -c -MM $$fil.c >> deps.depssmcalc.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSSMCALC) > $(PREPSMCALC) )
	#






########################
#
# COMMANDS to compile codes:
#
########################


$(FINISH):
#	rm -rf objsonepoint objsmanypoints objspuredouble

#$(CMD):   	GROUP1 GROUP2 makefile freespace.o tau_neededbyharm.o testfpp.o
#	echo objsmanypoints: $(OBJSMANYPOINTS)
#	echo objsonepoint: $(OBJSONEPOINT)
#	echo objspuredouble: $(OBJSPUREDOUBLE)
#	$(MYCOMMANDS3)
#	$(MYCOMMANDS2)
#	$(MYCOMMANDS1)


$(CMD):	$(PREP) $(OBJSMANYPOINTS) $(OBJSONEPOINT) $(OBJSPUREDOUBLE) makefile testfpp.P testfpp.c 
	$(CC)  $(CFLAGSNONPRECISE) -o $(CMD) $(OBJSMANYPOINTS) $(OBJSONEPOINT) $(OBJSPUREDOUBLE) $(LDFLAGS)


$(CMD2):       $(PREP) $(PREPLD) $(OBJSMANYPOINTS) $(OBJSONEPOINT) $(OBJSPUREDOUBLE) makefile libml.a $(OBJSOTHER)
	       $(CC2) $(GCCCFLAGS) -c polevll.c $(LDFLAGS) $(LDFLAGSOTHER)
	       $(CC)  $(CFLAGS) -o $(CMD2) $(OBJSMANYPOINTS) $(OBJSONEPOINT) $(OBJSPUREDOUBLE) polevll.o libml.a $(LDFLAGS) $(LDFLAGSOTHER)

$(LIAISONCMD): $(PREPLIAISON) $(OBJSLIAISON) makefile 
	       $(CC)  $(CFLAGS) -DUSINGLIAISON=1 -o $(LIAISONCMD) $(OBJSLIAISON) $(LDFLAGS) $(LDFLAGSOTHER)

$(JONINTERPCMD): $(PREPJONINTERP) $(OBJSJONINTERP) makefile
	       $(CC)  $(CFLAGS) -DUSINGJONINTERP=1 -o $(JONINTERPCMD) $(OBJSJONINTERP) $(LDFLAGS) $(LDFLAGSOTHER)


$(BIN2TXTCMD): $(PREPBIN2TXT) $(OBJSBIN2TXT) makefile
	       $(CC)  $(CFLAGS) -o $(BIN2TXTCMD) $(OBJSBIN2TXT) $(LDFLAGS) $(LDFLAGSOTHER)  $(BIN2TXTLIBS2)


$(SMCALCCMD): $(PREPSMCALC) $(OBJSSMCALC) makefile
	       $(CC)  $(CFLAGS)  -o $(SMCALCCMD) $(OBJSSMCALC) $(LDFLAGS) $(LDFLAGSOTHER) 

 
########################
#
# Includes of object lists and dependencies autogenerated during prepare step
# Note that if one does just "make" after changing the makefile, the below are included *before* the automatically-called prep command changes them. So not consistent.  Have to do "make prep ; make" if one changes the makefile itself.
#
########################


ifeq ($(MYMAKECMDGOALS),$(filter $(MYMAKECMDGOALS),$(BINARYTARGETS)))
-include objsonepoint.cmd.list
-include objsmanypoints.cmd.list
-include objspuredouble.cmd.list
-include deps.manypoints.list
-include deps.onepoint.list
-include deps.puredouble.list
endif

ifeq ($(MYMAKECMDGOALS),$(CMD2))
-include objslong.cmd.list
-include deps.long.list
endif

ifeq ($(MYMAKECMDGOALS),$(LIAISONCMD))
-include objsdepsliaison.cmd.list
-include deps.depsliaison.list
endif

ifeq ($(MYMAKECMDGOALS),$(JONINTERPCMD))
-include objsdepsjoninterp.cmd.list
-include deps.depsjoninterp.list
endif

ifeq ($(MYMAKECMDGOALS),$(BIN2TXTCMD))
-include objsdepsbin2txt.cmd.list
-include deps.depsbin2txt.list
endif

ifeq ($(MYMAKECMDGOALS),$(SMCALCCMD))
-include objsdepssmcalc.cmd.list
-include deps.depssmcalc.list
endif



# below automatic now (original has special code for assembly but doesn't work anymore on ki-rh42)
#-include maketail.ldouble.deps.inc


########################
#
# How to compile other non-standard parts of code.
#
# Old Intel compiler wouldn't compile freespace properly due to sys/vfs.h stuff.
# Old Intel compiler wouldn't compile tau_needbyharm properly due to f2c conversion.  Now tau_neededbyharm doesn't contain anything that requires f2c libraries.
#
# So, now (06/04/09) seems fine.
#
# This all was spurred by trying to avoid gcc since on TACC's Ranger, their gcc changed with mpicc -cc=gcc error due to -i-dynamic not recognized by gcc.  Stupid Intel specific commands in their mpicc.
#
########################

#freespace.o : freespace.c $(PREP) makefile maketail.inc defs.h defs.general.h
#	    $(CC2) $(GCCCFLAGS) -c freespace.c $(LDFLAGS)
#
#tau_neededbyharm.o : tau_neededbyharm.c $(PREP) makefile maketail.inc defs.h defs.general.h
##	    $(CC2) $(GCCCFLAGS) -c tau_neededbyharm.c -lm
##	    $(CC2) $(GCCCFLAGS) -c tau_neededbyharm.c -lf2c -lm
#testfpp.o : testfpp.c $(PREP) makefile maketail.inc defs.h defs.general.h
##	    $(CC2) $(GCCCFLAGS) -c testfpp.c -lm
##	    $(CC2) $(GCCCFLAGS) -c testfpp.c -lf2c -lm

# Might generally want to compile with -lf2c, unless sure have all needed functions.
#If getting "undefined reference to MAIN__", can fix by either 
#    making sure that you're linking to libf2c.a and not libf2c.so, 
#    or adding option "-u MAIN__" to linker.



########################
#
# how to compile long double library
#
########################


# Use $(OBJS) in ar command for libml.a if possible; else *.o
libml.a: $(OBJSOTHER) mconf.h objslong.cmd.list deps.long.list
	 ar -rv libml.a $(OBJSOTHER)
	 ranlib libml.a
#	 ar -rv libml.a *.o







########################
#
# END of maketail.inc
#
########################
