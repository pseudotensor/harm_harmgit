

ifeq ($(USEMPI),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) -cc=gcc $(GCCCFLAGS)
endif # endif usempich==1

ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=gcc  $(GCCCFLAGS)
endif # endif usempich==0

include maketail.harm.inc

# "make prep" must be called first so these are set correctly
MYCOMMANDS1 := $(shell cat objsonepoint.cmd.list ) ;
MYCOMMANDS2 := $(shell cat objsmanypoints.cmd.list ) ; 


#all:	$(PREP) $(CMD) $(FINISH)
all:	$(CMD) $(FINISH)



$(PREP):
	( sh ./makedecs.h.sh )
	( sh ./makempidecs.h.sh )
	( echo $(OBJSONEPOINT) > objsonepoint )
	( echo $(OBJSMANYPOINTS) > objsmanypoints )
#	( for fil in `cat objsonepoint` ; do echo $$fil ; done )
	( rm -rf objsonepoint.cmd.list )
	( for fil in `cat objsonepoint` ; do echo $(CC) $(CFLAGSPRENONPRECISE) -c $$fil $(LDFLAGS) >> objsonepoint.cmd.list ; done )
	( rm -rf objsmanypoints.cmd.list )
	( for fil in `cat objsmanypoints` ; do echo $(CC) $(CFLAGS) -c	$$fil $(LDFLAGS) >> objsmanypoints.cmd.list ; done )
	
$(FINISH):
	rm -rf objsonepoint objsmanypoints

#$(CMD):	$(OBJSMANYPOINTS) makefile $(OBJSONEPOINT) freespace.o
#	$(CC) $(CFLAGS) -c $(SRCSONEPOINT) $(LDFLAGS)
#	$(CC) $(CFLAGSPRENONPRECISE) -c $(SRCSONEPOINT) $(LDFLAGS)

$(CMD):	GROUP1 GROUP2 makefile freespace.o
	$(CC2) $(GCCCFLAGS) -c freespace.c $(LDFLAGS)
	$(CC)  $(CFLAGS) -o $(CMD) $(OBJSMANYPOINTS) $(OBJSONEPOINT) freespace.o $(LDFLAGS)

GROUP2:	$(OBJSMANYPOINTS) makefile
#	$(CC) $(CFLAGS) -c $(SRCSONEPOINT) $(LDFLAGS)
#	CFLAGS=$(CFLAGS)
	$(MYCOMMANDS2)

GROUP1:	$(OBJSONEPOINT) makefile
#	$(CC) $(CFLAGSPRENONPRECISE) -c $(SRCSONEPOINT) $(LDFLAGS)
#	CFLAGS=$(CFLAGSNONPRECISE)
	$(MYCOMMANDS1)

clean:
	rm *.o

cleani:
	rm *.o *.il

cleanall:
	rm *.o *.il *~

cleanbackup:
	rm *~



include maketail.harm.deps.inc

