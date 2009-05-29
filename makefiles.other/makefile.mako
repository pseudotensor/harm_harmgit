#
#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
USEMPI=1

USEICC=0
USEGCC=1
USECCC=0

ifeq ($(USEMPI),1)
# MCC = mpicc
MCC=/usr/local/mpich-1.2/bin/mpicc
#MCC=/usr/local/p4mpich-1.2.5-icc-noshmem/bin/mpicc
endif    

LAPACKLIB=lapack
BLASLIB=blas
F2CLIB=g2c  #-l$(F2CLIB)

PREP = prep
CMD=grmhd


#
# Define preprocessor and compile flags, and linked libraries

ifeq ($(USEGCC),1)
#COMP=gcc
COMP=/usr/local/bin/gcc3.2

# gcc type flags
CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX
#CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -pg -g -a
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -ffast-math -finline-functions -funroll-loops
#CFLAGS=-O0 -g -Wall -wunused-label -wunused-parameter
#CFLAGS=-O0 -g -Wall
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -g
#-pg
#-pg -g  source lines
#-pg -g -a   line count
# gprof -l <file> > out.txt
# gprof -A -I<sourcedir>
# gprof -l -A -x s

#below is typical flags for double precision...can take -pg off for no profile
#add -mstack-align-double if using pgcc
#CFLAGS = -Wall -O0 -g
#  -fomit-frame-pointer



#CFLAGS = -Wall -O0
#CFLAGS = -O6 -g
#CFLAGS = -O0 -pg -g
LDFLAGS = -lm
# -l$(LAPACKLIB) -l$(BLASLIB)  -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -l$(F2CLIB) 

#CC = cc
#AR	=	ar r
#RANLIB	=	ranlib
endif

ifeq ($(USEICC),1)
COMP=icc
# -ipo is not same for multiple names in c files, unless static names
# -ipo_wp is generally not safe, and not for HARM
# -ip is safe, always
CFLAGS=-O3 -tpp7 -axiMKW -ipo -unroll -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310
GCCCFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math
#CFLAGS = -O3 -g -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310

# P4 (don't use -wp_ipo -- incompat code)
#CFLAGS=-O3 -tpp7 -axiMKW -wp_ipo -unroll -w1
#CFLAGS=-O3 -tpp7 -axiMKW -ipo -unroll -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310 -pg
#CFLAGS=-O3 -unroll -axiMKW -unroll -pg
# -rcd causes problems, like asymetries in interp program I've noticed
#CFLAGS = -Wall -w2 -O3 -axiMKW -unroll -ipo -tpp7 -march=pentium4 -mcpu=pentium4# -p\
#arallel
# -w2 displays more info for warnings
#CFLAGS=-O3 -axiMKW -g -ipo -pg
# P3
# CFLAGS=-O3 -tpp6 -axiMK -ipo -unroll -w1
# GAMMIE
# CFLAGS = -O3 -ipo


LDFLAGS = -lm
endif


ifeq ($(USECCC),1)
COMP=ccc
LDFLAGS =  -lm -lcxml

#CDEBUG = -g3 # -g3 for higher opts than -O0
#CDEBUG = -g
# do profile (profile screws up speed for loops, doesn't unroll them, etc.)
#CDEBUG = -pg -g3
# production level
CFLAGS3 = -Wall -O4 -fast -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2 -msg_disable longlongtype -finline-functions -funroll-loops
#CFLAGS3 = -Wall -O0
# super annoying develop level
#CFLAGS3 = -Wall -O2 -fast
#CFLAGS3 = -fast -arch ev67
# debug level
#CFLAGS3 = -Wall -O0 -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2
#CFLAGS3 = -Wall -O0
#CFLAGS3 = -Wall -O2

endif



ifeq ($(USEMPI),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) -cc=gcc $(GCCCFLAGS)
endif # endif usempich==1

ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=gcc  $(GCCCFLAGS)
endif # endif usempich==0



SRCS = \
bounds.c boundsint.c coord.c diag.c dudp_calc_3vel.c dudp_calc.c dump.c dump_ener.c fixup.c gaussj.c \
image.c initbase.c init.c interp.c lubksb.c ludcmp.c main.c metric.c \
mnewt.c nrutil.c phys.c ranc.c restart.c set_arrays.c set_grid.c \
step_ch.c tensor.c utoprim.c utoprim_2d.c utoprim_1d.c utoprim_1d_opt.c vchar.c transforms.c fail.c utoprim_ldz.c \
init_mpi.c boundmpi.c boundmpiint.c sources.c
 
OBJS = \
bounds.o boundsint.o coord.o diag.o dudp_calc_3vel.o dudp_calc.o dump.o dump_ener.o fixup.o gaussj.o \
image.o initbase.o init.o interp.o lubksb.o ludcmp.o main.o metric.o \
mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
step_ch.o tensor.o utoprim.o utoprim_2d.o utoprim_1d.o utoprim_1d_opt.o vchar.o transforms.o fail.o utoprim_ldz.o \
init_mpi.o boundmpi.o boundmpiint.o sources.o

all:	$(PREP) $(CMD)


$(PREP):
	( sh ./makedecs.h.sh )
	( sh ./makempidecs.h.sh )

$(CMD):	$(OBJS) makefile
	$(CC2) $(GCCCFLAGS) -c freespace.c $(LDFLAGS)
	$(CC)  $(CFLAGS) -o $(CMD) $(OBJS) freespace.o $(LDFLAGS)

clean:
	rm *.o

cleani:
	rm *.o *.il

cleanall:
	rm *.o *.il *~

cleanbackup:
	rm *~

# dependencies
$(OBJS) : global.h defs.h mympi.h mpidefs.h makefile
phys.o : metric.h
metric.o : metric.h
transforms.o : metric.h
utoprim_2d.o : utoprim_2d.h
utoprim_1d.o : utoprim_1d.h
utoprim_1d_opt.o : utoprim_1d_opt.h
step_ch.o : step_ch.h
interp.o : interp.h



