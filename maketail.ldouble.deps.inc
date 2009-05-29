
# Use $(OBJS) in ar command for libml.a if possible; else *.o
libml.a: $(OBJS) mconf.h
	ar -rv libml.a $(OBJS)
	ranlib libml.a
#	ar -rv libml.a *.o

# Use assembly language sqrt and floor if available.
# sqrtl.c does not give strictly rounded results.
# The floor.c routines are slow.
#floorl.o: floorl.387
#	as -o floorl.o floorl.387
#sqrtl.o: sqrtl.387
#	as -o sqrtl.o sqrtl.387
#setprec.o: setprec.387
#	as -o setprec.o setprec.387

# ELF versions for linux
#floorl.o: floorlelf.387
#	as -o floorl.o floorlelf.387
sqrtl.o: sqrtlelf.387
	as -o sqrtl.o sqrtlelf.387
setprec.o: setprelf.387
	as -o setprec.o setprelf.387

# 68K routines
# For Sun 3
#floorl.o: floorl.sun
#	as -o floorl.o floorl.sun
#setprec.o: setprec.688
#	as -o setprec.o setprec.688
#sqrtl.o: sqrtl.sun
#	as -o sqrtl.o sqrtl.sun

# For Motorola 68k sysv
# Thanks to Philippe De Muyter <phdm@info.ucl.ac.be>
#floorl.o: floorl.mot
#	as -o floorl.o floorl.mot
#setprec.o: setprec.mot
#	as -o setprec.o setprec.mot
#sqrtl.o: sqrtl.mot
#	as -o sqrtl.o sqrtl.mot

fltestl: fltestl.c libml.a
	$(CC) $(CFLAGSOTHER) -o fltestl fltestl.c libml.a

fltestl.o: fltestl.c

flrtstl: flrtstl.c libml.a
	$(CC) $(CFLAGSOTHER) -o flrtstl flrtstl.c libml.a

flrtstl.o: flrtstl.c

nantst: nantst.c libml.a
	$(CC) $(CFLAGSOTHER) -o nantst nantst.c libml.a

nantst.o: nantst.c

testvect: testvect.o libml.a
	$(CC) $(CFLAGSOTHER) -o testvect testvect.o libml.a

testvect.o: testvect.c
	$(CC) -g -c -o testvect.o testvect.c

monotl: monotl.o libml.a
	$(CC) $(CFLAGSOTHER) -o monotl monotl.o libml.a

monotl.o: monotl.c
	$(CC) -g -c -o monotl.o monotl.c

acoshl.o: acoshl.c
asinhl.o: asinhl.c
asinl.o: asinl.c
atanhl.o: atanhl.c
atanl.o: atanl.c
bdtrl.o: bdtrl.c
btdtrl.o: btdtrl.c
cbrtl.o: cbrtl.c
chdtrl.o: chdtrl.c
coshl.o: coshl.c
ellpel.o: ellpel.c
ellpkl.o: ellpkl.c
elliel.o: elliel.c
ellikl.o: ellikl.c
ellpjl.o: ellpjl.c
exp10l.o: exp10l.c
exp2l.o: exp2l.c
expl.o: expl.c
expx2l.o: expx2l.c
fdtrl.o: fdtrl.c
gammal.o: gammal.c
gdtrl.o: gdtrl.c
igamil.o: igamil.c
igaml.o: igaml.c
incbetl.o: incbetl.c
incbil.o: incbil.c
isnanl.o: isnanl.c
j0l.o: j0l.c
j1l.o: j1l.c
jnl.o: jnl.c
ldrand.o: ldrand.c
log10l.o: log10l.c
log2l.o: log2l.c
logl.o: logl.c
nbdtrl.o: nbdtrl.c
ndtril.o: ndtril.c
ndtrl.o: ndtrl.c
pdtrl.o: pdtrl.c
powl.o: powl.c
powil.o: powil.c
sinhl.o: sinhl.c
sinl.o: sinl.c
sqrtl.o: sqrtl.c
stdtrl.o: stdtrl.c
tanhl.o: tanhl.c
tanl.o: tanl.c
unityl.o: unityl.c
ynl.o: ynl.c
floorl.o: floorl.c
polevll.o: polevll.c
unityl.o: unityl.c
mtherr.o: mtherr.c
