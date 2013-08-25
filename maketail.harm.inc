
SRCSMANYPOINTS=\
init.c init.tools.c initbase.c initbase.enerregions.c initbase.tools.c initbase.gridsectioning.c restart.c restart.checks.c set_arrays.c set_arrays_perpoint_perline.c set_arrays_multidimen.c mpi_set_arrays.c set_grid.c \
advance.c copyandinit_functions.c reconstructeno.c paraenohybrid.c flux.c flux.mergedc2ea2cmethod.c higherorder_pointavg.c wavespeeds.c boundsflux.c boundsvpot.c \
fluxct.c fluxctstag.c fluxvpot.c \
bounds.c bounds.tools.c boundsint.c \
diag.c dump.c dumpgen.c dump_ener.c image.c \
fail.c fixup.c \
mpi_init.c mpi_init_grmhd.c mpi_fileio.c mpi_fprintfs.c boundmpi.c boundmpiint.c mpi_grmhd_grray_liaison.c \
interppoint.c interppoint.para.c interpline.c interpline.para.c interpline.mono.c interpline.jmono.c interpline.smono.c \
geodesic.c \
set_arrays_multidimen.rad.c
#END


OBJSMANYPOINTS=\
init.o init.tools.o initbase.o initbase.enerregions.o initbase.tools.o initbase.gridsectioning.o restart.o restart.checks.o set_arrays.o set_arrays_perpoint_perline.o set_arrays_multidimen.o mpi_set_arrays.o set_grid.o \
advance.o copyandinit_functions.o reconstructeno.o paraenohybrid.o flux.o  flux.mergedc2ea2cmethod.o higherorder_pointavg.o wavespeeds.o boundsflux.o boundsvpot.o \
fluxct.o fluxctstag.o fluxvpot.o \
bounds.o bounds.tools.o boundsint.o \
diag.o dump.o dumpgen.o dump_ener.o image.o \
fail.o fixup.o \
mpi_init.o mpi_init_grmhd.o  mpi_fileio.o mpi_fprintfs.o boundmpi.o boundmpiint.o mpi_grmhd_grray_liaison.o \
interppoint.o interppoint.para.o interpline.o interpline.para.o interpline.mono.o interpline.jmono.o interpline.smono.o \
geodesic.o \
set_arrays_multidimen.rad.o
#END

# phys.c does currents across multiple points, but otherwise one pointed
# wavespeeds kinda at more than one point

# for some reason main.c really matters alot.  Needs to be compiled
# withOUT -mp -pc64 for precision in inversion to be ok  ODD!  GODMARK

SRCSONEPOINT=\
initbase.boundloop.c \
main.c mytime.c metric.c metric.tools.c math.tools.c metric_selfgravity_or_evolvemetric.c coord.c phys.c phys.tools.c step_ch.c fluxcompute.c \
phys.ffde.c phys.coldgrmhd.c eos.c \
vchar.c transforms.c sources.c rescale_interp.c \
gaussj.c lubksb.c ludcmp.c mnewt.c nrutil.c ranc.c tensor.c tetrad.c \
utoprimgen.c dudp_calc_3vel.c dudp_calc.c \
utoprim.orig.c utoprim_2d.c utoprim_1d.c utoprim_1d_opt.c utoprim_ldz.c \
utoprim_jon.c utoprim_1d_final.c utoprim_2d_final.c utoprim_5d2_final.c \
freespace.c tau_neededbyharm.c testfpp.c \
phys.tools.rad.c
#END

OBJSONEPOINT=\
initbase.boundloop.o \
main.o mytime.o metric.o metric.tools.o math.tools.o metric_selfgravity_or_evolvemetric.o coord.o phys.o phys.tools.o step_ch.o fluxcompute.o \
phys.ffde.o phys.coldgrmhd.o eos.o \
vchar.o transforms.o sources.o rescale_interp.o \
gaussj.o lubksb.o ludcmp.o mnewt.o nrutil.o ranc.o tensor.o tetrad.o \
utoprimgen.o dudp_calc_3vel.o dudp_calc.o \
utoprim.orig.o utoprim_2d.o utoprim_1d.o utoprim_1d_opt.o utoprim_ldz.o \
utoprim_jon.o utoprim_1d_final.o utoprim_2d_final.o utoprim_5d2_final.o \
freespace.o tau_neededbyharm.o testfpp.o \
phys.tools.rad.o
#END




SRCPUREDOUBLE=\
tetlapack_func_doubleonly.c
#END

OBJSPUREDOUBLE=\
tetlapack_func_doubleonly.o
#END



