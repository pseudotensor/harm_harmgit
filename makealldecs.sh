#!/bin/bash
sh ./makedecs.h.sh defs.h &
sh ./makedecs.h.sh mpidefs.h &
sh ./makedecs.h.sh rancdefs.h &
sh ./makedecs.h.sh kazfulleos.defsglobalprivate.h &
sh ./makedecs.h.sh kazfulleos.superdefs.h &
sh ./makedecs.h.sh superdefs.h &
sh ./makedecs.h.sh superdefs.pointers.h &
sh ./makedecs.h.sh superdefs.rad.h &
sh ./makedecs.h.sh superdefs.pointers.rad.h &
sh ./makedecs.h.sh defs.user.h &
sh ./makedecs.h.sh defs.general.h &
sh ./makedecs.h.sh supermpidefs.h &
sh ./makedecs.h.sh mpidefs.mpi_grmhd_grray_liaison.h &
wait
