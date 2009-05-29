#ifndef RECONSTRUCTENO_DEFS_H
#define RECONSTRUCTENO_DEFS_H

// This file is only included in reconstructeno.c, so no need to convert defs.h -> decs.h unless to be used for other files.

#include "reconstructeno_superdefs.h" // OPENMPNOTE: only global arrays with per i,j,k allowed, so applicable to OpenMP

int do_weno_lower_order_fraction; // OPENMPNOTE: ok global variable since dont for all i,j,k


// Below is just a check, not written to:
// this is defined NON-static to ensure this file is only included once
// Otherwise multiply defined functions and memory would be created per file and not conflict but also not share as probably desired
// This should generate a multiply-defined error during linking
FTYPE a_reconstruct_static_h_can_only_be_included_in_reconstructeno_dot_c;
FTYPE (*reconstruct_static_h_can_only_be_included_in_reconstructeno_dot_c) = (FTYPE (*)) (&a_reconstruct_static_h_can_only_be_included_in_reconstructeno_dot_c);






#endif
