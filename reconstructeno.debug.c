//#include "decs.h"
//#include "reconstructeno.h"
//#include "reconstructeno_static.h"



// this file is included in reconstructeno.c




int check_symmetry_of_weno_input( int cvt_type, int pl, int bs, int bf, FTYPE *yin )
{
  int i;
  int is_not_symmetric = 0;

  //if( iterglobal != 2 ) return 0; //only do comparisons for 2nd direction //SASMARKx
  for( i = 0; i <= bf-bs; i++ ) {
    if( (pl < U1 && yin[i+bs] != yin[bf-i]) || (pl == U1 && yin[i+bs] != - yin[bf-i]) ) {
      dualfprintf( fail_file, "Asymmetry in the input to cvt_type = %d:  nstep = %ld, steppart = %d, i = %d\n", cvt_type, nstep, steppart, i + bs );
      is_not_symmetric = 1;
    }
  }
  return( is_not_symmetric );
}






int check_symmetry_of_weno_output( int cvt_type, FTYPE (*monoindicator)[NBIGM], int pl, int ps, int pf, int maxorder, FTYPE *yout[MAX_NO_OF_WENO_OUTPUTS], weno_weights_t *stencil_weights_array )
{
  int i, j;
  int is_not_symmetric = 0;

  //differentiate reconstructions according to the number of outputs that they provide, this is not general
  //assumes that if there is one output, it is centered, if there are two outputs -- they are at the edges
  //if( iterglobal != 2 ) return 0; //only do comparisons for 2nd direction //SASMARKx

  if( weno_outputs[cvt_type].len == 1 ) {
    ////
    //check the symmetry of the reconstruction that has 1 output (a2c or c2a)
    ////

    for( i = 0; i <= pf-ps; i++ ) {
      for( j = 0; j < maxorder; j++ ) {
        if( stencil_weights_array[i+ps].weights[j] != stencil_weights_array[pf-i].weights[maxorder-1-j] ) {
          dualfprintf( fail_file, "Asymmetry in the weights of a2c/c2a:  nstep = %ld, steppart = %d\n", nstep, steppart );
          is_not_symmetric = 1;
        }
      }
      for( j = 0; j < maxorder; j++ ) {
        if( stencil_weights_array[i+ps].weights[j+maxorder] != stencil_weights_array[pf-i].weights[maxorder + maxorder-1-j] ) {
          dualfprintf( fail_file, "Asymmetry in the optimal weights of a2c/c2a:  nstep = %ld, steppart = %d\n", nstep, steppart );
          is_not_symmetric = 1;
        }
      }
    }

    //check symmetry of output solution atch symmetrize
    for( i = 0; i <= pf-ps; i++ ) {
      if( (pl < U1 && yout[0][i+ps] != yout[0][pf-i]) || (pl == U1 && yout[0][i+ps] != - yout[0][pf-i]) ) {
        dualfprintf( fail_file, "Asymmetry in the output of a2c/c2a (%d):  nstep = %ld, steppart = %d, i = %d\n", cvt_type, nstep, steppart, i + ps );
        is_not_symmetric = 1;
      }
    }

  }
  else if( weno_outputs[cvt_type].len == 2 ) {
    ////
    //check the symmetry of the reconstruction that has 2 outputs (c2e)
    ////
    for( i = 0; i <= pf-ps; i++ ) {
      if( 
         (pl < U1 && yout[0][i+ps] != yout[1][pf-i]) || 
         (pl == U1 && yout[0][i+ps] != - yout[1][pf-i]) ) {
        dualfprintf( fail_file, "Asymmetry in the output of c2e:  nstep = %ld, steppart = %d, i = %d\n", nstep, steppart, i + ps );
        is_not_symmetric = 1;
        for( j = 0; j < maxorder; j++ ) {
          if( stencil_weights_array[i+ps].weights[j] != stencil_weights_array[pf-i].weights[maxorder-1-j] ) {
            dualfprintf( fail_file, "Asymmetry in the weights of c2e:  nstep = %ld, steppart = %d, i = %d\n", nstep, steppart, i + ps );
            is_not_symmetric = 1;
          }
        }
        for( j = 0; j < 2 * maxorder; j++ ) {
          if( stencil_weights_array[i+ps].weights[j+maxorder] != stencil_weights_array[pf-i].weights[maxorder + 2*maxorder-1-j] ) {
            dualfprintf( fail_file, "Asymmetry in the optimal weights of c2e:  nstep = %ld, steppart = %d, i = %d\n", nstep, steppart, i + ps );
            is_not_symmetric = 1;
          }
        }
      }
    }


    for( i = 0; i <= pf-ps; i++ ) {
      for( j = 0; j < maxorder; j++ ) {
        if( stencil_weights_array[i+ps].weights[j] != stencil_weights_array[pf-i].weights[maxorder-1-j] || monoindicator[MONOYIN][i+ps] != monoindicator[MONOYIN][pf-i] ) {
          dualfprintf( fail_file, "Asymmetry in the weights of c2e:  nstep = %ld, steppart = %d\n", nstep, steppart );
          is_not_symmetric = 1;
        }
      }
      for( j = 0; j < 2 * maxorder; j++ ) {
        if( stencil_weights_array[i+ps].weights[j+maxorder] != stencil_weights_array[pf-i].weights[maxorder + 2*maxorder-1-j] && monoindicator[MONOYIN][i+ps] != 1.0 ) {
          dualfprintf( fail_file, "Asymmetry in the optimal weights of c2e:  nstep = %ld, steppart = %d\n", nstep, steppart );
          is_not_symmetric = 1;
        }
      }
    }
  }
  //atch symmetrize

  return( is_not_symmetric );

}



/////////////////////////
//
// UNCORRECTED CODE BELOW
//
/////////////////////////



    
 
#if(0) // sasha has to fix all this // JONMARK GODMARK SUPERGODMARK

void debugeno_compute(FTYPE (*p)[NSTORE2][NSTORE3][NPR],FTYPE (*U)[NSTORE2][NSTORE3][NPR],FTYPE (*debugvars)[NSTORE2][NSTORE3][NUMENODEBUGS])
// NPR*3 order*3*NPR   NPR
// Uavg Upoint  order*2 per point reduce
// NPR*2 order*2*NPR NPR
{
  void para4(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout, FTYPE *dqleft, FTYPE *dqcenter, FTYPE *dqright);
  extern int get_shock_indicator(int whichprimtype, int interporflux, int dir, int bs, int ps, int pe, int be,  int i, int j, int k, int idel, int jdel, int kdel, FTYPE (*yin)[2][NBIGM], FTYPE *V, FTYPE *P, FTYPE *shockindicator);

  int a_minorderit[NBIGM];  //should actually be arrays but currently weno routines use only the 1st element of those arrays
  int a_maxorderit[NBIGM];
  int *minorderit;
  int *maxorderit;
  int preforder = 5;

  FTYPE dqleft; 
  FTYPE dqcenter; 
  FTYPE dqright;

  int shiftit = 0;  //not used; in the future would determine the shift of the weno stencil when weno stencil shifting is impelemented
  FTYPE a_yin[NBIGM]; 
  FTYPE a_yout_left[NBIGM], a_yout_right[NBIGM];
  FTYPE *yin = &a_yin[NBIGBND];
  FTYPE *yin1;
  FTYPE *yout_left = &a_yout_left[NBIGBND];
  FTYPE *yout_right = &a_yout_right[NBIGBND];
  int counter; 

  FTYPE a_allpl_yin[NPR][2][NBIGM];
  FTYPE (*allpl_yin)[2][NBIGM] = (FTYPE (*)[2][NBIGM]) (&(a_allpl_yin[0][0][NBIGM]));

  const int start_index = 0; //i, j, k
  const int start_prim = start_index + 3;  //NPR primitives
  const int start_p_l = start_prim + NPR;  //NPR left primitives
  const int start_p_r = start_p_l + NPR;   //NPR right primitives
  const int start_para_p_l = start_p_r + NPR;
  const int start_para_p_r = start_para_p_l + NPR;
  const int start_weights_c2e = start_para_p_r + NPR;  //NPR * 3 * ( 2 *order - 1 )

  //c2e finished
  //a2c starts
  const int order = 3;
  const int maxorder = 3;
  const int minorder = 2;
  const int start_Uavg = start_weights_c2e + NPR * 3 * ( 2 *order - 1 );
  const int start_Upoint = start_Uavg + NPR;
  const int start_weights_a2c = start_Upoint + NPR;  //NPR * 2 * ( 2 *order - 1 )
  const int start_shockindicator = start_weights_a2c + NPR * 2 * ( 2 *order - 1 );
 
  //a2c finished

  int i, j, k, pl, pliter;
  int bs, ps, pf, bf;
  const int whichreduce = 0;
  int idel, jdel, kdel;

  int num_orders;
  int index;

  const int num_stencil_sets_c2e = 3;
  const int num_stencil_sets_a2c = 2;
  int interporflux1;

  //para4
  FTYPE left, right;

  weno_weights_t *p_stencil_weights_array_to_be_used;

  FTYPE a_shockindicator[NBIGM], *shockindicator = a_shockindicator + NBIGBND;

  assert( N2 != 1 || N3 != 1, "Debug output can only work in 1D so far; please make sure that N2 and N3 are equal to unity.\n" ) ;

  minorderit = a_minorderit + NBIGBND;
  maxorderit = a_maxorderit + NBIGBND;

  bs = -NBIGBND;
  bf = N1 - 1 + NBIGBND;
  ps = bs + 3;
  pf = bf - 3;

  for( i = bs; i <= bf; i++ ) {
    minorderit[i] = 2 * minorder - 1;
    maxorderit[i] = 2 * maxorder - 1;
    shockindicator[i] = 0.0;  //hard-core no reduction to 1st order since don't know it
  }

  //loop over primitive quantities
  PLOOP(pliter,pl) {
    yin1 = allpl_yin[pl][0];

    //extract a line with primitives
    for( i = bs, j = 0, k = 0; i <= bf; i++ ) {
      //copy centered primitives to the debugvars array
      MACP0A1(debugvars,i,j,k,pl+start_prim) = MACP0A1(p,i,j,k,pl);

      //copy centered primitives to a line to feed to reconstruction
      yin1[i] = MACP0A1(p,i,j,k,pl);
    }

    //call the reconstruction; the weights are in stencil_weights_array[] -- but we do not used them; calculate them separately
    //because need to know if the stencil was reduced
    eno_line_c2e( whichreduce, dir, preforder, pl, bs, ps, pf, bf, minorderit, maxorderit, &shiftit, shockindicator, yin1, yout_left, yout_right ); 
  
    //copy p_l and p_r to the debug array for that pl
    for( i = ps, j = 0, k = 0; i <= pf; i++ ) {
      MACP0A1(debugvars,i,j,k,start_index) = (FTYPE)i;
      MACP0A1(debugvars,i,j,k,start_index+1) = (FTYPE)j;
      MACP0A1(debugvars,i,j,k,start_index+2) = (FTYPE)k;

      //generate para4() left and right states
      para4( pl, &yin1[i], &left, &right, &dqleft, &dqcenter, &dqright );
      MACP0A1(debugvars,i  ,j,k,pl+start_para_p_r) = left;
      MACP0A1(debugvars,i+1,j,k,pl+start_para_p_l) = right;
   
      //atch modify indices in such a way that p_l and p_r sit on the same interface for the same i
      MACP0A1(debugvars,i  ,j,k,pl+start_p_r) = yout_left[i];
      MACP0A1(debugvars,i+1,j,k,pl+start_p_l) = yout_right[i];
   
#if( DO_ENO_STENCIL_REDUCTION )
      //for each point: choose the order of interpolation
      //choose whether to reduce the stencil and compute the weights
      num_orders = choose_weno_order( CVT_C2E, whichreduce, maxorder, minorder, i, pl, bs, bf, shockindicator, P, yin1, stencil_weights_array, &p_stencil_weights_array_to_be_used );
#else
      num_orders = 1;
      p_stencil_weights_array_to_be_used = &stencil_weights_array[i];
#endif
  
      //cycle over the c2e weights, 3 groups of them: roughness, optimized left, optimized right for 3d order
      for( counter = 0; counter < num_stencil_sets_c2e * order; counter++ ) {
        index = start_weights_c2e 
          + num_stencil_sets_c2e * (order + order - 1) * pl  //shift to the start of high-order weights
          + counter;                     //shift to the current weight

        if( counter < p_stencil_weights_array_to_be_used->len && p_stencil_weights_array_to_be_used->order == 3 ) {
          MACP0A1(debugvars,i,j,k,index) = p_stencil_weights_array_to_be_used->weights[counter];
        }
        else {
          //if the weights do not exist, feed zeroes so that SM does not crash
          MACP0A1(debugvars,i,j,k,index) = 0.0;
        }
      }

      //if there are two orders whose linear combination is to be taken, retrieve the lower order weights by shifting a pointer to point to them
      if( num_orders == 2 && p_stencil_weights_array_to_be_used->order == 3 ) {
        p_stencil_weights_array_to_be_used += 1;
      }

      //lower order weights
      for( counter = 0; counter < num_stencil_sets_c2e * (order - 1); counter++ ) {
        index = start_weights_c2e 
          + num_stencil_sets_c2e * (order + order - 1) * pl  //shift to the start of high-order weights
          + num_stencil_sets_c2e * order                     //shift to the start fo lower-order weights
          + counter;                     //shift to the current weight

        if( counter < p_stencil_weights_array_to_be_used->len /*&& p_stencil_weights_array_to_be_used->order == 2*/ ) {
          MACP0A1(debugvars,i,j,k,index) = p_stencil_weights_array_to_be_used->weights[counter];
        }
        else {
          //if the weights do not exist, feed zeroes so that SM does not crash
          MACP0A1(debugvars,i,j,k,index) = 0.0;
        }
      }

      //c2e stuff finished.
    }
 
    //a2c stuff starts
  
    //extract a line with primitives
    for( i = bs, j = 0, k = 0; i <= bf; i++ ) {
      //copy centered primitives to the debugvars array
      MACP0A1(debugvars,i,j,k,start_Uavg+pl) = MACP0A1(U,i,j,k,pl);

      //copy averaged conserved quantities to a line to feed to reconstruction
      yin[i] = MACP0A1(U,i,j,k,pl);
    }

    //call the reconstruction; the weights are in stencil_weights_array[] -- but we do not used them; calculate them separately
    //because need to know if the stencil was reduced
    //yout_left is actually not "left", just a bad name for a temp variable 
    eno_line_a2c(1, 5, pl, -2, 0, N1-1, N1 + 1, minorderit, maxorderit, &shiftit, shockindicator, yin, yout_left ); 
  
    bs = -2;
    ps = 0;
    pf = N1 - 1;
    bf = N1 + 1;

    //copy Upoint to the debug array for that pl
    for( i = ps, j = 0, k = 0; i <= pf; i++ ) {
      //atch modify indices in such a way that p_l and p_r sit on the same interface for the same i
      MACP0A1(debugvars,i,j,k,start_Upoint+pl) = yout_left[i];
   
      //choose whether to reduce the stencil and compute the weights
#if( DO_ENO_STENCIL_REDUCTION )
      num_orders = choose_weno_order( CVT_A2C, whichreduce, maxorder, minorder, i, pl, bs, bf, shockindicator, yin, stencil_weights_array, &p_stencil_weights_array_to_be_used );
#else
      num_orders = 1;
      p_stencil_weights_array_to_be_used = &stencil_weights_array[i];
#endif

      //cycle over the a2c weights, 2 groups of them: unoptimal, optimal
      for( counter = 0; counter < num_stencil_sets_a2c * order; counter++ ) {
        if( counter < p_stencil_weights_array_to_be_used->len && p_stencil_weights_array_to_be_used->order == 3 ) {
          MACP0A1(debugvars,i,j,k,start_weights_a2c + 2 * (order + order - 1) * pl + counter) = p_stencil_weights_array_to_be_used->weights[counter];
        }
        else {
          //if the weights do not exist, feed zeroes so that SM does not crash
          MACP0A1(debugvars,i,j,k,start_weights_a2c + 2 * (order + order - 1) * pl + counter) = 0.0;
        }
      }

      //if there are two orders whose linear combination is to be taken, retrieve the lower order weights by shifting a pointer to point to them
      if( num_orders == 2 && p_stencil_weights_array_to_be_used->order == 3 ) {
        p_stencil_weights_array_to_be_used += 1;
      }
    
      //lower order weights, therefore re
      for( counter = 0; counter < num_stencil_sets_a2c * (order - 1); counter++ ) {
        index = start_weights_a2c 
          + num_stencil_sets_a2c * (order + order - 1) * pl  //shift to the start of high-order weights
          + num_stencil_sets_a2c * order                     //shift to the start fo lower-order weights
          + counter;                     //shift to the current weight

        if( counter < p_stencil_weights_array_to_be_used->len && p_stencil_weights_array_to_be_used->order == 2 ) {
          MACP0A1(debugvars,i,j,k,index) = p_stencil_weights_array_to_be_used->weights[counter];
        }
        else {
          //if the weights do not exist, feed zeroes so that SM does not crash
          MACP0A1(debugvars,i,j,k,index) = 0.0;
        }
      }

      //c2e stuff finished.
    }
  }
 
  idel = jdel = kdel = 1;
  i = 0; 
  j = 0;
  k = 0;
  interporflux1 = ENOINTERPTYPE;
  get_shock_indicator(1, interporflux1, dir, bs, ps, pf, bf,  i, j, k, idel, jdel, kdel, allpl_yin, V, P, shockindicator);

  index = start_shockindicator;
  for( i = ps, j = 0, k = 0; i <= pf; i++ ) {
    MACP0A1(debugvars,i,j,k,index) = shockindicator[i];
  }

} 

#else

void debugeno_compute(FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR],FTYPE (*debugvars)[NSTORE2][NSTORE3][NUMENODEBUGS])
{
  // does nothing
}

#endif


