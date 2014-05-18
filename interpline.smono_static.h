
/*! \file interpline.smono_static.h
     \brief Sasha Mono static function declarations to be included elsewhere

*/

static FTYPE transition_function( FTYPE x, FTYPE x1, FTYPE x2 ) ;
static FTYPE compute_mono_indicator_point_eno5( FTYPE *yin, FTYPE epsilon );
static FTYPE compute_mono_indicator_average_eno5( FTYPE *yin, FTYPE epsilon );
static void set_as_rough(int recontype, int i, FTYPE *yin, FTYPE (*yout)[NBIGM],FTYPE (*youtpolycoef)[NBIGM], FTYPE (*monoindicator)[NBIGM]);

////SMONO settings 
//#define DO_SMONO_C2A 1
//#define DO_SMONO_A2C 1
//#define DO_SMONO_C2E 1
//
//#define DO_MONO_1ST_DERIVATIVE 1
//#define DO_3RD_DER 1
//#define DO_SMONO_CUSP_INDICATOR 0
