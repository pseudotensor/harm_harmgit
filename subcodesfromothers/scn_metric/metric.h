#define MCOORD 1		// coordinates for Kerr metric
// 0 : Boyer-Lindquist (based on r theta)
// 1 : Kerr-Schild (based on bl coords)
// contains metric definitions
// Kerr-Schild: future regularized, ep=-1, k=1

#define ks_gcov00 (-1. + 2.*r/rho2)
#define ks_gcov01 (2.*r/rho2)
#define ks_gcov02 (0)
#define ks_gcov03 (-2.*a*r*s2/rho2)
#define ks_gcov10 (ks_gcov01)
#define ks_gcov11 (1. + 2.*r/rho2)
#define ks_gcov12 (0)
#define ks_gcov13 (-a*s2*(1. + 2.*r/rho2))
#define ks_gcov20 (0)
#define ks_gcov21 (0)
#define ks_gcov22 (rho2)
#define ks_gcov23 (0)
#define ks_gcov30 (ks_gcov03)
#define ks_gcov31 (ks_gcov13)
#define ks_gcov32 (0)
#define ks_gcov33 (s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)))

#define ks_gcon00 (-(1.+2.*r/rho2))
#define ks_gcon01 (2.*r/rho2)
#define ks_gcon02 (0)
#define ks_gcon03 (0)
#define ks_gcon10 (ks_gcon01)
#define ks_gcon11 ((r*(r-2.)+a*a)/rho2)
#define ks_gcon12 (0)
#define ks_gcon13 (a/rho2)
#define ks_gcon20 (ks_gcon02)
#define ks_gcon21 (ks_gcon12)
#define ks_gcon22 (1./rho2)
#define ks_gcon23 (0)
#define ks_gcon30 (ks_gcon03)
#define ks_gcon31 (ks_gcon13)
#define ks_gcon32 (ks_gcon23)
#define ks_gcon33 (1./(rho2*s2))

#define bl_gdet   (r*r*fabs(sin(th))*(1. + 0.5*(a2/r2)*(1. + cos(2.*th))))

#define bl_gcov00 (-(1. - 2./(r*mu)))
#define bl_gcov01 (0)
#define bl_gcov02 (0)
#define bl_gcov03 (-2.*a*s2/(r*mu))
#define bl_gcov10 (0)
#define bl_gcov11 (mu/DD)
#define bl_gcov12 (0)
#define bl_gcov13 (0)
#define bl_gcov20 (0)
#define bl_gcov21 (0)
#define bl_gcov22 (r2*mu)
#define bl_gcov23 (0)
#define bl_gcov30 (bl_gcov03)
#define bl_gcov31 (0)
#define bl_gcov32 (0)
#define bl_gcov33 (r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu)))

#define bl_gcon00 (-1. - 2.*(1. + a2/r2)/(r*DD*mu))
#define bl_gcon01 (0)
#define bl_gcon02 (0)
#define bl_gcon03 (-2.*a/(r3*DD*mu))
#define bl_gcon10 (0)
#define bl_gcon11 (DD/mu)
#define bl_gcon12 (0)
#define bl_gcon13 (0)
#define bl_gcon20 (0)
#define bl_gcon21 (0)
#define bl_gcon22 (1./(r2*mu))
#define bl_gcon23 (0)
#define bl_gcon30 (bl_gcon03)
#define bl_gcon31 (0)
#define bl_gcon32 (0)
#define bl_gcon33 ((1. - 2./(r*mu))/(r2*sth*sth*DD))


// transformation from bl only
// identity is no transformation
#define bl_trans00   (1)
#define bl_trans01   (0)
#define bl_trans02   (0)
#define bl_trans03   (0)
#define bl_trans10   (0)
#define bl_trans11   (1)
#define bl_trans12   (0)
#define bl_trans13   (0)
#define bl_trans20   (0)
#define bl_trans21   (0)
#define bl_trans22   (1)
#define bl_trans23   (0)
#define bl_trans30   (0)
#define bl_trans31   (0)
#define bl_trans32   (0)
#define bl_trans33   (1)


#define ks_trans00   (1)
#define ks_trans01   (2.*r/(r*r - 2.*r + a*a))
#define ks_trans02   (0)
#define ks_trans03   (0)
#define ks_trans10   (0)
#define ks_trans11   (1)
#define ks_trans12   (0)
#define ks_trans13   (0)
#define ks_trans20   (0)
#define ks_trans21   (0)
#define ks_trans22   (1)
#define ks_trans23   (0)
#define ks_trans30   (0)
#define ks_trans31   (a/(r*r - 2.*r + a*a))
#define ks_trans32   (0)
#define ks_trans33   (1)

#define ks2bl_trans00   (1)
#define ks2bl_trans01   (-2.*r/(r*r - 2.*r + a*a))
#define ks2bl_trans02   (0)
#define ks2bl_trans03   (0)
#define ks2bl_trans10   (0)
#define ks2bl_trans11   (1)
#define ks2bl_trans12   (0)
#define ks2bl_trans13   (0)
#define ks2bl_trans20   (0)
#define ks2bl_trans21   (0)
#define ks2bl_trans22   (1)
#define ks2bl_trans23   (0)
#define ks2bl_trans30   (0)
#define ks2bl_trans31   (-a/(r*r - 2.*r + a*a))
#define ks2bl_trans32   (0)
#define ks2bl_trans33   (1)




#if(MCOORD==0)

#define gcov00 bl_gcov00
#define gcov01 bl_gcov01
#define gcov02 bl_gcov02
#define gcov03 bl_gcov03
#define gcov10 bl_gcov10
#define gcov11 bl_gcov11
#define gcov12 bl_gcov12
#define gcov13 bl_gcov13
#define gcov20 bl_gcov20
#define gcov21 bl_gcov21
#define gcov22 bl_gcov22
#define gcov23 bl_gcov23
#define gcov30 bl_gcov30
#define gcov31 bl_gcov31
#define gcov32 bl_gcov32
#define gcov33 bl_gcov33


#define trans00 bl_trans00
#define trans01 bl_trans01
#define trans02 bl_trans02
#define trans03 bl_trans03
#define trans10 bl_trans10
#define trans11 bl_trans11
#define trans12 bl_trans12
#define trans13 bl_trans13
#define trans20 bl_trans20
#define trans21 bl_trans21
#define trans22 bl_trans22
#define trans23 bl_trans23
#define trans30 bl_trans30
#define trans31 bl_trans31
#define trans32 bl_trans32
#define trans33 bl_trans33

#elif(MCOORD==1)

#define gcov00 ks_gcov00
#define gcov01 ks_gcov01
#define gcov02 ks_gcov02
#define gcov03 ks_gcov03
#define gcov10 ks_gcov10
#define gcov11 ks_gcov11
#define gcov12 ks_gcov12
#define gcov13 ks_gcov13
#define gcov20 ks_gcov20
#define gcov21 ks_gcov21
#define gcov22 ks_gcov22
#define gcov23 ks_gcov23
#define gcov30 ks_gcov30
#define gcov31 ks_gcov31
#define gcov32 ks_gcov32
#define gcov33 ks_gcov33

#define gcon00 ks_gcon00
#define gcon01 ks_gcon01
#define gcon02 ks_gcon02
#define gcon03 ks_gcon03
#define gcon10 ks_gcon10
#define gcon11 ks_gcon11
#define gcon12 ks_gcon12
#define gcon13 ks_gcon13
#define gcon20 ks_gcon20
#define gcon21 ks_gcon21
#define gcon22 ks_gcon22
#define gcon23 ks_gcon23
#define gcon30 ks_gcon30
#define gcon31 ks_gcon31
#define gcon32 ks_gcon32
#define gcon33 ks_gcon33


#define trans00 ks_trans00
#define trans01 ks_trans01
#define trans02 ks_trans02
#define trans03 ks_trans03
#define trans10 ks_trans10
#define trans11 ks_trans11
#define trans12 ks_trans12
#define trans13 ks_trans13
#define trans20 ks_trans20
#define trans21 ks_trans21
#define trans22 ks_trans22
#define trans23 ks_trans23
#define trans30 ks_trans30
#define trans31 ks_trans31
#define trans32 ks_trans32
#define trans33 ks_trans33

#define met2bl_trans00 ks2bl_trans00
#define met2bl_trans01 ks2bl_trans01
#define met2bl_trans02 ks2bl_trans02
#define met2bl_trans03 ks2bl_trans03
#define met2bl_trans10 ks2bl_trans10
#define met2bl_trans11 ks2bl_trans11
#define met2bl_trans12 ks2bl_trans12
#define met2bl_trans13 ks2bl_trans13
#define met2bl_trans20 ks2bl_trans20
#define met2bl_trans21 ks2bl_trans21
#define met2bl_trans22 ks2bl_trans22
#define met2bl_trans23 ks2bl_trans23
#define met2bl_trans30 ks2bl_trans30
#define met2bl_trans31 ks2bl_trans31
#define met2bl_trans32 ks2bl_trans32
#define met2bl_trans33 ks2bl_trans33

#endif
