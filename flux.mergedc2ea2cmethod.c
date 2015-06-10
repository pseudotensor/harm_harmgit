#include "decs.h"


/* So good news is that if I just replace FullSimplify by Simplify, then even for 7 products of 2nd order polynomials it can get the result after about 30 minutes that is actually simplified.  The optimization format package for C output then helps reduce the redundancy. */

/* FYI, the total number of multiplications is roughly, for each number of products is *twice* the below numbers (one for left and one for right): */

/* 1) 2 */
/* 2) 12 - 15 */
/* 3) 46 - 50 */
/* 4) 92 - 120 */
/* 5) 229 - 250 */
/* 6) 381 - 400 */
/* 7) 628 - 700 */



/* So things scale roughly as a power of 2 at large product number, which make sense. */

/* I have to do this for both matter and electromagnetic terms, that involve *each* of the 3 dimensions doing: */

/* Matter terms seem quite simple.  In the most general case one has for *each* left and right fluxes */

/* RHO: 3-product (\detg \rho_0 u^i) */
/* UU,U1,U2,U3 each: 4-product (\detg [ (\rho_0 + u + p)u^\mu u_\nu] ) */
/* U^{dir} 2-product (\detg [ p \delta^\mu_\nu ] ) */

/* For the EM terms it's much more involved: */

/* UU-U3 each: This one is crazy!  One has most generally for each of the 4 EOMs: */

/* 3X 7-product (\detg/(u^t)^2) (B^2 u^\mu u_\nu) */
/* 3X 5-product (\detg/(u^t)^2) ( (1/2)B^2 \delta^\mu_\nu) */
/* 6X 7-product (\detg/(u^t)^2) ( (1/2)(u_i B^i)^2 \delta^\mu_\nu) */
/* 1X 5-product (\detg/(u^t)^2) (B^\mu B_\nu) */
/* 3X 7-product (\detg/(u^t)^2) (B^\mu u_\nu (u_i B^i)) */
/* 3X 7-product (\detg/(u^t)^2) (u^\mu B_\nu (u_i B^i)) */

/* The explosion in the number of terms is because of the vector nature of the field. */

/* B1,B2,B3 each: 3-product but has to be done in 2D not 1D (i.e. simultaneous deconvolution in both directions in plane).  This I've not done yet.  I expect it to be similar to a 6-product in complexity: \detg (B^i v^j - B^j v^i) */

/* Now, an important point is that for the EMF I interpolation the 3-velocity since it would require many more interpolations to interpolate the 4-velocity to the corner.  So I suppose it wouldn't be crazy to also interpolate the 3-velocity in the stress-energy tensor.  For the 7-products, this reduces it to a 5-product.  It doesn't help the 5-products that have no 4-velocity to absorb into them.  This would be overall more similar to using the 4-field but where I've analytically canceled that bad \gamma^2 term (or after division by \gamma^2 is a bad order unity term where other terms are divided by \gamma of some power). */

/* My question is how bad do you think this is to treat the 3-velocity instead of the 4-velocity?  As you know the 3-velocity is not used in HARM for the *primitive* interpolation because at the interface the metric has changed and the 3-velocity there can easily imply no solution (i.e. \gamma is complex).  That's certainly bad. */

/* However, the flux correction is local, so the metric is unchanged.  That is, even though we use a distribution in all quantities, we only actually use their dependence at the location we get the correction.  The 0th order flux will be no different since the face values are as from interpolating the primitives as however we choose them and normally involve good 4-velocities. */

/* I think it would be bad to group, e.g., u.B or B^2 together for the same reason that directly interpolating u.B or B^2 is a bad idea as we know.  At least for B^2, when the sign changes in B but B^2 is constant, then the interpolation wouldn't be aware of this. */

/* However, at the same time, if B^2 is directly manifested in the EOM, then any result of it won't be affected by that sign issue.  The only concern then is perhaps that B^2 is treated to a lower accuracy as a parabola.  It may be that this B^2 is exactly the thing that creates the error when B^2/\rho_0\gg 1, so I'm quite weary of grouping B's together. */

/* Any thoughts?  Perhaps we can chat for less than an hour some time if you have any thoughts since typing is not a good idea. */







// local declarations
static int deconvolve_flux(int dir, int odir1, int odir2, FTYPE *EOSextra, struct of_state *stateleft, struct of_state *statecent, struct of_state *stateright, FTYPE *Fleft, FTYPE *Fright);
static int deconvolve_flux_ma(int dir, int odir1, int odir2, FTYPE *EOSextra, struct of_state *stateleft, struct of_state *statecent, struct of_state *stateright, FTYPE *Fleft, FTYPE *Fright);
static int deconvolve_flux_em(int dir, int odir1, int odir2, FTYPE *EOSextra, struct of_state *stateleft, struct of_state *statecent, struct of_state *stateright, FTYPE *Fleft, FTYPE *Fright);

static int onetermdeconvolution(
                                FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                ,FTYPE *Fleft, FTYPE *Fright
                                );

static int twotermdeconvolution(
                                FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                ,FTYPE udirl, FTYPE udirc, FTYPE udirr
                                ,FTYPE *Fleft, FTYPE *Fright
                                );

static int threetermdeconvolution(
                                  FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                  ,FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                  ,FTYPE udirl, FTYPE udirc, FTYPE udirr
                                  ,FTYPE *Fleft, FTYPE *Fright
                                  );

static int fourtermdeconvolution(
                                 FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                 ,FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                 ,FTYPE udirl, FTYPE udirc, FTYPE udirr
                                 ,FTYPE udnul, FTYPE udnuc, FTYPE udnur
                                 ,FTYPE *Fleft, FTYPE *Fright
                                 );

static int fivetermdeconvolution(
                                 FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                 ,FTYPE Bconl, FTYPE Bconc, FTYPE Bconr
                                 ,FTYPE Bcovl, FTYPE Bcovc, FTYPE Bcovr
                                 ,FTYPE uu0l, FTYPE uu0c, FTYPE uu0r
                                 ,FTYPE ud0l, FTYPE ud0c, FTYPE ud0r
                                 ,FTYPE *Fleft, FTYPE *Fright
                                 );

static int sixtermdeconvolution(
                                FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                ,FTYPE Bconl, FTYPE Bconc, FTYPE Bconr
                                ,FTYPE Bcovl, FTYPE Bcovc, FTYPE Bcovr
                                ,FTYPE uu0l, FTYPE uu0c, FTYPE uu0r
                                ,FTYPE ud0l, FTYPE ud0c, FTYPE ud0r
                                ,FTYPE udnul, FTYPE udnuc, FTYPE udnur
                                ,FTYPE *Fleft, FTYPE *Fright
                                );


static int seventermdeconvolution(
                                  FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                  ,FTYPE Bconl, FTYPE Bconc, FTYPE Bconr
                                  ,FTYPE Bcovl, FTYPE Bcovc, FTYPE Bcovr
                                  ,FTYPE uu0l, FTYPE uu0c, FTYPE uu0r
                                  ,FTYPE ud0l, FTYPE ud0c, FTYPE ud0r
                                  ,FTYPE udnul, FTYPE udnuc, FTYPE udnur
                                  ,FTYPE uunul, FTYPE uunuc, FTYPE uunur
                                  ,FTYPE *Fleft, FTYPE *Fright
                                  );

static int twoterm2Ddeconvolution(
                                  FTYPE Bconc,FTYPE Bconld, FTYPE Bconrd, FTYPE Bconlu, FTYPE Bconru, FTYPE Bconl, FTYPE Bconr, FTYPE Bcond, FTYPE Bconu
                                  ,FTYPE vconc,FTYPE vconld, FTYPE vconrd, FTYPE vconlu, FTYPE vconru, FTYPE vconl, FTYPE vconr, FTYPE vcond, FTYPE vconu
                                  ,FTYPE *Fld, FTYPE *Frd, FTYPE *Flu, FTYPE *Fru
                                  );


static int threeterm2Ddeconvolution(
                                    FTYPE gdetc,FTYPE gdetld, FTYPE gdetrd, FTYPE gdetlu, FTYPE gdetru, FTYPE gdetl, FTYPE gdetr, FTYPE gdetd, FTYPE gdetu
                                    ,FTYPE Bconc,FTYPE Bconld, FTYPE Bconrd, FTYPE Bconlu, FTYPE Bconru, FTYPE Bconl, FTYPE Bconr, FTYPE Bcond, FTYPE Bconu
                                    ,FTYPE vconc,FTYPE vconld, FTYPE vconrd, FTYPE vconlu, FTYPE vconru, FTYPE vconl, FTYPE vconr, FTYPE vcond, FTYPE vconu
                                    ,FTYPE *Fld, FTYPE *Frd, FTYPE *Flu, FTYPE *Fru
                                    );



static int a2cflux_from_prim(int dir, FTYPE (*prim_coef_list)[MAXSPACEORDER]);

static int deconvolve_emf_1d(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
static int deconvolve_emf_2d(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);




//TODO:
//1) Setup fully 3D SPC setup:
//a) Just average (somehow) flux, using \delta term that is larger, not just average?
//b) Test in 1D with (e.g.) shock.


// NOTE: When called with merged method, the flux loop over which compute_and_store_fluxstate() is called should have been over slightly larger range to store the full cell's left-right values (e.g. for i=-1 we compute interpolation and normally only need pright -> p_l[0] -> F[0], but for merged method we also need pleft->p_r[-1] to do a2c -> F[0] as well)


// TODO:
// 1) Worry about boundary edges and that don't store entire cell interpolation at very edges where only previously needed flux at face.  Might need to store more p_l, p_r, pl_ct, pr_ct, and state info
//    a) For stag method also need to *compute* additional interpolation to other corners at outer-most grid points not where EMF is computed but part of same cell
//    b) For non-stag method, already compute interpolation to outermost face not where computing flux, but don't store or compute state for it yet
// 2) Work out rest of EM stress-energy tensor
// 3) Work out EMF sorting issue
// 4) 



// get corner \detg for EMF calculation when doing merged method
// Only needed when CORNGDETVERSION==1, since otherwise \detg already put into field that goes into EMF
void store_geomcorn(int corner, int odir1, int odir2,FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int is,ie,js,je,ks,ke,di,dj,dk;
  int i,j,k;
  struct of_gdetgeom geomcorndontuse;
  struct of_gdetgeom *ptrgeomcorn=&geomcorndontuse;



#if(MERGEDC2EA2CMETHODEM&&FIELDSTAGMEM)
  // then supposed to be here
#else
  dualfprintf(fail_file,"Why doing store_geomcorn()?\n");
  myexit(3468346);
#endif





  // define loop range
  set_interppoint_loop_ranges_geomcorn_formerged(ENOINTERPTYPE, corner, odir1, odir2, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);

  // loop over all i,j,k
  COMPZSLOOP(is,ie,js,je,ks,ke){
      
    // assume ptrgeom->g all that's needed, not ptrgeom->EOMFUNCMAC(pl) since for field ptrgeom->EOMFUNCMAC's must be consistent!
    get_geometry_gdetonly(i, j, k, CORN1-1+corner, ptrgeomcorn); // at CORN[dir]
      
    // then need to store geometry for merged method
    MACP1A0(geomcorn,corner,i,j,k)=ptrgeomcorn->gdet; // SUPERGODMARK: Should be avoided since already stored metric if doing new method
  }// end loop over i,j,k


}








// used for flux.mergedc2aa2cmethod.c to obtain values necessary for 2D deconvolution
// return \detg B^i always, instead of ever B^i alone
// Note velocity needs to be 3-velocity corresponding to 3-field
// This means v^i = u^i/u^t in lab-frame and B^i = *F^{it}
// vcon1/2 and Bcon1/2 correspond to vcon[odir1/odir2] and Bcon[odir1/odir2] referenced from dir=corner
// below left/right are as for emf with NUMCS data size.  Here added CENT4EMF for centered position
//
// NEWMARK: Need to work out signature issues
//
int setup_9value_vB(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k,
                    //      FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                    FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                    struct of_state (*fluxstatecent)[NSTORE2][NSTORE3],
                    struct of_state (*fluxstate)[NSTORE1][NSTORE2][NSTORE3][NUMLEFTRIGHT],
                    FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                    FTYPE (*vcon)[NUMCS+1][NUMCS+1],
                    FTYPE (*gdetBcon)[NUMCS+1][NUMCS+1]
                    )
{
  int ileft, jleft, kleft;
  int irightodir1, jrightodir1, krightodir1;
  int irightodir2, jrightodir2, krightodir2;
  int ijkdimen,ijkcorn[NDIM][NUMCS+1][NUMCS+1];
  FTYPE localgeomcorn;
  int m,l;
  int ii,jj,kk;


#if(MERGEDC2EA2CMETHOD)

  // face and cent values
  // MACP1A1(fluxstate,dimen,i,j,k,ISLEFT)
  // MAC(fluxstatecent,i,j,k)
  // MACP1A1(fluxstate,dimen,i,j,k,ISRIGHT)


  // OPTMARK: Should compute these things outside ijk loop and pass to this function
  // left-flux position (centered or left in odir1 and/or odir2)
  ileft=i;
  jleft=j;
  kleft=k;
  
  // right-flux position for odir1 direction and centered in odir2 dir
  irightodir1=i+(1==odir1)*N1NOT1;
  jrightodir1=j+(2==odir1)*N2NOT1;
  krightodir1=k+(3==odir1)*N3NOT1;

  // right-flux position for odir2 direction and centered in odir1 dir
  irightodir2=i+(1==odir2)*N1NOT1;
  jrightodir2=j+(2==odir2)*N2NOT1;
  krightodir2=k+(3==odir2)*N3NOT1;

  // +-odir1 and +-odir2 w.r.t. cell center
  // set default i,j,k position
  //  DIMENLOOP(ijkdimen){
  for(m=0;m<NUMCS;m++){
    for(l=0;l<NUMCS;l++){
      ijkcorn[1][m][l]=i;
      ijkcorn[2][m][l]=j;
      ijkcorn[3][m][l]=k;
    }
  }
  //}

  // shift i,j,k up only as required
  // shift each
  ijkcorn[odir1][RIGHT4EMF][LEFT4EMF]+=NNOT1vec[odir1];
  ijkcorn[odir2][LEFT4EMF][RIGHT4EMF]+=NNOT1vec[odir2];
  // shift both
  ijkcorn[odir1][RIGHT4EMF][RIGHT4EMF]+=NNOT1vec[odir1];
  ijkcorn[odir2][RIGHT4EMF][RIGHT4EMF]+=NNOT1vec[odir2];


  /////////////////////
  //
  // lab-frame 3-velocity: v^i
  // lab-frame 3-field: B^i
  //
  /////////////////////

  // cent-cent
  vcon[odir1][CENT4EMF][CENT4EMF]=MAC(fluxstatecent,i,j,k).vcon[odir1]; // v^{odir1} = u^{odir1}/u^t
  vcon[odir2][CENT4EMF][CENT4EMF]=MAC(fluxstatecent,i,j,k).vcon[odir2]; // v^{odir2} = u^{odir2}/u^t

  gdetBcon[odir1][CENT4EMF][CENT4EMF]=MAC(fluxstatecent,i,j,k).gdetBcon[odir1]; // \detg B^{odir1}
  gdetBcon[odir2][CENT4EMF][CENT4EMF]=MAC(fluxstatecent,i,j,k).gdetBcon[odir2]; // \detg B^{odir2}

  // left-center (note: As in flux.mergedc2ea2cmethod.c, position of left value is stored in local cell's ISRIGHT)
  vcon[odir1][LEFT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,ileft,jleft,kleft,ISRIGHT).vcon[odir1]; // v^{odir1} = u^{odir1}/u^t
  vcon[odir2][LEFT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,ileft,jleft,kleft,ISRIGHT).vcon[odir2]; // v^{odir2} = u^{odir2}/u^t

  gdetBcon[odir1][LEFT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,ileft,jleft,kleft,ISRIGHT).gdetBcon[odir1]; // // \detg B^{odir1}
  gdetBcon[odir2][LEFT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,ileft,jleft,kleft,ISRIGHT).gdetBcon[odir2]; // // \detg B^{odir2}

  // right-center
  vcon[odir1][RIGHT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,irightodir1,jrightodir1,krightodir1,ISLEFT).vcon[odir1]; // v^{odir1} = u^{odir1}/u^t
  vcon[odir2][RIGHT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,irightodir1,jrightodir1,krightodir1,ISLEFT).vcon[odir2]; // v^{odir2} = u^{odir2}/u^t

  gdetBcon[odir1][RIGHT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,irightodir1,jrightodir1,krightodir1,ISLEFT).gdetBcon[odir1];
  gdetBcon[odir2][RIGHT4EMF][CENT4EMF]=MACP1A1(fluxstate,odir1,irightodir1,jrightodir1,krightodir1,ISLEFT).gdetBcon[odir2];

  // center-left
  vcon[odir1][CENT4EMF][LEFT4EMF]=MACP1A1(fluxstate,odir2,ileft,jleft,kleft,ISRIGHT).vcon[odir1];
  vcon[odir2][CENT4EMF][LEFT4EMF]=MACP1A1(fluxstate,odir2,ileft,jleft,kleft,ISRIGHT).vcon[odir2];

  gdetBcon[odir1][CENT4EMF][LEFT4EMF]=MACP1A1(fluxstate,odir2,ileft,jleft,kleft,ISRIGHT).gdetBcon[odir1];
  gdetBcon[odir2][CENT4EMF][LEFT4EMF]=MACP1A1(fluxstate,odir2,ileft,jleft,kleft,ISRIGHT).gdetBcon[odir2];

  // right-center
  vcon[odir1][CENT4EMF][RIGHT4EMF]=MACP1A1(fluxstate,odir2,irightodir2,jrightodir2,krightodir2,ISLEFT).vcon[odir1];
  vcon[odir2][CENT4EMF][RIGHT4EMF]=MACP1A1(fluxstate,odir2,irightodir2,jrightodir2,krightodir2,ISLEFT).vcon[odir2];

  gdetBcon[odir1][CENT4EMF][RIGHT4EMF]=MACP1A1(fluxstate,odir2,irightodir2,jrightodir2,krightodir2,ISLEFT).gdetBcon[odir1];
  gdetBcon[odir2][CENT4EMF][RIGHT4EMF]=MACP1A1(fluxstate,odir2,irightodir2,jrightodir2,krightodir2,ISLEFT).gdetBcon[odir2];



  //////////////////
  //
  // CORNER lab-frame 3-velocity: v^i and B^i
  //
  //////////////////
  // pvcorn[which corner][which component in pl form][+-odir1][+-odir2]
  // pbcorn[which corner][which component in pl form][+-remaining direction that is not corn nor pl-dir]
  // for example, emf3[+-x][+-y] = By[+-x]*vx[+-x][+-y] - Bx[+-y]*vy[+-x][+-y]


  // -odir1 -odir2 w.r.t. cell center i,j,k
  ii=ijkcorn[1][LEFT4EMF][LEFT4EMF];
  jj=ijkcorn[2][LEFT4EMF][LEFT4EMF];
  kk=ijkcorn[3][LEFT4EMF][LEFT4EMF];
  vcon[odir1][LEFT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,RIGHT4EMF,RIGHT4EMF);
  vcon[odir2][LEFT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,RIGHT4EMF,RIGHT4EMF);

  // +odir1 -odir2
  ii=ijkcorn[1][RIGHT4EMF][LEFT4EMF];
  jj=ijkcorn[2][RIGHT4EMF][LEFT4EMF];
  kk=ijkcorn[3][RIGHT4EMF][LEFT4EMF];
  vcon[odir1][RIGHT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,LEFT4EMF,RIGHT4EMF);
  vcon[odir2][RIGHT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,LEFT4EMF,RIGHT4EMF);

  // -odir1 +odir2
  ii=ijkcorn[1][LEFT4EMF][RIGHT4EMF];
  jj=ijkcorn[2][LEFT4EMF][RIGHT4EMF];
  kk=ijkcorn[3][LEFT4EMF][RIGHT4EMF];
  vcon[odir1][LEFT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,RIGHT4EMF,LEFT4EMF);
  vcon[odir2][LEFT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,RIGHT4EMF,LEFT4EMF);

  // +odir1 +odir2
  ii=ijkcorn[1][RIGHT4EMF][RIGHT4EMF];
  jj=ijkcorn[2][RIGHT4EMF][RIGHT4EMF];
  kk=ijkcorn[3][RIGHT4EMF][RIGHT4EMF];
  vcon[odir1][RIGHT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,LEFT4EMF,LEFT4EMF);
  vcon[odir2][RIGHT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,LEFT4EMF,LEFT4EMF);

  //////////////////
  //
  // CORNER lab-frame 3-field: B^i
  //
  // Note pbcorn[corner][B1-1+odir1][LEFT4EMF][iodir1] = pbcorn[corner][B1-1+odir1][RIGHT4EMF][iodir1], so within an i,j,k cell, field along itself on upper edge is obtained from LEFT @ iodir+1, not right that doesn't exist
  // That is, pbcorn[LEFT4EMF/RIGHT4EMF] always refers to transverse direction to field pl-direction
  //
  //////////////////


  ////////////////////////
  //
  // first deal with geometry to be consistent with how stored pbcorn[] when doing face2corn()
  //
  ////////////////////////


  // -odir1 -odir2 w.r.t. cell center i,j,k
  ii=ijkcorn[1][LEFT4EMF][LEFT4EMF];
  jj=ijkcorn[2][LEFT4EMF][LEFT4EMF];
  kk=ijkcorn[3][LEFT4EMF][LEFT4EMF];
#if(CORNGDETVERSION==1)
  // means pbcorn has no \detg, so must put it back
  localgeomcorn=MACP1A0(geomcorn,corner,ii,jj,kk); // SUPERGODMARK: Should use stored metric at whatever location required
#else
  // means pbcorn has \detg already
  localgeomcorn=1.0;
#endif
  gdetBcon[odir1][LEFT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,NUMCS,RIGHT4EMF)*localgeomcorn;
  gdetBcon[odir2][LEFT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,NUMCS,RIGHT4EMF)*localgeomcorn;

  // +odir1 -odir2 w.r.t. cell center i,j,k
  ii=ijkcorn[1][RIGHT4EMF][LEFT4EMF];
  jj=ijkcorn[2][RIGHT4EMF][LEFT4EMF];
  kk=ijkcorn[3][RIGHT4EMF][LEFT4EMF];
#if(CORNGDETVERSION==1)
  // means pbcorn has no \detg, so must put it back
  localgeomcorn=MACP1A0(geomcorn,corner,ii,jj,kk);
#else
  // means pbcorn has \detg already
  localgeomcorn=1.0;
#endif
  gdetBcon[odir1][RIGHT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,NUMCS,RIGHT4EMF)*localgeomcorn;
  gdetBcon[odir2][RIGHT4EMF][LEFT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,NUMCS,LEFT4EMF)*localgeomcorn;

  // -odir1 +odir2 w.r.t. cell center i,j,k
  ii=ijkcorn[1][LEFT4EMF][RIGHT4EMF];
  jj=ijkcorn[2][LEFT4EMF][RIGHT4EMF];
  kk=ijkcorn[3][LEFT4EMF][RIGHT4EMF];
#if(CORNGDETVERSION==1)
  // means pbcorn has no \detg, so must put it back
  localgeomcorn=MACP1A0(geomcorn,corner,ii,jj,kk);
#else
  // means pbcorn has \detg already
  localgeomcorn=1.0;
#endif
  gdetBcon[odir1][LEFT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,NUMCS,LEFT4EMF)*localgeomcorn;
  gdetBcon[odir2][LEFT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,NUMCS,RIGHT4EMF)*localgeomcorn;

  // +odir1 +odir2 w.r.t. cell center i,j,k
  ii=ijkcorn[1][RIGHT4EMF][RIGHT4EMF];
  jj=ijkcorn[2][RIGHT4EMF][RIGHT4EMF];
  kk=ijkcorn[3][RIGHT4EMF][RIGHT4EMF];
#if(CORNGDETVERSION==1)
  // means pbcorn has no \detg, so must put it back
  localgeomcorn=MACP1A0(geomcorn,corner,ii,jj,kk);
#else
  // means pbcorn has \detg already
  localgeomcorn=1.0;
#endif
  gdetBcon[odir1][RIGHT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir1,NUMCS,LEFT4EMF)*localgeomcorn;
  gdetBcon[odir2][RIGHT4EMF][RIGHT4EMF]=MACP1A3(pvbcorn,corner,ii,jj,kk,odir2,NUMCS,LEFT4EMF)*localgeomcorn;



#endif// end if merged method


  return(0);

}




// used to obtain 1-D line in case where EMF only varies in one transverse direction rather than 2 transverse directions
// OPTMARK: for now don't optimize, just get full 9 value and return along needed line
int setup_3value_vB(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k,
                    //      FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                    FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                    struct of_state (*fluxstatecent)[NSTORE2][NSTORE3],
                    struct of_state (*fluxstate)[NSTORE1][NSTORE2][NSTORE3][NUMLEFTRIGHT],
                    FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                    FTYPE (*vcon)[NUMCS+1],
                    FTYPE (*gdetBcon)[NUMCS+1]
                    )
{
  int setup_9value_vB(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k,
                      //        FTYPE (*pbcorn)[COMPDIM][NUMCS][NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                      FTYPE (*pvbcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3][COMPDIM][NUMCS+1][NUMCS],
                      struct of_state (*fluxstatecent)[NSTORE2][NSTORE3],
                      struct of_state (*fluxstate)[NSTORE1][NSTORE2][NSTORE3][NUMLEFTRIGHT],
                      FTYPE (*geomcorn)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],
                      FTYPE (*vcon)[NUMCS+1][NUMCS+1],
                      FTYPE (*gdetBcon)[NUMCS+1][NUMCS+1]
                      );
  FTYPE vcon9[NDIM][NUMCS+1][NUMCS+1];
  FTYPE gdetBcon9[NDIM][NUMCS+1][NUMCS+1];
  int jj;
  int positer;
  int refodir1,refodir2;
  int odir1pos,odir2pos;


  setup_9value_vB(corner, odir1, odir2 , Nvec, NNOT1vec, i,j,k, pvbcorn,fluxstatecent,fluxstate,geomcorn, vcon9 ,gdetBcon9 );
        
  // choose starting position for positer loop
  //  This chooses relevant odir direction for positer loop below
  // use of CENT4EMF as opposed to other positions is not relevant if using setup_9value_vB()
  refodir1=NNOT1vec[odir1]*START4EMF + (1-NNOT1vec[odir1])*CENT4EMF;
  refodir2=NNOT1vec[odir2]*START4EMF + (1-NNOT1vec[odir2])*CENT4EMF;


  SLOOPA(jj){ // OPTMARK: not all jj are needed, just the odir1 and odir2 components, so could somehow skip jj==corner if it avoided conditional and that might be better
    for(positer=0;positer<NUMPOS4EMF;positer++){
      
      odir1pos=refodir1+positer*NNOT1vec[odir1];
      odir2pos=refodir2+positer*NNOT1vec[odir2];
      
      vcon[jj][positer]=vcon9[jj][odir1pos][odir2pos];
      gdetBcon[jj][positer]=gdetBcon9[jj][odir1pos][odir2pos];
    }
  }



  return(0);

}








// called from flux.c from fluxcalc() after all dimensions fluxes are computed and need to correct
void mergedc2ea2cmethod_compute(int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int dimen;
  int corner;
  FTYPE Fleft[NPR],Fright[NPR];
  int is,ie,js,je,ks,ke,di,dj,dk;
  int i,j,k;
  int ileft,jleft,kleft;
  int iright,jright,kright;
  int pl,pliter;
  int odir1,odir2;
  int odimen1,odimen2;
  int NNOT1vec[NDIM];



  // setup which directions are not 1 in size
  NNOT1vec[0]=0;
  NNOT1vec[1]=N1NOT1;
  NNOT1vec[2]=N2NOT1;
  NNOT1vec[3]=N3NOT1;


  // avoid changing EMFs if doing FLUXB==FLUXCTTOTH since this merged call is after the FLUXCTTOTH correction in flux.c
  // No correction if Nvec[dimen]==1

  // We loop over i,j,k corresponding to center of cell that has primitive distribution within it that will be used to update the conserved quantity at CENT
  // fluxes are corrected on edges of this box



  // NOTE: Uses global variables fluxstate and fluxstatecent
  DIMENLOOP(dimen){
    if(Nvec[dimen]>1){

      get_odirs(dimen,&odimen1,&odimen2);

      // define loop range
      // loop range should be same as where computed fluxes, or where formed primitive interpolation, as in interppoint.c:

      // set range of positions interpolated to (di,dj,dk useless)

      int loc=CENT;
      int continuous=0;
      set_interppoint_loop_ranges(ENOINTERPTYPE, dimen, loc, continuous, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);

      // loop over all i,j,k
      COMPZSLOOP(is,ie,js,je,ks,ke){

        // initialize Fleft and Fright updates
        PLOOP(pliter,pl){
          Fleft[pl]=Fright[pl]=0.0;
        }
 
        // left-flux position
        ileft=i;
        jleft=j;
        kleft=k;
  
        // right-flux position
        iright=i+(1==dimen)*N1NOT1;
        jright=j+(2==dimen)*N2NOT1;
        kright=k+(3==dimen)*N3NOT1;



        // deconvolve flux
        // Note meaning of Fleft and Fright is not same as F_l and F_r:
        //
        // |                   F_l|F_r                                |
        // |         fs[ISLEFT][i]|fs[ISRIGHT][i]      fs[ISLEFT][i+1]|
        // |                      |                i                  |
        // |                      |                                   |
        // |                      |Fleft                        Fright|
        // |                      |fluxvec                            |
        // |                      |              EOSextra             |
        //
        // Assume EOSextra can be DONOR accurate
        // And note that F_l=fluxstate[ISLEFT] and F_r=fluxstate[ISRIGHT]
        // so single cell needs Fleft=state[ISRIGHT][i] and Fright=state[ISLEFT][i+1]
        // And note that fluxvec is single-valued at face and so naturally identified by [i] for cell [i]
        // So Fleft modifies fluxvec[i] and Fright modifies fluxvec[i+1]
        //
        deconvolve_flux(dimen, odimen1, odimen2, GLOBALMAC(EOSextraglobal,i,j,k), &GLOBALMACP1A1(fluxstate,dimen,ileft,jleft,kleft,ISRIGHT), &GLOBALMAC(fluxstatecent,i,j,k), &GLOBALMACP1A1(fluxstate,dimen,iright,jright,kright,ISLEFT), Fleft, Fright);

        // Now realize that at faces each correction is applied ultimately as average of 2 fluxes for non-EMF terms
        // ultimately if smooth, then is actually full correction since both sides of face will contribute same correction
        // Note that using interpolation loop means apply correction beyond where needed, but this is ok as long as 2 boundary condition cells (for upper edge) and 1 for lower edge

        // correct flux on left part of cell
        PLOOP(pliter,pl) MACP1A1(fluxvec,dimen,ileft,jleft,kleft,pl) += 0.5*Fleft[pl];

        // correct flux on right part of cell
        PLOOP(pliter,pl) MACP1A1(fluxvec,dimen,iright,jright,kright,pl) += 0.5*Fright[pl];

      }// end over i,j,k
    }// end if dimension exists
  }// end over dimensions











  // TODO: GODMARK: Not done --
  if(MERGEDC2EA2CMETHODEM&&(FLUXB==FLUXCTSTAG)){ // STAG has information at corners, while FLUXCTTOTH and FLUXCTHLL methods don't
    // Don't apply EMF correct if doing Toth method since inconsistent with Toth averaging



    // NOTE: Uses global variables fluxstate and fluxstatecent
    DIMENLOOP(corner){

      get_odirs(corner,&odir1,&odir2);
      if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){
        // only here if doing correction with non-zero (non-cancelling) EMF


        if(CORNGDETVERSION==1){
          // need to get corner geometry over all domain before processing EMF
          // OPTMARK: if metric is changing, then need to do this every time, else could avoid
          store_geomcorn(corner,odir1,odir2,GLOBALPOINT(geomcornglobal));
        }

        // set range of positions interpolated to (di,dj,dk useless)
        // loop range should be over all cell centers that have a final flux on one face (upper or lower) that needs to be corrected
        set_interppoint_loop_ranges_2D_EMF_formerged(ENOINTERPTYPE, corner, odir1, odir2, &is, &ie, &js, &je, &ks, &ke, &di, &dj, &dk);

  

        // deconvolve flux
        if(Nvec[odir1]==1 || Nvec[odir2]==1){
          // then only need to do 1D deconvolution in odir2 direction
          // e.g. if corner=CORN2, then odir1=3 and odir2=1 and in 2D(N1>1 and N2>1 and N3==1), then correction of E2 in odir2=1 needs to be done.
          // e.g. cont. In this limit, deconvolution is 1D and takes place using values at FACE1(i), CENT(i), FACE1(i+1) since no offset in odir1=3 direction as would be true if N3>1
          //
          // loop over all i,j,k
          COMPZSLOOP(is,ie,js,je,ks,ke){
            deconvolve_emf_1d(corner, odir1, odir2, Nvec, NNOT1vec, i,j,k, fluxvec);
          }
        }
        else{
          // full 2D deconvolution
          // correction done inside
          // loop over all i,j,k
          COMPZSLOOP(is,ie,js,je,ks,ke){
            deconvolve_emf_2d(corner, odir1, odir2, Nvec, NNOT1vec, i,j,k, fluxvec);
          }
        }
 
      }// end if emf has transverse direction to correct
    }// end over corners
  }

}


// should operate on a single dimension
// can turn on/off MA or EM terms separately
static int deconvolve_flux(int dir, int odir1, int odir2, FTYPE *EOSextra, struct of_state *stateleft, struct of_state *statecent, struct of_state *stateright, FTYPE *Fleft, FTYPE *Fright)
{

#if(MERGEDC2EA2CMETHODMA)
  // do matter parts first
  deconvolve_flux_ma(dir, odir1, odir2, EOSextra, stateleft, statecent, stateright, Fleft, Fright);
#endif

#if(MERGEDC2EA2CMETHODEM)
  // then do EM parts
  deconvolve_flux_em(dir, odir1, odir2, EOSextra, stateleft, statecent, stateright, Fleft, Fright);
#endif


  return(0);
}


// should operate on a single dimension
static int deconvolve_flux_ma(int dir, int odir1, int odir2, FTYPE *EOSextra, struct of_state *stateleft, struct of_state *statecent, struct of_state *stateright, FTYPE *Fleft, FTYPE *Fright)
{
  FTYPE gdetl,gdetc,gdetr;
  FTYPE rhol,rhoc,rhor;
  FTYPE iel,iec,ier;
  FTYPE pressl,pressc,pressr;
  FTYPE udirl,udirc,udirr;
  FTYPE udnul,udnuc,udnur;
  FTYPE Flefttemp,Frighttemp;
  FTYPE Fpressleft,Fpressright;
  int jj;
  FTYPE genucovl[NDIM];
  FTYPE genucovc[NDIM];
  FTYPE genucovr[NDIM];

  // use fluxstate, fluxstatecent per cell per dimension to get corrections to flux for each dimension and each MA-type of flux
  // use e[pl] for \gdet for each equation of motion

  // F^i_\rho = \rho u^i
  // F_E^i = (\rho+u+p)u^i u_t
  // F_{\rho v^j}^i = (\rho+u+p)u^i u_j + p\delta^i_j  (note F_{\rho v^i}^i = F_E^i u_i/u_t + p
  //   generally, F_{\rho v^j}^i = F_E^i u_i/u_t + p \delta^i_j


  ////////////////////////
  //
  // Setup things common to all EOMs
  //
  ////////////////////////

#if(MERGEDC2EA2CMETHOD)
  // \detg for \rho_0
  // overridden by below if WHICHEOM!=WITHGDET
  gdetl=stateleft->gdet;
  gdetc=statecent->gdet;
  gdetr=stateright->gdet;


  // \rho_0
  rhol=stateleft->prim[RHO];
  rhoc=statecent->prim[RHO];
  rhor=stateright->prim[RHO];

  // ie (internal energy density)
  iel=stateleft->prim[UU];
  iec=statecent->prim[UU];
  ier=stateright->prim[UU];
#endif
  // pressure
  pressl=stateleft->pressure;
  pressc=statecent->pressure;
  pressr=stateright->pressure;

  // 4-velocity in coordinate frame
  udirl=stateleft->ucon[dir];
  udirc=statecent->ucon[dir];
  udirr=stateright->ucon[dir];



  ///////////////////
  //
  // REST-MASS FLUX: \detg \rho_0 u^{dir}
  //
  ///////////////////

#if(WHICHEOM!=WITHGDET&&MERGEDC2EA2CMETHOD)
  // need gdet for each variable
  gdetl=stateleft->EOMFUNCMAC(RHO);
  gdetc=statecent->EOMFUNCMAC(RHO);
  gdetr=stateright->EOMFUNCMAC(RHO);
#endif

#if(MERGEDC2EA2CMETHOD)

#if(MCOORD==CARTMINKMETRIC)
  twotermdeconvolution(rhol, rhoc, rhor  ,udirl, udirc, udirr  ,&Flefttemp, &Frighttemp);
  Fleft[RHO]  = gdetc*Flefttemp;
  Fright[RHO] = gdetc*Frighttemp;
#else
  threetermdeconvolution(gdetl, gdetc, gdetr,   rhol, rhoc, rhor  ,udirl, udirc, udirr  ,&Fleft[RHO], &Fright[RHO]);
#endif

#endif


  ///////////////////
  //
  // Energy-Momentum FLUX: \detg [ (\rho_0 + u + P)u^i u_\nu + \delta^i_\nu P ]
  //
  // Note: \rho_0+u+P term is linear so no deconvolution correction, so can just sum those up *before* deconvolution and result will be same
  //
  //
  // when REMOVERESTMASSFROMUU==1,2 then really have for TT term:
  //
  // \detg [ (\rho_0 [1+u_t] + (u+p)u_t)u^i + \delta^i_\nu P ]
  // so in this case have to break \rho away from u+p term
  //
  //
  //
  ///////////////////


  /////////////
  //
  // setup 1+u_t or u_t depending upon REMOVERESTMASSFROMUU==0 or 1,2
  // by the time we reach here, UtoU() has operated on flux and if REMOVERESTMASSFROMUU>0 then final flux in UEVOLVE form has no rest-mass term
  //
  //////////////
  SLOOPA(jj){
    genucovl[jj]=stateleft->ucov[jj];
    genucovc[jj]=statecent->ucov[jj];
    genucovr[jj]=stateright->ucov[jj];
  }
  jj=TT;
  genucovl[jj]=stateleft->ifremoverestplus1ud0elseud0;
  genucovc[jj]=statecent->ifremoverestplus1ud0elseud0;
  genucovr[jj]=stateright->ifremoverestplus1ud0elseud0;



  ///////////////
  //
  // Loop over T^{dir}_\mu type term
  //
  ///////////////

  DLOOPA(jj){

#if(WHICHEOM!=WITHGDET&&MERGEDC2EA2CMETHOD)
    // need gdet for each variable
    gdetl=stateleft->EOMFUNCMAC(UU+jj);
    gdetc=statecent->EOMFUNCMAC(UU+jj);
    gdetr=stateright->EOMFUNCMAC(UU+jj);
#endif



    // lower 4-velocity in coordinate frame
    udnul=stateleft->ucov[jj];
    udnuc=statecent->ucov[jj];
    udnur=stateright->ucov[jj];

    
    // initialize as we loop over terms that are simply part of sum that has no extra coupling for deconvolution
    Fleft[UU+jj]  = 0.0;
    Fright[UU+jj] = 0.0;

#if(MCOORD==CARTMINKMETRIC)


#if(REMOVERESTMASSFROMUU==1 || REMOVERESTMASSFROMUU==2)
    threetermdeconvolution(rhol, rhoc, rhor  ,udirl, udirc, udirr,  genucovl[jj], genucovc[jj], genucovr[jj]   ,&Flefttemp, &Frighttemp);
    Fleft[UU+jj] += Flefttemp;
    Fright[UU+jj] += Frighttemp;

    threetermdeconvolution(iel+pressl, iec+pressc, ier+pressr  ,udirl, udirc, udirr,  udnul, udnuc, udnur   ,&Flefttemp, &Frighttemp);
    Fleft[UU+jj] += Flefttemp;
    Fright[UU+jj] += Frighttemp;
#else
    // the sum on rho+ie+P is allowed because single term with otherwise same multiplications, each term is itself rather than some function, so deconvolution integral is same for those terms together or apart
    threetermdeconvolution(rhol+iel+pressl, rhoc+iec+pressc, rhor+ier+pressr  ,udirl, udirc, udirr,  udnul, udnuc, udnur   ,&Flefttemp, &Frighttemp);
    Fleft[UU+jj] += Flefttemp;
    Fright[UU+jj] += Frighttemp;
#endif

    if(jj==dir){
      onetermdeconvolution(pressl, pressc, pressr   ,&Flefttemp, &Frighttemp);
      Fleft[UU+jj] += Flefttemp;
      Fright[UU+jj] += Frighttemp;      
    }

    // finally apply constant \detg factor
    Fleft[UU+jj]*=gdetc;
    Fright[UU+jj]*=gdetc;

#else


#if(REMOVERESTMASSFROMUU==1 || REMOVERESTMASSFROMUU==2)

    fourtermdeconvolution(gdetl, gdetc, gdetr,  rhol, rhoc, rhor  ,udirl, udirc, udirr,  genucovl[jj], genucovc[jj], genucovr[jj]   ,&Flefttemp, &Frighttemp);
    Fleft[UU+jj] += Flefttemp;
    Fright[UU+jj] += Frighttemp;

    fourtermdeconvolution(gdetl, gdetc, gdetr,  iel+pressl, iec+pressc, ier+pressr  ,udirl, udirc, udirr,  udnul, udnuc, udnur   ,&Flefttemp, &Frighttemp);
    Fleft[UU+jj] += Flefttemp;
    Fright[UU+jj] += Frighttemp;

#else
    // the sum on rho+ie+P is allowed because single term with otherwise same multiplications, each term is itself rather than some function, so deconvolution integral is same for those terms together or apart
    fourtermdeconvolution(gdetl, gdetc, gdetr,  rhol+iel+pressl, rhoc+iec+pressc, rhor+ier+pressr  ,udirl, udirc, udirr,  udnul, udnuc, udnur   ,&Flefttemp, &Frighttemp);
    Fleft[UU+jj] += Flefttemp;
    Fright[UU+jj] += Frighttemp;
#endif

    if(jj==dir){
      twotermdeconvolution(gdetl, gdetc, gdetr,  pressl, pressc, pressr   ,&Flefttemp, &Frighttemp);
      Fleft[UU+jj] += Flefttemp;
      Fright[UU+jj] += Frighttemp;      
    }

#endif
  }// end loop over jj = \nu={0,1,2,3}



  //////////
  //
  // Setup radiation terms
  //
  //////////


  // KORALTODO

  //////////
  //
  // setup Y_L and Y_\nu terms
  //
  //////////

#if(MERGEDC2EA2CMETHOD)
#if(DOYNU!=DONOYNU)

  FTYPE YNUl,YNUc,YNUr;

  // Y_\nu:
  YNUl=stateleft->prim[YNU];
  YNUc=statecent->prim[YNU];
  YNUr=stateright->prim[YNU];
#endif
#endif



#if(MERGEDC2EA2CMETHOD)
#if(DOYL!=DONOYL)

  FTYPE YLl,YLc,YLr;

  // Y_l:
  YLl=stateleft->prim[YL];
  YLc=statecent->prim[YL];
  YLr=stateright->prim[YL];
#endif
#endif

#if(MERGEDC2EA2CMETHOD) // SUPERGODMARK: none done for YFLx
#if(DOYFL!=0)

  FTYPE YFLl,YFLc,YFLr;

  // Y_l:
  YFLl=stateleft->prim[YFL1];
  YFLc=statecent->prim[YFL1];
  YFLr=stateright->prim[YFL1];
#endif
#endif



  ///////////////////
  //
  // Flux of Y_fl : \detg \rho_0 Y_fl u^i
  //
  ///////////////////

#if(MERGEDC2EA2CMETHOD)
#if(DOYFL==1)

  FTYPE yfl2advectl,yfl2advectc,yfl2advectr;

#if(WHICHEOS==KAZFULL)
  yfl2advect_kazfull(EOSextra,YFLl,YNUl,&yfl2advectl);
  yfl2advect_kazfull(EOSextra,YFLc,YNUc,&yfl2advectc);
  yfl2advect_kazfull(EOSextra,YFLr,YNUr,&yfl2advectr);
#else
  yfl2advectl=YFLl;
  yfl2advectc=YFLc;
  yfl2advectr=YFLr;
#endif
  

#if(MCOORD==CARTMINKMETRIC)
  threetermdeconvolution(rhol, rhoc, rhor  ,udirl, udirc, udirr,  yfl2advectl, yfl2advectc, yfl2advectr   ,&Fleft[YFL1], &Fright[YFL1]);

  // Apply constant \detg factor
  Fleft[YFL1]*=gdetc;
  Fright[YFL1]*=gdetc;

#else
  fourtermdeconvolution(gdetl, gdetc, gdetr,  rhol, rhoc, rhor  ,udirl, udirc, udirr,  yfl2advectl, yfl2advectc, yfl2advectr   ,&Fleft[YFL1], &Fright[YFL1]);
#endif

#endif // end if doing YFL
#endif


  ///////////////////
  //
  // Flux of Y_l : \detg \rho_0 Y_l u^i
  //
  ///////////////////

#if(MERGEDC2EA2CMETHOD)
#if(DOYL!=DONOYL)

  FTYPE yl2advectl,yl2advectc,yl2advectr;

#if(WHICHEOS==KAZFULL)
  yl2advect_kazfull(EOSextra,YLl,YNUl,&yl2advectl);
  yl2advect_kazfull(EOSextra,YLc,YNUc,&yl2advectc);
  yl2advect_kazfull(EOSextra,YLr,YNUr,&yl2advectr);
#else
  yl2advectl=YLl;
  yl2advectc=YLc;
  yl2advectr=YLr;
#endif
  

#if(MCOORD==CARTMINKMETRIC)
  threetermdeconvolution(rhol, rhoc, rhor  ,udirl, udirc, udirr,  yl2advectl, yl2advectc, yl2advectr   ,&Fleft[YL], &Fright[YL]);

  // Apply constant \detg factor
  Fleft[YL]*=gdetc;
  Fright[YL]*=gdetc;

#else
  fourtermdeconvolution(gdetl, gdetc, gdetr,  rhol, rhoc, rhor  ,udirl, udirc, udirr,  yl2advectl, yl2advectc, yl2advectr   ,&Fleft[YL], &Fright[YL]);
#endif

#endif // end if doing YL
#endif



  ///////////////////
  //
  // Flux of Y_\nu : \detg \rho_0 Y_\nu u^i
  //
  ///////////////////

#if(MERGEDC2EA2CMETHOD)
#if(DOYNU!=DONOYNU)

  FTYPE ynu2advectl,ynu2advectc,ynu2advectr;

#if(WHICHEOS==KAZFULL)
  ynu2advect_kazfull(EOSextra,YLl,YNUl,&ynu2advectl);
  ynu2advect_kazfull(EOSextra,YLc,YNUc,&ynu2advectc);
  ynu2advect_kazfull(EOSextra,YLr,YNUr,&ynu2advectr);
#else
  ynu2advectl=YNUl;
  ynu2advectc=YNUc;
  ynu2advectr=YNUr;
#endif

#if(MCOORD==CARTMINKMETRIC)
  threetermdeconvolution(rhol, rhoc, rhor  ,udirl, udirc, udirr,  ynu2advectl, ynu2advectc, ynu2advectr   ,&Fleft[YNU], &Fright[YNU]);

  // Apply constant \detg factor
  Fleft[YNU]*=gdetc;
  Fright[YNU]*=gdetc;

#else
  fourtermdeconvolution(gdetl, gdetc, gdetr,  rhol, rhoc, rhor  ,udirl, udirc, udirr,  ynu2advectl, ynu2advectc, ynu2advectr   ,&Fleft[YNU], &Fright[YNU]);
#endif

#endif // end if doing YNU
#endif


  

  
  ///////////////////
  //
  // Flux of entropy : \detg (entropy) u^i
  //
  ///////////////////

#if(DOENTROPY!=DONOENTROPY)

  FTYPE entropyl,entropyc,entropyr;

  // entropy:
  entropyl=stateleft->entropy;
  entropyc=statecent->entropy;
  entropyr=stateright->entropy;

#if(MCOORD==CARTMINKMETRIC)
  twotermdeconvolution(entropyl, entropyc, entropyr   ,udirl, udirc, udirr  ,&Fleft[ENTROPY], &Fright[ENTROPY]);

  // Apply constant \detg factor
  Fleft[ENTROPY]*=gdetc;
  Fright[ENTROPY]*=gdetc;

#else
  threetermdeconvolution(gdetl, gdetc, gdetr,  entropyl, entropyc, entropyr  ,udirl, udirc, udirr,  &Fleft[ENTROPY], &Fright[ENTROPY]);
#endif
    
#endif // end if doing ENTROPY



  return(0);
}


// should operate on a single dimension
static int deconvolve_flux_em(int dir, int odir1, int odir2, FTYPE *EOSextra, struct of_state *stateleft, struct of_state *statecent, struct of_state *stateright, FTYPE *Fleft, FTYPE *Fright)
{
  FTYPE gdetl,gdetc,gdetr;
  //  FTYPE uu0l,uu0c,uu0r;
  FTYPE overutl,overutc,overutr;
  FTYPE uul[NDIM],uuc[NDIM],uur[NDIM];
  FTYPE udl[NDIM],udc[NDIM],udr[NDIM];
  FTYPE Bul[NDIM],Buc[NDIM],Bur[NDIM];
  FTYPE Bdl[NDIM],Bdc[NDIM],Bdr[NDIM];
  FTYPE Flefttemp,Frighttemp;
  int jj,kk;
  int iter,myloop[NDIM];
  int otherdiriter;
  int firstodir,secondodir;
  FTYPE fluxsign;



  // use fluxstate, fluxstatecent
  // Note that will be using pr[B1,B2,B3], and need to ensure primitive-type for magnetic field is as expected (e.g. could be \alpha B^i instead)
  // recall, use \detg B^{dir}

  // 

  //c) Can form a2c using \tilde{B}^i "directly", no need to use B^i at face or center.  No, can't use \tilde{B}^i directly since need to interpolate dB/dx and have reduction.
  //d) Using 3 dB/dx^i's -> 4 \tilde{B}^i's for parabolic interpolation, so use 1 extra \tilde{B} always
  //e) Do each cell at a time using l c r values per cell.  Pull out values and then operate.  Probably as expensive as para, which is ok
  //f) FLUXCTSTAG: Form full array \delta\tilde B^i[COMPDIM] that has to include \detg.  Store p_left, p_c, p_right for each dir=1,2,3 corresponding to each B1,B2,B3.
  //   FLUXCTSTAG: Feed pl, pc, pr into global array for use on next step (field p's since used reduction to get p at faces and cent as point value)
  //g) FLUXCTSTAG: EMF: In 2D plane have all data necessary without extra interpolations.
  //   a) That is, Bx, v2, v1, By all exist in fluxstag, fluxstagcent, or at CORN3, within cell so that should have enough info to form 2D deconvolution at all CORN3's
  // INTERPOLATING FIELDS ALONG SELF:
  // a) form full array \delta{Bhat}[COMPDIM][i][j][k] with \detg included
  // b) Interpolate this \delta like normal CENT quantity
  // c) compute Bc, Bleft, Bright and store so can use on next timestep for merged method -- only need this extra step for higher-than-para methods since do already store in "state" the primitive field, so have it at left,right faces and CENT as required/needed -- instead use pstagglobal for FLUXCTSTAG for field along itself.  As discussed in step_ch.c in RK methods, pstagglobal contains last pf-related staggered updated field that will always (currently) be used as pb on this step
  //... other notes related to this.
  // See also Londrillo & Del Zanna (2004) JCP 195, 17-48 around equations 22.


  // for higher than para methods:
  // a) obtain result at i-1_CENT, i_FACE1, i_CENT, i+1_FACE1, i+1_CENT  -- OR even better: i_FACE1, i+1/4 i_CENT, i+3/4, i+1_FACE1 using Sasha's coefficients result
  // b) geometry: use 5th order ENO for geometry using i_FACE1, i_CENT, i+1_FACE1 as base points to keep compact since have geometry at these compact locations
  //    worry about dx/2 vs. dx in using Sasha's coefficient extraction
  // c) aux functions (e.g. like pressure or other desired aux functions), form at 2 more points inside (i.e. i+1/4 and i+3/4) from higher-order primitives.  Then just use 5th order ENO without reduction to obtain coefficients (or use data directly!).  The resulting aux functions will of course be perfectly consistent with primitives, but not perfectly conservative-primitive consistent since expanded number of independent funtions that come into deconvolving the flux.  However, if start with compactly interpolated primitives, then should be quite good still.

  // Note that in the limit that the aux function is the flux, if one is only using 3-points per flux to do the deconvolution, then no help since originally flux error is 3rd order anyways.  If wanted to try this, then should store Fleft, Fright when computed in fluxcompute.c (store into state?) and then have to compute flux at center (and store into state?)  GODMARK: CHECK THIS IN MATHEMATICA -- checked, true: check_centered_diff_forflux.nb
  // Well, while true that F\propto O(dx^3), but a2c correction force error in flux to be exactly vanishing.  This a2c produces a correction order dx^2 due to second derivative of flux that is what allows for the error in dF/dx to vanish. It's this vanishing that is important for stiff problems where terms within the flux can be much larger(smaller) than others and be dominating(dominated by) other terms so that the smaller scale is forced to have huge error given by larger term.

  // Energy-Momentum:
  // Note \gamma^2 type terms cancel from stress if one uses B^i instead of b^\mu
  // Still end up needing B_\mu for B^2 term, so now that is stored in state

  // EMF: \detg *F^{ij} = \detg [b^i u^j - b^j u^i] = \detg [B^i v^j - B^j v^i], where j=flux dir and i is which field to flux





  ///////////////////
  //
  // Energy-Momentum FLUX: \detg [ ( b^2 )u^j u_\nu + \delta^j_\nu b^2/2 ]
  //
  // Note that using 4-field b^\mu causes term to be present that is \propto \gamma^2 that is larger than all other terms.  This can swap the other terms due to machine errors in the larger term.  If one uses 3-field in the lab-frame B^i, then this term can be analytically cancelled.  So this is what we do here.
  //
  // Final result is:
  //
  // \detg T^j_\nu[EM] = \detg (u^t)^{-2} [ (B^i B_i u^j u_\nu + (1/2)(B^i B_i + (u_i B^i)^2) \delta^j_\nu  - B^j B_\nu - (u_i B^i)(u^j B_\nu + B^j u_\nu)  ]
  //
  // in keeping with merged method only dealing with multiplicative factors, treat (1/u^t) as independent variable since we don't want to mix with B or u^i or u_i since then like interpolating 3-velocity that corresponds to no solution within the domain
  //
  // Compromised:  Absorb u^t into u^\mu or u_\mu if possible since for EMF do interpolate 3-velocity.  GODMARK: In case where perfect balance between MA and EM terms is required, interpolating same quantities is desired (i.e. all 4-vel)
  // Also: Abosrbed \detg into \detg B^{dir} if in term.  In one case two B^{dir}'s exist, but one is related to u.B, so assume scalar so ignore geometry.  GODMARK: In perfect case, should convolve exactly those things interpolated, but for this term one would have division by \detg to make this happen!
  // 
  //
  ///////////////////


  ////////////////////////
  //
  // Setup things common to all EOMs
  //
  ////////////////////////

  // \detg for \rho_0
#if(MERGEDC2EA2CMETHOD)
  // overridden by below if WHICHEOM!=WITHGDET
  gdetl=stateleft->gdet;
  gdetc=statecent->gdet;
  gdetr=stateright->gdet;
#endif

  // upper 4-velocity in coordinate frame
  DLOOPA(jj){
    uul[jj]=stateleft->ucon[jj];
    uuc[jj]=statecent->ucon[jj];
    uur[jj]=stateright->ucon[jj];

    udl[jj]=stateleft->ucov[jj];
    udc[jj]=statecent->ucov[jj];
    udr[jj]=stateright->ucov[jj];

#if(MERGEDC2EA2CMETHOD)
    // lower 3-field in coordinate frame
    Bdl[jj]=stateleft->Blower[jj];
    Bdc[jj]=statecent->Blower[jj];
    Bdr[jj]=stateright->Blower[jj];
#endif
  }

#if(MERGEDC2EA2CMETHOD)
  Bul[TT]=Buc[TT]=Bur[TT]=0.0;
  SLOOPA(jj){
    // upper 3-field in coordinate frame
    Bul[jj]=stateleft->prim[B1-1+jj];
    Buc[jj]=statecent->prim[B1-1+jj];
    Bur[jj]=stateright->prim[B1-1+jj];
  }
#endif

  // upper 4-velocity in coordinate frame
  overutl=1.0/(stateleft->ucon[TT]);
  overutc=1.0/(statecent->ucon[TT]);
  overutr=1.0/(stateright->ucon[TT]);


  // setup loop over non-diagonal components
  myloop[0]=TT;
  if(dir==1){ 
    myloop[1]=2;
    myloop[2]=3;
  }
  else if(dir==2){
    myloop[1]=1;
    myloop[2]=3;
  }
  else if(dir==3){
    myloop[1]=1;
    myloop[2]=2;
  }

  ///////////////////////////
  //
  // loop over \nu=0,1,2,3 for EM EOMs except diagonal (nu=dir) term
  //
  ///////////////////////////
  //DLOOPA(jj){
  // only loop over non-diagonal terms, doing diagonal term last outside loop (done to avoid conditional inside loop)
  for(iter=0;iter<NDIM-1;iter++){
    jj=myloop[iter];


#if(WHICHEOM!=WITHGDET&&MERGEDC2EA2CMETHOD)
    // need gdet for each variable
    gdetl=stateleft->EOMFUNCMAC(UU+jj);
    gdetc=statecent->EOMFUNCMAC(UU+jj);
    gdetr=stateright->EOMFUNCMAC(UU+jj);
#endif
    
    // initialize as we loop over terms that are simply part of sum that has no extra coupling for deconvolution
    Fleft[UU+jj]  = 0.0;
    Fright[UU+jj] = 0.0;


#if(0&&MCOORD==CARTMINKMETRIC)
    //
    // GODMARK: Not setup yet, using full \detg method
    //
#else

    
    // See check_tmunu_em_cancellations.nb (cancellations after 3-field is used):
    // applies to only dir=1,2,3 and \nu=0,1,2,3 (i.e. reduced so no longer symmetric for time-spatial terms! and not valid for time-time term)

    // non-daginal terms:
    //        t[1 ]=+(-gdet*Bu[dir])*(Bd[jj])*(1/uu[TT])*(1/uu[TT]);
    fourtermdeconvolution(-gdetl*Bul[dir],-gdetc*Buc[dir],-gdetr*Bur[dir],  Bdl[jj],Bdc[jj],Bdr[jj],  overutl,overutc,overutr,   overutl,overutc,overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[2 ]=+(-gdet*Bu[dir])*(Bu[odir1])*(ud[jj]/uu[TT])*(ud[odir1]/uu[TT]);
    fourtermdeconvolution(-gdetl*Bul[dir],-gdetc*Buc[dir],-gdetr*Bur[dir],  Bul[odir1],Buc[odir1],Bur[odir1],  udl[jj]*overutl,udc[jj]*overutc,udr[jj]*overutr,   udl[odir1]*overutl,udc[odir1]*overutc,udr[odir1]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[3 ]=+(-gdet*Bu[dir])*(Bu[odir2])*(ud[jj]/uu[TT])*(ud[odir2]/uu[TT]);
    fourtermdeconvolution(-gdetl*Bul[dir],-gdetc*Buc[dir],-gdetr*Bur[dir],  Bul[odir2],Buc[odir2],Bur[odir2],  udl[jj]*overutl,udc[jj]*overutc,udr[jj]*overutr,   udl[odir2]*overutl,udc[odir2]*overutc,udr[odir2]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[4 ]=+(gdet*Bu[dir])*(Bd[dir])*(uu[dir]/uu[TT])*(ud[jj]/uu[TT]);
    fourtermdeconvolution(gdetl*Bul[dir],gdetc*Buc[dir],gdetr*Bur[dir],  Bdl[dir],Bdc[dir],Bdr[dir],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[jj]*overutl,udc[jj]*overutc,udr[jj]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[5 ]=+(gdet)*(Bu[odir1])*(Bd[odir1])*(uu[dir]/uu[TT])*(ud[jj]/uu[TT]);
    fivetermdeconvolution(gdetl,gdetc,gdetr,   Bul[odir1],Buc[odir1],Bur[odir1],   Bdl[odir1],Bdc[odir1],Bdr[odir1],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[jj]*overutl,udc[jj]*overutc,udr[jj]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[6 ]=+(gdet)*(Bu[odir2])*(Bd[odir2])*(uu[dir]/uu[TT])*(ud[jj]/uu[TT]);
    fivetermdeconvolution(gdetl,gdetc,gdetr,   Bul[odir2],Buc[odir2],Bur[odir2],   Bdl[odir2],Bdc[odir2],Bdr[odir2],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[jj]*overutl,udc[jj]*overutc,udr[jj]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[7 ]=+(-gdet*Bu[dir])*(Bd[jj])*(uu[dir]/uu[TT])*(ud[dir]/uu[TT]);
    fourtermdeconvolution(-gdetl*Bul[dir],-gdetc*Buc[dir],-gdetr*Bur[dir],  Bdl[jj],Bdc[jj],Bdr[jj],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[dir]*overutl,udc[dir]*overutc,udr[dir]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[8 ]=+(-gdet)*(Bu[odir1])*(Bd[jj])*(uu[dir]/uu[TT])*(ud[odir1]/uu[TT]);
    fivetermdeconvolution(-gdetl,-gdetc,-gdetr,   Bul[odir1],Buc[odir1],Bur[odir1],   Bdl[jj],Bdc[jj],Bdr[jj],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[odir1]*overutl,udc[odir1]*overutc,udr[odir1]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[9 ]=+(-gdet)*(Bu[odir2])*(Bd[jj])*(uu[dir]/uu[TT])*(ud[odir2]/uu[TT]);
    fivetermdeconvolution(-gdetl,-gdetc,-gdetr,   Bul[odir2],Buc[odir2],Bur[odir2],   Bdl[jj],Bdc[jj],Bdr[jj],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[odir2]*overutl,udc[odir2]*overutc,udr[odir2]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

    // t[10]=+(-gdet*Bu[dir])*(Bu[dir])*(ud[jj]/uu[TT])*(ud[dir]/uu[TT]);
    fourtermdeconvolution(-gdetl*Bul[dir],-gdetc*Buc[dir],-gdetr*Bur[dir],  Bul[dir],Buc[dir],Bur[dir],  udl[jj]*overutl,udc[jj]*overutc,udr[jj]*overutr,   udl[dir]*overutl,udc[dir]*overutc,udr[dir]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;
#endif

  }// end loop over jj = \nu={0,1,2,3} except jj=dir




  /////////////////////////
  //
  // now do diagonal term
  //
  ////////////////////////
  jj=dir;

#if(WHICHEOM!=WITHGDET&&MERGEDC2EA2CMETHOD)
  // need gdet for each variable
  gdetl=stateleft->EOMFUNCMAC(UU+jj);
  gdetc=statecent->EOMFUNCMAC(UU+jj);
  gdetr=stateright->EOMFUNCMAC(UU+jj);
#endif
   
  // initialize as we loop over terms that are simply part of sum that has no extra coupling for deconvolution
  Fleft[UU+jj]  = 0.0;
  Fright[UU+jj] = 0.0;

#if(0&&MCOORD==CARTMINKMETRIC)
  //
  // GODMARK: Not setup yet, using full \detg method
  //
#else
  // diagonal spatial terms
  // t[1 ]=+(-0.5*gdet*Bu[dir])*(Bd[dir])*(1/uu[TT])*(1/uu[TT]);
  fourtermdeconvolution(-0.5*gdetl*Bul[dir],-0.5*gdetc*Buc[dir],-0.5*gdetr*Bur[dir],  Bdl[dir],Bdc[dir],Bdr[dir],  overutl,overutc,overutr,   overutl,overutc,overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[2 ]=+(0.5*gdet)*(Bu[odir1])*(Bd[odir1])*(1/uu[TT])*(1/uu[TT]);
  fivetermdeconvolution(0.5*gdetl,0.5*gdetc,0.5*gdetr,  Bul[odir1],Buc[odir1],Bur[odir1],  Bdl[odir1],Bdc[odir1],Bdr[odir1],  overutl,overutc,overutr,   overutl,overutc,overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[3 ]=+(0.5*gdet)*(Bu[odir2])*(Bd[odir2])*(1/uu[TT])*(1/uu[TT]);
  fivetermdeconvolution(0.5*gdetl,0.5*gdetc,0.5*gdetr,  Bul[odir2],Buc[odir2],Bur[odir2],  Bdl[odir2],Bdc[odir2],Bdr[odir2],  overutl,overutc,overutr,   overutl,overutc,overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[4 ]=+(gdet)*(Bu[odir1])*(Bu[odir2])*(ud[odir1]/uu[TT])*(ud[odir2]/uu[TT]);
  fivetermdeconvolution(gdetl,gdetc,gdetr,   Bul[odir1],Buc[odir1],Bur[odir1],   Bul[odir2],Buc[odir2],Bur[odir2],  udl[odir1]*overutl,udc[odir1]*overutc,udr[odir1]*overutr,   udl[odir2]*overutl,udc[odir2]*overutc,udr[odir2]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[5 ]=+(-gdet)*(Bu[odir1])*(Bd[dir])*(uu[dir]/uu[TT])*(ud[odir1]/uu[TT]);
  fivetermdeconvolution(-gdetl,-gdetc,-gdetr,   Bul[odir1],Buc[odir1],Bur[odir1],   Bdl[dir],Bdc[dir],Bdr[dir],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[odir1]*overutl,udc[odir1]*overutc,udr[odir1]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[6 ]=+(-gdet)*(Bu[odir2])*(Bd[dir])*(uu[dir]/uu[TT])*(ud[odir2]/uu[TT]);
  fivetermdeconvolution(-gdetl,-gdetc,-gdetr,   Bul[odir2],Buc[odir2],Bur[odir2],   Bdl[dir],Bdc[dir],Bdr[dir],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[odir2]*overutl,udc[odir2]*overutc,udr[odir2]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[7 ]=+(gdet)*(Bu[odir1])*(Bd[odir1])*(uu[dir]/uu[TT])*(ud[dir]/uu[TT]);
  fivetermdeconvolution(gdetl,gdetc,gdetr,   Bul[odir1],Buc[odir1],Bur[odir1],   Bdl[odir1],Bdc[odir1],Bdr[odir1],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[dir]*overutl,udc[dir]*overutc,udr[dir]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[8 ]=+(gdet)*(Bu[odir2])*(Bd[odir2])*(uu[dir]/uu[TT])*(ud[dir]/uu[TT]);
  fivetermdeconvolution(gdetl,gdetc,gdetr,   Bul[odir2],Buc[odir2],Bur[odir2],   Bdl[odir2],Bdc[odir2],Bdr[odir2],  uul[dir]*overutl,uuc[dir]*overutc,uur[dir]*overutr,   udl[dir]*overutl,udc[dir]*overutc,udr[dir]*overutr, &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[9 ]=+(-0.5*gdet*Bu[dir])*(Bu[dir])*(ud[dir]/uu[TT])*(ud[dir]/uu[TT]);
  fourtermdeconvolution(-0.5*gdetl*Bul[dir],-0.5*gdetc*Buc[dir],-0.5*gdetr*Bur[dir],  Bul[dir],Buc[dir],Bur[dir],  udl[dir]*overutl,udc[dir]*overutc,udr[dir]*overutr,   udl[dir]*overutl,udc[dir]*overutc,udr[dir]*overutr,  &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[10]=+(0.5*gdet)*(Bu[odir1])*(Bu[odir1])*(ud[odir1]/uu[TT])*(ud[odir1]/uu[TT]);
  fivetermdeconvolution(0.5*gdetl,0.5*gdetc,0.5*gdetr,   Bul[odir1],Buc[odir1],Bur[odir1],   Bul[odir1],Buc[odir1],Bur[odir1],  udl[odir1]*overutl,udc[odir1]*overutc,udr[odir1]*overutr,   udl[odir1]*overutl,udc[odir1]*overutc,udr[odir1]*overutr,  &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;

  // t[11]=+(0.5*gdet)*(Bu[odir2])*(Bu[odir2])*(ud[odir2]/uu[TT])*(ud[odir2]/uu[TT]);
  fivetermdeconvolution(0.5*gdetl,0.5*gdetc,0.5*gdetr,   Bul[odir2],Buc[odir2],Bur[odir2],   Bul[odir2],Buc[odir2],Bur[odir2],  udl[odir2]*overutl,udc[odir2]*overutc,udr[odir2]*overutr,   udl[odir2]*overutl,udc[odir2]*overutc,udr[odir2]*overutr,  &Flefttemp, &Frighttemp);   Fleft[UU+jj] += Flefttemp; Fright[UU+jj] += Frighttemp;
#endif











  // FLUXCTHLL method that operates along dimension lines as normal non-EMF flux:
  // This makes the method higher-order but this method does not preserve divB=0
  // Here use 3-term deconvolution since not expensive compared to above and no worry about preserving divB=0 and grouping of \detg with B^{dir}
  if(FLUXB==FLUXCTHLL){

    // see fluxct.c where below was copied from:
    // emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]
    // emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]
    // emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]

    // F^dir_{Bnu} = B^{nu} v^{dir} - B^{dir} v^{nu}

    // iterate only over spatial EOMs that are not along dir (\nu==dir has flux=0 identically)
    // loop here is differently used than above since here start with iter=1
    for(iter=1;iter<NDIM-1;iter++){
      jj=myloop[iter]; // term is jj-EOM, not flux in jj direction

#if(WHICHEOM!=WITHGDET&&MERGEDC2EA2CMETHOD)
      // need gdet for each variable
      gdetl=stateleft->EOMFUNCMAC(B1-1+jj);
      gdetc=statecent->EOMFUNCMAC(B1-1+jj);
      gdetr=stateright->EOMFUNCMAC(B1-1+jj);
#endif


      // get corrections for both terms in EMF
      for(otherdiriter=0;otherdiriter<=1;otherdiriter++){

        // choose sign: otherdiriter==0 has positive sign
        fluxsign=otherdiriter*(-1.0) + (1-otherdiriter)*(1.0);

        // these odir choices are not directions of interpolation, but instead of components of vectors
        // choose which odir is first
        firstodir=jj*(1-otherdiriter) + dir*otherdiriter;
        // second is always other odir
        secondodir=dir*(1-otherdiriter) + jj*otherdiriter;

        // NEWMARK: Signature!
        // EMF[corner with odir1=dir and odir2=jj] = +\detg B^{jj} v^{dir} -\detg B^{dir} v^{jj}
        threetermdeconvolution(gdetl, gdetc, gdetr,   Bul[firstodir], Buc[firstodir], Bur[firstodir]   ,uul[secondodir]*overutl, uuc[secondodir]*overutc, uur[secondodir]*overutr   ,&Flefttemp, &Frighttemp);
        //  note that unlike FLUXCTSTAG or FLUXCTTOTH, FLUXCTHLL has different fluxes for (e.g.) F2[B1] and F1[B2] rather than just sign change since at different locations
        // 0.5 is since 2 flux corrections will be made for each face, and assumes LAXF or limit of HLL fluxmethod
        Fleft[B1-1+jj]  += 0.5*fluxsign*Flefttemp;
        Fright[B1-1+jj] += 0.5*fluxsign*Frighttemp;
      }
    }// end over flux directions (not corners)

  } // end if FLUXB==FLUXCTHLL





  return(0);
}




// Note that always \hat{odir1}\times\hat{odir2}=\hat{corner,dir}
// F[+-odir1 l/r][+-odir2 d/u]
static int deconvolve_emf_2d(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  FTYPE vcon[NDIM][NUMCS+1][NUMCS+1];
  FTYPE gdetBcon[NDIM][NUMCS+1][NUMCS+1];
  int ii,jj,kk;
  int otherdiriter,firstodir,secondodir;
  FTYPE Fld,Frd,Flu,Fru;


  // get 3x3 matrix of values for 2D deconvolution
  setup_9value_vB(corner, odir1, odir2, Nvec, NNOT1vec, i, j, k, GLOBALPOINT(pvbcorninterp), GLOBALPOINT(fluxstatecent), GLOBALPOINT(fluxstate), GLOBALPOINT(geomcornglobal), vcon, gdetBcon);


  // get corrections for both terms in EMF
  for(otherdiriter=0;otherdiriter<=1;otherdiriter++){

    // choose which odir is first
    firstodir=odir1*(1-otherdiriter) + odir2*otherdiriter;
    // second is always other odir
    secondodir=odir2*(1-otherdiriter) + odir1*otherdiriter;

    // NEWMARK: signature!  (normally sign flips between each, but right now don't know which is negative or positive
    twoterm2Ddeconvolution(
                           gdetBcon[firstodir][CENT4EMF][CENT4EMF],gdetBcon[firstodir][LEFT4EMF][LEFT4EMF],gdetBcon[firstodir][RIGHT4EMF][LEFT4EMF],gdetBcon[firstodir][LEFT4EMF][RIGHT4EMF],gdetBcon[firstodir][RIGHT4EMF][RIGHT4EMF],gdetBcon[firstodir][LEFT4EMF][CENT4EMF],gdetBcon[firstodir][RIGHT4EMF][CENT4EMF],gdetBcon[firstodir][CENT4EMF][LEFT4EMF],gdetBcon[firstodir][CENT4EMF][RIGHT4EMF]
                           ,vcon[secondodir][CENT4EMF][CENT4EMF],vcon[secondodir][LEFT4EMF][LEFT4EMF],vcon[secondodir][RIGHT4EMF][LEFT4EMF],vcon[secondodir][LEFT4EMF][RIGHT4EMF],vcon[secondodir][RIGHT4EMF][RIGHT4EMF],vcon[secondodir][LEFT4EMF][CENT4EMF],vcon[secondodir][RIGHT4EMF][CENT4EMF],vcon[secondodir][CENT4EMF][LEFT4EMF],vcon[secondodir][CENT4EMF][RIGHT4EMF]
                           ,&Fld, &Frd, &Flu, &Fru);

    // correct left-down corner
    // e.g., if corner==3, then odir1=1 and odir2=2, then EMF at left-down CORN3 at i,j
    // NEWMARK: signature!
    if(NNOT1vec[odir1]) MACP1A1(fluxvec,odir1,i,j,k,B1-1+odir2) += 0.25*Fld;
    if(NNOT1vec[odir2]) MACP1A1(fluxvec,odir2,i,j,k,B1-1+odir1) += -0.25*Fld;
  
    // correct right-down corner
    // e.g., if corner==3, then odir1=1 and odir2=2, then EMF at right-down CORN3 at i+1,j
    ii=i+(1==odir1)*NNOT1vec[1];
    jj=j+(2==odir1)*NNOT1vec[2];
    kk=k+(3==odir1)*NNOT1vec[3];
    // NEWMARK: signature!
    if(NNOT1vec[odir1]) MACP1A1(fluxvec,odir1,ii,jj,kk,B1-1+odir2) += 0.25*Frd;
    if(NNOT1vec[odir2]) MACP1A1(fluxvec,odir2,ii,jj,kk,B1-1+odir1) += -0.25*Frd;

    // correct left-up corner
    // e.g., if corner==3, then odir1=1 and odir2=2, then EMF at left-up CORN3 at i,j+1
    ii=i+(1==odir2)*NNOT1vec[1];
    jj=j+(2==odir2)*NNOT1vec[2];
    kk=k+(3==odir2)*NNOT1vec[3];
    // NEWMARK: signature!
    if(NNOT1vec[odir1]) MACP1A1(fluxvec,odir1,ii,jj,kk,B1-1+odir2) += 0.25*Flu;
    if(NNOT1vec[odir2]) MACP1A1(fluxvec,odir2,ii,jj,kk,B1-1+odir1) += -0.25*Flu;

    // correct right-up corner
    // e.g., if corner==3, then odir1=1 and odir2=2, then EMF at right-up CORN3 at i+1,j+1
    // odir1 and odir2 will never be equal, so this sum is good
    ii=i+(1==odir1+(1==odir2))*NNOT1vec[1];
    jj=j+(2==odir1+(2==odir2))*NNOT1vec[2];
    kk=k+(3==odir1+(3==odir2))*NNOT1vec[3];
    // NEWMARK: signature!
    if(NNOT1vec[odir1]) MACP1A1(fluxvec,odir1,ii,jj,kk,B1-1+odir2) += 0.25*Fru;
    if(NNOT1vec[odir2]) MACP1A1(fluxvec,odir2,ii,jj,kk,B1-1+odir1) += -0.25*Fru;

  }  
  


  return(0);
}




static int deconvolve_emf_1d(int corner, int odir1, int odir2, int *Nvec, int *NNOT1vec, int i, int j, int k, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  FTYPE Fleft[NPR],Fright[NPR];
  FTYPE Flefttemp,Frighttemp;
  int doingodir, notdoingodir;
  int firstodir, secondodir;
  int ileft,jleft,kleft;
  int iright,jright,kright;
  FTYPE vcon[NDIM][NUMCS+1];
  FTYPE gdetBcon[NDIM][NUMCS+1];
  int otherdiriter;


  // set 3 values in correct non-trivial transverse direction
  setup_3value_vB(corner, odir1, odir2, Nvec, NNOT1vec, i, j, k, GLOBALPOINT(pvbcorninterp), GLOBALPOINT(fluxstatecent), GLOBALPOINT(fluxstate), GLOBALPOINT(geomcornglobal), vcon, gdetBcon);

  // perform 1D deconvolution using 2-product version twice

  // choose which odir we are doing (only 1 of 2 is chosen)
  doingodir=NNOT1vec[odir1]*odir1 + NNOT1vec[odir2]*odir2;
  // choose other dir we are not doing (only 1 of 2 is chosen)
  notdoingodir=(NNOT1vec[odir1]-1)*odir1 + (1-NNOT1vec[odir2])*odir2;

  Fleft[B1-1+doingodir]=Fright[B1-1+doingodir]=0.0;

  // get corrections for both terms in EMF
  for(otherdiriter=0;otherdiriter<=1;otherdiriter++){

    // choose which odir is first
    firstodir=doingodir*(1-otherdiriter) + notdoingodir*otherdiriter;
    // second is always other odir
    secondodir=notdoingodir*(1-otherdiriter) + doingodir*otherdiriter;


    // NEWMARK: Need to work out signature issue
    twotermdeconvolution(vcon[firstodir][LEFT4EMF], vcon[firstodir][CENT4EMF], vcon[firstodir][RIGHT4EMF], gdetBcon[secondodir][LEFT4EMF], gdetBcon[secondodir][CENT4EMF], gdetBcon[secondodir][RIGHT4EMF]  ,&Flefttemp, &Frighttemp);     Fleft[B1-1+doingodir] += Flefttemp; Fright[B1-1+doingodir] += Frighttemp;

  }

  // left-flux position (centered or left in odir1 and/or odir2)
  ileft=i;
  jleft=j;
  kleft=k;
  
  // correct the flux that is the emf
  // in 1D, factor is 1/2 since implicitly sum of 0.25+0.25 for LAXF (HLL more complicated in general)
  
  // correct left (-doingodir direction) flux
  MACP1A1(fluxvec,doingodir,ileft,jleft,kleft,B1-1+notdoingodir) += 0.5*Fleft[B1-1+doingodir];

  // right-flux position for doingodir direction
  iright=i+(1==doingodir)*NNOT1vec[1];
  jright=j+(2==doingodir)*NNOT1vec[2];
  kright=k+(3==doingodir)*NNOT1vec[3];
  
  // correct left (+doingodir direction) flux
  MACP1A1(fluxvec,doingodir,iright,jright,kright,B1-1+notdoingodir) += 0.5*Fright[B1-1+doingodir];

  // Do not correct fluxvec[notdoingodir][B1-1+doingodir] since that flux doesn't exist (or implicitly the EMF cancels itself across the box in the dimension that is not studied)


  return(0);

}







// Some WENO5BND-related things not setup yet
static int a2cflux_from_prim(int dir, FTYPE (*prim_coef_list)[MAXSPACEORDER])
{


  return(0);





}



 


// 1-term trivial deconvolution
static int onetermdeconvolution(
                                FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                ,FTYPE *Fleft, FTYPE *Fright
                                )
{

  //LEFT:
  *Fleft = SIXTH*(2.*rhoc - rhol - rhor);

  //RIGHT:
  *Fright = SIXTH*(2.*rhoc - rhol - rhor);

  return(0);

}


// 2-terms multiplying deconvolution
// inputs are named after \rho_0 u^i (without \detg) as could be useful when \detg is constant
// See: a2confluxes_eom1_restmassflux.nb
static int twotermdeconvolution(
                                FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                ,FTYPE udirl, FTYPE udirc, FTYPE udirr
                                ,FTYPE *Fleft, FTYPE *Fright
                                )
{

  //Left:
  *Fleft = (1./30.)*(rhol*(33.*udirc - 29.*udirl - 9.*udirr) + rhor*(3.*udirc - 9.*udirl + udirr) + rhoc*(-26.*udirc + 33.*udirl + 3.*udirr));
  
  //Right:
  *Fright = (1./30.)*(rhor*(33.*udirc - 9.*udirl - 29.*udirr) + rhol*(3.*udirc + udirl - 9.*udirr) + rhoc*(-26.*udirc + 3.*udirl + 33.*udirr));

  return(0);

}


// 3-terms multiplying deconvolution
// inputs are named after rest-mass flux term; \detg \rho_0 u^i, but can be used for any 3-term multiplication
// See: a2confluxes_eom1_restmassflux.nb
static int threetermdeconvolution(
                                  FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                  ,FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                  ,FTYPE udirl, FTYPE udirc, FTYPE udirr
                                  ,FTYPE *Fleft, FTYPE *Fright
                                  )
{
  FTYPE t[42]; // 46 total temp vars (42 in t[] for left, and 23 in t[] for right)

  // 3 terms for gdet*\rho_0 udir^i


  FTYPE vartemp63310=-45.*udirc;
  FTYPE vartemp63315=-19.*udirc;
  FTYPE vartemp63316=udirl+udirr;
  FTYPE vartemp63317=6.*vartemp63316;
  FTYPE vartemp63318=vartemp63315+vartemp63317;
  t[1]=-866.*rhoc*udirc;
  t[2]=447.*rhol*udirc;
  t[3]=237.*rhor*udirc;
  t[4]=447.*rhoc*udirl;
  t[5]=-45.*rhol*udirl;
  t[6]=-171.*rhor*udirl;
  t[7]=79.*rhoc;
  t[8]=-57.*rhol;
  t[9]=-15.*rhor;
  t[10]=udirr*(t[7]+t[8]+t[9]);
  t[11]=t[6];
  t[12]=3.*t[10];
  t[13]=t[1]+t[2]+t[3];
  t[14]=t[4]+t[5]+t[11]+t[12];
  t[15]=79.*udirc;
  t[16]=-57.*udirl;
  t[17]=-15.*udirr;
  t[18]=rhoc*(t[15]+t[16]+t[17]);
  t[19]=54.*udirl;
  t[20]=-2.*udirr+vartemp63310;
  t[21]=9.*rhol*vartemp63318;
  t[22]=rhor*(t[19]+t[20]);
  t[23]=t[21];
  t[24]=3.*t[18];
  t[25]=t[22]+t[23];
  t[26]=149.*udirc;
  t[27]=-15.*udirl;
  t[28]=-57.*udirr;
  t[29]=rhoc*(t[26]+t[27]+t[28]);
  t[30]=-212.*udirl;
  t[31]=54.*udirr+vartemp63310;
  t[32]=9.*rhor*vartemp63318;
  t[33]=rhol*(t[30]+t[31]);
  t[34]=t[32];
  t[35]=3.*t[29];
  t[36]=t[33]+t[34];
  t[37]=gdetr*(t[24]+t[25]);
  t[38]=gdetl*(t[35]+t[36]);
  t[39]=gdetc*(t[13]+t[14]);
  t[40]=t[37]+t[38];
  t[41]=t[39]+t[40];
  *Fleft=0.0047619047619047619047619047619048*t[41];







  FTYPE vartemp63400=149.*udirc;
  FTYPE vartemp63401=-57.*udirl;
  FTYPE vartemp63402=-15.*udirr;
  FTYPE vartemp63403=vartemp63400+vartemp63401+vartemp63402;
  FTYPE vartemp63395=79.*udirc;
  FTYPE vartemp63396=-15.*udirl;
  FTYPE vartemp63397=-57.*udirr;
  FTYPE vartemp63398=vartemp63395+vartemp63396+vartemp63397;
  FTYPE vartemp63414=-45.*udirc;
  FTYPE vartemp63420=-19.*udirc;
  FTYPE vartemp63421=udirl+udirr;
  FTYPE vartemp63422=6.*vartemp63421;
  FTYPE vartemp63423=vartemp63420+vartemp63422;
  t[1]=-866.*udirc;
  t[2]=237.*udirl;
  t[3]=447.*udirr;
  t[4]=rhol*vartemp63398+rhor*vartemp63403;
  t[5]=rhoc*(t[1]+t[2]+t[3]);
  t[6]=3.*t[4];
  t[7]=3.*rhoc*vartemp63403;
  t[8]=54.*udirl;
  t[9]=-212.*udirr+vartemp63414;
  t[10]=9.*rhol*vartemp63423;
  t[11]=rhor*(t[8]+t[9]);
  t[12]=t[10];
  t[13]=3.*rhoc*vartemp63398;
  t[14]=-2.*udirl;
  t[15]=54.*udirr+vartemp63414;
  t[16]=9.*rhor*vartemp63423;
  t[17]=rhol*(t[14]+t[15]);
  t[18]=t[16];
  t[19]=gdetr*(t[7]+t[11]+t[12]);
  t[20]=gdetl*(t[13]+t[17]+t[18]);
  t[21]=gdetc*(t[5]+t[6]);
  t[22]=t[19]+t[20];
  t[23]=t[21]+t[22];
  *Fright=0.0047619047619047619047619047619048*t[23];


  return(0);

}


// 4-terms multiplying deconvolution
// inputs are named after T^i_\nu term; \detg \rho_0 u^i u_\nu, but can be used for any 4-term multiplication
// See: a2confluxes_eom2_5_energymomentumflux_MAtermsonly.nb
static int fourtermdeconvolution(
                                 FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                 ,FTYPE rhol, FTYPE rhoc, FTYPE rhor
                                 ,FTYPE udirl, FTYPE udirc, FTYPE udirr
                                 ,FTYPE udnul, FTYPE udnuc, FTYPE udnur
                                 ,FTYPE *Fleft, FTYPE *Fright
                                 )
{
  FTYPE t[49]; // 92 total temp vars, 49 in t[].


  FTYPE vartemp331067=-106.*udnuc;
  FTYPE vartemp331055=106.*udnuc;
  FTYPE vartemp331050=102.*udnuc;
  FTYPE vartemp331060=-155.*udnuc;
  FTYPE vartemp331068=53.*udnul;
  FTYPE vartemp331069=35.*udnur;
  FTYPE vartemp331070=vartemp331067+vartemp331068+vartemp331069;
  FTYPE vartemp331071=udirr*vartemp331070;
  FTYPE vartemp331072=35.*udnul;
  FTYPE vartemp331073=53.*udnur;
  FTYPE vartemp331074=vartemp331067+vartemp331072+vartemp331073;
  FTYPE vartemp331075=udirl*vartemp331074;
  FTYPE vartemp331076=269.*udnuc;
  FTYPE vartemp331077=udnul+udnur;
  FTYPE vartemp331078=-106.*vartemp331077;
  FTYPE vartemp331079=vartemp331076+vartemp331078;
  FTYPE vartemp331080=udirc*vartemp331079;
  FTYPE vartemp331081=vartemp331071+vartemp331075+vartemp331080;
  FTYPE vartemp331092=106.*udnul;
  FTYPE vartemp331093=34.*udnur;
  FTYPE vartemp331094=vartemp331060+vartemp331092+vartemp331093;
  FTYPE vartemp331037=-155.*udirl*udnuc;
  FTYPE vartemp331038=-269.*udirr*udnuc;
  FTYPE vartemp331039=34.*udirl*udnul;
  FTYPE vartemp331040=106.*udirr*udnul;
  FTYPE vartemp331041=573.*udnuc;
  FTYPE vartemp331042=-155.*udnul;
  FTYPE vartemp331043=-269.*udnur;
  FTYPE vartemp331044=vartemp331041+vartemp331042+vartemp331043;
  FTYPE vartemp331045=udirc*vartemp331044;
  FTYPE vartemp331046=udirl+udirr;
  FTYPE vartemp331047=106.*udnur*vartemp331046;
  FTYPE vartemp331048=vartemp331037+vartemp331038+vartemp331039+vartemp331040+vartemp331045+vartemp331047;
  FTYPE vartemp331099=503.*udnuc;
  FTYPE vartemp331100=-269.*udnul;
  FTYPE vartemp331101=-155.*udnur;
  FTYPE vartemp331102=vartemp331099+vartemp331100+vartemp331101;
  FTYPE vartemp331103=udirc*vartemp331102;
  FTYPE vartemp331104=udirr*vartemp331094;
  FTYPE vartemp331105=-269.*udnuc;
  FTYPE vartemp331106=106.*vartemp331077;
  FTYPE vartemp331107=vartemp331105+vartemp331106;
  FTYPE vartemp331108=udirl*vartemp331107;
  FTYPE vartemp331109=vartemp331103+vartemp331104+vartemp331108;
  t[1]=3.*rhoc*vartemp331048;
  t[2]=-209.*udnul;
  t[3]=-105.*udnur+vartemp331050;
  t[4]=-35.*udnul;
  t[5]=-53.*udnur+vartemp331055;
  t[6]=udirr*(t[4]+t[5]);
  t[7]=34.*udnul;
  t[8]=106.*udnur+vartemp331060;
  t[9]=udirc*(t[7]+t[8]);
  t[10]=3.*t[6];
  t[11]=3.*t[9];
  t[12]=udirl*(t[2]+t[3]);
  t[13]=t[10]+t[11];
  t[14]=-3.*rhor*vartemp331081;
  t[15]=rhol*(t[12]+t[13]);
  t[16]=t[14];
  t[17]=-3.*rhol*vartemp331081;
  t[18]=udirr*(-105.*udnul+udnur+vartemp331050);
  t[19]=-53.*udnul;
  t[20]=-35.*udnur+vartemp331055;
  t[21]=udirl*(t[19]+t[20]);
  t[22]=3.*udirc*vartemp331094;
  t[23]=3.*t[21];
  t[24]=t[22];
  t[25]=3.*rhoc*vartemp331109;
  t[26]=rhor*(t[18]+t[23]+t[24]);
  t[27]=t[25];
  t[28]=-4094.*udirc*udnuc;
  t[29]=1719.*udirl*udnuc;
  t[30]=1509.*udirr*udnuc;
  t[31]=1719.*udirc*udnul;
  t[32]=-465.*udirl*udnul;
  t[33]=-807.*udirr*udnul;
  t[34]=1509.*udirc*udnur;
  t[35]=-807.*udirl*udnur;
  t[36]=-465.*udirr*udnur;
  t[37]=t[32]+t[33];
  t[38]=t[34]+t[35]+t[36];
  t[39]=t[28]+t[29]+t[30];
  t[40]=t[31]+t[37]+t[38];
  t[41]=rhol*vartemp331048+rhor*vartemp331109;
  t[42]=rhoc*(t[39]+t[40]);
  t[43]=3.*t[41];
  t[44]=gdetr*(t[17]+t[26]+t[27]);
  t[45]=gdetc*(t[42]+t[43]);
  t[46]=gdetl*(t[1]+t[15]+t[16]);
  t[47]=t[44]+t[45];
  t[48]=t[46]+t[47];
  *Fleft=0.0047619047619047619047619047619048*t[48];







  FTYPE vartemp353378=-106.*udnuc;
  FTYPE vartemp353362=102.*udnuc;
  FTYPE vartemp353366=106.*udnuc;
  FTYPE vartemp353371=-155.*udnuc;
  FTYPE vartemp353379=53.*udnul;
  FTYPE vartemp353380=35.*udnur;
  FTYPE vartemp353381=vartemp353378+vartemp353379+vartemp353380;
  FTYPE vartemp353382=udirr*vartemp353381;
  FTYPE vartemp353383=35.*udnul;
  FTYPE vartemp353384=53.*udnur;
  FTYPE vartemp353385=vartemp353378+vartemp353383+vartemp353384;
  FTYPE vartemp353386=udirl*vartemp353385;
  FTYPE vartemp353387=269.*udnuc;
  FTYPE vartemp353388=udnul+udnur;
  FTYPE vartemp353389=-106.*vartemp353388;
  FTYPE vartemp353390=vartemp353387+vartemp353389;
  FTYPE vartemp353391=udirc*vartemp353390;
  FTYPE vartemp353392=vartemp353382+vartemp353386+vartemp353391;
  FTYPE vartemp353404=106.*udnul;
  FTYPE vartemp353405=34.*udnur;
  FTYPE vartemp353406=vartemp353371+vartemp353404+vartemp353405;
  FTYPE vartemp353349=-155.*udirl*udnuc;
  FTYPE vartemp353350=-269.*udirr*udnuc;
  FTYPE vartemp353351=34.*udirl*udnul;
  FTYPE vartemp353352=106.*udirr*udnul;
  FTYPE vartemp353353=503.*udnuc;
  FTYPE vartemp353354=-155.*udnul;
  FTYPE vartemp353355=-269.*udnur;
  FTYPE vartemp353356=vartemp353353+vartemp353354+vartemp353355;
  FTYPE vartemp353357=udirc*vartemp353356;
  FTYPE vartemp353358=udirl+udirr;
  FTYPE vartemp353359=106.*udnur*vartemp353358;
  FTYPE vartemp353360=vartemp353349+vartemp353350+vartemp353351+vartemp353352+vartemp353357+vartemp353359;
  FTYPE vartemp353411=573.*udnuc;
  FTYPE vartemp353412=-269.*udnul;
  FTYPE vartemp353413=-155.*udnur;
  FTYPE vartemp353414=vartemp353411+vartemp353412+vartemp353413;
  FTYPE vartemp353415=udirc*vartemp353414;
  FTYPE vartemp353416=udirr*vartemp353406;
  FTYPE vartemp353417=-269.*udnuc;
  FTYPE vartemp353418=106.*vartemp353388;
  FTYPE vartemp353419=vartemp353417+vartemp353418;
  FTYPE vartemp353420=udirl*vartemp353419;
  FTYPE vartemp353421=vartemp353415+vartemp353416+vartemp353420;
  t[1]=3.*rhoc*vartemp353360;
  t[2]=udirl*(udnul-105.*udnur+vartemp353362);
  t[3]=-35.*udnul;
  t[4]=-53.*udnur+vartemp353366;
  t[5]=udirr*(t[3]+t[4]);
  t[6]=34.*udnul;
  t[7]=106.*udnur+vartemp353371;
  t[8]=udirc*(t[6]+t[7]);
  t[9]=3.*t[5];
  t[10]=3.*t[8];
  t[11]=-3.*rhor*vartemp353392;
  t[12]=rhol*(t[2]+t[9]+t[10]);
  t[13]=t[11];
  t[14]=-3.*rhol*vartemp353392;
  t[15]=-105.*udnul;
  t[16]=-209.*udnur+vartemp353362;
  t[17]=-53.*udnul;
  t[18]=-35.*udnur+vartemp353366;
  t[19]=udirl*(t[17]+t[18]);
  t[20]=3.*udirc*vartemp353406;
  t[21]=3.*t[19];
  t[22]=t[20];
  t[23]=udirr*(t[15]+t[16]);
  t[24]=t[21]+t[22];
  t[25]=3.*rhoc*vartemp353421;
  t[26]=rhor*(t[23]+t[24]);
  t[27]=t[25];
  t[28]=-4094.*udirc*udnuc;
  t[29]=1509.*udirl*udnuc;
  t[30]=1719.*udirr*udnuc;
  t[31]=1509.*udirc*udnul;
  t[32]=-465.*udirl*udnul;
  t[33]=-807.*udirr*udnul;
  t[34]=1719.*udirc*udnur;
  t[35]=-807.*udirl*udnur;
  t[36]=-465.*udirr*udnur;
  t[37]=t[32]+t[33];
  t[38]=t[34]+t[35]+t[36];
  t[39]=t[28]+t[29]+t[30];
  t[40]=t[31]+t[37]+t[38];
  t[41]=rhol*vartemp353360+rhor*vartemp353421;
  t[42]=rhoc*(t[39]+t[40]);
  t[43]=3.*t[41];
  t[44]=gdetr*(t[14]+t[26]+t[27]);
  t[45]=gdetc*(t[42]+t[43]);
  t[46]=gdetl*(t[1]+t[12]+t[13]);
  t[47]=t[44]+t[45];
  t[48]=t[46]+t[47];
  *Fright=0.0047619047619047619047619047619048*t[48];



  return(0);
}







// 5-terms multiplying deconvolution
// See: a2confluxes_5terms.nb
static int fivetermdeconvolution(
                                 FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                 ,FTYPE Bconl, FTYPE Bconc, FTYPE Bconr
                                 ,FTYPE Bcovl, FTYPE Bcovc, FTYPE Bcovr
                                 ,FTYPE uu0l, FTYPE uu0c, FTYPE uu0r
                                 ,FTYPE ud0l, FTYPE ud0c, FTYPE ud0r
                                 ,FTYPE *Fleft, FTYPE *Fright
                                 )
{
  FTYPE t[124]; // 229 total temp variables, 124 t[]'s.  True for both left and right.




  FTYPE vartemp3181=-2314.*uu0l;
  FTYPE vartemp3180=5211.*uu0c;
  FTYPE vartemp3182=-2314.*uu0r;
  FTYPE vartemp3183=vartemp3180+vartemp3181+vartemp3182;
  FTYPE vartemp3185=4045.*uu0c;
  FTYPE vartemp3189=-10422.*uu0c;
  FTYPE vartemp3190=5211.*uu0l;
  FTYPE vartemp3191=4045.*uu0r;
  FTYPE vartemp3192=vartemp3189+vartemp3190+vartemp3191;
  FTYPE vartemp3200=4045.*uu0l;
  FTYPE vartemp3201=5211.*uu0r;
  FTYPE vartemp3202=vartemp3189+vartemp3200+vartemp3201;
  FTYPE vartemp3227=-17885.*uu0c;
  FTYPE vartemp3228=10422.*uu0l;
  FTYPE vartemp3229=5758.*uu0r;
  FTYPE vartemp3230=vartemp3227+vartemp3228+vartemp3229;
  FTYPE vartemp3217=10422.*uu0c;
  FTYPE vartemp3223=-4045.*uu0l;
  FTYPE vartemp3218=-5211.*uu0l;
  FTYPE vartemp3219=-4045.*uu0r;
  FTYPE vartemp3220=vartemp3217+vartemp3218+vartemp3219;
  FTYPE vartemp3240=-23803.*uu0c;
  FTYPE vartemp3209=uu0l+uu0r;
  FTYPE vartemp3143=-17274.*ud0c*uu0c;
  FTYPE vartemp3148=-4038.*ud0r*uu0l;
  FTYPE vartemp3150=-4038.*ud0l*uu0r;
  FTYPE vartemp3167=10422.*ud0c*uu0c;
  FTYPE vartemp3172=2314.*ud0r*uu0l;
  FTYPE vartemp3174=2314.*ud0l*uu0r;
  FTYPE vartemp3157=4045.*ud0r*uu0l;
  FTYPE vartemp3158=17885.*uu0c;
  FTYPE vartemp3163=4045.*ud0l*uu0r;
  FTYPE vartemp3184=ud0l*vartemp3183;
  FTYPE vartemp3186=-1346.*uu0r;
  FTYPE vartemp3187=vartemp3181+vartemp3185+vartemp3186;
  FTYPE vartemp3188=ud0r*vartemp3187;
  FTYPE vartemp3193=ud0c*vartemp3192;
  FTYPE vartemp3194=vartemp3184+vartemp3188+vartemp3193;
  FTYPE vartemp3195=gdetr*vartemp3194;
  FTYPE vartemp3196=ud0r*vartemp3183;
  FTYPE vartemp3197=-1346.*uu0l;
  FTYPE vartemp3198=vartemp3182+vartemp3185+vartemp3197;
  FTYPE vartemp3199=ud0l*vartemp3198;
  FTYPE vartemp3203=ud0c*vartemp3202;
  FTYPE vartemp3204=vartemp3196+vartemp3199+vartemp3203;
  FTYPE vartemp3205=gdetl*vartemp3204;
  FTYPE vartemp3206=ud0r*vartemp3192;
  FTYPE vartemp3207=ud0l*vartemp3202;
  FTYPE vartemp3208=23803.*uu0c;
  FTYPE vartemp3210=-10422.*vartemp3209;
  FTYPE vartemp3211=vartemp3208+vartemp3210;
  FTYPE vartemp3212=ud0c*vartemp3211;
  FTYPE vartemp3213=vartemp3206+vartemp3207+vartemp3212;
  FTYPE vartemp3214=gdetc*vartemp3213;
  FTYPE vartemp3215=vartemp3195+vartemp3205+vartemp3214;
  FTYPE vartemp3246=-5211.*uu0r;
  FTYPE vartemp3247=vartemp3217+vartemp3223+vartemp3246;
  FTYPE vartemp3222=5758.*uu0c;
  FTYPE vartemp3241=10422.*uu0r;
  FTYPE vartemp3292=5758.*uu0l;
  FTYPE vartemp3293=vartemp3227+vartemp3241+vartemp3292;
  FTYPE vartemp3242=vartemp3228+vartemp3240+vartemp3241;
  FTYPE vartemp3248=ud0l*vartemp3247;
  FTYPE vartemp3249=ud0r*vartemp3220;
  FTYPE vartemp3250=10422.*vartemp3209;
  FTYPE vartemp3251=vartemp3240+vartemp3250;
  FTYPE vartemp3252=ud0c*vartemp3251;
  FTYPE vartemp3253=vartemp3248+vartemp3249+vartemp3252;
  FTYPE vartemp3329=10422.*ud0r*uu0l;
  FTYPE vartemp3331=10422.*ud0l*uu0r;
  FTYPE vartemp3221=ud0l*vartemp3220;
  FTYPE vartemp3224=-1339.*uu0r;
  FTYPE vartemp3225=vartemp3222+vartemp3223+vartemp3224;
  FTYPE vartemp3226=ud0r*vartemp3225;
  FTYPE vartemp3231=ud0c*vartemp3230;
  FTYPE vartemp3232=vartemp3221+vartemp3226+vartemp3231;
  FTYPE vartemp3233=gdetr*vartemp3232;
  FTYPE vartemp3234=47221.*uu0c;
  FTYPE vartemp3235=-23803.*uu0l;
  FTYPE vartemp3236=-17885.*uu0r;
  FTYPE vartemp3237=vartemp3234+vartemp3235+vartemp3236;
  FTYPE vartemp3238=ud0c*vartemp3237;
  FTYPE vartemp3239=ud0r*vartemp3230;
  FTYPE vartemp3243=ud0l*vartemp3242;
  FTYPE vartemp3244=vartemp3238+vartemp3239+vartemp3243;
  FTYPE vartemp3245=gdetc*vartemp3244;
  FTYPE vartemp3254=gdetl*vartemp3253;
  FTYPE vartemp3255=vartemp3233+vartemp3245+vartemp3254;
  FTYPE vartemp3288=ud0r*vartemp3247;
  FTYPE vartemp3289=-1339.*uu0l;
  FTYPE vartemp3290=vartemp3219+vartemp3222+vartemp3289;
  FTYPE vartemp3291=ud0l*vartemp3290;
  FTYPE vartemp3294=ud0c*vartemp3293;
  FTYPE vartemp3295=vartemp3288+vartemp3291+vartemp3294;
  FTYPE vartemp3296=gdetl*vartemp3295;
  FTYPE vartemp3297=47991.*uu0c;
  FTYPE vartemp3298=-17885.*uu0l;
  FTYPE vartemp3299=-23803.*uu0r;
  FTYPE vartemp3300=vartemp3297+vartemp3298+vartemp3299;
  FTYPE vartemp3301=ud0c*vartemp3300;
  FTYPE vartemp3302=ud0l*vartemp3293;
  FTYPE vartemp3303=ud0r*vartemp3242;
  FTYPE vartemp3304=vartemp3301+vartemp3302+vartemp3303;
  FTYPE vartemp3305=gdetc*vartemp3304;
  FTYPE vartemp3306=gdetr*vartemp3253;
  FTYPE vartemp3307=vartemp3296+vartemp3305+vartemp3306;
  t[1]=12135.*ud0l*uu0c;
  t[2]=4017.*ud0r*uu0c;
  t[3]=12135.*ud0c*uu0l;
  t[4]=-6942.*ud0l*uu0l;
  t[5]=4017.*ud0c*uu0r;
  t[6]=10.*ud0r*uu0r;
  t[7]=vartemp3143+vartemp3148+vartemp3150;
  t[8]=t[1]+t[2]+t[3];
  t[9]=t[4]+t[5]+t[6]+t[7];
  t[10]=-10422.*ud0l*uu0c;
  t[11]=-5758.*ud0r*uu0c;
  t[12]=5211.*ud0l*uu0l;
  t[13]=1339.*ud0r*uu0r+vartemp3157;
  t[14]=-10422.*uu0l;
  t[15]=-5758.*uu0r+vartemp3158;
  t[16]=t[13];
  t[17]=ud0c*(t[14]+t[15]);
  t[18]=vartemp3163+t[10]+t[11];
  t[19]=t[12]+t[16]+t[17];
  t[20]=gdetc*(t[18]+t[19]);
  t[21]=-5211.*ud0l*uu0c;
  t[22]=-4045.*ud0r*uu0c;
  t[23]=-5211.*ud0c*uu0l;
  t[24]=2314.*ud0l*uu0l;
  t[25]=-4045.*ud0c*uu0r;
  t[26]=1346.*ud0r*uu0r;
  t[27]=vartemp3167+vartemp3172+vartemp3174;
  t[28]=t[21]+t[22]+t[23];
  t[29]=t[24]+t[25]+t[26]+t[27];
  t[30]=gdetl*(t[28]+t[29]);
  t[31]=3.*t[20];
  t[32]=-3.*t[30];
  t[33]=gdetr*(t[8]+t[9]);
  t[34]=t[31]+t[32];
  t[35]=3.*Bcovl*vartemp3215;
  t[36]=-3.*Bcovc*vartemp3255;
  t[37]=Bcovr*(t[33]+t[34]);
  t[38]=t[35]+t[36];
  t[39]=Bconr*(t[37]+t[38]);
  t[40]=4017.*ud0l*uu0c;
  t[41]=12135.*ud0r*uu0c;
  t[42]=4017.*ud0c*uu0l;
  t[43]=2320.*ud0l*uu0l;
  t[44]=12135.*ud0c*uu0r;
  t[45]=-6942.*ud0r*uu0r;
  t[46]=vartemp3143+vartemp3148+vartemp3150;
  t[47]=t[40]+t[41]+t[42];
  t[48]=t[43]+t[44]+t[45]+t[46];
  t[49]=-5758.*ud0l*uu0c;
  t[50]=-10422.*ud0r*uu0c;
  t[51]=1339.*ud0l*uu0l;
  t[52]=5211.*ud0r*uu0r+vartemp3157;
  t[53]=-5758.*uu0l;
  t[54]=-10422.*uu0r+vartemp3158;
  t[55]=t[52];
  t[56]=ud0c*(t[53]+t[54]);
  t[57]=vartemp3163+t[49]+t[50];
  t[58]=t[51]+t[55]+t[56];
  t[59]=gdetc*(t[57]+t[58]);
  t[60]=-4045.*ud0l*uu0c;
  t[61]=-5211.*ud0r*uu0c;
  t[62]=-4045.*ud0c*uu0l;
  t[63]=1346.*ud0l*uu0l;
  t[64]=-5211.*ud0c*uu0r;
  t[65]=2314.*ud0r*uu0r;
  t[66]=vartemp3167+vartemp3172+vartemp3174;
  t[67]=t[60]+t[61]+t[62];
  t[68]=t[63]+t[64]+t[65]+t[66];
  t[69]=gdetr*(t[67]+t[68]);
  t[70]=3.*t[59];
  t[71]=-3.*t[69];
  t[72]=gdetl*(t[47]+t[48]);
  t[73]=t[70]+t[71];
  t[74]=3.*Bcovr*vartemp3215;
  t[75]=-3.*Bcovc*vartemp3307;
  t[76]=Bcovl*(t[72]+t[73]);
  t[77]=t[74]+t[75];
  t[78]=Bconl*(t[76]+t[77]);
  t[79]=Bcovr*vartemp3255+Bcovl*vartemp3307;
  t[80]=-143973.*ud0l*uu0c;
  t[81]=-141663.*ud0r*uu0c;
  t[82]=53655.*ud0l*uu0l;
  t[83]=71409.*ud0r*uu0l;
  t[84]=330670.*uu0c;
  t[85]=-143973.*uu0l;
  t[86]=-141663.*uu0r;
  t[87]=t[83];
  t[88]=ud0c*(t[84]+t[85]+t[86]);
  t[89]=71409.*ud0l*uu0r;
  t[90]=53655.*ud0r*uu0r;
  t[91]=t[80]+t[81]+t[82];
  t[92]=t[87]+t[88]+t[89]+t[90];
  t[93]=47221.*ud0c*uu0c;
  t[94]=-23803.*ud0l*uu0c;
  t[95]=-17885.*ud0r*uu0c;
  t[96]=-23803.*ud0c*uu0l;
  t[97]=10422.*ud0l*uu0l;
  t[98]=-17885.*ud0c*uu0r;
  t[99]=5758.*ud0r*uu0r+vartemp3329+vartemp3331;
  t[100]=t[93]+t[94]+t[95];
  t[101]=t[96]+t[97]+t[98]+t[99];
  t[102]=47991.*ud0c*uu0c;
  t[103]=-17885.*ud0l*uu0c;
  t[104]=-23803.*ud0r*uu0c;
  t[105]=-17885.*ud0c*uu0l;
  t[106]=5758.*ud0l*uu0l;
  t[107]=-23803.*ud0c*uu0r;
  t[108]=10422.*ud0r*uu0r+vartemp3329+vartemp3331;
  t[109]=t[102]+t[103]+t[104];
  t[110]=t[105]+t[106]+t[107]+t[108];
  t[111]=gdetr*(t[100]+t[101]);
  t[112]=gdetl*(t[109]+t[110]);
  t[113]=t[111]+t[112];
  t[114]=gdetc*(t[91]+t[92]);
  t[115]=-3.*t[113];
  t[116]=-3.*t[79];
  t[117]=Bcovc*(t[114]+t[115]);
  t[118]=Bconc*(t[116]+t[117]);
  t[119]=-1.*t[78];
  t[120]=-1.*t[118];
  t[121]=-1.*t[39];
  t[122]=t[119]+t[120];
  t[123]=t[121]+t[122];
  *Fleft=0.00043290043290043290043290043290043*t[123];


  FTYPE vartemp3736=-2314.*uu0l;
  FTYPE vartemp3735=5211.*uu0c;
  FTYPE vartemp3737=-2314.*uu0r;
  FTYPE vartemp3738=vartemp3735+vartemp3736+vartemp3737;
  FTYPE vartemp3740=4045.*uu0c;
  FTYPE vartemp3744=-10422.*uu0c;
  FTYPE vartemp3745=5211.*uu0l;
  FTYPE vartemp3746=4045.*uu0r;
  FTYPE vartemp3747=vartemp3744+vartemp3745+vartemp3746;
  FTYPE vartemp3755=4045.*uu0l;
  FTYPE vartemp3756=5211.*uu0r;
  FTYPE vartemp3757=vartemp3744+vartemp3755+vartemp3756;
  FTYPE vartemp3782=-17885.*uu0c;
  FTYPE vartemp3783=10422.*uu0l;
  FTYPE vartemp3784=5758.*uu0r;
  FTYPE vartemp3785=vartemp3782+vartemp3783+vartemp3784;
  FTYPE vartemp3772=10422.*uu0c;
  FTYPE vartemp3778=-4045.*uu0l;
  FTYPE vartemp3773=-5211.*uu0l;
  FTYPE vartemp3774=-4045.*uu0r;
  FTYPE vartemp3775=vartemp3772+vartemp3773+vartemp3774;
  FTYPE vartemp3795=-23803.*uu0c;
  FTYPE vartemp3764=uu0l+uu0r;
  FTYPE vartemp3722=-17274.*ud0c*uu0c;
  FTYPE vartemp3727=-4038.*ud0r*uu0l;
  FTYPE vartemp3729=-4038.*ud0l*uu0r;
  FTYPE vartemp3711=10422.*ud0c*uu0c;
  FTYPE vartemp3716=2314.*ud0r*uu0l;
  FTYPE vartemp3718=2314.*ud0l*uu0r;
  FTYPE vartemp3701=4045.*ud0r*uu0l;
  FTYPE vartemp3702=17885.*uu0c;
  FTYPE vartemp3707=4045.*ud0l*uu0r;
  FTYPE vartemp3739=ud0l*vartemp3738;
  FTYPE vartemp3741=-1346.*uu0r;
  FTYPE vartemp3742=vartemp3736+vartemp3740+vartemp3741;
  FTYPE vartemp3743=ud0r*vartemp3742;
  FTYPE vartemp3748=ud0c*vartemp3747;
  FTYPE vartemp3749=vartemp3739+vartemp3743+vartemp3748;
  FTYPE vartemp3750=gdetr*vartemp3749;
  FTYPE vartemp3751=ud0r*vartemp3738;
  FTYPE vartemp3752=-1346.*uu0l;
  FTYPE vartemp3753=vartemp3737+vartemp3740+vartemp3752;
  FTYPE vartemp3754=ud0l*vartemp3753;
  FTYPE vartemp3758=ud0c*vartemp3757;
  FTYPE vartemp3759=vartemp3751+vartemp3754+vartemp3758;
  FTYPE vartemp3760=gdetl*vartemp3759;
  FTYPE vartemp3761=ud0r*vartemp3747;
  FTYPE vartemp3762=ud0l*vartemp3757;
  FTYPE vartemp3763=23803.*uu0c;
  FTYPE vartemp3765=-10422.*vartemp3764;
  FTYPE vartemp3766=vartemp3763+vartemp3765;
  FTYPE vartemp3767=ud0c*vartemp3766;
  FTYPE vartemp3768=vartemp3761+vartemp3762+vartemp3767;
  FTYPE vartemp3769=gdetc*vartemp3768;
  FTYPE vartemp3770=vartemp3750+vartemp3760+vartemp3769;
  FTYPE vartemp3801=-5211.*uu0r;
  FTYPE vartemp3802=vartemp3772+vartemp3778+vartemp3801;
  FTYPE vartemp3777=5758.*uu0c;
  FTYPE vartemp3796=10422.*uu0r;
  FTYPE vartemp3847=5758.*uu0l;
  FTYPE vartemp3848=vartemp3782+vartemp3796+vartemp3847;
  FTYPE vartemp3797=vartemp3783+vartemp3795+vartemp3796;
  FTYPE vartemp3803=ud0l*vartemp3802;
  FTYPE vartemp3804=ud0r*vartemp3775;
  FTYPE vartemp3805=10422.*vartemp3764;
  FTYPE vartemp3806=vartemp3795+vartemp3805;
  FTYPE vartemp3807=ud0c*vartemp3806;
  FTYPE vartemp3808=vartemp3803+vartemp3804+vartemp3807;
  FTYPE vartemp3884=10422.*ud0r*uu0l;
  FTYPE vartemp3886=10422.*ud0l*uu0r;
  FTYPE vartemp3776=ud0l*vartemp3775;
  FTYPE vartemp3779=-1339.*uu0r;
  FTYPE vartemp3780=vartemp3777+vartemp3778+vartemp3779;
  FTYPE vartemp3781=ud0r*vartemp3780;
  FTYPE vartemp3786=ud0c*vartemp3785;
  FTYPE vartemp3787=vartemp3776+vartemp3781+vartemp3786;
  FTYPE vartemp3788=gdetr*vartemp3787;
  FTYPE vartemp3789=47991.*uu0c;
  FTYPE vartemp3790=-23803.*uu0l;
  FTYPE vartemp3791=-17885.*uu0r;
  FTYPE vartemp3792=vartemp3789+vartemp3790+vartemp3791;
  FTYPE vartemp3793=ud0c*vartemp3792;
  FTYPE vartemp3794=ud0r*vartemp3785;
  FTYPE vartemp3798=ud0l*vartemp3797;
  FTYPE vartemp3799=vartemp3793+vartemp3794+vartemp3798;
  FTYPE vartemp3800=gdetc*vartemp3799;
  FTYPE vartemp3809=gdetl*vartemp3808;
  FTYPE vartemp3810=vartemp3788+vartemp3800+vartemp3809;
  FTYPE vartemp3843=ud0r*vartemp3802;
  FTYPE vartemp3844=-1339.*uu0l;
  FTYPE vartemp3845=vartemp3774+vartemp3777+vartemp3844;
  FTYPE vartemp3846=ud0l*vartemp3845;
  FTYPE vartemp3849=ud0c*vartemp3848;
  FTYPE vartemp3850=vartemp3843+vartemp3846+vartemp3849;
  FTYPE vartemp3851=gdetl*vartemp3850;
  FTYPE vartemp3852=47221.*uu0c;
  FTYPE vartemp3853=-17885.*uu0l;
  FTYPE vartemp3854=-23803.*uu0r;
  FTYPE vartemp3855=vartemp3852+vartemp3853+vartemp3854;
  FTYPE vartemp3856=ud0c*vartemp3855;
  FTYPE vartemp3857=ud0l*vartemp3848;
  FTYPE vartemp3858=ud0r*vartemp3797;
  FTYPE vartemp3859=vartemp3856+vartemp3857+vartemp3858;
  FTYPE vartemp3860=gdetc*vartemp3859;
  FTYPE vartemp3861=gdetr*vartemp3808;
  FTYPE vartemp3862=vartemp3851+vartemp3860+vartemp3861;
  t[1]=-10422.*ud0l*uu0c;
  t[2]=-5758.*ud0r*uu0c;
  t[3]=5211.*ud0l*uu0l;
  t[4]=1339.*ud0r*uu0r+vartemp3701;
  t[5]=-10422.*uu0l;
  t[6]=-5758.*uu0r+vartemp3702;
  t[7]=t[4];
  t[8]=ud0c*(t[5]+t[6]);
  t[9]=vartemp3707+t[1]+t[2];
  t[10]=t[3]+t[7]+t[8];
  t[11]=gdetc*(t[9]+t[10]);
  t[12]=-5211.*ud0l*uu0c;
  t[13]=-4045.*ud0r*uu0c;
  t[14]=-5211.*ud0c*uu0l;
  t[15]=2314.*ud0l*uu0l;
  t[16]=-4045.*ud0c*uu0r;
  t[17]=1346.*ud0r*uu0r;
  t[18]=vartemp3711+vartemp3716+vartemp3718;
  t[19]=t[12]+t[13]+t[14];
  t[20]=t[15]+t[16]+t[17]+t[18];
  t[21]=gdetl*(t[19]+t[20]);
  t[22]=12135.*ud0l*uu0c;
  t[23]=4017.*ud0r*uu0c;
  t[24]=12135.*ud0c*uu0l;
  t[25]=-6942.*ud0l*uu0l;
  t[26]=4017.*ud0c*uu0r;
  t[27]=2320.*ud0r*uu0r;
  t[28]=vartemp3722+vartemp3727+vartemp3729;
  t[29]=t[22]+t[23]+t[24];
  t[30]=t[25]+t[26]+t[27]+t[28];
  t[31]=-3.*t[21];
  t[32]=gdetr*(t[29]+t[30]);
  t[33]=3.*t[11];
  t[34]=t[31]+t[32];
  t[35]=3.*Bcovl*vartemp3770;
  t[36]=-3.*Bcovc*vartemp3810;
  t[37]=Bcovr*(t[33]+t[34]);
  t[38]=t[35]+t[36];
  t[39]=Bconr*(t[37]+t[38]);
  t[40]=-5758.*ud0l*uu0c;
  t[41]=-10422.*ud0r*uu0c;
  t[42]=1339.*ud0l*uu0l;
  t[43]=5211.*ud0r*uu0r+vartemp3701;
  t[44]=-5758.*uu0l;
  t[45]=-10422.*uu0r+vartemp3702;
  t[46]=t[43];
  t[47]=ud0c*(t[44]+t[45]);
  t[48]=vartemp3707+t[40]+t[41];
  t[49]=t[42]+t[46]+t[47];
  t[50]=gdetc*(t[48]+t[49]);
  t[51]=-4045.*ud0l*uu0c;
  t[52]=-5211.*ud0r*uu0c;
  t[53]=-4045.*ud0c*uu0l;
  t[54]=1346.*ud0l*uu0l;
  t[55]=-5211.*ud0c*uu0r;
  t[56]=2314.*ud0r*uu0r;
  t[57]=vartemp3711+vartemp3716+vartemp3718;
  t[58]=t[51]+t[52]+t[53];
  t[59]=t[54]+t[55]+t[56]+t[57];
  t[60]=gdetr*(t[58]+t[59]);
  t[61]=4017.*ud0l*uu0c;
  t[62]=12135.*ud0r*uu0c;
  t[63]=4017.*ud0c*uu0l;
  t[64]=10.*ud0l*uu0l;
  t[65]=12135.*ud0c*uu0r;
  t[66]=-6942.*ud0r*uu0r;
  t[67]=vartemp3722+vartemp3727+vartemp3729;
  t[68]=t[61]+t[62]+t[63];
  t[69]=t[64]+t[65]+t[66]+t[67];
  t[70]=-3.*t[60];
  t[71]=gdetl*(t[68]+t[69]);
  t[72]=3.*t[50];
  t[73]=t[70]+t[71];
  t[74]=3.*Bcovr*vartemp3770;
  t[75]=-3.*Bcovc*vartemp3862;
  t[76]=Bcovl*(t[72]+t[73]);
  t[77]=t[74]+t[75];
  t[78]=Bconl*(t[76]+t[77]);
  t[79]=Bcovr*vartemp3810+Bcovl*vartemp3862;
  t[80]=-141663.*ud0l*uu0c;
  t[81]=-143973.*ud0r*uu0c;
  t[82]=53655.*ud0l*uu0l;
  t[83]=71409.*ud0r*uu0l;
  t[84]=330670.*uu0c;
  t[85]=-141663.*uu0l;
  t[86]=-143973.*uu0r;
  t[87]=t[83];
  t[88]=ud0c*(t[84]+t[85]+t[86]);
  t[89]=71409.*ud0l*uu0r;
  t[90]=53655.*ud0r*uu0r;
  t[91]=t[80]+t[81]+t[82];
  t[92]=t[87]+t[88]+t[89]+t[90];
  t[93]=47991.*ud0c*uu0c;
  t[94]=-23803.*ud0l*uu0c;
  t[95]=-17885.*ud0r*uu0c;
  t[96]=-23803.*ud0c*uu0l;
  t[97]=10422.*ud0l*uu0l;
  t[98]=-17885.*ud0c*uu0r;
  t[99]=5758.*ud0r*uu0r+vartemp3884+vartemp3886;
  t[100]=t[93]+t[94]+t[95];
  t[101]=t[96]+t[97]+t[98]+t[99];
  t[102]=47221.*ud0c*uu0c;
  t[103]=-17885.*ud0l*uu0c;
  t[104]=-23803.*ud0r*uu0c;
  t[105]=-17885.*ud0c*uu0l;
  t[106]=5758.*ud0l*uu0l;
  t[107]=-23803.*ud0c*uu0r;
  t[108]=10422.*ud0r*uu0r+vartemp3884+vartemp3886;
  t[109]=t[102]+t[103]+t[104];
  t[110]=t[105]+t[106]+t[107]+t[108];
  t[111]=gdetr*(t[100]+t[101]);
  t[112]=gdetl*(t[109]+t[110]);
  t[113]=t[111]+t[112];
  t[114]=gdetc*(t[91]+t[92]);
  t[115]=-3.*t[113];
  t[116]=-3.*t[79];
  t[117]=Bcovc*(t[114]+t[115]);
  t[118]=Bconc*(t[116]+t[117]);
  t[119]=-1.*t[78];
  t[120]=-1.*t[118];
  t[121]=-1.*t[39];
  t[122]=t[119]+t[120];
  t[123]=t[121]+t[122];
  *Fright=0.00043290043290043290043290043290043*t[123];


  return(0);


}

// 6-terms multiplying deconvolution
// See: a2confluxes_6terms.nb
static int sixtermdeconvolution(
                                FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                ,FTYPE Bconl, FTYPE Bconc, FTYPE Bconr
                                ,FTYPE Bcovl, FTYPE Bcovc, FTYPE Bcovr
                                ,FTYPE uu0l, FTYPE uu0c, FTYPE uu0r
                                ,FTYPE ud0l, FTYPE ud0c, FTYPE ud0r
                                ,FTYPE udnul, FTYPE udnuc, FTYPE udnur
                                ,FTYPE *Fleft, FTYPE *Fright
                                )
{
  FTYPE t[159]; // 381 total temps vars, 159 t[]'s for both left/right


 
  FTYPE vartemp6857=-405594.*uu0c;
  FTYPE vartemp6858=202797.*uu0l;
  FTYPE vartemp6859=172715.*uu0r;
  FTYPE vartemp6860=vartemp6857+vartemp6858+vartemp6859;
  FTYPE vartemp6853=172715.*uu0l;
  FTYPE vartemp6847=743445.*uu0c;
  FTYPE vartemp6848=-405594.*uu0l;
  FTYPE vartemp6849=-285266.*uu0r;
  FTYPE vartemp6850=vartemp6847+vartemp6848+vartemp6849;
  FTYPE vartemp6868=32553.*uu0c;
  FTYPE vartemp6869=uu0l+uu0r;
  FTYPE vartemp6870=-15022.*vartemp6869;
  FTYPE vartemp6871=vartemp6868+vartemp6870;
  FTYPE vartemp6852=-285266.*uu0c;
  FTYPE vartemp6865=202797.*uu0r;
  FTYPE vartemp6866=vartemp6853+vartemp6857+vartemp6865;
  FTYPE vartemp6864=udnur*vartemp6860;
  FTYPE vartemp6867=udnul*vartemp6866;
  FTYPE vartemp6872=27.*udnuc*vartemp6871;
  FTYPE vartemp6873=vartemp6864+vartemp6867+vartemp6872;
  FTYPE vartemp6886=-285266.*uu0l;
  FTYPE vartemp6887=-405594.*uu0r;
  FTYPE vartemp6888=vartemp6847+vartemp6886+vartemp6887;
  FTYPE vartemp6876=65106.*uu0c;
  FTYPE vartemp6898=-27535.*uu0l;
  FTYPE vartemp6899=-32553.*uu0r;
  FTYPE vartemp6900=vartemp6876+vartemp6898+vartemp6899;
  FTYPE vartemp6877=-32553.*uu0l;
  FTYPE vartemp6878=-27535.*uu0r;
  FTYPE vartemp6879=vartemp6876+vartemp6877+vartemp6878;
  FTYPE vartemp6875=udnur*vartemp6850;
  FTYPE vartemp6880=-27.*udnuc*vartemp6879;
  FTYPE vartemp6881=27.*udnul*vartemp6871;
  FTYPE vartemp6882=vartemp6875+vartemp6880+vartemp6881;
  FTYPE vartemp6897=udnul*vartemp6888;
  FTYPE vartemp6901=-27.*udnuc*vartemp6900;
  FTYPE vartemp6902=27.*udnur*vartemp6871;
  FTYPE vartemp6903=vartemp6897+vartemp6901+vartemp6902;
  FTYPE vartemp6960=405594.*uu0c;
  FTYPE vartemp6966=-172715.*uu0l;
  FTYPE vartemp6961=-202797.*uu0l;
  FTYPE vartemp6962=-172715.*uu0r;
  FTYPE vartemp6963=vartemp6960+vartemp6961+vartemp6962;
  FTYPE vartemp6970=-743445.*uu0c;
  FTYPE vartemp6971=405594.*uu0l;
  FTYPE vartemp6972=285266.*uu0r;
  FTYPE vartemp6973=vartemp6970+vartemp6971+vartemp6972;
  FTYPE vartemp6981=-878931.*uu0c;
  FTYPE vartemp6982=405594.*vartemp6869;
  FTYPE vartemp6983=vartemp6981+vartemp6982;
  FTYPE vartemp6999=-1.215918e6*uu0c;
  FTYPE vartemp7000=743445.*uu0l;
  FTYPE vartemp7001=397619.*uu0r;
  FTYPE vartemp7002=vartemp6999+vartemp7000+vartemp7001;
  FTYPE vartemp6987=27.*udnuc*vartemp6879;
  FTYPE vartemp6988=udnur*vartemp6973;
  FTYPE vartemp6989=udnul*vartemp6983;
  FTYPE vartemp6990=vartemp6987+vartemp6988+vartemp6989;
  FTYPE vartemp6907=udnul*vartemp6900;
  FTYPE vartemp6908=udnur*vartemp6879;
  FTYPE vartemp6909=vartemp6907+vartemp6908;
  FTYPE vartemp6938=1.215918e6*uu0c;
  FTYPE vartemp6939=-743445.*uu0l;
  FTYPE vartemp6940=-397619.*uu0r;
  FTYPE vartemp6941=vartemp6938+vartemp6939+vartemp6940;
  FTYPE vartemp6942=udnuc*vartemp6941;
  FTYPE vartemp7023=285266.*uu0l;
  FTYPE vartemp7010=32553.*uu0l;
  FTYPE vartemp7011=22517.*uu0r;
  FTYPE vartemp7012=vartemp7010+vartemp7011;
  FTYPE vartemp7053=-3.206285e6*uu0c;
  FTYPE vartemp7035=-3.825163e6*uu0c;
  FTYPE vartemp7066=1.757862e6*uu0l;
  FTYPE vartemp7034=27.*vartemp6909;
  FTYPE vartemp7036=1.757862e6*vartemp6869;
  FTYPE vartemp7037=vartemp7035+vartemp7036;
  FTYPE vartemp7038=udnuc*vartemp7037;
  FTYPE vartemp7039=vartemp7034+vartemp7038;
  FTYPE vartemp7051=udnur*vartemp6941;
  FTYPE vartemp7052=27.*udnul*vartemp6879;
  FTYPE vartemp7054=54.*vartemp7012;
  FTYPE vartemp7055=vartemp7053+vartemp7054;
  FTYPE vartemp7056=udnuc*vartemp7055;
  FTYPE vartemp7057=vartemp7051+vartemp7052+vartemp7056;
  FTYPE vartemp6851=udnuc*vartemp6850;
  FTYPE vartemp6854=95053.*uu0r;
  FTYPE vartemp6855=vartemp6852+vartemp6853+vartemp6854;
  FTYPE vartemp6856=udnur*vartemp6855;
  FTYPE vartemp6861=udnul*vartemp6860;
  FTYPE vartemp6862=vartemp6851+vartemp6856+vartemp6861;
  FTYPE vartemp6863=ud0r*vartemp6862;
  FTYPE vartemp6874=ud0l*vartemp6873;
  FTYPE vartemp6883=ud0c*vartemp6882;
  FTYPE vartemp6884=vartemp6863+vartemp6874+vartemp6883;
  FTYPE vartemp6885=gdetr*vartemp6884;
  FTYPE vartemp6889=udnuc*vartemp6888;
  FTYPE vartemp6890=95053.*uu0l;
  FTYPE vartemp6891=vartemp6852+vartemp6859+vartemp6890;
  FTYPE vartemp6892=udnul*vartemp6891;
  FTYPE vartemp6893=udnur*vartemp6866;
  FTYPE vartemp6894=vartemp6889+vartemp6892+vartemp6893;
  FTYPE vartemp6895=ud0l*vartemp6894;
  FTYPE vartemp6896=ud0r*vartemp6873;
  FTYPE vartemp6904=ud0c*vartemp6903;
  FTYPE vartemp6905=vartemp6895+vartemp6896+vartemp6904;
  FTYPE vartemp6906=gdetl*vartemp6905;
  FTYPE vartemp6910=-27.*vartemp6909;
  FTYPE vartemp6911=3.825163e6*uu0c;
  FTYPE vartemp6912=-1.757862e6*vartemp6869;
  FTYPE vartemp6913=vartemp6911+vartemp6912;
  FTYPE vartemp6914=udnuc*vartemp6913;
  FTYPE vartemp6915=vartemp6910+vartemp6914;
  FTYPE vartemp6916=ud0c*vartemp6915;
  FTYPE vartemp6917=ud0r*vartemp6882;
  FTYPE vartemp6918=ud0l*vartemp6903;
  FTYPE vartemp6919=vartemp6916+vartemp6917+vartemp6918;
  FTYPE vartemp6920=gdetc*vartemp6919;
  FTYPE vartemp6921=vartemp6885+vartemp6906+vartemp6920;
  FTYPE vartemp6923=-1.192857e6*udnuc*uu0c;
  FTYPE vartemp6928=-285159.*udnur*uu0l;
  FTYPE vartemp6930=-285159.*udnul*uu0r;
  FTYPE vartemp6947=743445.*udnuc*uu0c;
  FTYPE vartemp6952=172715.*udnur*uu0l;
  FTYPE vartemp6954=172715.*udnul*uu0r;
  FTYPE vartemp6937=285266.*udnur*uu0l;
  FTYPE vartemp6943=285266.*udnul*uu0r;
  FTYPE vartemp6977=-202797.*uu0r;
  FTYPE vartemp6978=vartemp6960+vartemp6966+vartemp6977;
  FTYPE vartemp6965=285266.*uu0c;
  FTYPE vartemp7024=405594.*uu0r;
  FTYPE vartemp7025=vartemp6970+vartemp7023+vartemp7024;
  FTYPE vartemp7022=27.*udnuc*vartemp6900;
  FTYPE vartemp7026=udnul*vartemp7025;
  FTYPE vartemp7027=-32553.*uu0c;
  FTYPE vartemp7028=15022.*vartemp6869;
  FTYPE vartemp7029=vartemp7027+vartemp7028;
  FTYPE vartemp7030=27.*udnur*vartemp7029;
  FTYPE vartemp7031=vartemp7022+vartemp7026+vartemp7030;
  FTYPE vartemp6979=udnul*vartemp6978;
  FTYPE vartemp6980=udnur*vartemp6963;
  FTYPE vartemp6984=udnuc*vartemp6983;
  FTYPE vartemp6985=vartemp6979+vartemp6980+vartemp6984;
  FTYPE vartemp6995=397619.*uu0c;
  FTYPE vartemp7127=397619.*uu0l;
  FTYPE vartemp7128=743445.*uu0r;
  FTYPE vartemp7129=vartemp6999+vartemp7127+vartemp7128;
  FTYPE vartemp7009=3.206285e6*uu0c;
  FTYPE vartemp7136=32553.*uu0r;
  FTYPE vartemp7032=ud0l*vartemp7031;
  FTYPE vartemp7033=ud0r*vartemp6990;
  FTYPE vartemp7040=ud0c*vartemp7039;
  FTYPE vartemp7041=vartemp7032+vartemp7033+vartemp7040;
  FTYPE vartemp7103=-397619.*uu0l;
  FTYPE vartemp7104=-743445.*uu0r;
  FTYPE vartemp7105=vartemp6938+vartemp7103+vartemp7104;
  FTYPE vartemp7106=udnuc*vartemp7105;
  FTYPE vartemp7043=-397619.*uu0c;
  FTYPE vartemp7140=22517.*uu0l;
  FTYPE vartemp7141=vartemp7136+vartemp7140;
  FTYPE vartemp7070=1.757862e6*uu0r;
  FTYPE vartemp7071=vartemp7035+vartemp7066+vartemp7070;
  FTYPE vartemp7159=udnul*vartemp7105;
  FTYPE vartemp7160=27.*udnur*vartemp6900;
  FTYPE vartemp7161=54.*vartemp7141;
  FTYPE vartemp7162=vartemp7053+vartemp7161;
  FTYPE vartemp7163=udnuc*vartemp7162;
  FTYPE vartemp7164=vartemp7159+vartemp7160+vartemp7163;
  FTYPE vartemp7042=gdetl*vartemp7041;
  FTYPE vartemp7044=94946.*uu0r;
  FTYPE vartemp7045=vartemp7023+vartemp7043+vartemp7044;
  FTYPE vartemp7046=udnur*vartemp7045;
  FTYPE vartemp7047=udnul*vartemp6973;
  FTYPE vartemp7048=vartemp6942+vartemp7046+vartemp7047;
  FTYPE vartemp7049=ud0r*vartemp7048;
  FTYPE vartemp7050=ud0l*vartemp6990;
  FTYPE vartemp7058=ud0c*vartemp7057;
  FTYPE vartemp7059=vartemp7049+vartemp7050+vartemp7058;
  FTYPE vartemp7060=gdetr*vartemp7059;
  FTYPE vartemp7061=7.645321e6*uu0c;
  FTYPE vartemp7062=-3.825163e6*uu0l;
  FTYPE vartemp7063=-3.206285e6*uu0r;
  FTYPE vartemp7064=vartemp7061+vartemp7062+vartemp7063;
  FTYPE vartemp7065=udnuc*vartemp7064;
  FTYPE vartemp7067=1.215918e6*uu0r;
  FTYPE vartemp7068=vartemp7053+vartemp7066+vartemp7067;
  FTYPE vartemp7069=udnur*vartemp7068;
  FTYPE vartemp7072=udnul*vartemp7071;
  FTYPE vartemp7073=vartemp7065+vartemp7069+vartemp7072;
  FTYPE vartemp7074=ud0c*vartemp7073;
  FTYPE vartemp7075=ud0l*vartemp7039;
  FTYPE vartemp7076=ud0r*vartemp7057;
  FTYPE vartemp7077=vartemp7074+vartemp7075+vartemp7076;
  FTYPE vartemp7078=gdetc*vartemp7077;
  FTYPE vartemp7079=vartemp7042+vartemp7060+vartemp7078;
  FTYPE vartemp7151=gdetr*vartemp7041;
  FTYPE vartemp7152=94946.*uu0l;
  FTYPE vartemp7153=vartemp6972+vartemp7043+vartemp7152;
  FTYPE vartemp7154=udnul*vartemp7153;
  FTYPE vartemp7155=udnur*vartemp7025;
  FTYPE vartemp7156=vartemp7106+vartemp7154+vartemp7155;
  FTYPE vartemp7157=ud0l*vartemp7156;
  FTYPE vartemp7158=ud0r*vartemp7031;
  FTYPE vartemp7165=ud0c*vartemp7164;
  FTYPE vartemp7166=vartemp7157+vartemp7158+vartemp7165;
  FTYPE vartemp7167=gdetl*vartemp7166;
  FTYPE vartemp7168=7.655331e6*uu0c;
  FTYPE vartemp7169=-3.206285e6*uu0l;
  FTYPE vartemp7170=-3.825163e6*uu0r;
  FTYPE vartemp7171=vartemp7168+vartemp7169+vartemp7170;
  FTYPE vartemp7172=udnuc*vartemp7171;
  FTYPE vartemp7173=1.215918e6*uu0l;
  FTYPE vartemp7174=vartemp7053+vartemp7070+vartemp7173;
  FTYPE vartemp7175=udnul*vartemp7174;
  FTYPE vartemp7176=udnur*vartemp7071;
  FTYPE vartemp7177=vartemp7172+vartemp7175+vartemp7176;
  FTYPE vartemp7178=ud0c*vartemp7177;
  FTYPE vartemp7179=ud0r*vartemp7039;
  FTYPE vartemp7180=ud0l*vartemp7164;
  FTYPE vartemp7181=vartemp7178+vartemp7179+vartemp7180;
  FTYPE vartemp7182=gdetc*vartemp7181;
  FTYPE vartemp7183=vartemp7151+vartemp7167+vartemp7182;
  FTYPE vartemp7211=1.757862e6*udnur*uu0l;
  FTYPE vartemp7213=1.757862e6*udnul*uu0r;
  t[1]=-3.*Bcovl*vartemp6921;
  t[2]=855798.*udnul*uu0c;
  t[3]=284838.*udnur*uu0c;
  t[4]=855798.*udnuc*uu0l;
  t[5]=-518145.*udnul*uu0l;
  t[6]=284838.*udnuc*uu0r;
  t[7]=191.*udnur*uu0r;
  t[8]=vartemp6923+vartemp6928+vartemp6930;
  t[9]=t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=-743445.*udnul*uu0c;
  t[12]=-397619.*udnur*uu0c;
  t[13]=405594.*udnul*uu0l;
  t[14]=94946.*udnur*uu0r+vartemp6937+vartemp6942+vartemp6943;
  t[15]=t[11]+t[12]+t[13]+t[14];
  t[16]=-405594.*udnul*uu0c;
  t[17]=-285266.*udnur*uu0c;
  t[18]=-405594.*udnuc*uu0l;
  t[19]=202797.*udnul*uu0l;
  t[20]=-285266.*udnuc*uu0r;
  t[21]=95053.*udnur*uu0r;
  t[22]=vartemp6947+vartemp6952+vartemp6954;
  t[23]=t[16]+t[17]+t[18];
  t[24]=t[19]+t[20]+t[21]+t[22];
  t[25]=ud0l*(t[23]+t[24]);
  t[26]=3.*ud0c*t[15];
  t[27]=-3.*t[25];
  t[28]=ud0r*(t[9]+t[10]);
  t[29]=t[26]+t[27];
  t[30]=udnul*vartemp6963;
  t[31]=udnur*(-95053.*uu0r+vartemp6965+vartemp6966);
  t[32]=udnuc*vartemp6973;
  t[33]=ud0l*vartemp6985+ud0c*vartemp6990;
  t[34]=ud0r*(t[30]+t[31]+t[32]);
  t[35]=t[33];
  t[36]=gdetl*(t[34]+t[35]);
  t[37]=ud0l*vartemp6882;
  t[38]=udnul*vartemp6850;
  t[39]=udnur*(-94946.*uu0r+vartemp6886+vartemp6995);
  t[40]=udnuc*vartemp7002;
  t[41]=-27.*udnul*vartemp6879;
  t[42]=udnur*vartemp7002;
  t[43]=udnuc*(vartemp7009-54.*vartemp7012);
  t[44]=ud0r*(t[38]+t[39]+t[40]);
  t[45]=ud0c*(t[41]+t[42]+t[43]);
  t[46]=gdetc*(t[37]+t[44]+t[45]);
  t[47]=3.*t[36];
  t[48]=-3.*t[46];
  t[49]=gdetr*(t[28]+t[29]);
  t[50]=t[47]+t[48];
  t[51]=3.*Bcovc*vartemp7079;
  t[52]=Bcovr*(t[49]+t[50]);
  t[53]=t[51];
  t[54]=-3.*Bcovr*vartemp6921;
  t[55]=ud0r*vartemp6985;
  t[56]=udnul*(-95053.*uu0l+vartemp6962+vartemp6965);
  t[57]=udnur*vartemp6978+udnuc*vartemp7025;
  t[58]=ud0c*vartemp7031;
  t[59]=ud0l*(t[56]+t[57]);
  t[60]=t[58];
  t[61]=gdetr*(t[55]+t[59]+t[60]);
  t[62]=284838.*udnul*uu0c;
  t[63]=855798.*udnur*uu0c;
  t[64]=284838.*udnuc*uu0l;
  t[65]=-29839.*udnul*uu0l;
  t[66]=855798.*udnuc*uu0r;
  t[67]=-518145.*udnur*uu0r;
  t[68]=vartemp6923+vartemp6928+vartemp6930;
  t[69]=t[62]+t[63]+t[64];
  t[70]=t[65]+t[66]+t[67]+t[68];
  t[71]=-285266.*udnul*uu0c;
  t[72]=-405594.*udnur*uu0c;
  t[73]=-285266.*udnuc*uu0l;
  t[74]=95053.*udnul*uu0l;
  t[75]=-405594.*udnuc*uu0r;
  t[76]=202797.*udnur*uu0r;
  t[77]=vartemp6947+vartemp6952+vartemp6954;
  t[78]=t[71]+t[72]+t[73];
  t[79]=t[74]+t[75]+t[76]+t[77];
  t[80]=ud0r*(t[78]+t[79]);
  t[81]=-397619.*udnul*uu0c;
  t[82]=-743445.*udnur*uu0c;
  t[83]=94946.*udnul*uu0l;
  t[84]=405594.*udnur*uu0r+vartemp6937+vartemp6943+vartemp7106;
  t[85]=t[81]+t[82]+t[83]+t[84];
  t[86]=-3.*t[80];
  t[87]=3.*ud0c*t[85];
  t[88]=ud0l*(t[69]+t[70]);
  t[89]=t[86]+t[87];
  t[90]=ud0r*vartemp6903;
  t[91]=udnur*vartemp6888;
  t[92]=udnul*(-94946.*uu0l+vartemp6849+vartemp6995);
  t[93]=udnuc*vartemp7129;
  t[94]=udnul*vartemp7129;
  t[95]=-65106.*uu0c;
  t[96]=27535.*uu0l+vartemp7136;
  t[97]=udnur*(t[95]+t[96]);
  t[98]=udnuc*(vartemp7009-54.*vartemp7141);
  t[99]=27.*t[97];
  t[100]=t[98];
  t[101]=ud0l*(t[91]+t[92]+t[93]);
  t[102]=ud0c*(t[94]+t[99]+t[100]);
  t[103]=gdetc*(t[90]+t[101]+t[102]);
  t[104]=gdetl*(t[88]+t[89]);
  t[105]=-3.*t[103];
  t[106]=3.*t[61];
  t[107]=t[104]+t[105];
  t[108]=3.*Bcovc*vartemp7183;
  t[109]=Bcovl*(t[106]+t[107]);
  t[110]=t[108];
  t[111]=Bcovr*vartemp7079+Bcovl*vartemp7183;
  t[112]=gdetr*vartemp7077+gdetl*vartemp7181;
  t[113]=-5.0200666e7*uu0c;
  t[114]=2.2965993e7*uu0l;
  t[115]=2.2935963e7*uu0r;
  t[116]=7.655331e6*udnul*uu0c;
  t[117]=7.645321e6*udnur*uu0c;
  t[118]=-3.206285e6*udnul*uu0l;
  t[119]=-3.825163e6*udnur*uu0l;
  t[120]=-3.825163e6*udnul*uu0r;
  t[121]=-3.206285e6*udnur*uu0r;
  t[122]=t[116]+t[117]+t[118];
  t[123]=t[119]+t[120]+t[121];
  t[124]=t[122]+t[123];
  t[125]=udnuc*(t[113]+t[114]+t[115]);
  t[126]=3.*t[124];
  t[127]=7.645321e6*udnuc*uu0c;
  t[128]=-3.825163e6*udnul*uu0c;
  t[129]=-3.206285e6*udnur*uu0c;
  t[130]=-3.825163e6*udnuc*uu0l;
  t[131]=1.757862e6*udnul*uu0l;
  t[132]=-3.206285e6*udnuc*uu0r;
  t[133]=1.215918e6*udnur*uu0r+vartemp7211+vartemp7213;
  t[134]=t[127]+t[128]+t[129];
  t[135]=t[130]+t[131]+t[132]+t[133];
  t[136]=7.655331e6*udnuc*uu0c;
  t[137]=-3.206285e6*udnul*uu0c;
  t[138]=-3.825163e6*udnur*uu0c;
  t[139]=-3.206285e6*udnuc*uu0l;
  t[140]=1.215918e6*udnul*uu0l;
  t[141]=-3.825163e6*udnuc*uu0r;
  t[142]=1.757862e6*udnur*uu0r+vartemp7211+vartemp7213;
  t[143]=t[136]+t[137]+t[138];
  t[144]=t[139]+t[140]+t[141]+t[142];
  t[145]=ud0r*(t[134]+t[135]);
  t[146]=ud0l*(t[143]+t[144]);
  t[147]=t[145]+t[146];
  t[148]=ud0c*(t[125]+t[126]);
  t[149]=3.*t[147];
  t[150]=3.*t[112];
  t[151]=gdetc*(t[148]+t[149]);
  t[152]=3.*t[111];
  t[153]=Bcovc*(t[150]+t[151]);
  t[154]=Bconl*(t[54]+t[109]+t[110]);
  t[155]=Bconc*(t[152]+t[153]);
  t[156]=Bconr*(t[1]+t[52]+t[53]);
  t[157]=t[154]+t[155];
  t[158]=t[156]+t[157];
  *Fleft=0.000033300033300033300033300033300033*t[158];




  FTYPE vartemp7861=-405594.*uu0c;
  FTYPE vartemp7862=202797.*uu0l;
  FTYPE vartemp7863=172715.*uu0r;
  FTYPE vartemp7864=vartemp7861+vartemp7862+vartemp7863;
  FTYPE vartemp7857=172715.*uu0l;
  FTYPE vartemp7851=743445.*uu0c;
  FTYPE vartemp7852=-405594.*uu0l;
  FTYPE vartemp7853=-285266.*uu0r;
  FTYPE vartemp7854=vartemp7851+vartemp7852+vartemp7853;
  FTYPE vartemp7872=32553.*uu0c;
  FTYPE vartemp7873=uu0l+uu0r;
  FTYPE vartemp7874=-15022.*vartemp7873;
  FTYPE vartemp7875=vartemp7872+vartemp7874;
  FTYPE vartemp7856=-285266.*uu0c;
  FTYPE vartemp7869=202797.*uu0r;
  FTYPE vartemp7870=vartemp7857+vartemp7861+vartemp7869;
  FTYPE vartemp7868=udnur*vartemp7864;
  FTYPE vartemp7871=udnul*vartemp7870;
  FTYPE vartemp7876=27.*udnuc*vartemp7875;
  FTYPE vartemp7877=vartemp7868+vartemp7871+vartemp7876;
  FTYPE vartemp7890=-285266.*uu0l;
  FTYPE vartemp7891=-405594.*uu0r;
  FTYPE vartemp7892=vartemp7851+vartemp7890+vartemp7891;
  FTYPE vartemp7880=65106.*uu0c;
  FTYPE vartemp7902=-27535.*uu0l;
  FTYPE vartemp7903=-32553.*uu0r;
  FTYPE vartemp7904=vartemp7880+vartemp7902+vartemp7903;
  FTYPE vartemp7881=-32553.*uu0l;
  FTYPE vartemp7882=-27535.*uu0r;
  FTYPE vartemp7883=vartemp7880+vartemp7881+vartemp7882;
  FTYPE vartemp7879=udnur*vartemp7854;
  FTYPE vartemp7884=-27.*udnuc*vartemp7883;
  FTYPE vartemp7885=27.*udnul*vartemp7875;
  FTYPE vartemp7886=vartemp7879+vartemp7884+vartemp7885;
  FTYPE vartemp7901=udnul*vartemp7892;
  FTYPE vartemp7905=-27.*udnuc*vartemp7904;
  FTYPE vartemp7906=27.*udnur*vartemp7875;
  FTYPE vartemp7907=vartemp7901+vartemp7905+vartemp7906;
  FTYPE vartemp7964=405594.*uu0c;
  FTYPE vartemp7970=-172715.*uu0l;
  FTYPE vartemp7965=-202797.*uu0l;
  FTYPE vartemp7966=-172715.*uu0r;
  FTYPE vartemp7967=vartemp7964+vartemp7965+vartemp7966;
  FTYPE vartemp7974=-743445.*uu0c;
  FTYPE vartemp7975=405594.*uu0l;
  FTYPE vartemp7976=285266.*uu0r;
  FTYPE vartemp7977=vartemp7974+vartemp7975+vartemp7976;
  FTYPE vartemp7985=-878931.*uu0c;
  FTYPE vartemp7986=405594.*vartemp7873;
  FTYPE vartemp7987=vartemp7985+vartemp7986;
  FTYPE vartemp8003=-1.215918e6*uu0c;
  FTYPE vartemp8004=743445.*uu0l;
  FTYPE vartemp8005=397619.*uu0r;
  FTYPE vartemp8006=vartemp8003+vartemp8004+vartemp8005;
  FTYPE vartemp7991=27.*udnuc*vartemp7883;
  FTYPE vartemp7992=udnur*vartemp7977;
  FTYPE vartemp7993=udnul*vartemp7987;
  FTYPE vartemp7994=vartemp7991+vartemp7992+vartemp7993;
  FTYPE vartemp7911=udnul*vartemp7904;
  FTYPE vartemp7912=udnur*vartemp7883;
  FTYPE vartemp7913=vartemp7911+vartemp7912;
  FTYPE vartemp7942=1.215918e6*uu0c;
  FTYPE vartemp7943=-743445.*uu0l;
  FTYPE vartemp7944=-397619.*uu0r;
  FTYPE vartemp7945=vartemp7942+vartemp7943+vartemp7944;
  FTYPE vartemp7946=udnuc*vartemp7945;
  FTYPE vartemp8027=285266.*uu0l;
  FTYPE vartemp8014=32553.*uu0l;
  FTYPE vartemp8015=22517.*uu0r;
  FTYPE vartemp8016=vartemp8014+vartemp8015;
  FTYPE vartemp8057=-3.206285e6*uu0c;
  FTYPE vartemp8039=-3.825163e6*uu0c;
  FTYPE vartemp8070=1.757862e6*uu0l;
  FTYPE vartemp8038=27.*vartemp7913;
  FTYPE vartemp8040=1.757862e6*vartemp7873;
  FTYPE vartemp8041=vartemp8039+vartemp8040;
  FTYPE vartemp8042=udnuc*vartemp8041;
  FTYPE vartemp8043=vartemp8038+vartemp8042;
  FTYPE vartemp8055=udnur*vartemp7945;
  FTYPE vartemp8056=27.*udnul*vartemp7883;
  FTYPE vartemp8058=54.*vartemp8016;
  FTYPE vartemp8059=vartemp8057+vartemp8058;
  FTYPE vartemp8060=udnuc*vartemp8059;
  FTYPE vartemp8061=vartemp8055+vartemp8056+vartemp8060;
  FTYPE vartemp7855=udnuc*vartemp7854;
  FTYPE vartemp7858=95053.*uu0r;
  FTYPE vartemp7859=vartemp7856+vartemp7857+vartemp7858;
  FTYPE vartemp7860=udnur*vartemp7859;
  FTYPE vartemp7865=udnul*vartemp7864;
  FTYPE vartemp7866=vartemp7855+vartemp7860+vartemp7865;
  FTYPE vartemp7867=ud0r*vartemp7866;
  FTYPE vartemp7878=ud0l*vartemp7877;
  FTYPE vartemp7887=ud0c*vartemp7886;
  FTYPE vartemp7888=vartemp7867+vartemp7878+vartemp7887;
  FTYPE vartemp7889=gdetr*vartemp7888;
  FTYPE vartemp7893=udnuc*vartemp7892;
  FTYPE vartemp7894=95053.*uu0l;
  FTYPE vartemp7895=vartemp7856+vartemp7863+vartemp7894;
  FTYPE vartemp7896=udnul*vartemp7895;
  FTYPE vartemp7897=udnur*vartemp7870;
  FTYPE vartemp7898=vartemp7893+vartemp7896+vartemp7897;
  FTYPE vartemp7899=ud0l*vartemp7898;
  FTYPE vartemp7900=ud0r*vartemp7877;
  FTYPE vartemp7908=ud0c*vartemp7907;
  FTYPE vartemp7909=vartemp7899+vartemp7900+vartemp7908;
  FTYPE vartemp7910=gdetl*vartemp7909;
  FTYPE vartemp7914=-27.*vartemp7913;
  FTYPE vartemp7915=3.825163e6*uu0c;
  FTYPE vartemp7916=-1.757862e6*vartemp7873;
  FTYPE vartemp7917=vartemp7915+vartemp7916;
  FTYPE vartemp7918=udnuc*vartemp7917;
  FTYPE vartemp7919=vartemp7914+vartemp7918;
  FTYPE vartemp7920=ud0c*vartemp7919;
  FTYPE vartemp7921=ud0r*vartemp7886;
  FTYPE vartemp7922=ud0l*vartemp7907;
  FTYPE vartemp7923=vartemp7920+vartemp7921+vartemp7922;
  FTYPE vartemp7924=gdetc*vartemp7923;
  FTYPE vartemp7925=vartemp7889+vartemp7910+vartemp7924;
  FTYPE vartemp7927=-1.192857e6*udnuc*uu0c;
  FTYPE vartemp7932=-285159.*udnur*uu0l;
  FTYPE vartemp7934=-285159.*udnul*uu0r;
  FTYPE vartemp7951=743445.*udnuc*uu0c;
  FTYPE vartemp7956=172715.*udnur*uu0l;
  FTYPE vartemp7958=172715.*udnul*uu0r;
  FTYPE vartemp7941=285266.*udnur*uu0l;
  FTYPE vartemp7947=285266.*udnul*uu0r;
  FTYPE vartemp7981=-202797.*uu0r;
  FTYPE vartemp7982=vartemp7964+vartemp7970+vartemp7981;
  FTYPE vartemp7969=285266.*uu0c;
  FTYPE vartemp8028=405594.*uu0r;
  FTYPE vartemp8029=vartemp7974+vartemp8027+vartemp8028;
  FTYPE vartemp8026=27.*udnuc*vartemp7904;
  FTYPE vartemp8030=udnul*vartemp8029;
  FTYPE vartemp8031=-32553.*uu0c;
  FTYPE vartemp8032=15022.*vartemp7873;
  FTYPE vartemp8033=vartemp8031+vartemp8032;
  FTYPE vartemp8034=27.*udnur*vartemp8033;
  FTYPE vartemp8035=vartemp8026+vartemp8030+vartemp8034;
  FTYPE vartemp7983=udnul*vartemp7982;
  FTYPE vartemp7984=udnur*vartemp7967;
  FTYPE vartemp7988=udnuc*vartemp7987;
  FTYPE vartemp7989=vartemp7983+vartemp7984+vartemp7988;
  FTYPE vartemp7999=397619.*uu0c;
  FTYPE vartemp8131=397619.*uu0l;
  FTYPE vartemp8132=743445.*uu0r;
  FTYPE vartemp8133=vartemp8003+vartemp8131+vartemp8132;
  FTYPE vartemp8013=3.206285e6*uu0c;
  FTYPE vartemp8140=32553.*uu0r;
  FTYPE vartemp8036=ud0l*vartemp8035;
  FTYPE vartemp8037=ud0r*vartemp7994;
  FTYPE vartemp8044=ud0c*vartemp8043;
  FTYPE vartemp8045=vartemp8036+vartemp8037+vartemp8044;
  FTYPE vartemp8107=-397619.*uu0l;
  FTYPE vartemp8108=-743445.*uu0r;
  FTYPE vartemp8109=vartemp7942+vartemp8107+vartemp8108;
  FTYPE vartemp8110=udnuc*vartemp8109;
  FTYPE vartemp8047=-397619.*uu0c;
  FTYPE vartemp8144=22517.*uu0l;
  FTYPE vartemp8145=vartemp8140+vartemp8144;
  FTYPE vartemp8074=1.757862e6*uu0r;
  FTYPE vartemp8075=vartemp8039+vartemp8070+vartemp8074;
  FTYPE vartemp8163=udnul*vartemp8109;
  FTYPE vartemp8164=27.*udnur*vartemp7904;
  FTYPE vartemp8165=54.*vartemp8145;
  FTYPE vartemp8166=vartemp8057+vartemp8165;
  FTYPE vartemp8167=udnuc*vartemp8166;
  FTYPE vartemp8168=vartemp8163+vartemp8164+vartemp8167;
  FTYPE vartemp8046=gdetl*vartemp8045;
  FTYPE vartemp8048=94946.*uu0r;
  FTYPE vartemp8049=vartemp8027+vartemp8047+vartemp8048;
  FTYPE vartemp8050=udnur*vartemp8049;
  FTYPE vartemp8051=udnul*vartemp7977;
  FTYPE vartemp8052=vartemp7946+vartemp8050+vartemp8051;
  FTYPE vartemp8053=ud0r*vartemp8052;
  FTYPE vartemp8054=ud0l*vartemp7994;
  FTYPE vartemp8062=ud0c*vartemp8061;
  FTYPE vartemp8063=vartemp8053+vartemp8054+vartemp8062;
  FTYPE vartemp8064=gdetr*vartemp8063;
  FTYPE vartemp8065=7.655331e6*uu0c;
  FTYPE vartemp8066=-3.825163e6*uu0l;
  FTYPE vartemp8067=-3.206285e6*uu0r;
  FTYPE vartemp8068=vartemp8065+vartemp8066+vartemp8067;
  FTYPE vartemp8069=udnuc*vartemp8068;
  FTYPE vartemp8071=1.215918e6*uu0r;
  FTYPE vartemp8072=vartemp8057+vartemp8070+vartemp8071;
  FTYPE vartemp8073=udnur*vartemp8072;
  FTYPE vartemp8076=udnul*vartemp8075;
  FTYPE vartemp8077=vartemp8069+vartemp8073+vartemp8076;
  FTYPE vartemp8078=ud0c*vartemp8077;
  FTYPE vartemp8079=ud0l*vartemp8043;
  FTYPE vartemp8080=ud0r*vartemp8061;
  FTYPE vartemp8081=vartemp8078+vartemp8079+vartemp8080;
  FTYPE vartemp8082=gdetc*vartemp8081;
  FTYPE vartemp8083=vartemp8046+vartemp8064+vartemp8082;
  FTYPE vartemp8155=gdetr*vartemp8045;
  FTYPE vartemp8156=94946.*uu0l;
  FTYPE vartemp8157=vartemp7976+vartemp8047+vartemp8156;
  FTYPE vartemp8158=udnul*vartemp8157;
  FTYPE vartemp8159=udnur*vartemp8029;
  FTYPE vartemp8160=vartemp8110+vartemp8158+vartemp8159;
  FTYPE vartemp8161=ud0l*vartemp8160;
  FTYPE vartemp8162=ud0r*vartemp8035;
  FTYPE vartemp8169=ud0c*vartemp8168;
  FTYPE vartemp8170=vartemp8161+vartemp8162+vartemp8169;
  FTYPE vartemp8171=gdetl*vartemp8170;
  FTYPE vartemp8172=7.645321e6*uu0c;
  FTYPE vartemp8173=-3.206285e6*uu0l;
  FTYPE vartemp8174=-3.825163e6*uu0r;
  FTYPE vartemp8175=vartemp8172+vartemp8173+vartemp8174;
  FTYPE vartemp8176=udnuc*vartemp8175;
  FTYPE vartemp8177=1.215918e6*uu0l;
  FTYPE vartemp8178=vartemp8057+vartemp8074+vartemp8177;
  FTYPE vartemp8179=udnul*vartemp8178;
  FTYPE vartemp8180=udnur*vartemp8075;
  FTYPE vartemp8181=vartemp8176+vartemp8179+vartemp8180;
  FTYPE vartemp8182=ud0c*vartemp8181;
  FTYPE vartemp8183=ud0r*vartemp8043;
  FTYPE vartemp8184=ud0l*vartemp8168;
  FTYPE vartemp8185=vartemp8182+vartemp8183+vartemp8184;
  FTYPE vartemp8186=gdetc*vartemp8185;
  FTYPE vartemp8187=vartemp8155+vartemp8171+vartemp8186;
  FTYPE vartemp8215=1.757862e6*udnur*uu0l;
  FTYPE vartemp8217=1.757862e6*udnul*uu0r;
  t[1]=-3.*Bcovl*vartemp7925;
  t[2]=855798.*udnul*uu0c;
  t[3]=284838.*udnur*uu0c;
  t[4]=855798.*udnuc*uu0l;
  t[5]=-518145.*udnul*uu0l;
  t[6]=284838.*udnuc*uu0r;
  t[7]=-29839.*udnur*uu0r;
  t[8]=vartemp7927+vartemp7932+vartemp7934;
  t[9]=t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=-743445.*udnul*uu0c;
  t[12]=-397619.*udnur*uu0c;
  t[13]=405594.*udnul*uu0l;
  t[14]=94946.*udnur*uu0r+vartemp7941+vartemp7946+vartemp7947;
  t[15]=t[11]+t[12]+t[13]+t[14];
  t[16]=-405594.*udnul*uu0c;
  t[17]=-285266.*udnur*uu0c;
  t[18]=-405594.*udnuc*uu0l;
  t[19]=202797.*udnul*uu0l;
  t[20]=-285266.*udnuc*uu0r;
  t[21]=95053.*udnur*uu0r;
  t[22]=vartemp7951+vartemp7956+vartemp7958;
  t[23]=t[16]+t[17]+t[18];
  t[24]=t[19]+t[20]+t[21]+t[22];
  t[25]=ud0l*(t[23]+t[24]);
  t[26]=3.*ud0c*t[15];
  t[27]=-3.*t[25];
  t[28]=ud0r*(t[9]+t[10]);
  t[29]=t[26]+t[27];
  t[30]=udnul*vartemp7967;
  t[31]=udnur*(-95053.*uu0r+vartemp7969+vartemp7970);
  t[32]=udnuc*vartemp7977;
  t[33]=ud0l*vartemp7989+ud0c*vartemp7994;
  t[34]=ud0r*(t[30]+t[31]+t[32]);
  t[35]=t[33];
  t[36]=gdetl*(t[34]+t[35]);
  t[37]=ud0l*vartemp7886;
  t[38]=udnul*vartemp7854;
  t[39]=udnur*(-94946.*uu0r+vartemp7890+vartemp7999);
  t[40]=udnuc*vartemp8006;
  t[41]=-27.*udnul*vartemp7883;
  t[42]=udnur*vartemp8006;
  t[43]=udnuc*(vartemp8013-54.*vartemp8016);
  t[44]=ud0r*(t[38]+t[39]+t[40]);
  t[45]=ud0c*(t[41]+t[42]+t[43]);
  t[46]=gdetc*(t[37]+t[44]+t[45]);
  t[47]=3.*t[36];
  t[48]=-3.*t[46];
  t[49]=gdetr*(t[28]+t[29]);
  t[50]=t[47]+t[48];
  t[51]=3.*Bcovc*vartemp8083;
  t[52]=Bcovr*(t[49]+t[50]);
  t[53]=t[51];
  t[54]=-3.*Bcovr*vartemp7925;
  t[55]=ud0r*vartemp7989;
  t[56]=udnul*(-95053.*uu0l+vartemp7966+vartemp7969);
  t[57]=udnur*vartemp7982+udnuc*vartemp8029;
  t[58]=ud0c*vartemp8035;
  t[59]=ud0l*(t[56]+t[57]);
  t[60]=t[58];
  t[61]=gdetr*(t[55]+t[59]+t[60]);
  t[62]=284838.*udnul*uu0c;
  t[63]=855798.*udnur*uu0c;
  t[64]=284838.*udnuc*uu0l;
  t[65]=191.*udnul*uu0l;
  t[66]=855798.*udnuc*uu0r;
  t[67]=-518145.*udnur*uu0r;
  t[68]=vartemp7927+vartemp7932+vartemp7934;
  t[69]=t[62]+t[63]+t[64];
  t[70]=t[65]+t[66]+t[67]+t[68];
  t[71]=-285266.*udnul*uu0c;
  t[72]=-405594.*udnur*uu0c;
  t[73]=-285266.*udnuc*uu0l;
  t[74]=95053.*udnul*uu0l;
  t[75]=-405594.*udnuc*uu0r;
  t[76]=202797.*udnur*uu0r;
  t[77]=vartemp7951+vartemp7956+vartemp7958;
  t[78]=t[71]+t[72]+t[73];
  t[79]=t[74]+t[75]+t[76]+t[77];
  t[80]=ud0r*(t[78]+t[79]);
  t[81]=-397619.*udnul*uu0c;
  t[82]=-743445.*udnur*uu0c;
  t[83]=94946.*udnul*uu0l;
  t[84]=405594.*udnur*uu0r+vartemp7941+vartemp7947+vartemp8110;
  t[85]=t[81]+t[82]+t[83]+t[84];
  t[86]=-3.*t[80];
  t[87]=3.*ud0c*t[85];
  t[88]=ud0l*(t[69]+t[70]);
  t[89]=t[86]+t[87];
  t[90]=ud0r*vartemp7907;
  t[91]=udnur*vartemp7892;
  t[92]=udnul*(-94946.*uu0l+vartemp7853+vartemp7999);
  t[93]=udnuc*vartemp8133;
  t[94]=udnul*vartemp8133;
  t[95]=-65106.*uu0c;
  t[96]=27535.*uu0l+vartemp8140;
  t[97]=udnur*(t[95]+t[96]);
  t[98]=udnuc*(vartemp8013-54.*vartemp8145);
  t[99]=27.*t[97];
  t[100]=t[98];
  t[101]=ud0l*(t[91]+t[92]+t[93]);
  t[102]=ud0c*(t[94]+t[99]+t[100]);
  t[103]=gdetc*(t[90]+t[101]+t[102]);
  t[104]=gdetl*(t[88]+t[89]);
  t[105]=-3.*t[103];
  t[106]=3.*t[61];
  t[107]=t[104]+t[105];
  t[108]=3.*Bcovc*vartemp8187;
  t[109]=Bcovl*(t[106]+t[107]);
  t[110]=t[108];
  t[111]=Bcovr*vartemp8083+Bcovl*vartemp8187;
  t[112]=gdetr*vartemp8081+gdetl*vartemp8185;
  t[113]=-5.0200666e7*uu0c;
  t[114]=2.2935963e7*uu0l;
  t[115]=2.2965993e7*uu0r;
  t[116]=7.645321e6*udnul*uu0c;
  t[117]=7.655331e6*udnur*uu0c;
  t[118]=-3.206285e6*udnul*uu0l;
  t[119]=-3.825163e6*udnur*uu0l;
  t[120]=-3.825163e6*udnul*uu0r;
  t[121]=-3.206285e6*udnur*uu0r;
  t[122]=t[116]+t[117]+t[118];
  t[123]=t[119]+t[120]+t[121];
  t[124]=t[122]+t[123];
  t[125]=udnuc*(t[113]+t[114]+t[115]);
  t[126]=3.*t[124];
  t[127]=7.655331e6*udnuc*uu0c;
  t[128]=-3.825163e6*udnul*uu0c;
  t[129]=-3.206285e6*udnur*uu0c;
  t[130]=-3.825163e6*udnuc*uu0l;
  t[131]=1.757862e6*udnul*uu0l;
  t[132]=-3.206285e6*udnuc*uu0r;
  t[133]=1.215918e6*udnur*uu0r+vartemp8215+vartemp8217;
  t[134]=t[127]+t[128]+t[129];
  t[135]=t[130]+t[131]+t[132]+t[133];
  t[136]=7.645321e6*udnuc*uu0c;
  t[137]=-3.206285e6*udnul*uu0c;
  t[138]=-3.825163e6*udnur*uu0c;
  t[139]=-3.206285e6*udnuc*uu0l;
  t[140]=1.215918e6*udnul*uu0l;
  t[141]=-3.825163e6*udnuc*uu0r;
  t[142]=1.757862e6*udnur*uu0r+vartemp8215+vartemp8217;
  t[143]=t[136]+t[137]+t[138];
  t[144]=t[139]+t[140]+t[141]+t[142];
  t[145]=ud0r*(t[134]+t[135]);
  t[146]=ud0l*(t[143]+t[144]);
  t[147]=t[145]+t[146];
  t[148]=ud0c*(t[125]+t[126]);
  t[149]=3.*t[147];
  t[150]=3.*t[112];
  t[151]=gdetc*(t[148]+t[149]);
  t[152]=3.*t[111];
  t[153]=Bcovc*(t[150]+t[151]);
  t[154]=Bconl*(t[54]+t[109]+t[110]);
  t[155]=Bconc*(t[152]+t[153]);
  t[156]=Bconr*(t[1]+t[52]+t[53]);
  t[157]=t[154]+t[155];
  t[158]=t[156]+t[157];
  *Fright=0.000033300033300033300033300033300033*t[158];




  return(0);



}


// 7-terms multiplying deconvolution
// See: a2confluxes_7terms.nb
static int seventermdeconvolution(
                                  FTYPE gdetl, FTYPE gdetc, FTYPE gdetr
                                  ,FTYPE Bconl, FTYPE Bconc, FTYPE Bconr
                                  ,FTYPE Bcovl, FTYPE Bcovc, FTYPE Bcovr
                                  ,FTYPE uu0l, FTYPE uu0c, FTYPE uu0r
                                  ,FTYPE ud0l, FTYPE ud0c, FTYPE ud0r
                                  ,FTYPE udnul, FTYPE udnuc, FTYPE udnur
                                  ,FTYPE uunul, FTYPE uunuc, FTYPE uunur
                                  ,FTYPE *Fleft, FTYPE *Fright
                                  )
{
  FTYPE t[191]; // 628 total, 191 t[]


  

  FTYPE vartemp26178=5.810922e6*uunuc;
  FTYPE vartemp26179=-3.311055e6*uunul;
  FTYPE vartemp26180=-2.214601e6*uunur;
  FTYPE vartemp26181=vartemp26178+vartemp26179+vartemp26180;
  FTYPE vartemp26183=1.554185e6*uunuc;
  FTYPE vartemp26184=-825922.*uunul;
  FTYPE vartemp26185=-645658.*uunur;
  FTYPE vartemp26186=vartemp26183+vartemp26184+vartemp26185;
  FTYPE vartemp26201=-9.256023e6*uunuc;
  FTYPE vartemp26202=5.810922e6*uunul;
  FTYPE vartemp26203=3.047482e6*uunur;
  FTYPE vartemp26204=vartemp26201+vartemp26202+vartemp26203;
  FTYPE vartemp26208=2.717734e6*uunuc;
  FTYPE vartemp26209=-1.554185e6*uunul;
  FTYPE vartemp26210=-1.028447e6*uunur;
  FTYPE vartemp26211=vartemp26208+vartemp26209+vartemp26210;
  FTYPE vartemp26212=9.*uu0c*vartemp26211;
  FTYPE vartemp26182=uu0r*vartemp26181;
  FTYPE vartemp26187=-9.*uu0c*vartemp26186;
  FTYPE vartemp26188=825922.*uunuc;
  FTYPE vartemp26189=-412961.*uunul;
  FTYPE vartemp26190=-367895.*uunur;
  FTYPE vartemp26191=vartemp26188+vartemp26189+vartemp26190;
  FTYPE vartemp26192=9.*uu0l*vartemp26191;
  FTYPE vartemp26193=vartemp26182+vartemp26187+vartemp26192;
  FTYPE vartemp26223=3.499006e6*uunuc;
  FTYPE vartemp26224=-1.749503e6*uunul;
  FTYPE vartemp26225=-1.554185e6*uunur;
  FTYPE vartemp26226=vartemp26223+vartemp26224+vartemp26225;
  FTYPE vartemp26227=uu0c*vartemp26226;
  FTYPE vartemp26259=825922.*uunul;
  FTYPE vartemp26269=-825922.*uunuc;
  FTYPE vartemp26257=udnur*vartemp26193;
  FTYPE vartemp26258=-1.554185e6*uunuc;
  FTYPE vartemp26260=645658.*uunur;
  FTYPE vartemp26261=vartemp26258+vartemp26259+vartemp26260;
  FTYPE vartemp26262=uu0r*vartemp26261;
  FTYPE vartemp26263=-1.749503e6*uunuc;
  FTYPE vartemp26264=825922.*uunur;
  FTYPE vartemp26265=vartemp26259+vartemp26263+vartemp26264;
  FTYPE vartemp26266=uu0l*vartemp26265;
  FTYPE vartemp26267=vartemp26227+vartemp26262+vartemp26266;
  FTYPE vartemp26268=9.*udnuc*vartemp26267;
  FTYPE vartemp26270=412961.*uunul;
  FTYPE vartemp26271=367895.*uunur;
  FTYPE vartemp26272=vartemp26269+vartemp26270+vartemp26271;
  FTYPE vartemp26273=uu0r*vartemp26272;
  FTYPE vartemp26274=367895.*uunul;
  FTYPE vartemp26275=412961.*uunur;
  FTYPE vartemp26276=vartemp26269+vartemp26274+vartemp26275;
  FTYPE vartemp26277=uu0l*vartemp26276;
  FTYPE vartemp26278=1.749503e6*uunuc;
  FTYPE vartemp26279=uunul+uunur;
  FTYPE vartemp26280=-825922.*vartemp26279;
  FTYPE vartemp26281=vartemp26278+vartemp26280;
  FTYPE vartemp26282=uu0c*vartemp26281;
  FTYPE vartemp26283=vartemp26273+vartemp26277+vartemp26282;
  FTYPE vartemp26284=-9.*udnul*vartemp26283;
  FTYPE vartemp26285=vartemp26257+vartemp26268+vartemp26284;
  FTYPE vartemp26197=-2.214601e6*uunul;
  FTYPE vartemp26219=-1.749503e6*uu0l*uunuc;
  FTYPE vartemp26220=-1.554185e6*uu0r*uunuc;
  FTYPE vartemp26221=825922.*uu0l*uunul;
  FTYPE vartemp26222=825922.*uu0r*uunul;
  FTYPE vartemp26228=825922.*uu0l*uunur;
  FTYPE vartemp26229=645658.*uu0r*uunur;
  FTYPE vartemp26315=3.499006e6*uu0c*uunuc;
  FTYPE vartemp26303=-1.749503e6*uunur;
  FTYPE vartemp26304=vartemp26209+vartemp26223+vartemp26303;
  FTYPE vartemp26230=vartemp26219+vartemp26220+vartemp26221+vartemp26222+vartemp26227+vartemp26228+vartemp26229;
  FTYPE vartemp26231=9.*udnul*vartemp26230;
  FTYPE vartemp26232=-1.3987665e7*uu0l*uunuc;
  FTYPE vartemp26233=-9.256023e6*uu0r*uunuc;
  FTYPE vartemp26234=7.433298e6*uu0l*uunul;
  FTYPE vartemp26235=5.810922e6*uu0r*uunul;
  FTYPE vartemp26236=5.810922e6*uu0l*uunur;
  FTYPE vartemp26237=3.047482e6*uu0r*uunur;
  FTYPE vartemp26238=vartemp26212+vartemp26232+vartemp26233+vartemp26234+vartemp26235+vartemp26236+vartemp26237;
  FTYPE vartemp26239=udnur*vartemp26238;
  FTYPE vartemp26240=-5.9156945e7*uunuc;
  FTYPE vartemp26241=3.1491054e7*uunul;
  FTYPE vartemp26242=2.4459606e7*uunur;
  FTYPE vartemp26243=vartemp26240+vartemp26241+vartemp26242;
  FTYPE vartemp26244=uu0c*vartemp26243;
  FTYPE vartemp26245=3.499006e6*uu0l*uunuc;
  FTYPE vartemp26246=2.717734e6*uu0r*uunuc;
  FTYPE vartemp26247=-1.749503e6*uu0l*uunul;
  FTYPE vartemp26248=-1.554185e6*uu0r*uunul;
  FTYPE vartemp26249=-1.554185e6*uu0l*uunur;
  FTYPE vartemp26250=-1.028447e6*uu0r*uunur;
  FTYPE vartemp26251=vartemp26245+vartemp26246+vartemp26247+vartemp26248+vartemp26249+vartemp26250;
  FTYPE vartemp26252=9.*vartemp26251;
  FTYPE vartemp26253=vartemp26244+vartemp26252;
  FTYPE vartemp26254=udnuc*vartemp26253;
  FTYPE vartemp26255=vartemp26231+vartemp26239+vartemp26254;
  FTYPE vartemp26316=-1.749503e6*uu0c*uunul;
  FTYPE vartemp26317=-1.554185e6*uu0c*uunur;
  FTYPE vartemp26318=vartemp26219+vartemp26220+vartemp26221+vartemp26222+vartemp26228+vartemp26229+vartemp26315+vartemp26316+vartemp26317;
  FTYPE vartemp26319=udnur*vartemp26318;
  FTYPE vartemp26320=-1.554185e6*uu0l*uunuc;
  FTYPE vartemp26321=-1.749503e6*uu0r*uunuc;
  FTYPE vartemp26322=-1.554185e6*uu0c*uunul;
  FTYPE vartemp26323=645658.*uu0l*uunul;
  FTYPE vartemp26324=-1.749503e6*uu0c*uunur;
  FTYPE vartemp26325=825922.*uu0r*uunur;
  FTYPE vartemp26326=vartemp26222+vartemp26228+vartemp26315+vartemp26320+vartemp26321+vartemp26322+vartemp26323+vartemp26324+vartemp26325;
  FTYPE vartemp26327=udnul*vartemp26326;
  FTYPE vartemp26328=vartemp26319+vartemp26327;
  FTYPE vartemp26329=9.*vartemp26328;
  FTYPE vartemp26330=uu0l*vartemp26304;
  FTYPE vartemp26331=uu0r*vartemp26226;
  FTYPE vartemp26332=vartemp26330+vartemp26331;
  FTYPE vartemp26333=9.*vartemp26332;
  FTYPE vartemp26334=-6.6807271e7*uunuc;
  FTYPE vartemp26335=3.1491054e7*vartemp26279;
  FTYPE vartemp26336=vartemp26334+vartemp26335;
  FTYPE vartemp26337=uu0c*vartemp26336;
  FTYPE vartemp26338=vartemp26333+vartemp26337;
  FTYPE vartemp26339=udnuc*vartemp26338;
  FTYPE vartemp26340=vartemp26329+vartemp26339;
  FTYPE vartemp26305=uu0c*vartemp26304;
  FTYPE vartemp26294=-825922.*uunur;
  FTYPE vartemp26381=vartemp26184+vartemp26278+vartemp26294;
  FTYPE vartemp26293=-645658.*uunul;
  FTYPE vartemp26295=vartemp26183+vartemp26293+vartemp26294;
  FTYPE vartemp26384=-3.499006e6*uunuc;
  FTYPE vartemp26393=-5.810922e6*uunuc;
  FTYPE vartemp26382=uu0l*vartemp26381;
  FTYPE vartemp26383=uu0r*vartemp26186;
  FTYPE vartemp26385=1.749503e6*uunul;
  FTYPE vartemp26386=1.554185e6*uunur;
  FTYPE vartemp26387=vartemp26384+vartemp26385+vartemp26386;
  FTYPE vartemp26388=uu0c*vartemp26387;
  FTYPE vartemp26389=vartemp26382+vartemp26383+vartemp26388;
  FTYPE vartemp26390=9.*udnuc*vartemp26389;
  FTYPE vartemp26391=9.*uu0c*vartemp26186;
  FTYPE vartemp26392=-9.*uu0l*vartemp26191;
  FTYPE vartemp26394=3.311055e6*uunul;
  FTYPE vartemp26395=2.214601e6*uunur;
  FTYPE vartemp26396=vartemp26393+vartemp26394+vartemp26395;
  FTYPE vartemp26397=uu0r*vartemp26396;
  FTYPE vartemp26398=vartemp26391+vartemp26392+vartemp26397;
  FTYPE vartemp26399=udnur*vartemp26398;
  FTYPE vartemp26400=9.*udnul*vartemp26283;
  FTYPE vartemp26401=vartemp26390+vartemp26399+vartemp26400;
  FTYPE vartemp26425=3.311055e6*uunuc;
  FTYPE vartemp26439=-9.*uu0c*vartemp26191;
  FTYPE vartemp26440=412961.*uunuc;
  FTYPE vartemp26441=-195214.*vartemp26279;
  FTYPE vartemp26442=vartemp26440+vartemp26441;
  FTYPE vartemp26443=9.*uu0l*vartemp26442;
  FTYPE vartemp26444=878463.*uunul;
  FTYPE vartemp26445=690707.*uunur;
  FTYPE vartemp26446=vartemp26444+vartemp26445;
  FTYPE vartemp26447=-2.*vartemp26446;
  FTYPE vartemp26448=vartemp26425+vartemp26447;
  FTYPE vartemp26449=uu0r*vartemp26448;
  FTYPE vartemp26450=vartemp26439+vartemp26443+vartemp26449;
  FTYPE vartemp26297=-367895.*uunul;
  FTYPE vartemp26298=-412961.*uunur;
  FTYPE vartemp26299=vartemp26188+vartemp26297+vartemp26298;
  FTYPE vartemp26403=uu0r*vartemp26381;
  FTYPE vartemp26404=uu0l*vartemp26295;
  FTYPE vartemp26405=1.554185e6*uunul;
  FTYPE vartemp26406=1.749503e6*uunur;
  FTYPE vartemp26407=vartemp26384+vartemp26405+vartemp26406;
  FTYPE vartemp26408=uu0c*vartemp26407;
  FTYPE vartemp26409=vartemp26403+vartemp26404+vartemp26408;
  FTYPE vartemp26410=9.*udnuc*vartemp26409;
  FTYPE vartemp26411=9.*uu0c*vartemp26295;
  FTYPE vartemp26412=9.*uu0r*vartemp26276;
  FTYPE vartemp26413=2.214601e6*uunul;
  FTYPE vartemp26414=3.311055e6*uunur;
  FTYPE vartemp26415=vartemp26393+vartemp26413+vartemp26414;
  FTYPE vartemp26416=uu0l*vartemp26415;
  FTYPE vartemp26417=vartemp26411+vartemp26412+vartemp26416;
  FTYPE vartemp26418=udnul*vartemp26417;
  FTYPE vartemp26419=9.*udnur*vartemp26283;
  FTYPE vartemp26420=vartemp26410+vartemp26418+vartemp26419;
  FTYPE vartemp26454=9.*udnuc*vartemp26283;
  FTYPE vartemp26455=udnur*vartemp26450;
  FTYPE vartemp26456=-9.*uu0c*vartemp26299;
  FTYPE vartemp26457=9.*uu0r*vartemp26442;
  FTYPE vartemp26458=690707.*uunul;
  FTYPE vartemp26459=878463.*uunur;
  FTYPE vartemp26460=vartemp26458+vartemp26459;
  FTYPE vartemp26461=-2.*vartemp26460;
  FTYPE vartemp26462=vartemp26425+vartemp26461;
  FTYPE vartemp26463=uu0l*vartemp26462;
  FTYPE vartemp26464=vartemp26456+vartemp26457+vartemp26463;
  FTYPE vartemp26465=udnul*vartemp26464;
  FTYPE vartemp26466=vartemp26454+vartemp26455+vartemp26465;
  FTYPE vartemp26431=-1.381414e6*uunul;
  FTYPE vartemp26430=2.214601e6*uunuc;
  FTYPE vartemp26427=-1.381414e6*uunur;
  FTYPE vartemp26497=-3.311055e6*uunuc;
  FTYPE vartemp26498=1.756926e6*uunul;
  FTYPE vartemp26499=1.381414e6*uunur;
  FTYPE vartemp26500=vartemp26497+vartemp26498+vartemp26499;
  FTYPE vartemp26504=9.*uu0c*vartemp26191;
  FTYPE vartemp26505=uu0r*vartemp26500;
  FTYPE vartemp26506=-9.*uu0l*vartemp26442;
  FTYPE vartemp26507=vartemp26504+vartemp26505+vartemp26506;
  FTYPE vartemp26493=1.381414e6*uunul;
  FTYPE vartemp26194=udnul*vartemp26193;
  FTYPE vartemp26195=uu0l*vartemp26181;
  FTYPE vartemp26196=3.047482e6*uunuc;
  FTYPE vartemp26198=-737935.*uunur;
  FTYPE vartemp26199=vartemp26196+vartemp26197+vartemp26198;
  FTYPE vartemp26200=uu0r*vartemp26199;
  FTYPE vartemp26205=uu0c*vartemp26204;
  FTYPE vartemp26206=vartemp26195+vartemp26200+vartemp26205;
  FTYPE vartemp26207=udnur*vartemp26206;
  FTYPE vartemp26213=-9.*uu0l*vartemp26186;
  FTYPE vartemp26214=uu0r*vartemp26204;
  FTYPE vartemp26215=vartemp26212+vartemp26213+vartemp26214;
  FTYPE vartemp26216=udnuc*vartemp26215;
  FTYPE vartemp26217=vartemp26194+vartemp26207+vartemp26216;
  FTYPE vartemp26426=-1.756926e6*uunul;
  FTYPE vartemp26428=vartemp26425+vartemp26426+vartemp26427;
  FTYPE vartemp26429=uu0l*vartemp26428;
  FTYPE vartemp26432=-738134.*uunur;
  FTYPE vartemp26433=vartemp26430+vartemp26431+vartemp26432;
  FTYPE vartemp26434=uu0r*vartemp26433;
  FTYPE vartemp26435=uu0c*vartemp26396;
  FTYPE vartemp26436=vartemp26429+vartemp26434+vartemp26435;
  FTYPE vartemp26437=udnur*vartemp26436;
  FTYPE vartemp26438=udnuc*vartemp26398;
  FTYPE vartemp26451=udnul*vartemp26450;
  FTYPE vartemp26452=vartemp26437+vartemp26438+vartemp26451;
  FTYPE vartemp26548=9.256023e6*uunuc;
  FTYPE vartemp26549=-5.810922e6*uunul;
  FTYPE vartemp26550=-3.047482e6*uunur;
  FTYPE vartemp26551=vartemp26548+vartemp26549+vartemp26550;
  FTYPE vartemp26552=uu0c*vartemp26551;
  FTYPE vartemp26582=uu0r*vartemp26551;
  FTYPE vartemp26583=9.*uu0l*vartemp26186;
  FTYPE vartemp26584=-2.4459606e7*uunuc;
  FTYPE vartemp26585=1.3987665e7*uunul;
  FTYPE vartemp26586=9.256023e6*uunur;
  FTYPE vartemp26587=vartemp26584+vartemp26585+vartemp26586;
  FTYPE vartemp26588=uu0c*vartemp26587;
  FTYPE vartemp26589=vartemp26582+vartemp26583+vartemp26588;
  FTYPE vartemp26290=-3.311055e6*uunur;
  FTYPE vartemp26291=vartemp26178+vartemp26197+vartemp26290;
  FTYPE vartemp26292=uu0l*vartemp26291;
  FTYPE vartemp26296=-9.*uu0c*vartemp26295;
  FTYPE vartemp26300=9.*uu0r*vartemp26299;
  FTYPE vartemp26301=vartemp26292+vartemp26296+vartemp26300;
  FTYPE vartemp26306=645658.*uunul;
  FTYPE vartemp26307=vartemp26258+vartemp26264+vartemp26306;
  FTYPE vartemp26623=3.047482e6*uunul;
  FTYPE vartemp26624=5.810922e6*uunur;
  FTYPE vartemp26625=vartemp26201+vartemp26623+vartemp26624;
  FTYPE vartemp26308=uu0l*vartemp26307;
  FTYPE vartemp26309=uu0r*vartemp26265;
  FTYPE vartemp26310=vartemp26305+vartemp26308+vartemp26309;
  FTYPE vartemp26629=-1.028447e6*uunul;
  FTYPE vartemp26630=vartemp26208+vartemp26225+vartemp26629;
  FTYPE vartemp26631=9.*uu0c*vartemp26630;
  FTYPE vartemp26632=9.*uu0r*vartemp26307;
  FTYPE vartemp26633=uu0l*vartemp26625;
  FTYPE vartemp26634=vartemp26631+vartemp26632+vartemp26633;
  FTYPE vartemp26302=udnul*vartemp26301;
  FTYPE vartemp26311=9.*udnuc*vartemp26310;
  FTYPE vartemp26312=-9.*udnur*vartemp26283;
  FTYPE vartemp26313=vartemp26302+vartemp26311+vartemp26312;
  FTYPE vartemp26289=ud0r*vartemp26285;
  FTYPE vartemp26314=ud0l*vartemp26313;
  FTYPE vartemp26341=ud0c*vartemp26340;
  FTYPE vartemp26342=vartemp26289+vartemp26314+vartemp26341;
  FTYPE vartemp26638=9.*udnur*vartemp26310;
  FTYPE vartemp26639=udnul*vartemp26634;
  FTYPE vartemp26640=uu0r*vartemp26304;
  FTYPE vartemp26641=uu0l*vartemp26630;
  FTYPE vartemp26642=vartemp26640+vartemp26641;
  FTYPE vartemp26643=9.*vartemp26642;
  FTYPE vartemp26644=2.4459606e7*uunul;
  FTYPE vartemp26645=3.1491054e7*uunur;
  FTYPE vartemp26646=vartemp26240+vartemp26644+vartemp26645;
  FTYPE vartemp26647=uu0c*vartemp26646;
  FTYPE vartemp26648=vartemp26643+vartemp26647;
  FTYPE vartemp26649=udnuc*vartemp26648;
  FTYPE vartemp26650=vartemp26638+vartemp26639+vartemp26649;
  FTYPE vartemp26349=3.1491054e7*uu0r*uunul;
  FTYPE vartemp26355=3.1491054e7*uu0l*uunur;
  FTYPE vartemp26367=udnur*vartemp26230;
  FTYPE vartemp26368=vartemp26222+vartemp26228+vartemp26305+vartemp26320+vartemp26321+vartemp26323+vartemp26325;
  FTYPE vartemp26369=udnul*vartemp26368;
  FTYPE vartemp26370=vartemp26367+vartemp26369;
  FTYPE vartemp26371=-9.*vartemp26370;
  FTYPE vartemp26372=-9.*vartemp26332;
  FTYPE vartemp26373=6.6807271e7*uunuc;
  FTYPE vartemp26374=-3.1491054e7*vartemp26279;
  FTYPE vartemp26375=vartemp26373+vartemp26374;
  FTYPE vartemp26376=uu0c*vartemp26375;
  FTYPE vartemp26377=vartemp26372+vartemp26376;
  FTYPE vartemp26378=udnuc*vartemp26377;
  FTYPE vartemp26379=vartemp26371+vartemp26378;
  FTYPE vartemp26380=ud0c*vartemp26379;
  FTYPE vartemp26402=ud0r*vartemp26401;
  FTYPE vartemp26421=ud0l*vartemp26420;
  FTYPE vartemp26422=vartemp26380+vartemp26402+vartemp26421;
  FTYPE vartemp26423=gdetc*vartemp26422;
  FTYPE vartemp26424=ud0c*vartemp26401;
  FTYPE vartemp26453=ud0r*vartemp26452;
  FTYPE vartemp26467=ud0l*vartemp26466;
  FTYPE vartemp26468=vartemp26424+vartemp26453+vartemp26467;
  FTYPE vartemp26469=gdetr*vartemp26468;
  FTYPE vartemp26470=ud0c*vartemp26420;
  FTYPE vartemp26471=ud0r*vartemp26466;
  FTYPE vartemp26472=-1.756926e6*uunur;
  FTYPE vartemp26473=vartemp26425+vartemp26431+vartemp26472;
  FTYPE vartemp26474=uu0r*vartemp26473;
  FTYPE vartemp26475=-738134.*uunul;
  FTYPE vartemp26476=vartemp26427+vartemp26430+vartemp26475;
  FTYPE vartemp26477=uu0l*vartemp26476;
  FTYPE vartemp26478=uu0c*vartemp26415;
  FTYPE vartemp26479=vartemp26474+vartemp26477+vartemp26478;
  FTYPE vartemp26480=udnul*vartemp26479;
  FTYPE vartemp26481=udnuc*vartemp26417;
  FTYPE vartemp26482=udnur*vartemp26464;
  FTYPE vartemp26483=vartemp26480+vartemp26481+vartemp26482;
  FTYPE vartemp26484=ud0l*vartemp26483;
  FTYPE vartemp26485=vartemp26470+vartemp26471+vartemp26484;
  FTYPE vartemp26486=gdetl*vartemp26485;
  FTYPE vartemp26487=vartemp26423+vartemp26469+vartemp26486;
  FTYPE vartemp26492=-2.214601e6*uunuc;
  FTYPE vartemp26513=1.756926e6*uunur;
  FTYPE vartemp26514=vartemp26493+vartemp26497+vartemp26513;
  FTYPE vartemp26512=9.*uu0c*vartemp26299;
  FTYPE vartemp26515=uu0l*vartemp26514;
  FTYPE vartemp26516=-412961.*uunuc;
  FTYPE vartemp26517=195214.*vartemp26279;
  FTYPE vartemp26518=vartemp26516+vartemp26517;
  FTYPE vartemp26519=9.*uu0r*vartemp26518;
  FTYPE vartemp26520=vartemp26512+vartemp26515+vartemp26519;
  FTYPE vartemp26511=udnur*vartemp26507;
  FTYPE vartemp26521=udnul*vartemp26520;
  FTYPE vartemp26522=uu0l*vartemp26299;
  FTYPE vartemp26523=uu0r*vartemp26191;
  FTYPE vartemp26524=825922.*vartemp26279;
  FTYPE vartemp26525=vartemp26263+vartemp26524;
  FTYPE vartemp26526=uu0c*vartemp26525;
  FTYPE vartemp26527=vartemp26522+vartemp26523+vartemp26526;
  FTYPE vartemp26528=9.*udnuc*vartemp26527;
  FTYPE vartemp26529=vartemp26511+vartemp26521+vartemp26528;
  FTYPE vartemp26533=-9.142446e6*uu0c*uunuc;
  FTYPE vartemp26538=-2.214402e6*uu0r*uunul;
  FTYPE vartemp26540=-2.214402e6*uu0l*uunur;
  FTYPE vartemp26557=5.810922e6*uu0c*uunuc;
  FTYPE vartemp26562=1.381414e6*uu0r*uunul;
  FTYPE vartemp26564=1.381414e6*uu0l*uunur;
  FTYPE vartemp26547=2.214601e6*uu0r*uunul;
  FTYPE vartemp26553=2.214601e6*uu0l*uunur;
  FTYPE vartemp26618=udnur*vartemp26301;
  FTYPE vartemp26619=uu0r*vartemp26291;
  FTYPE vartemp26620=-737935.*uunul;
  FTYPE vartemp26621=vartemp26180+vartemp26196+vartemp26620;
  FTYPE vartemp26622=uu0l*vartemp26621;
  FTYPE vartemp26626=uu0c*vartemp26625;
  FTYPE vartemp26627=vartemp26619+vartemp26622+vartemp26626;
  FTYPE vartemp26628=udnul*vartemp26627;
  FTYPE vartemp26635=udnuc*vartemp26634;
  FTYPE vartemp26636=vartemp26618+vartemp26628+vartemp26635;
  FTYPE vartemp26719=-3.047482e6*uunul;
  FTYPE vartemp26720=-5.810922e6*uunur;
  FTYPE vartemp26721=vartemp26548+vartemp26719+vartemp26720;
  FTYPE vartemp26722=uu0c*vartemp26721;
  FTYPE vartemp26574=-3.047482e6*uunuc;
  FTYPE vartemp26739=uu0l*vartemp26721;
  FTYPE vartemp26740=9.*uu0r*vartemp26295;
  FTYPE vartemp26741=9.256023e6*uunul;
  FTYPE vartemp26742=1.3987665e7*uunur;
  FTYPE vartemp26743=vartemp26584+vartemp26741+vartemp26742;
  FTYPE vartemp26744=uu0c*vartemp26743;
  FTYPE vartemp26745=vartemp26739+vartemp26740+vartemp26744;
  FTYPE vartemp26596=5.9156945e7*uunuc;
  FTYPE vartemp26601=-2.717734e6*uunuc;
  FTYPE vartemp26218=ud0r*vartemp26217;
  FTYPE vartemp26256=ud0c*vartemp26255;
  FTYPE vartemp26286=ud0l*vartemp26285;
  FTYPE vartemp26287=vartemp26218+vartemp26256+vartemp26286;
  FTYPE vartemp26288=gdetr*vartemp26287;
  FTYPE vartemp26343=gdetl*vartemp26342;
  FTYPE vartemp26344=ud0r*vartemp26255;
  FTYPE vartemp26345=ud0l*vartemp26340;
  FTYPE vartemp26346=-6.6807271e7*uu0l*uunuc;
  FTYPE vartemp26347=-5.9156945e7*uu0r*uunuc;
  FTYPE vartemp26348=3.1491054e7*uu0l*uunul;
  FTYPE vartemp26350=1.33609537e8*uunuc;
  FTYPE vartemp26351=-6.6807271e7*uunul;
  FTYPE vartemp26352=-5.9156945e7*uunur;
  FTYPE vartemp26353=vartemp26350+vartemp26351+vartemp26352;
  FTYPE vartemp26354=uu0c*vartemp26353;
  FTYPE vartemp26356=2.4459606e7*uu0r*uunur;
  FTYPE vartemp26357=vartemp26346+vartemp26347+vartemp26348+vartemp26349+vartemp26354+vartemp26355+vartemp26356;
  FTYPE vartemp26358=udnuc*vartemp26357;
  FTYPE vartemp26359=udnur*vartemp26253;
  FTYPE vartemp26360=udnul*vartemp26338;
  FTYPE vartemp26361=vartemp26358+vartemp26359+vartemp26360;
  FTYPE vartemp26362=ud0c*vartemp26361;
  FTYPE vartemp26363=vartemp26344+vartemp26345+vartemp26362;
  FTYPE vartemp26364=gdetc*vartemp26363;
  FTYPE vartemp26365=vartemp26288+vartemp26343+vartemp26364;
  FTYPE vartemp26637=ud0l*vartemp26636;
  FTYPE vartemp26651=ud0c*vartemp26650;
  FTYPE vartemp26652=ud0r*vartemp26313;
  FTYPE vartemp26653=vartemp26637+vartemp26651+vartemp26652;
  FTYPE vartemp26654=gdetl*vartemp26653;
  FTYPE vartemp26655=gdetr*vartemp26342;
  FTYPE vartemp26656=ud0l*vartemp26650;
  FTYPE vartemp26657=ud0r*vartemp26340;
  FTYPE vartemp26658=-5.9156945e7*uu0l*uunuc;
  FTYPE vartemp26659=-6.6807271e7*uu0r*uunuc;
  FTYPE vartemp26660=2.4459606e7*uu0l*uunul;
  FTYPE vartemp26661=1.33619547e8*uunuc;
  FTYPE vartemp26662=-5.9156945e7*uunul;
  FTYPE vartemp26663=-6.6807271e7*uunur;
  FTYPE vartemp26664=vartemp26661+vartemp26662+vartemp26663;
  FTYPE vartemp26665=uu0c*vartemp26664;
  FTYPE vartemp26666=3.1491054e7*uu0r*uunur;
  FTYPE vartemp26667=vartemp26349+vartemp26355+vartemp26658+vartemp26659+vartemp26660+vartemp26665+vartemp26666;
  FTYPE vartemp26668=udnuc*vartemp26667;
  FTYPE vartemp26669=2.717734e6*uu0l*uunuc;
  FTYPE vartemp26670=3.499006e6*uu0r*uunuc;
  FTYPE vartemp26671=-1.028447e6*uu0l*uunul;
  FTYPE vartemp26672=-1.749503e6*uu0r*uunur;
  FTYPE vartemp26673=vartemp26248+vartemp26249+vartemp26669+vartemp26670+vartemp26671+vartemp26672;
  FTYPE vartemp26674=9.*vartemp26673;
  FTYPE vartemp26675=vartemp26647+vartemp26674;
  FTYPE vartemp26676=udnul*vartemp26675;
  FTYPE vartemp26677=udnur*vartemp26338;
  FTYPE vartemp26678=vartemp26668+vartemp26676+vartemp26677;
  FTYPE vartemp26679=ud0c*vartemp26678;
  FTYPE vartemp26680=vartemp26656+vartemp26657+vartemp26679;
  FTYPE vartemp26681=gdetc*vartemp26680;
  FTYPE vartemp26682=vartemp26654+vartemp26655+vartemp26681;
  t[1]=-3.*Bcovc*vartemp26365;
  t[2]=3.*Bcovl*vartemp26487;
  t[3]=ud0c*vartemp26285;
  t[4]=udnuc*vartemp26193;
  t[5]=uu0c*vartemp26181;
  t[6]=uu0r*(738134.*uunur+vartemp26492+vartemp26493);
  t[7]=uu0l*vartemp26500;
  t[8]=udnul*vartemp26507;
  t[9]=udnur*(t[5]+t[6]+t[7]);
  t[10]=t[8];
  t[11]=ud0l*vartemp26529;
  t[12]=ud0r*(t[4]+t[9]+t[10]);
  t[13]=t[11];
  t[14]=gdetl*(t[3]+t[12]+t[13]);
  t[15]=-3.*ud0c*vartemp26217;
  t[16]=3.*ud0l*vartemp26452;
  t[17]=6.643803e6*uu0l*uunuc;
  t[18]=2.213805e6*uu0r*uunuc;
  t[19]=6.643803e6*uu0c*uunul;
  t[20]=-4.144242e6*uu0l*uunul;
  t[21]=2.213805e6*uu0c*uunur;
  t[22]=406.*uu0r*uunur;
  t[23]=vartemp26533+vartemp26538+vartemp26540;
  t[24]=t[17]+t[18]+t[19];
  t[25]=t[20]+t[21]+t[22]+t[23];
  t[26]=-5.810922e6*uu0l*uunuc;
  t[27]=-3.047482e6*uu0r*uunuc;
  t[28]=3.311055e6*uu0l*uunul;
  t[29]=737935.*uu0r*uunur+vartemp26547+vartemp26552+vartemp26553;
  t[30]=t[26]+t[27]+t[28]+t[29];
  t[31]=-3.311055e6*uu0l*uunuc;
  t[32]=-2.214601e6*uu0r*uunuc;
  t[33]=-3.311055e6*uu0c*uunul;
  t[34]=1.756926e6*uu0l*uunul;
  t[35]=-2.214601e6*uu0c*uunur;
  t[36]=738134.*uu0r*uunur;
  t[37]=vartemp26557+vartemp26562+vartemp26564;
  t[38]=t[31]+t[32]+t[33];
  t[39]=t[34]+t[35]+t[36]+t[37];
  t[40]=udnul*(t[38]+t[39]);
  t[41]=3.*udnuc*t[30];
  t[42]=-3.*t[40];
  t[43]=udnur*(t[24]+t[25]);
  t[44]=t[41]+t[42];
  t[45]=t[16];
  t[46]=ud0r*(t[43]+t[44]);
  t[47]=ud0l*vartemp26401;
  t[48]=udnul*vartemp26398;
  t[49]=uu0l*vartemp26396;
  t[50]=vartemp26552+uu0r*(737935.*uunur+vartemp26413+vartemp26574);
  t[51]=udnuc*vartemp26589;
  t[52]=udnur*(t[49]+t[50]);
  t[53]=t[51];
  t[54]=-9.*udnul*vartemp26267;
  t[55]=udnur*vartemp26589;
  t[56]=-3.1491054e7*uunul;
  t[57]=-2.4459606e7*uunur+vartemp26596;
  t[58]=uu0l*vartemp26387;
  t[59]=uu0r*(1.028447e6*uunur+vartemp26405+vartemp26601);
  t[60]=t[58]+t[59];
  t[61]=uu0c*(t[56]+t[57]);
  t[62]=9.*t[60];
  t[63]=t[55];
  t[64]=udnuc*(t[61]+t[62]);
  t[65]=ud0r*(t[48]+t[52]+t[53]);
  t[66]=ud0c*(t[54]+t[63]+t[64]);
  t[67]=gdetc*(t[47]+t[65]+t[66]);
  t[68]=gdetr*(t[15]+t[45]+t[46]);
  t[69]=3.*t[67];
  t[70]=-3.*t[14];
  t[71]=t[68]+t[69];
  t[72]=t[2];
  t[73]=Bcovr*(t[70]+t[71]);
  t[74]=Bconr*(t[1]+t[72]+t[73]);
  t[75]=8.51887918e8*uunuc;
  t[76]=-4.00858641e8*uunul;
  t[77]=-4.00828611e8*uunur;
  t[78]=-1.33619547e8*uu0l*uunuc;
  t[79]=-1.33609537e8*uu0r*uunuc;
  t[80]=5.9156945e7*uu0l*uunul;
  t[81]=6.6807271e7*uu0r*uunul;
  t[82]=6.6807271e7*uu0l*uunur;
  t[83]=5.9156945e7*uu0r*uunur;
  t[84]=t[78]+t[79]+t[80];
  t[85]=t[81]+t[82]+t[83];
  t[86]=t[84]+t[85];
  t[87]=uu0c*(t[75]+t[76]+t[77]);
  t[88]=3.*t[86];
  t[89]=1.33609537e8*uu0c*uunuc;
  t[90]=-6.6807271e7*uu0c*uunul;
  t[91]=-5.9156945e7*uu0c*uunur+vartemp26346;
  t[92]=vartemp26347+vartemp26348+vartemp26349+vartemp26355+vartemp26356;
  t[93]=t[89]+t[90]+t[91]+t[92];
  t[94]=1.33619547e8*uu0c*uunuc;
  t[95]=-5.9156945e7*uu0c*uunul;
  t[96]=-6.6807271e7*uu0c*uunur+vartemp26349;
  t[97]=vartemp26355+vartemp26658+vartemp26659+vartemp26660+vartemp26666;
  t[98]=t[94]+t[95]+t[96]+t[97];
  t[99]=udnur*t[93]+udnul*t[98];
  t[100]=udnuc*(t[87]+t[88]);
  t[101]=-3.*t[99];
  t[102]=ud0r*vartemp26361+ud0l*vartemp26678;
  t[103]=ud0c*(t[100]+t[101]);
  t[104]=-3.*t[102];
  t[105]=gdetr*vartemp26363+gdetl*vartemp26680;
  t[106]=gdetc*(t[103]+t[104]);
  t[107]=-3.*t[105];
  t[108]=Bcovr*vartemp26365+Bcovl*vartemp26682;
  t[109]=Bcovc*(t[106]+t[107]);
  t[110]=-3.*t[108];
  t[111]=Bconc*(t[109]+t[110]);
  t[112]=3.*Bcovr*vartemp26487;
  t[113]=-3.*Bcovc*vartemp26682;
  t[114]=ud0c*vartemp26313;
  t[115]=udnuc*vartemp26301;
  t[116]=uu0c*vartemp26291;
  t[117]=uu0l*(738134.*uunul+vartemp26492+vartemp26499);
  t[118]=uu0r*vartemp26514;
  t[119]=udnur*vartemp26520;
  t[120]=udnul*(t[116]+t[117]+t[118]);
  t[121]=t[119];
  t[122]=ud0r*vartemp26529;
  t[123]=ud0l*(t[115]+t[120]+t[121]);
  t[124]=t[122];
  t[125]=gdetr*(t[114]+t[123]+t[124]);
  t[126]=3.*ud0r*vartemp26483;
  t[127]=-3.*ud0c*vartemp26636;
  t[128]=2.213805e6*uu0l*uunuc;
  t[129]=6.643803e6*uu0r*uunuc;
  t[130]=2.213805e6*uu0c*uunul;
  t[131]=30436.*uu0l*uunul;
  t[132]=6.643803e6*uu0c*uunur;
  t[133]=-4.144242e6*uu0r*uunur;
  t[134]=vartemp26533+vartemp26538+vartemp26540;
  t[135]=t[128]+t[129]+t[130];
  t[136]=t[131]+t[132]+t[133]+t[134];
  t[137]=-2.214601e6*uu0l*uunuc;
  t[138]=-3.311055e6*uu0r*uunuc;
  t[139]=-2.214601e6*uu0c*uunul;
  t[140]=738134.*uu0l*uunul;
  t[141]=-3.311055e6*uu0c*uunur;
  t[142]=1.756926e6*uu0r*uunur;
  t[143]=vartemp26557+vartemp26562+vartemp26564;
  t[144]=t[137]+t[138]+t[139];
  t[145]=t[140]+t[141]+t[142]+t[143];
  t[146]=udnur*(t[144]+t[145]);
  t[147]=-3.047482e6*uu0l*uunuc;
  t[148]=-5.810922e6*uu0r*uunuc;
  t[149]=737935.*uu0l*uunul;
  t[150]=3.311055e6*uu0r*uunur+vartemp26547+vartemp26553+vartemp26722;
  t[151]=t[147]+t[148]+t[149]+t[150];
  t[152]=-3.*t[146];
  t[153]=3.*udnuc*t[151];
  t[154]=udnul*(t[135]+t[136]);
  t[155]=t[152]+t[153];
  t[156]=t[127];
  t[157]=ud0l*(t[154]+t[155]);
  t[158]=ud0r*vartemp26420;
  t[159]=udnur*vartemp26417;
  t[160]=uu0r*vartemp26415;
  t[161]=uu0l*(737935.*uunul+vartemp26395+vartemp26574)+vartemp26722;
  t[162]=udnuc*vartemp26745;
  t[163]=udnul*(t[160]+t[161]);
  t[164]=t[162];
  t[165]=-9.*udnur*vartemp26310;
  t[166]=-2.4459606e7*uunul;
  t[167]=-3.1491054e7*uunur+vartemp26596;
  t[168]=uu0r*vartemp26407;
  t[169]=uu0l*(1.028447e6*uunul+vartemp26386+vartemp26601);
  t[170]=t[168]+t[169];
  t[171]=uu0c*(t[166]+t[167]);
  t[172]=9.*t[170];
  t[173]=udnul*vartemp26745;
  t[174]=udnuc*(t[171]+t[172]);
  t[175]=t[173];
  t[176]=ud0l*(t[159]+t[163]+t[164]);
  t[177]=ud0c*(t[165]+t[174]+t[175]);
  t[178]=gdetc*(t[158]+t[176]+t[177]);
  t[179]=gdetl*(t[126]+t[156]+t[157]);
  t[180]=3.*t[178];
  t[181]=-3.*t[125];
  t[182]=t[179]+t[180];
  t[183]=t[113];
  t[184]=Bcovl*(t[181]+t[182]);
  t[185]=Bconl*(t[112]+t[183]+t[184]);
  t[186]=-1.*t[111];
  t[187]=-1.*t[185];
  t[188]=-1.*t[74];
  t[189]=t[186]+t[187];
  t[190]=t[188]+t[189];
  *Fleft=0.000033300033300033300033300033300033*t[190];

  FTYPE vartemp27894=5.810922e6*uunuc;
  FTYPE vartemp27895=-3.311055e6*uunul;
  FTYPE vartemp27896=-2.214601e6*uunur;
  FTYPE vartemp27897=vartemp27894+vartemp27895+vartemp27896;
  FTYPE vartemp27899=1.554185e6*uunuc;
  FTYPE vartemp27900=-825922.*uunul;
  FTYPE vartemp27901=-645658.*uunur;
  FTYPE vartemp27902=vartemp27899+vartemp27900+vartemp27901;
  FTYPE vartemp27917=-9.256023e6*uunuc;
  FTYPE vartemp27918=5.810922e6*uunul;
  FTYPE vartemp27919=3.047482e6*uunur;
  FTYPE vartemp27920=vartemp27917+vartemp27918+vartemp27919;
  FTYPE vartemp27924=2.717734e6*uunuc;
  FTYPE vartemp27925=-1.554185e6*uunul;
  FTYPE vartemp27926=-1.028447e6*uunur;
  FTYPE vartemp27927=vartemp27924+vartemp27925+vartemp27926;
  FTYPE vartemp27928=9.*uu0c*vartemp27927;
  FTYPE vartemp27898=uu0r*vartemp27897;
  FTYPE vartemp27903=-9.*uu0c*vartemp27902;
  FTYPE vartemp27904=825922.*uunuc;
  FTYPE vartemp27905=-412961.*uunul;
  FTYPE vartemp27906=-367895.*uunur;
  FTYPE vartemp27907=vartemp27904+vartemp27905+vartemp27906;
  FTYPE vartemp27908=9.*uu0l*vartemp27907;
  FTYPE vartemp27909=vartemp27898+vartemp27903+vartemp27908;
  FTYPE vartemp27939=3.499006e6*uunuc;
  FTYPE vartemp27940=-1.749503e6*uunul;
  FTYPE vartemp27941=-1.554185e6*uunur;
  FTYPE vartemp27942=vartemp27939+vartemp27940+vartemp27941;
  FTYPE vartemp27943=uu0c*vartemp27942;
  FTYPE vartemp27975=825922.*uunul;
  FTYPE vartemp27985=-825922.*uunuc;
  FTYPE vartemp27973=udnur*vartemp27909;
  FTYPE vartemp27974=-1.554185e6*uunuc;
  FTYPE vartemp27976=645658.*uunur;
  FTYPE vartemp27977=vartemp27974+vartemp27975+vartemp27976;
  FTYPE vartemp27978=uu0r*vartemp27977;
  FTYPE vartemp27979=-1.749503e6*uunuc;
  FTYPE vartemp27980=825922.*uunur;
  FTYPE vartemp27981=vartemp27975+vartemp27979+vartemp27980;
  FTYPE vartemp27982=uu0l*vartemp27981;
  FTYPE vartemp27983=vartemp27943+vartemp27978+vartemp27982;
  FTYPE vartemp27984=9.*udnuc*vartemp27983;
  FTYPE vartemp27986=412961.*uunul;
  FTYPE vartemp27987=367895.*uunur;
  FTYPE vartemp27988=vartemp27985+vartemp27986+vartemp27987;
  FTYPE vartemp27989=uu0r*vartemp27988;
  FTYPE vartemp27990=367895.*uunul;
  FTYPE vartemp27991=412961.*uunur;
  FTYPE vartemp27992=vartemp27985+vartemp27990+vartemp27991;
  FTYPE vartemp27993=uu0l*vartemp27992;
  FTYPE vartemp27994=1.749503e6*uunuc;
  FTYPE vartemp27995=uunul+uunur;
  FTYPE vartemp27996=-825922.*vartemp27995;
  FTYPE vartemp27997=vartemp27994+vartemp27996;
  FTYPE vartemp27998=uu0c*vartemp27997;
  FTYPE vartemp27999=vartemp27989+vartemp27993+vartemp27998;
  FTYPE vartemp28000=-9.*udnul*vartemp27999;
  FTYPE vartemp28001=vartemp27973+vartemp27984+vartemp28000;
  FTYPE vartemp27913=-2.214601e6*uunul;
  FTYPE vartemp27935=-1.749503e6*uu0l*uunuc;
  FTYPE vartemp27936=-1.554185e6*uu0r*uunuc;
  FTYPE vartemp27937=825922.*uu0l*uunul;
  FTYPE vartemp27938=825922.*uu0r*uunul;
  FTYPE vartemp27944=825922.*uu0l*uunur;
  FTYPE vartemp27945=645658.*uu0r*uunur;
  FTYPE vartemp28031=3.499006e6*uu0c*uunuc;
  FTYPE vartemp28019=-1.749503e6*uunur;
  FTYPE vartemp28020=vartemp27925+vartemp27939+vartemp28019;
  FTYPE vartemp27946=vartemp27935+vartemp27936+vartemp27937+vartemp27938+vartemp27943+vartemp27944+vartemp27945;
  FTYPE vartemp27947=9.*udnul*vartemp27946;
  FTYPE vartemp27948=-1.3987665e7*uu0l*uunuc;
  FTYPE vartemp27949=-9.256023e6*uu0r*uunuc;
  FTYPE vartemp27950=7.433298e6*uu0l*uunul;
  FTYPE vartemp27951=5.810922e6*uu0r*uunul;
  FTYPE vartemp27952=5.810922e6*uu0l*uunur;
  FTYPE vartemp27953=3.047482e6*uu0r*uunur;
  FTYPE vartemp27954=vartemp27928+vartemp27948+vartemp27949+vartemp27950+vartemp27951+vartemp27952+vartemp27953;
  FTYPE vartemp27955=udnur*vartemp27954;
  FTYPE vartemp27956=-5.9156945e7*uunuc;
  FTYPE vartemp27957=3.1491054e7*uunul;
  FTYPE vartemp27958=2.4459606e7*uunur;
  FTYPE vartemp27959=vartemp27956+vartemp27957+vartemp27958;
  FTYPE vartemp27960=uu0c*vartemp27959;
  FTYPE vartemp27961=3.499006e6*uu0l*uunuc;
  FTYPE vartemp27962=2.717734e6*uu0r*uunuc;
  FTYPE vartemp27963=-1.749503e6*uu0l*uunul;
  FTYPE vartemp27964=-1.554185e6*uu0r*uunul;
  FTYPE vartemp27965=-1.554185e6*uu0l*uunur;
  FTYPE vartemp27966=-1.028447e6*uu0r*uunur;
  FTYPE vartemp27967=vartemp27961+vartemp27962+vartemp27963+vartemp27964+vartemp27965+vartemp27966;
  FTYPE vartemp27968=9.*vartemp27967;
  FTYPE vartemp27969=vartemp27960+vartemp27968;
  FTYPE vartemp27970=udnuc*vartemp27969;
  FTYPE vartemp27971=vartemp27947+vartemp27955+vartemp27970;
  FTYPE vartemp28032=-1.749503e6*uu0c*uunul;
  FTYPE vartemp28033=-1.554185e6*uu0c*uunur;
  FTYPE vartemp28034=vartemp27935+vartemp27936+vartemp27937+vartemp27938+vartemp27944+vartemp27945+vartemp28031+vartemp28032+vartemp28033;
  FTYPE vartemp28035=udnur*vartemp28034;
  FTYPE vartemp28036=-1.554185e6*uu0l*uunuc;
  FTYPE vartemp28037=-1.749503e6*uu0r*uunuc;
  FTYPE vartemp28038=-1.554185e6*uu0c*uunul;
  FTYPE vartemp28039=645658.*uu0l*uunul;
  FTYPE vartemp28040=-1.749503e6*uu0c*uunur;
  FTYPE vartemp28041=825922.*uu0r*uunur;
  FTYPE vartemp28042=vartemp27938+vartemp27944+vartemp28031+vartemp28036+vartemp28037+vartemp28038+vartemp28039+vartemp28040+vartemp28041;
  FTYPE vartemp28043=udnul*vartemp28042;
  FTYPE vartemp28044=vartemp28035+vartemp28043;
  FTYPE vartemp28045=9.*vartemp28044;
  FTYPE vartemp28046=uu0l*vartemp28020;
  FTYPE vartemp28047=uu0r*vartemp27942;
  FTYPE vartemp28048=vartemp28046+vartemp28047;
  FTYPE vartemp28049=9.*vartemp28048;
  FTYPE vartemp28050=-6.6807271e7*uunuc;
  FTYPE vartemp28051=3.1491054e7*vartemp27995;
  FTYPE vartemp28052=vartemp28050+vartemp28051;
  FTYPE vartemp28053=uu0c*vartemp28052;
  FTYPE vartemp28054=vartemp28049+vartemp28053;
  FTYPE vartemp28055=udnuc*vartemp28054;
  FTYPE vartemp28056=vartemp28045+vartemp28055;
  FTYPE vartemp28021=uu0c*vartemp28020;
  FTYPE vartemp28010=-825922.*uunur;
  FTYPE vartemp28097=vartemp27900+vartemp27994+vartemp28010;
  FTYPE vartemp28009=-645658.*uunul;
  FTYPE vartemp28011=vartemp27899+vartemp28009+vartemp28010;
  FTYPE vartemp28100=-3.499006e6*uunuc;
  FTYPE vartemp28109=-5.810922e6*uunuc;
  FTYPE vartemp28098=uu0l*vartemp28097;
  FTYPE vartemp28099=uu0r*vartemp27902;
  FTYPE vartemp28101=1.749503e6*uunul;
  FTYPE vartemp28102=1.554185e6*uunur;
  FTYPE vartemp28103=vartemp28100+vartemp28101+vartemp28102;
  FTYPE vartemp28104=uu0c*vartemp28103;
  FTYPE vartemp28105=vartemp28098+vartemp28099+vartemp28104;
  FTYPE vartemp28106=9.*udnuc*vartemp28105;
  FTYPE vartemp28107=9.*uu0c*vartemp27902;
  FTYPE vartemp28108=-9.*uu0l*vartemp27907;
  FTYPE vartemp28110=3.311055e6*uunul;
  FTYPE vartemp28111=2.214601e6*uunur;
  FTYPE vartemp28112=vartemp28109+vartemp28110+vartemp28111;
  FTYPE vartemp28113=uu0r*vartemp28112;
  FTYPE vartemp28114=vartemp28107+vartemp28108+vartemp28113;
  FTYPE vartemp28115=udnur*vartemp28114;
  FTYPE vartemp28116=9.*udnul*vartemp27999;
  FTYPE vartemp28117=vartemp28106+vartemp28115+vartemp28116;
  FTYPE vartemp28141=3.311055e6*uunuc;
  FTYPE vartemp28155=-9.*uu0c*vartemp27907;
  FTYPE vartemp28156=412961.*uunuc;
  FTYPE vartemp28157=-195214.*vartemp27995;
  FTYPE vartemp28158=vartemp28156+vartemp28157;
  FTYPE vartemp28159=9.*uu0l*vartemp28158;
  FTYPE vartemp28160=878463.*uunul;
  FTYPE vartemp28161=690707.*uunur;
  FTYPE vartemp28162=vartemp28160+vartemp28161;
  FTYPE vartemp28163=-2.*vartemp28162;
  FTYPE vartemp28164=vartemp28141+vartemp28163;
  FTYPE vartemp28165=uu0r*vartemp28164;
  FTYPE vartemp28166=vartemp28155+vartemp28159+vartemp28165;
  FTYPE vartemp28013=-367895.*uunul;
  FTYPE vartemp28014=-412961.*uunur;
  FTYPE vartemp28015=vartemp27904+vartemp28013+vartemp28014;
  FTYPE vartemp28119=uu0r*vartemp28097;
  FTYPE vartemp28120=uu0l*vartemp28011;
  FTYPE vartemp28121=1.554185e6*uunul;
  FTYPE vartemp28122=1.749503e6*uunur;
  FTYPE vartemp28123=vartemp28100+vartemp28121+vartemp28122;
  FTYPE vartemp28124=uu0c*vartemp28123;
  FTYPE vartemp28125=vartemp28119+vartemp28120+vartemp28124;
  FTYPE vartemp28126=9.*udnuc*vartemp28125;
  FTYPE vartemp28127=9.*uu0c*vartemp28011;
  FTYPE vartemp28128=9.*uu0r*vartemp27992;
  FTYPE vartemp28129=2.214601e6*uunul;
  FTYPE vartemp28130=3.311055e6*uunur;
  FTYPE vartemp28131=vartemp28109+vartemp28129+vartemp28130;
  FTYPE vartemp28132=uu0l*vartemp28131;
  FTYPE vartemp28133=vartemp28127+vartemp28128+vartemp28132;
  FTYPE vartemp28134=udnul*vartemp28133;
  FTYPE vartemp28135=9.*udnur*vartemp27999;
  FTYPE vartemp28136=vartemp28126+vartemp28134+vartemp28135;
  FTYPE vartemp28170=9.*udnuc*vartemp27999;
  FTYPE vartemp28171=udnur*vartemp28166;
  FTYPE vartemp28172=-9.*uu0c*vartemp28015;
  FTYPE vartemp28173=9.*uu0r*vartemp28158;
  FTYPE vartemp28174=690707.*uunul;
  FTYPE vartemp28175=878463.*uunur;
  FTYPE vartemp28176=vartemp28174+vartemp28175;
  FTYPE vartemp28177=-2.*vartemp28176;
  FTYPE vartemp28178=vartemp28141+vartemp28177;
  FTYPE vartemp28179=uu0l*vartemp28178;
  FTYPE vartemp28180=vartemp28172+vartemp28173+vartemp28179;
  FTYPE vartemp28181=udnul*vartemp28180;
  FTYPE vartemp28182=vartemp28170+vartemp28171+vartemp28181;
  FTYPE vartemp28147=-1.381414e6*uunul;
  FTYPE vartemp28146=2.214601e6*uunuc;
  FTYPE vartemp28143=-1.381414e6*uunur;
  FTYPE vartemp28213=-3.311055e6*uunuc;
  FTYPE vartemp28214=1.756926e6*uunul;
  FTYPE vartemp28215=1.381414e6*uunur;
  FTYPE vartemp28216=vartemp28213+vartemp28214+vartemp28215;
  FTYPE vartemp28220=9.*uu0c*vartemp27907;
  FTYPE vartemp28221=uu0r*vartemp28216;
  FTYPE vartemp28222=-9.*uu0l*vartemp28158;
  FTYPE vartemp28223=vartemp28220+vartemp28221+vartemp28222;
  FTYPE vartemp28209=1.381414e6*uunul;
  FTYPE vartemp27910=udnul*vartemp27909;
  FTYPE vartemp27911=uu0l*vartemp27897;
  FTYPE vartemp27912=3.047482e6*uunuc;
  FTYPE vartemp27914=-737935.*uunur;
  FTYPE vartemp27915=vartemp27912+vartemp27913+vartemp27914;
  FTYPE vartemp27916=uu0r*vartemp27915;
  FTYPE vartemp27921=uu0c*vartemp27920;
  FTYPE vartemp27922=vartemp27911+vartemp27916+vartemp27921;
  FTYPE vartemp27923=udnur*vartemp27922;
  FTYPE vartemp27929=-9.*uu0l*vartemp27902;
  FTYPE vartemp27930=uu0r*vartemp27920;
  FTYPE vartemp27931=vartemp27928+vartemp27929+vartemp27930;
  FTYPE vartemp27932=udnuc*vartemp27931;
  FTYPE vartemp27933=vartemp27910+vartemp27923+vartemp27932;
  FTYPE vartemp28142=-1.756926e6*uunul;
  FTYPE vartemp28144=vartemp28141+vartemp28142+vartemp28143;
  FTYPE vartemp28145=uu0l*vartemp28144;
  FTYPE vartemp28148=-738134.*uunur;
  FTYPE vartemp28149=vartemp28146+vartemp28147+vartemp28148;
  FTYPE vartemp28150=uu0r*vartemp28149;
  FTYPE vartemp28151=uu0c*vartemp28112;
  FTYPE vartemp28152=vartemp28145+vartemp28150+vartemp28151;
  FTYPE vartemp28153=udnur*vartemp28152;
  FTYPE vartemp28154=udnuc*vartemp28114;
  FTYPE vartemp28167=udnul*vartemp28166;
  FTYPE vartemp28168=vartemp28153+vartemp28154+vartemp28167;
  FTYPE vartemp28264=9.256023e6*uunuc;
  FTYPE vartemp28265=-5.810922e6*uunul;
  FTYPE vartemp28266=-3.047482e6*uunur;
  FTYPE vartemp28267=vartemp28264+vartemp28265+vartemp28266;
  FTYPE vartemp28268=uu0c*vartemp28267;
  FTYPE vartemp28298=uu0r*vartemp28267;
  FTYPE vartemp28299=9.*uu0l*vartemp27902;
  FTYPE vartemp28300=-2.4459606e7*uunuc;
  FTYPE vartemp28301=1.3987665e7*uunul;
  FTYPE vartemp28302=9.256023e6*uunur;
  FTYPE vartemp28303=vartemp28300+vartemp28301+vartemp28302;
  FTYPE vartemp28304=uu0c*vartemp28303;
  FTYPE vartemp28305=vartemp28298+vartemp28299+vartemp28304;
  FTYPE vartemp28006=-3.311055e6*uunur;
  FTYPE vartemp28007=vartemp27894+vartemp27913+vartemp28006;
  FTYPE vartemp28008=uu0l*vartemp28007;
  FTYPE vartemp28012=-9.*uu0c*vartemp28011;
  FTYPE vartemp28016=9.*uu0r*vartemp28015;
  FTYPE vartemp28017=vartemp28008+vartemp28012+vartemp28016;
  FTYPE vartemp28022=645658.*uunul;
  FTYPE vartemp28023=vartemp27974+vartemp27980+vartemp28022;
  FTYPE vartemp28339=3.047482e6*uunul;
  FTYPE vartemp28340=5.810922e6*uunur;
  FTYPE vartemp28341=vartemp27917+vartemp28339+vartemp28340;
  FTYPE vartemp28024=uu0l*vartemp28023;
  FTYPE vartemp28025=uu0r*vartemp27981;
  FTYPE vartemp28026=vartemp28021+vartemp28024+vartemp28025;
  FTYPE vartemp28345=-1.028447e6*uunul;
  FTYPE vartemp28346=vartemp27924+vartemp27941+vartemp28345;
  FTYPE vartemp28347=9.*uu0c*vartemp28346;
  FTYPE vartemp28348=9.*uu0r*vartemp28023;
  FTYPE vartemp28349=uu0l*vartemp28341;
  FTYPE vartemp28350=vartemp28347+vartemp28348+vartemp28349;
  FTYPE vartemp28018=udnul*vartemp28017;
  FTYPE vartemp28027=9.*udnuc*vartemp28026;
  FTYPE vartemp28028=-9.*udnur*vartemp27999;
  FTYPE vartemp28029=vartemp28018+vartemp28027+vartemp28028;
  FTYPE vartemp28005=ud0r*vartemp28001;
  FTYPE vartemp28030=ud0l*vartemp28029;
  FTYPE vartemp28057=ud0c*vartemp28056;
  FTYPE vartemp28058=vartemp28005+vartemp28030+vartemp28057;
  FTYPE vartemp28354=9.*udnur*vartemp28026;
  FTYPE vartemp28355=udnul*vartemp28350;
  FTYPE vartemp28356=uu0r*vartemp28020;
  FTYPE vartemp28357=uu0l*vartemp28346;
  FTYPE vartemp28358=vartemp28356+vartemp28357;
  FTYPE vartemp28359=9.*vartemp28358;
  FTYPE vartemp28360=2.4459606e7*uunul;
  FTYPE vartemp28361=3.1491054e7*uunur;
  FTYPE vartemp28362=vartemp27956+vartemp28360+vartemp28361;
  FTYPE vartemp28363=uu0c*vartemp28362;
  FTYPE vartemp28364=vartemp28359+vartemp28363;
  FTYPE vartemp28365=udnuc*vartemp28364;
  FTYPE vartemp28366=vartemp28354+vartemp28355+vartemp28365;
  FTYPE vartemp28065=3.1491054e7*uu0r*uunul;
  FTYPE vartemp28071=3.1491054e7*uu0l*uunur;
  FTYPE vartemp28083=udnur*vartemp27946;
  FTYPE vartemp28084=vartemp27938+vartemp27944+vartemp28021+vartemp28036+vartemp28037+vartemp28039+vartemp28041;
  FTYPE vartemp28085=udnul*vartemp28084;
  FTYPE vartemp28086=vartemp28083+vartemp28085;
  FTYPE vartemp28087=-9.*vartemp28086;
  FTYPE vartemp28088=-9.*vartemp28048;
  FTYPE vartemp28089=6.6807271e7*uunuc;
  FTYPE vartemp28090=-3.1491054e7*vartemp27995;
  FTYPE vartemp28091=vartemp28089+vartemp28090;
  FTYPE vartemp28092=uu0c*vartemp28091;
  FTYPE vartemp28093=vartemp28088+vartemp28092;
  FTYPE vartemp28094=udnuc*vartemp28093;
  FTYPE vartemp28095=vartemp28087+vartemp28094;
  FTYPE vartemp28096=ud0c*vartemp28095;
  FTYPE vartemp28118=ud0r*vartemp28117;
  FTYPE vartemp28137=ud0l*vartemp28136;
  FTYPE vartemp28138=vartemp28096+vartemp28118+vartemp28137;
  FTYPE vartemp28139=gdetc*vartemp28138;
  FTYPE vartemp28140=ud0c*vartemp28117;
  FTYPE vartemp28169=ud0r*vartemp28168;
  FTYPE vartemp28183=ud0l*vartemp28182;
  FTYPE vartemp28184=vartemp28140+vartemp28169+vartemp28183;
  FTYPE vartemp28185=gdetr*vartemp28184;
  FTYPE vartemp28186=ud0c*vartemp28136;
  FTYPE vartemp28187=ud0r*vartemp28182;
  FTYPE vartemp28188=-1.756926e6*uunur;
  FTYPE vartemp28189=vartemp28141+vartemp28147+vartemp28188;
  FTYPE vartemp28190=uu0r*vartemp28189;
  FTYPE vartemp28191=-738134.*uunul;
  FTYPE vartemp28192=vartemp28143+vartemp28146+vartemp28191;
  FTYPE vartemp28193=uu0l*vartemp28192;
  FTYPE vartemp28194=uu0c*vartemp28131;
  FTYPE vartemp28195=vartemp28190+vartemp28193+vartemp28194;
  FTYPE vartemp28196=udnul*vartemp28195;
  FTYPE vartemp28197=udnuc*vartemp28133;
  FTYPE vartemp28198=udnur*vartemp28180;
  FTYPE vartemp28199=vartemp28196+vartemp28197+vartemp28198;
  FTYPE vartemp28200=ud0l*vartemp28199;
  FTYPE vartemp28201=vartemp28186+vartemp28187+vartemp28200;
  FTYPE vartemp28202=gdetl*vartemp28201;
  FTYPE vartemp28203=vartemp28139+vartemp28185+vartemp28202;
  FTYPE vartemp28208=-2.214601e6*uunuc;
  FTYPE vartemp28229=1.756926e6*uunur;
  FTYPE vartemp28230=vartemp28209+vartemp28213+vartemp28229;
  FTYPE vartemp28228=9.*uu0c*vartemp28015;
  FTYPE vartemp28231=uu0l*vartemp28230;
  FTYPE vartemp28232=-412961.*uunuc;
  FTYPE vartemp28233=195214.*vartemp27995;
  FTYPE vartemp28234=vartemp28232+vartemp28233;
  FTYPE vartemp28235=9.*uu0r*vartemp28234;
  FTYPE vartemp28236=vartemp28228+vartemp28231+vartemp28235;
  FTYPE vartemp28227=udnur*vartemp28223;
  FTYPE vartemp28237=udnul*vartemp28236;
  FTYPE vartemp28238=uu0l*vartemp28015;
  FTYPE vartemp28239=uu0r*vartemp27907;
  FTYPE vartemp28240=825922.*vartemp27995;
  FTYPE vartemp28241=vartemp27979+vartemp28240;
  FTYPE vartemp28242=uu0c*vartemp28241;
  FTYPE vartemp28243=vartemp28238+vartemp28239+vartemp28242;
  FTYPE vartemp28244=9.*udnuc*vartemp28243;
  FTYPE vartemp28245=vartemp28227+vartemp28237+vartemp28244;
  FTYPE vartemp28249=-9.142446e6*uu0c*uunuc;
  FTYPE vartemp28254=-2.214402e6*uu0r*uunul;
  FTYPE vartemp28256=-2.214402e6*uu0l*uunur;
  FTYPE vartemp28273=5.810922e6*uu0c*uunuc;
  FTYPE vartemp28278=1.381414e6*uu0r*uunul;
  FTYPE vartemp28280=1.381414e6*uu0l*uunur;
  FTYPE vartemp28263=2.214601e6*uu0r*uunul;
  FTYPE vartemp28269=2.214601e6*uu0l*uunur;
  FTYPE vartemp28334=udnur*vartemp28017;
  FTYPE vartemp28335=uu0r*vartemp28007;
  FTYPE vartemp28336=-737935.*uunul;
  FTYPE vartemp28337=vartemp27896+vartemp27912+vartemp28336;
  FTYPE vartemp28338=uu0l*vartemp28337;
  FTYPE vartemp28342=uu0c*vartemp28341;
  FTYPE vartemp28343=vartemp28335+vartemp28338+vartemp28342;
  FTYPE vartemp28344=udnul*vartemp28343;
  FTYPE vartemp28351=udnuc*vartemp28350;
  FTYPE vartemp28352=vartemp28334+vartemp28344+vartemp28351;
  FTYPE vartemp28435=-3.047482e6*uunul;
  FTYPE vartemp28436=-5.810922e6*uunur;
  FTYPE vartemp28437=vartemp28264+vartemp28435+vartemp28436;
  FTYPE vartemp28438=uu0c*vartemp28437;
  FTYPE vartemp28290=-3.047482e6*uunuc;
  FTYPE vartemp28455=uu0l*vartemp28437;
  FTYPE vartemp28456=9.*uu0r*vartemp28011;
  FTYPE vartemp28457=9.256023e6*uunul;
  FTYPE vartemp28458=1.3987665e7*uunur;
  FTYPE vartemp28459=vartemp28300+vartemp28457+vartemp28458;
  FTYPE vartemp28460=uu0c*vartemp28459;
  FTYPE vartemp28461=vartemp28455+vartemp28456+vartemp28460;
  FTYPE vartemp28312=5.9156945e7*uunuc;
  FTYPE vartemp28317=-2.717734e6*uunuc;
  FTYPE vartemp27934=ud0r*vartemp27933;
  FTYPE vartemp27972=ud0c*vartemp27971;
  FTYPE vartemp28002=ud0l*vartemp28001;
  FTYPE vartemp28003=vartemp27934+vartemp27972+vartemp28002;
  FTYPE vartemp28004=gdetr*vartemp28003;
  FTYPE vartemp28059=gdetl*vartemp28058;
  FTYPE vartemp28060=ud0r*vartemp27971;
  FTYPE vartemp28061=ud0l*vartemp28056;
  FTYPE vartemp28062=-6.6807271e7*uu0l*uunuc;
  FTYPE vartemp28063=-5.9156945e7*uu0r*uunuc;
  FTYPE vartemp28064=3.1491054e7*uu0l*uunul;
  FTYPE vartemp28066=1.33619547e8*uunuc;
  FTYPE vartemp28067=-6.6807271e7*uunul;
  FTYPE vartemp28068=-5.9156945e7*uunur;
  FTYPE vartemp28069=vartemp28066+vartemp28067+vartemp28068;
  FTYPE vartemp28070=uu0c*vartemp28069;
  FTYPE vartemp28072=2.4459606e7*uu0r*uunur;
  FTYPE vartemp28073=vartemp28062+vartemp28063+vartemp28064+vartemp28065+vartemp28070+vartemp28071+vartemp28072;
  FTYPE vartemp28074=udnuc*vartemp28073;
  FTYPE vartemp28075=udnur*vartemp27969;
  FTYPE vartemp28076=udnul*vartemp28054;
  FTYPE vartemp28077=vartemp28074+vartemp28075+vartemp28076;
  FTYPE vartemp28078=ud0c*vartemp28077;
  FTYPE vartemp28079=vartemp28060+vartemp28061+vartemp28078;
  FTYPE vartemp28080=gdetc*vartemp28079;
  FTYPE vartemp28081=vartemp28004+vartemp28059+vartemp28080;
  FTYPE vartemp28353=ud0l*vartemp28352;
  FTYPE vartemp28367=ud0c*vartemp28366;
  FTYPE vartemp28368=ud0r*vartemp28029;
  FTYPE vartemp28369=vartemp28353+vartemp28367+vartemp28368;
  FTYPE vartemp28370=gdetl*vartemp28369;
  FTYPE vartemp28371=gdetr*vartemp28058;
  FTYPE vartemp28372=ud0l*vartemp28366;
  FTYPE vartemp28373=ud0r*vartemp28056;
  FTYPE vartemp28374=-5.9156945e7*uu0l*uunuc;
  FTYPE vartemp28375=-6.6807271e7*uu0r*uunuc;
  FTYPE vartemp28376=2.4459606e7*uu0l*uunul;
  FTYPE vartemp28377=1.33609537e8*uunuc;
  FTYPE vartemp28378=-5.9156945e7*uunul;
  FTYPE vartemp28379=-6.6807271e7*uunur;
  FTYPE vartemp28380=vartemp28377+vartemp28378+vartemp28379;
  FTYPE vartemp28381=uu0c*vartemp28380;
  FTYPE vartemp28382=3.1491054e7*uu0r*uunur;
  FTYPE vartemp28383=vartemp28065+vartemp28071+vartemp28374+vartemp28375+vartemp28376+vartemp28381+vartemp28382;
  FTYPE vartemp28384=udnuc*vartemp28383;
  FTYPE vartemp28385=2.717734e6*uu0l*uunuc;
  FTYPE vartemp28386=3.499006e6*uu0r*uunuc;
  FTYPE vartemp28387=-1.028447e6*uu0l*uunul;
  FTYPE vartemp28388=-1.749503e6*uu0r*uunur;
  FTYPE vartemp28389=vartemp27964+vartemp27965+vartemp28385+vartemp28386+vartemp28387+vartemp28388;
  FTYPE vartemp28390=9.*vartemp28389;
  FTYPE vartemp28391=vartemp28363+vartemp28390;
  FTYPE vartemp28392=udnul*vartemp28391;
  FTYPE vartemp28393=udnur*vartemp28054;
  FTYPE vartemp28394=vartemp28384+vartemp28392+vartemp28393;
  FTYPE vartemp28395=ud0c*vartemp28394;
  FTYPE vartemp28396=vartemp28372+vartemp28373+vartemp28395;
  FTYPE vartemp28397=gdetc*vartemp28396;
  FTYPE vartemp28398=vartemp28370+vartemp28371+vartemp28397;
  t[1]=-3.*Bcovc*vartemp28081;
  t[2]=3.*Bcovl*vartemp28203;
  t[3]=ud0c*vartemp28001;
  t[4]=udnuc*vartemp27909;
  t[5]=uu0c*vartemp27897;
  t[6]=uu0r*(738134.*uunur+vartemp28208+vartemp28209);
  t[7]=uu0l*vartemp28216;
  t[8]=udnul*vartemp28223;
  t[9]=udnur*(t[5]+t[6]+t[7]);
  t[10]=t[8];
  t[11]=ud0l*vartemp28245;
  t[12]=ud0r*(t[4]+t[9]+t[10]);
  t[13]=t[11];
  t[14]=gdetl*(t[3]+t[12]+t[13]);
  t[15]=-3.*ud0c*vartemp27933;
  t[16]=3.*ud0l*vartemp28168;
  t[17]=6.643803e6*uu0l*uunuc;
  t[18]=2.213805e6*uu0r*uunuc;
  t[19]=6.643803e6*uu0c*uunul;
  t[20]=-4.144242e6*uu0l*uunul;
  t[21]=2.213805e6*uu0c*uunur;
  t[22]=30436.*uu0r*uunur;
  t[23]=vartemp28249+vartemp28254+vartemp28256;
  t[24]=t[17]+t[18]+t[19];
  t[25]=t[20]+t[21]+t[22]+t[23];
  t[26]=-5.810922e6*uu0l*uunuc;
  t[27]=-3.047482e6*uu0r*uunuc;
  t[28]=3.311055e6*uu0l*uunul;
  t[29]=737935.*uu0r*uunur+vartemp28263+vartemp28268+vartemp28269;
  t[30]=t[26]+t[27]+t[28]+t[29];
  t[31]=-3.311055e6*uu0l*uunuc;
  t[32]=-2.214601e6*uu0r*uunuc;
  t[33]=-3.311055e6*uu0c*uunul;
  t[34]=1.756926e6*uu0l*uunul;
  t[35]=-2.214601e6*uu0c*uunur;
  t[36]=738134.*uu0r*uunur;
  t[37]=vartemp28273+vartemp28278+vartemp28280;
  t[38]=t[31]+t[32]+t[33];
  t[39]=t[34]+t[35]+t[36]+t[37];
  t[40]=udnul*(t[38]+t[39]);
  t[41]=3.*udnuc*t[30];
  t[42]=-3.*t[40];
  t[43]=udnur*(t[24]+t[25]);
  t[44]=t[41]+t[42];
  t[45]=t[16];
  t[46]=ud0r*(t[43]+t[44]);
  t[47]=ud0l*vartemp28117;
  t[48]=udnul*vartemp28114;
  t[49]=uu0l*vartemp28112;
  t[50]=vartemp28268+uu0r*(737935.*uunur+vartemp28129+vartemp28290);
  t[51]=udnuc*vartemp28305;
  t[52]=udnur*(t[49]+t[50]);
  t[53]=t[51];
  t[54]=-9.*udnul*vartemp27983;
  t[55]=udnur*vartemp28305;
  t[56]=-3.1491054e7*uunul;
  t[57]=-2.4459606e7*uunur+vartemp28312;
  t[58]=uu0l*vartemp28103;
  t[59]=uu0r*(1.028447e6*uunur+vartemp28121+vartemp28317);
  t[60]=t[58]+t[59];
  t[61]=uu0c*(t[56]+t[57]);
  t[62]=9.*t[60];
  t[63]=t[55];
  t[64]=udnuc*(t[61]+t[62]);
  t[65]=ud0r*(t[48]+t[52]+t[53]);
  t[66]=ud0c*(t[54]+t[63]+t[64]);
  t[67]=gdetc*(t[47]+t[65]+t[66]);
  t[68]=gdetr*(t[15]+t[45]+t[46]);
  t[69]=3.*t[67];
  t[70]=-3.*t[14];
  t[71]=t[68]+t[69];
  t[72]=t[2];
  t[73]=Bcovr*(t[70]+t[71]);
  t[74]=Bconr*(t[1]+t[72]+t[73]);
  t[75]=8.51887918e8*uunuc;
  t[76]=-4.00828611e8*uunul;
  t[77]=-4.00858641e8*uunur;
  t[78]=-1.33609537e8*uu0l*uunuc;
  t[79]=-1.33619547e8*uu0r*uunuc;
  t[80]=5.9156945e7*uu0l*uunul;
  t[81]=6.6807271e7*uu0r*uunul;
  t[82]=6.6807271e7*uu0l*uunur;
  t[83]=5.9156945e7*uu0r*uunur;
  t[84]=t[78]+t[79]+t[80];
  t[85]=t[81]+t[82]+t[83];
  t[86]=t[84]+t[85];
  t[87]=uu0c*(t[75]+t[76]+t[77]);
  t[88]=3.*t[86];
  t[89]=1.33619547e8*uu0c*uunuc;
  t[90]=-6.6807271e7*uu0c*uunul;
  t[91]=-5.9156945e7*uu0c*uunur+vartemp28062;
  t[92]=vartemp28063+vartemp28064+vartemp28065+vartemp28071+vartemp28072;
  t[93]=t[89]+t[90]+t[91]+t[92];
  t[94]=1.33609537e8*uu0c*uunuc;
  t[95]=-5.9156945e7*uu0c*uunul;
  t[96]=-6.6807271e7*uu0c*uunur+vartemp28065;
  t[97]=vartemp28071+vartemp28374+vartemp28375+vartemp28376+vartemp28382;
  t[98]=t[94]+t[95]+t[96]+t[97];
  t[99]=udnur*t[93]+udnul*t[98];
  t[100]=udnuc*(t[87]+t[88]);
  t[101]=-3.*t[99];
  t[102]=ud0r*vartemp28077+ud0l*vartemp28394;
  t[103]=ud0c*(t[100]+t[101]);
  t[104]=-3.*t[102];
  t[105]=gdetr*vartemp28079+gdetl*vartemp28396;
  t[106]=gdetc*(t[103]+t[104]);
  t[107]=-3.*t[105];
  t[108]=Bcovr*vartemp28081+Bcovl*vartemp28398;
  t[109]=Bcovc*(t[106]+t[107]);
  t[110]=-3.*t[108];
  t[111]=Bconc*(t[109]+t[110]);
  t[112]=3.*Bcovr*vartemp28203;
  t[113]=-3.*Bcovc*vartemp28398;
  t[114]=ud0c*vartemp28029;
  t[115]=udnuc*vartemp28017;
  t[116]=uu0c*vartemp28007;
  t[117]=uu0l*(738134.*uunul+vartemp28208+vartemp28215);
  t[118]=uu0r*vartemp28230;
  t[119]=udnur*vartemp28236;
  t[120]=udnul*(t[116]+t[117]+t[118]);
  t[121]=t[119];
  t[122]=ud0r*vartemp28245;
  t[123]=ud0l*(t[115]+t[120]+t[121]);
  t[124]=t[122];
  t[125]=gdetr*(t[114]+t[123]+t[124]);
  t[126]=3.*ud0r*vartemp28199;
  t[127]=-3.*ud0c*vartemp28352;
  t[128]=2.213805e6*uu0l*uunuc;
  t[129]=6.643803e6*uu0r*uunuc;
  t[130]=2.213805e6*uu0c*uunul;
  t[131]=406.*uu0l*uunul;
  t[132]=6.643803e6*uu0c*uunur;
  t[133]=-4.144242e6*uu0r*uunur;
  t[134]=vartemp28249+vartemp28254+vartemp28256;
  t[135]=t[128]+t[129]+t[130];
  t[136]=t[131]+t[132]+t[133]+t[134];
  t[137]=-2.214601e6*uu0l*uunuc;
  t[138]=-3.311055e6*uu0r*uunuc;
  t[139]=-2.214601e6*uu0c*uunul;
  t[140]=738134.*uu0l*uunul;
  t[141]=-3.311055e6*uu0c*uunur;
  t[142]=1.756926e6*uu0r*uunur;
  t[143]=vartemp28273+vartemp28278+vartemp28280;
  t[144]=t[137]+t[138]+t[139];
  t[145]=t[140]+t[141]+t[142]+t[143];
  t[146]=udnur*(t[144]+t[145]);
  t[147]=-3.047482e6*uu0l*uunuc;
  t[148]=-5.810922e6*uu0r*uunuc;
  t[149]=737935.*uu0l*uunul;
  t[150]=3.311055e6*uu0r*uunur+vartemp28263+vartemp28269+vartemp28438;
  t[151]=t[147]+t[148]+t[149]+t[150];
  t[152]=-3.*t[146];
  t[153]=3.*udnuc*t[151];
  t[154]=udnul*(t[135]+t[136]);
  t[155]=t[152]+t[153];
  t[156]=t[127];
  t[157]=ud0l*(t[154]+t[155]);
  t[158]=ud0r*vartemp28136;
  t[159]=udnur*vartemp28133;
  t[160]=uu0r*vartemp28131;
  t[161]=uu0l*(737935.*uunul+vartemp28111+vartemp28290)+vartemp28438;
  t[162]=udnuc*vartemp28461;
  t[163]=udnul*(t[160]+t[161]);
  t[164]=t[162];
  t[165]=-9.*udnur*vartemp28026;
  t[166]=-2.4459606e7*uunul;
  t[167]=-3.1491054e7*uunur+vartemp28312;
  t[168]=uu0r*vartemp28123;
  t[169]=uu0l*(1.028447e6*uunul+vartemp28102+vartemp28317);
  t[170]=t[168]+t[169];
  t[171]=uu0c*(t[166]+t[167]);
  t[172]=9.*t[170];
  t[173]=udnul*vartemp28461;
  t[174]=udnuc*(t[171]+t[172]);
  t[175]=t[173];
  t[176]=ud0l*(t[159]+t[163]+t[164]);
  t[177]=ud0c*(t[165]+t[174]+t[175]);
  t[178]=gdetc*(t[158]+t[176]+t[177]);
  t[179]=gdetl*(t[126]+t[156]+t[157]);
  t[180]=3.*t[178];
  t[181]=-3.*t[125];
  t[182]=t[179]+t[180];
  t[183]=t[113];
  t[184]=Bcovl*(t[181]+t[182]);
  t[185]=Bconl*(t[112]+t[183]+t[184]);
  t[186]=-1.*t[111];
  t[187]=-1.*t[185];
  t[188]=-1.*t[74];
  t[189]=t[186]+t[187];
  t[190]=t[188]+t[189];
  *Fright=0.000033300033300033300033300033300033*t[190];


 





  return(0);


}















// 2-terms multiplying 2D deconvolution
// See: a2confluxes_bi3terms_2productonly.nb
static int twoterm2Ddeconvolution(
                                  FTYPE Bconc,FTYPE Bconld, FTYPE Bconrd, FTYPE Bconlu, FTYPE Bconru, FTYPE Bconl, FTYPE Bconr, FTYPE Bcond, FTYPE Bconu
                                  ,FTYPE vconc,FTYPE vconld, FTYPE vconrd, FTYPE vconlu, FTYPE vconru, FTYPE vconl, FTYPE vconr, FTYPE vcond, FTYPE vconu
                                  ,FTYPE *Fld, FTYPE *Frd, FTYPE *Flu, FTYPE *Fru
                                  )
{
  FTYPE t[114]; // 114 total in t[]

  t[1]=-858.*Bconl*vconc;
  t[2]=1089.*Bconld*vconc;
  t[3]=99.*Bconlu*vconc;
  t[4]=-78.*Bconr*vconc;
  t[5]=99.*Bconrd*vconc;
  t[6]=9.*Bconru*vconc;
  t[7]=-78.*Bconu*vconc;
  t[8]=1089.*Bconl*vcond;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=33.*Bconld*vcond;
  t[12]=-297.*Bconlu*vcond;
  t[13]=99.*Bconr*vcond;
  t[14]=3.*Bconrd*vcond;
  t[15]=-27.*Bconru*vcond;
  t[16]=234.*Bconu*vcond;
  t[17]=-26.*Bconl*vconl;
  t[18]=33.*Bconld*vconl;
  t[19]=t[11]+t[12]+t[13]+t[14];
  t[20]=t[15]+t[16]+t[17]+t[18];
  t[21]=3.*Bconlu*vconl;
  t[22]=234.*Bconr*vconl;
  t[23]=-297.*Bconrd*vconl;
  t[24]=-27.*Bconru*vconl;
  t[25]=99.*Bconu*vconl;
  t[26]=33.*Bconl*vconld;
  t[27]=-899.*Bconld*vconld;
  t[28]=-9.*Bconlu*vconld;
  t[29]=t[21]+t[22]+t[23]+t[24];
  t[30]=t[25]+t[26]+t[27]+t[28];
  t[31]=-297.*Bconr*vconld;
  t[32]=-9.*Bconrd*vconld;
  t[33]=81.*Bconru*vconld;
  t[34]=-297.*Bconu*vconld;
  t[35]=3.*Bconl*vconlu;
  t[36]=-9.*Bconld*vconlu;
  t[37]=Bconlu*vconlu-27.*Bconr*vconlu;
  t[38]=t[31]+t[32]+t[33];
  t[39]=t[34]+t[35]+t[36]+t[37];
  t[40]=t[9]+t[10]+t[19]+t[20];
  t[41]=t[29]+t[30]+t[38]+t[39];
  t[42]=81.*Bconrd*vconlu;
  t[43]=-9.*Bconru*vconlu;
  t[44]=33.*Bconu*vconlu;
  t[45]=234.*Bconl*vconr;
  t[46]=-297.*Bconld*vconr;
  t[47]=-27.*Bconlu*vconr;
  t[48]=-26.*Bconr*vconr;
  t[49]=33.*Bconrd*vconr;
  t[50]=t[42]+t[43]+t[44]+t[45];
  t[51]=t[46]+t[47]+t[48]+t[49];
  t[52]=3.*Bconru*vconr;
  t[53]=9.*Bconu*vconr;
  t[54]=-297.*Bconl*vconrd;
  t[55]=-9.*Bconld*vconrd;
  t[56]=81.*Bconlu*vconrd;
  t[57]=33.*Bconr*vconrd;
  t[58]=Bconrd*vconrd-9.*Bconru*vconrd;
  t[59]=t[52]+t[53]+t[54];
  t[60]=t[55]+t[56]+t[57]+t[58];
  t[61]=-27.*Bconu*vconrd;
  t[62]=-27.*Bconl*vconru;
  t[63]=81.*Bconld*vconru;
  t[64]=-9.*Bconlu*vconru;
  t[65]=3.*Bconr*vconru;
  t[66]=-9.*Bconrd*vconru;
  t[67]=Bconru*vconru+3.*Bconu*vconru;
  t[68]=t[61]+t[62]+t[63];
  t[69]=t[64]+t[65]+t[66]+t[67];
  t[70]=676.*vconc;
  t[71]=-858.*vcond;
  t[72]=-858.*vconl;
  t[73]=1089.*vconld;
  t[74]=99.*vconlu;
  t[75]=-78.*vconr;
  t[76]=99.*vconrd;
  t[77]=9.*vconru;
  t[78]=-78.*vconu;
  t[79]=t[74]+t[75];
  t[80]=t[76]+t[77]+t[78];
  t[81]=t[70]+t[71]+t[72];
  t[82]=t[73]+t[79]+t[80];
  t[83]=99.*Bconl*vconu;
  t[84]=Bconc*(t[81]+t[82]);
  t[85]=t[83];
  t[86]=-297.*Bconld*vconu;
  t[87]=33.*Bconlu*vconu;
  t[88]=9.*Bconr*vconu;
  t[89]=-27.*Bconrd*vconu;
  t[90]=3.*Bconru*vconu;
  t[91]=-26.*Bconu*vconu;
  t[92]=-858.*vconc;
  t[93]=-26.*vcond;
  t[94]=363.*vconl;
  t[95]=11.*vconld;
  t[96]=-99.*vconlu;
  t[97]=33.*vconr+vconrd;
  t[98]=-9.*vconru;
  t[99]=78.*vconu;
  t[100]=t[94]+t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[100]+t[101];
  t[103]=t[93];
  t[104]=3.*t[102];
  t[105]=t[91];
  t[106]=Bcond*(t[92]+t[103]+t[104]);
  t[107]=t[88]+t[89];
  t[108]=t[90]+t[105]+t[106];
  t[109]=t[84]+t[85]+t[86];
  t[110]=t[87]+t[107]+t[108];
  t[111]=t[50]+t[51]+t[59]+t[60];
  t[112]=t[68]+t[69]+t[109]+t[110];
  t[113]=t[40]+t[41]+t[111]+t[112];
  *Fld=0.0011111111111111111111111111111111*t[113];



  t[1]=-78.*Bconl*vconc;
  t[2]=99.*Bconld*vconc;
  t[3]=9.*Bconlu*vconc;
  t[4]=-858.*Bconr*vconc;
  t[5]=1089.*Bconrd*vconc;
  t[6]=99.*Bconru*vconc;
  t[7]=-78.*Bconu*vconc;
  t[8]=99.*Bconl*vcond;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=3.*Bconld*vcond;
  t[12]=-27.*Bconlu*vcond;
  t[13]=1089.*Bconr*vcond;
  t[14]=33.*Bconrd*vcond;
  t[15]=-297.*Bconru*vcond;
  t[16]=234.*Bconu*vcond;
  t[17]=-26.*Bconl*vconl;
  t[18]=33.*Bconld*vconl;
  t[19]=t[11]+t[12]+t[13]+t[14];
  t[20]=t[15]+t[16]+t[17]+t[18];
  t[21]=3.*Bconlu*vconl;
  t[22]=234.*Bconr*vconl;
  t[23]=-297.*Bconrd*vconl;
  t[24]=-27.*Bconru*vconl;
  t[25]=9.*Bconu*vconl;
  t[26]=33.*Bconl*vconld;
  t[27]=Bconld*vconld-9.*Bconlu*vconld;
  t[28]=t[21]+t[22]+t[23];
  t[29]=t[24]+t[25]+t[26]+t[27];
  t[30]=-297.*Bconr*vconld;
  t[31]=-9.*Bconrd*vconld;
  t[32]=81.*Bconru*vconld;
  t[33]=-27.*Bconu*vconld;
  t[34]=3.*Bconl*vconlu;
  t[35]=-9.*Bconld*vconlu;
  t[36]=Bconlu*vconlu-27.*Bconr*vconlu;
  t[37]=t[30]+t[31]+t[32];
  t[38]=t[33]+t[34]+t[35]+t[36];
  t[39]=t[9]+t[10]+t[19]+t[20];
  t[40]=t[28]+t[29]+t[37]+t[38];
  t[41]=81.*Bconrd*vconlu;
  t[42]=-9.*Bconru*vconlu;
  t[43]=3.*Bconu*vconlu;
  t[44]=234.*Bconl*vconr;
  t[45]=-297.*Bconld*vconr;
  t[46]=-27.*Bconlu*vconr;
  t[47]=-26.*Bconr*vconr;
  t[48]=33.*Bconrd*vconr;
  t[49]=t[41]+t[42]+t[43]+t[44];
  t[50]=t[45]+t[46]+t[47]+t[48];
  t[51]=3.*Bconru*vconr;
  t[52]=99.*Bconu*vconr;
  t[53]=-297.*Bconl*vconrd;
  t[54]=-9.*Bconld*vconrd;
  t[55]=81.*Bconlu*vconrd;
  t[56]=33.*Bconr*vconrd;
  t[57]=-899.*Bconrd*vconrd;
  t[58]=-9.*Bconru*vconrd;
  t[59]=t[51]+t[52]+t[53]+t[54];
  t[60]=t[55]+t[56]+t[57]+t[58];
  t[61]=-297.*Bconu*vconrd;
  t[62]=-27.*Bconl*vconru;
  t[63]=81.*Bconld*vconru;
  t[64]=-9.*Bconlu*vconru;
  t[65]=3.*Bconr*vconru;
  t[66]=-9.*Bconrd*vconru;
  t[67]=Bconru*vconru+33.*Bconu*vconru;
  t[68]=t[61]+t[62]+t[63];
  t[69]=t[64]+t[65]+t[66]+t[67];
  t[70]=676.*vconc;
  t[71]=-858.*vcond;
  t[72]=-78.*vconl;
  t[73]=99.*vconld;
  t[74]=9.*vconlu;
  t[75]=-858.*vconr;
  t[76]=1089.*vconrd;
  t[77]=99.*vconru;
  t[78]=-78.*vconu;
  t[79]=t[74]+t[75];
  t[80]=t[76]+t[77]+t[78];
  t[81]=t[70]+t[71]+t[72];
  t[82]=t[73]+t[79]+t[80];
  t[83]=9.*Bconl*vconu;
  t[84]=Bconc*(t[81]+t[82]);
  t[85]=t[83];
  t[86]=-27.*Bconld*vconu;
  t[87]=3.*Bconlu*vconu;
  t[88]=99.*Bconr*vconu;
  t[89]=-297.*Bconrd*vconu;
  t[90]=33.*Bconru*vconu;
  t[91]=-26.*Bconu*vconu;
  t[92]=-858.*vconc;
  t[93]=-26.*vcond;
  t[94]=33.*vconl;
  t[95]=vconld-9.*vconlu;
  t[96]=363.*vconr;
  t[97]=11.*vconrd;
  t[98]=-99.*vconru;
  t[99]=78.*vconu;
  t[100]=t[94]+t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[100]+t[101];
  t[103]=t[93];
  t[104]=3.*t[102];
  t[105]=t[91];
  t[106]=Bcond*(t[92]+t[103]+t[104]);
  t[107]=t[88]+t[89];
  t[108]=t[90]+t[105]+t[106];
  t[109]=t[84]+t[85]+t[86];
  t[110]=t[87]+t[107]+t[108];
  t[111]=t[49]+t[50]+t[59]+t[60];
  t[112]=t[68]+t[69]+t[109]+t[110];
  t[113]=t[39]+t[40]+t[111]+t[112];
  *Frd=0.0011111111111111111111111111111111*t[113];


  t[1]=-858.*Bconl*vconc;
  t[2]=99.*Bconld*vconc;
  t[3]=1089.*Bconlu*vconc;
  t[4]=-78.*Bconr*vconc;
  t[5]=9.*Bconrd*vconc;
  t[6]=99.*Bconru*vconc;
  t[7]=-858.*Bconu*vconc;
  t[8]=99.*Bconl*vcond;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=33.*Bconld*vcond;
  t[12]=-297.*Bconlu*vcond;
  t[13]=9.*Bconr*vcond;
  t[14]=3.*Bconrd*vcond;
  t[15]=-27.*Bconru*vcond;
  t[16]=234.*Bconu*vcond;
  t[17]=-26.*Bconl*vconl;
  t[18]=3.*Bconld*vconl;
  t[19]=t[11]+t[12]+t[13]+t[14];
  t[20]=t[15]+t[16]+t[17]+t[18];
  t[21]=33.*Bconlu*vconl;
  t[22]=234.*Bconr*vconl;
  t[23]=-27.*Bconrd*vconl;
  t[24]=-297.*Bconru*vconl;
  t[25]=1089.*Bconu*vconl;
  t[26]=3.*Bconl*vconld;
  t[27]=Bconld*vconld-9.*Bconlu*vconld;
  t[28]=t[21]+t[22]+t[23];
  t[29]=t[24]+t[25]+t[26]+t[27];
  t[30]=-27.*Bconr*vconld;
  t[31]=-9.*Bconrd*vconld;
  t[32]=81.*Bconru*vconld;
  t[33]=-297.*Bconu*vconld;
  t[34]=33.*Bconl*vconlu;
  t[35]=-9.*Bconld*vconlu;
  t[36]=-899.*Bconlu*vconlu;
  t[37]=-297.*Bconr*vconlu;
  t[38]=t[30]+t[31]+t[32]+t[33];
  t[39]=t[34]+t[35]+t[36]+t[37];
  t[40]=t[9]+t[10]+t[19]+t[20];
  t[41]=t[28]+t[29]+t[38]+t[39];
  t[42]=81.*Bconrd*vconlu;
  t[43]=-9.*Bconru*vconlu;
  t[44]=33.*Bconu*vconlu;
  t[45]=234.*Bconl*vconr;
  t[46]=-27.*Bconld*vconr;
  t[47]=-297.*Bconlu*vconr;
  t[48]=-26.*Bconr*vconr;
  t[49]=3.*Bconrd*vconr;
  t[50]=t[42]+t[43]+t[44]+t[45];
  t[51]=t[46]+t[47]+t[48]+t[49];
  t[52]=33.*Bconru*vconr;
  t[53]=99.*Bconu*vconr;
  t[54]=-27.*Bconl*vconrd;
  t[55]=-9.*Bconld*vconrd;
  t[56]=81.*Bconlu*vconrd;
  t[57]=3.*Bconr*vconrd;
  t[58]=Bconrd*vconrd-9.*Bconru*vconrd;
  t[59]=t[52]+t[53]+t[54];
  t[60]=t[55]+t[56]+t[57]+t[58];
  t[61]=-27.*Bconu*vconrd;
  t[62]=-297.*Bconl*vconru;
  t[63]=81.*Bconld*vconru;
  t[64]=-9.*Bconlu*vconru;
  t[65]=33.*Bconr*vconru;
  t[66]=-9.*Bconrd*vconru;
  t[67]=Bconru*vconru+3.*Bconu*vconru;
  t[68]=t[61]+t[62]+t[63];
  t[69]=t[64]+t[65]+t[66]+t[67];
  t[70]=676.*vconc;
  t[71]=-78.*vcond;
  t[72]=-858.*vconl;
  t[73]=99.*vconld;
  t[74]=1089.*vconlu;
  t[75]=-78.*vconr;
  t[76]=9.*vconrd;
  t[77]=99.*vconru;
  t[78]=-858.*vconu;
  t[79]=t[74]+t[75];
  t[80]=t[76]+t[77]+t[78];
  t[81]=t[70]+t[71]+t[72];
  t[82]=t[73]+t[79]+t[80];
  t[83]=1089.*Bconl*vconu;
  t[84]=Bconc*(t[81]+t[82]);
  t[85]=t[83];
  t[86]=-297.*Bconld*vconu;
  t[87]=33.*Bconlu*vconu;
  t[88]=99.*Bconr*vconu;
  t[89]=-27.*Bconrd*vconu;
  t[90]=3.*Bconru*vconu;
  t[91]=-26.*Bconu*vconu;
  t[92]=-78.*vconc;
  t[93]=-26.*vcond;
  t[94]=33.*vconl;
  t[95]=11.*vconld;
  t[96]=-99.*vconlu;
  t[97]=3.*vconr+vconrd;
  t[98]=-9.*vconru;
  t[99]=78.*vconu;
  t[100]=t[94]+t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[100]+t[101];
  t[103]=t[93];
  t[104]=3.*t[102];
  t[105]=t[91];
  t[106]=Bcond*(t[92]+t[103]+t[104]);
  t[107]=t[88]+t[89];
  t[108]=t[90]+t[105]+t[106];
  t[109]=t[84]+t[85]+t[86];
  t[110]=t[87]+t[107]+t[108];
  t[111]=t[50]+t[51]+t[59]+t[60];
  t[112]=t[68]+t[69]+t[109]+t[110];
  t[113]=t[40]+t[41]+t[111]+t[112];
  *Flu=0.0011111111111111111111111111111111*t[113];

  t[1]=-78.*Bconl*vconc;
  t[2]=9.*Bconld*vconc;
  t[3]=99.*Bconlu*vconc;
  t[4]=-858.*Bconr*vconc;
  t[5]=99.*Bconrd*vconc;
  t[6]=1089.*Bconru*vconc;
  t[7]=-858.*Bconu*vconc;
  t[8]=9.*Bconl*vcond;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=3.*Bconld*vcond;
  t[12]=-27.*Bconlu*vcond;
  t[13]=99.*Bconr*vcond;
  t[14]=33.*Bconrd*vcond;
  t[15]=-297.*Bconru*vcond;
  t[16]=234.*Bconu*vcond;
  t[17]=-26.*Bconl*vconl;
  t[18]=3.*Bconld*vconl;
  t[19]=t[11]+t[12]+t[13]+t[14];
  t[20]=t[15]+t[16]+t[17]+t[18];
  t[21]=33.*Bconlu*vconl;
  t[22]=234.*Bconr*vconl;
  t[23]=-27.*Bconrd*vconl;
  t[24]=-297.*Bconru*vconl;
  t[25]=99.*Bconu*vconl;
  t[26]=3.*Bconl*vconld;
  t[27]=Bconld*vconld-9.*Bconlu*vconld;
  t[28]=t[21]+t[22]+t[23];
  t[29]=t[24]+t[25]+t[26]+t[27];
  t[30]=-27.*Bconr*vconld;
  t[31]=-9.*Bconrd*vconld;
  t[32]=81.*Bconru*vconld;
  t[33]=-27.*Bconu*vconld;
  t[34]=33.*Bconl*vconlu;
  t[35]=-9.*Bconld*vconlu;
  t[36]=Bconlu*vconlu-297.*Bconr*vconlu;
  t[37]=t[30]+t[31]+t[32];
  t[38]=t[33]+t[34]+t[35]+t[36];
  t[39]=t[9]+t[10]+t[19]+t[20];
  t[40]=t[28]+t[29]+t[37]+t[38];
  t[41]=81.*Bconrd*vconlu;
  t[42]=-9.*Bconru*vconlu;
  t[43]=3.*Bconu*vconlu;
  t[44]=234.*Bconl*vconr;
  t[45]=-27.*Bconld*vconr;
  t[46]=-297.*Bconlu*vconr;
  t[47]=-26.*Bconr*vconr;
  t[48]=3.*Bconrd*vconr;
  t[49]=t[41]+t[42]+t[43]+t[44];
  t[50]=t[45]+t[46]+t[47]+t[48];
  t[51]=33.*Bconru*vconr;
  t[52]=1089.*Bconu*vconr;
  t[53]=-27.*Bconl*vconrd;
  t[54]=-9.*Bconld*vconrd;
  t[55]=81.*Bconlu*vconrd;
  t[56]=3.*Bconr*vconrd;
  t[57]=Bconrd*vconrd-9.*Bconru*vconrd;
  t[58]=t[51]+t[52]+t[53];
  t[59]=t[54]+t[55]+t[56]+t[57];
  t[60]=-297.*Bconu*vconrd;
  t[61]=-297.*Bconl*vconru;
  t[62]=81.*Bconld*vconru;
  t[63]=-9.*Bconlu*vconru;
  t[64]=33.*Bconr*vconru;
  t[65]=-9.*Bconrd*vconru;
  t[66]=-899.*Bconru*vconru;
  t[67]=33.*Bconu*vconru;
  t[68]=t[60]+t[61]+t[62]+t[63];
  t[69]=t[64]+t[65]+t[66]+t[67];
  t[70]=676.*vconc;
  t[71]=-78.*vcond;
  t[72]=-78.*vconl;
  t[73]=9.*vconld;
  t[74]=99.*vconlu;
  t[75]=-858.*vconr;
  t[76]=99.*vconrd;
  t[77]=1089.*vconru;
  t[78]=-858.*vconu;
  t[79]=t[74]+t[75];
  t[80]=t[76]+t[77]+t[78];
  t[81]=t[70]+t[71]+t[72];
  t[82]=t[73]+t[79]+t[80];
  t[83]=99.*Bconl*vconu;
  t[84]=Bconc*(t[81]+t[82]);
  t[85]=t[83];
  t[86]=-27.*Bconld*vconu;
  t[87]=3.*Bconlu*vconu;
  t[88]=1089.*Bconr*vconu;
  t[89]=-297.*Bconrd*vconu;
  t[90]=33.*Bconru*vconu;
  t[91]=-26.*Bconu*vconu;
  t[92]=-78.*vconc;
  t[93]=-26.*vcond;
  t[94]=3.*vconl;
  t[95]=vconld-9.*vconlu;
  t[96]=33.*vconr;
  t[97]=11.*vconrd;
  t[98]=-99.*vconru;
  t[99]=78.*vconu;
  t[100]=t[94]+t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[100]+t[101];
  t[103]=t[93];
  t[104]=3.*t[102];
  t[105]=t[91];
  t[106]=Bcond*(t[92]+t[103]+t[104]);
  t[107]=t[88]+t[89];
  t[108]=t[90]+t[105]+t[106];
  t[109]=t[84]+t[85]+t[86];
  t[110]=t[87]+t[107]+t[108];
  t[111]=t[49]+t[50]+t[58]+t[59];
  t[112]=t[68]+t[69]+t[109]+t[110];
  t[113]=t[39]+t[40]+t[111]+t[112];
  *Fru=0.0011111111111111111111111111111111*t[113];




  return(0);



}







// 3-terms multiplying 2D deconvolution
// See: a2confluxes_bi3terms.nb
static int threeterm2Ddeconvolution(
                                    FTYPE gdetc,FTYPE gdetld, FTYPE gdetrd, FTYPE gdetlu, FTYPE gdetru, FTYPE gdetl, FTYPE gdetr, FTYPE gdetd, FTYPE gdetu
                                    ,FTYPE Bconc,FTYPE Bconld, FTYPE Bconrd, FTYPE Bconlu, FTYPE Bconru, FTYPE Bconl, FTYPE Bconr, FTYPE Bcond, FTYPE Bconu
                                    ,FTYPE vconc,FTYPE vconld, FTYPE vconrd, FTYPE vconlu, FTYPE vconru, FTYPE vconl, FTYPE vconr, FTYPE vcond, FTYPE vconu
                                    ,FTYPE *Fld, FTYPE *Frd, FTYPE *Flu, FTYPE *Fru
                                    )
{
  FTYPE t[1078]; // 1078+12 total and 1078 in t[] and 12 in vartemp





  FTYPE vartemp8718=129034.*vconc;
  FTYPE vartemp8719=-4330.*vcond;
  FTYPE vartemp8720=-22201.*vconl;
  FTYPE vartemp8721=2235.*vconld;
  FTYPE vartemp8722=8493.*vconlu;
  FTYPE vartemp8723=-11771.*vconr;
  FTYPE vartemp8724=1185.*vconrd;
  FTYPE vartemp8725=4503.*vconru;
  FTYPE vartemp8726=-16454.*vconu;
  FTYPE vartemp8727=vartemp8719+vartemp8720+vartemp8721+vartemp8722+vartemp8723+vartemp8724+vartemp8725+vartemp8726;
  FTYPE vartemp8728=3.*vartemp8727;
  FTYPE vartemp8729=vartemp8718+vartemp8728;
  t[1]=-387102.*Bconl*gdetc*vconc;
  t[2]=199809.*Bconld*gdetc*vconc;
  t[3]=105939.*Bconlu*gdetc*vconc;
  t[4]=-205242.*Bconr*gdetc*vconc;
  t[5]=105939.*Bconrd*gdetc*vconc;
  t[6]=56169.*Bconru*gdetc*vconc;
  t[7]=-205242.*Bconu*gdetc*vconc;
  t[8]=199809.*Bconl*gdetd*vconc;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=-20115.*Bconld*gdetd*vconc;
  t[12]=-76437.*Bconlu*gdetd*vconc;
  t[13]=105939.*Bconr*gdetd*vconc;
  t[14]=-10665.*Bconrd*gdetd*vconc;
  t[15]=-40527.*Bconru*gdetd*vconc;
  t[16]=148086.*Bconu*gdetd*vconc;
  t[17]=38970.*Bconl*gdetl*vconc;
  t[18]=-20115.*Bconld*gdetl*vconc;
  t[19]=-10665.*Bconlu*gdetl*vconc;
  t[20]=t[15]+t[16];
  t[21]=t[17]+t[18]+t[19];
  t[22]=t[11]+t[12]+t[13];
  t[23]=t[14]+t[20]+t[21];
  t[24]=148086.*Bconr*gdetl*vconc;
  t[25]=-76437.*Bconrd*gdetl*vconc;
  t[26]=-40527.*Bconru*gdetl*vconc;
  t[27]=105939.*Bconu*gdetl*vconc;
  t[28]=-20115.*Bconl*gdetld*vconc;
  t[29]=2025.*Bconld*gdetld*vconc;
  t[30]=7695.*Bconlu*gdetld*vconc;
  t[31]=-76437.*Bconr*gdetld*vconc;
  t[32]=7695.*Bconrd*gdetld*vconc;
  t[33]=t[28]+t[29];
  t[34]=t[30]+t[31]+t[32];
  t[35]=t[24]+t[25]+t[26];
  t[36]=t[27]+t[33]+t[34];
  t[37]=29241.*Bconru*gdetld*vconc;
  t[38]=-76437.*Bconu*gdetld*vconc;
  t[39]=-10665.*Bconl*gdetlu*vconc;
  t[40]=7695.*Bconld*gdetlu*vconc;
  t[41]=2025.*Bconlu*gdetlu*vconc;
  t[42]=-40527.*Bconr*gdetlu*vconc;
  t[43]=29241.*Bconrd*gdetlu*vconc;
  t[44]=7695.*Bconru*gdetlu*vconc;
  t[45]=-20115.*Bconu*gdetlu*vconc;
  t[46]=t[41]+t[42];
  t[47]=t[43]+t[44]+t[45];
  t[48]=t[37]+t[38]+t[39];
  t[49]=t[40]+t[46]+t[47];
  t[50]=t[9]+t[10]+t[22]+t[23];
  t[51]=t[35]+t[36]+t[48]+t[49];
  t[52]=148086.*Bconl*gdetr*vconc;
  t[53]=-76437.*Bconld*gdetr*vconc;
  t[54]=-40527.*Bconlu*gdetr*vconc;
  t[55]=38970.*Bconr*gdetr*vconc;
  t[56]=-20115.*Bconrd*gdetr*vconc;
  t[57]=-10665.*Bconru*gdetr*vconc;
  t[58]=56169.*Bconu*gdetr*vconc;
  t[59]=-76437.*Bconl*gdetrd*vconc;
  t[60]=7695.*Bconld*gdetrd*vconc;
  t[61]=t[56]+t[57];
  t[62]=t[58]+t[59]+t[60];
  t[63]=t[52]+t[53]+t[54];
  t[64]=t[55]+t[61]+t[62];
  t[65]=29241.*Bconlu*gdetrd*vconc;
  t[66]=-20115.*Bconr*gdetrd*vconc;
  t[67]=2025.*Bconrd*gdetrd*vconc;
  t[68]=7695.*Bconru*gdetrd*vconc;
  t[69]=-40527.*Bconu*gdetrd*vconc;
  t[70]=-40527.*Bconl*gdetru*vconc;
  t[71]=29241.*Bconld*gdetru*vconc;
  t[72]=7695.*Bconlu*gdetru*vconc;
  t[73]=-10665.*Bconr*gdetru*vconc;
  t[74]=t[69]+t[70];
  t[75]=t[71]+t[72]+t[73];
  t[76]=t[65]+t[66]+t[67];
  t[77]=t[68]+t[74]+t[75];
  t[78]=7695.*Bconrd*gdetru*vconc;
  t[79]=2025.*Bconru*gdetru*vconc;
  t[80]=-10665.*Bconu*gdetru*vconc;
  t[81]=105939.*Bconl*gdetu*vconc;
  t[82]=-76437.*Bconld*gdetu*vconc;
  t[83]=-20115.*Bconlu*gdetu*vconc;
  t[84]=56169.*Bconr*gdetu*vconc;
  t[85]=-40527.*Bconrd*gdetu*vconc;
  t[86]=-10665.*Bconru*gdetu*vconc;
  t[87]=t[82]+t[83];
  t[88]=t[84]+t[85]+t[86];
  t[89]=t[78]+t[79]+t[80];
  t[90]=t[81]+t[87]+t[88];
  t[91]=38970.*Bconu*gdetu*vconc;
  t[92]=199809.*Bconl*gdetc*vcond;
  t[93]=-20115.*Bconld*gdetc*vcond;
  t[94]=-76437.*Bconlu*gdetc*vcond;
  t[95]=105939.*Bconr*gdetc*vcond;
  t[96]=-10665.*Bconrd*gdetc*vcond;
  t[97]=-40527.*Bconru*gdetc*vcond;
  t[98]=148086.*Bconu*gdetc*vcond;
  t[99]=-20115.*Bconl*gdetd*vcond;
  t[100]=t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[91]+t[92]+t[93];
  t[103]=t[94]+t[100]+t[101];
  t[104]=t[63]+t[64]+t[76]+t[77];
  t[105]=t[89]+t[90]+t[102]+t[103];
  t[106]=-894.*Bconld*gdetd*vcond;
  t[107]=24138.*Bconlu*gdetd*vcond;
  t[108]=-10665.*Bconr*gdetd*vcond;
  t[109]=-474.*Bconrd*gdetd*vcond;
  t[110]=12798.*Bconru*gdetd*vcond;
  t[111]=-46764.*Bconu*gdetd*vcond;
  t[112]=-20115.*Bconl*gdetl*vcond;
  t[113]=2025.*Bconld*gdetl*vcond;
  t[114]=t[106]+t[107]+t[108]+t[109];
  t[115]=t[110]+t[111]+t[112]+t[113];
  t[116]=7695.*Bconlu*gdetl*vcond;
  t[117]=-76437.*Bconr*gdetl*vcond;
  t[118]=7695.*Bconrd*gdetl*vcond;
  t[119]=29241.*Bconru*gdetl*vcond;
  t[120]=-76437.*Bconu*gdetl*vcond;
  t[121]=2025.*Bconl*gdetld*vcond;
  t[122]=90.*Bconld*gdetld*vcond;
  t[123]=-2430.*Bconlu*gdetld*vcond;
  t[124]=7695.*Bconr*gdetld*vcond;
  t[125]=t[120]+t[121];
  t[126]=t[122]+t[123]+t[124];
  t[127]=t[116]+t[117]+t[118];
  t[128]=t[119]+t[125]+t[126];
  t[129]=342.*Bconrd*gdetld*vcond;
  t[130]=-9234.*Bconru*gdetld*vcond;
  t[131]=24138.*Bconu*gdetld*vcond;
  t[132]=7695.*Bconl*gdetlu*vcond;
  t[133]=-2430.*Bconld*gdetlu*vcond;
  t[134]=-2430.*Bconlu*gdetlu*vcond;
  t[135]=29241.*Bconr*gdetlu*vcond;
  t[136]=-9234.*Bconrd*gdetlu*vcond;
  t[137]=-9234.*Bconru*gdetlu*vcond;
  t[138]=t[133]+t[134];
  t[139]=t[135]+t[136]+t[137];
  t[140]=t[129]+t[130]+t[131];
  t[141]=t[132]+t[138]+t[139];
  t[142]=24138.*Bconu*gdetlu*vcond;
  t[143]=-76437.*Bconl*gdetr*vcond;
  t[144]=7695.*Bconld*gdetr*vcond;
  t[145]=29241.*Bconlu*gdetr*vcond;
  t[146]=-20115.*Bconr*gdetr*vcond;
  t[147]=2025.*Bconrd*gdetr*vcond;
  t[148]=7695.*Bconru*gdetr*vcond;
  t[149]=-40527.*Bconu*gdetr*vcond;
  t[150]=7695.*Bconl*gdetrd*vcond;
  t[151]=t[146]+t[147];
  t[152]=t[148]+t[149]+t[150];
  t[153]=t[142]+t[143]+t[144];
  t[154]=t[145]+t[151]+t[152];
  t[155]=t[114]+t[115]+t[127]+t[128];
  t[156]=t[140]+t[141]+t[153]+t[154];
  t[157]=342.*Bconld*gdetrd*vcond;
  t[158]=-9234.*Bconlu*gdetrd*vcond;
  t[159]=2025.*Bconr*gdetrd*vcond;
  t[160]=90.*Bconrd*gdetrd*vcond;
  t[161]=-2430.*Bconru*gdetrd*vcond;
  t[162]=12798.*Bconu*gdetrd*vcond;
  t[163]=29241.*Bconl*gdetru*vcond;
  t[164]=-9234.*Bconld*gdetru*vcond;
  t[165]=-9234.*Bconlu*gdetru*vcond;
  t[166]=t[161]+t[162];
  t[167]=t[163]+t[164]+t[165];
  t[168]=t[157]+t[158]+t[159];
  t[169]=t[160]+t[166]+t[167];
  t[170]=7695.*Bconr*gdetru*vcond;
  t[171]=-2430.*Bconrd*gdetru*vcond;
  t[172]=-2430.*Bconru*gdetru*vcond;
  t[173]=12798.*Bconu*gdetru*vcond;
  t[174]=-76437.*Bconl*gdetu*vcond;
  t[175]=24138.*Bconld*gdetu*vcond;
  t[176]=24138.*Bconlu*gdetu*vcond;
  t[177]=-40527.*Bconr*gdetu*vcond;
  t[178]=12798.*Bconrd*gdetu*vcond;
  t[179]=t[174]+t[175];
  t[180]=t[176]+t[177]+t[178];
  t[181]=t[170]+t[171]+t[172];
  t[182]=t[173]+t[179]+t[180];
  t[183]=12798.*Bconru*gdetu*vcond;
  t[184]=-46764.*Bconu*gdetu*vcond;
  t[185]=38970.*Bconl*gdetc*vconl;
  t[186]=-20115.*Bconld*gdetc*vconl;
  t[187]=-10665.*Bconlu*gdetc*vconl;
  t[188]=148086.*Bconr*gdetc*vconl;
  t[189]=-76437.*Bconrd*gdetc*vconl;
  t[190]=-40527.*Bconru*gdetc*vconl;
  t[191]=105939.*Bconu*gdetc*vconl;
  t[192]=t[187]+t[188];
  t[193]=t[189]+t[190]+t[191];
  t[194]=t[183]+t[184]+t[185];
  t[195]=t[186]+t[192]+t[193];
  t[196]=-20115.*Bconl*gdetd*vconl;
  t[197]=2025.*Bconld*gdetd*vconl;
  t[198]=7695.*Bconlu*gdetd*vconl;
  t[199]=-76437.*Bconr*gdetd*vconl;
  t[200]=7695.*Bconrd*gdetd*vconl;
  t[201]=29241.*Bconru*gdetd*vconl;
  t[202]=-76437.*Bconu*gdetd*vconl;
  t[203]=1732.*Bconl*gdetl*vconl;
  t[204]=-894.*Bconld*gdetl*vconl;
  t[205]=t[200]+t[201];
  t[206]=t[202]+t[203]+t[204];
  t[207]=t[196]+t[197]+t[198];
  t[208]=t[199]+t[205]+t[206];
  t[209]=t[168]+t[169]+t[181]+t[182];
  t[210]=t[194]+t[195]+t[207]+t[208];
  t[211]=t[50]+t[51]+t[104]+t[105];
  t[212]=t[155]+t[156]+t[209]+t[210];
  t[213]=-474.*Bconlu*gdetl*vconl;
  t[214]=-46764.*Bconr*gdetl*vconl;
  t[215]=24138.*Bconrd*gdetl*vconl;
  t[216]=12798.*Bconru*gdetl*vconl;
  t[217]=-10665.*Bconu*gdetl*vconl;
  t[218]=-894.*Bconl*gdetld*vconl;
  t[219]=90.*Bconld*gdetld*vconl;
  t[220]=342.*Bconlu*gdetld*vconl;
  t[221]=t[213]+t[214]+t[215]+t[216];
  t[222]=t[217]+t[218]+t[219]+t[220];
  t[223]=24138.*Bconr*gdetld*vconl;
  t[224]=-2430.*Bconrd*gdetld*vconl;
  t[225]=-9234.*Bconru*gdetld*vconl;
  t[226]=7695.*Bconu*gdetld*vconl;
  t[227]=-474.*Bconl*gdetlu*vconl;
  t[228]=342.*Bconld*gdetlu*vconl;
  t[229]=90.*Bconlu*gdetlu*vconl;
  t[230]=12798.*Bconr*gdetlu*vconl;
  t[231]=-9234.*Bconrd*gdetlu*vconl;
  t[232]=t[227]+t[228];
  t[233]=t[229]+t[230]+t[231];
  t[234]=t[223]+t[224]+t[225];
  t[235]=t[226]+t[232]+t[233];
  t[236]=-2430.*Bconru*gdetlu*vconl;
  t[237]=2025.*Bconu*gdetlu*vconl;
  t[238]=-46764.*Bconl*gdetr*vconl;
  t[239]=24138.*Bconld*gdetr*vconl;
  t[240]=12798.*Bconlu*gdetr*vconl;
  t[241]=-46764.*Bconr*gdetr*vconl;
  t[242]=24138.*Bconrd*gdetr*vconl;
  t[243]=12798.*Bconru*gdetr*vconl;
  t[244]=-40527.*Bconu*gdetr*vconl;
  t[245]=t[240]+t[241];
  t[246]=t[242]+t[243]+t[244];
  t[247]=t[236]+t[237]+t[238];
  t[248]=t[239]+t[245]+t[246];
  t[249]=24138.*Bconl*gdetrd*vconl;
  t[250]=-2430.*Bconld*gdetrd*vconl;
  t[251]=-9234.*Bconlu*gdetrd*vconl;
  t[252]=24138.*Bconr*gdetrd*vconl;
  t[253]=-2430.*Bconrd*gdetrd*vconl;
  t[254]=-9234.*Bconru*gdetrd*vconl;
  t[255]=29241.*Bconu*gdetrd*vconl;
  t[256]=12798.*Bconl*gdetru*vconl;
  t[257]=-9234.*Bconld*gdetru*vconl;
  t[258]=t[253]+t[254];
  t[259]=t[255]+t[256]+t[257];
  t[260]=t[249]+t[250]+t[251];
  t[261]=t[252]+t[258]+t[259];
  t[262]=t[221]+t[222]+t[234]+t[235];
  t[263]=t[247]+t[248]+t[260]+t[261];
  t[264]=-2430.*Bconlu*gdetru*vconl;
  t[265]=12798.*Bconr*gdetru*vconl;
  t[266]=-9234.*Bconrd*gdetru*vconl;
  t[267]=-2430.*Bconru*gdetru*vconl;
  t[268]=7695.*Bconu*gdetru*vconl;
  t[269]=-10665.*Bconl*gdetu*vconl;
  t[270]=7695.*Bconld*gdetu*vconl;
  t[271]=2025.*Bconlu*gdetu*vconl;
  t[272]=-40527.*Bconr*gdetu*vconl;
  t[273]=t[268]+t[269];
  t[274]=t[270]+t[271]+t[272];
  t[275]=t[264]+t[265]+t[266];
  t[276]=t[267]+t[273]+t[274];
  t[277]=29241.*Bconrd*gdetu*vconl;
  t[278]=7695.*Bconru*gdetu*vconl;
  t[279]=-20115.*Bconu*gdetu*vconl;
  t[280]=-20115.*Bconl*gdetc*vconld;
  t[281]=2025.*Bconld*gdetc*vconld;
  t[282]=7695.*Bconlu*gdetc*vconld;
  t[283]=-76437.*Bconr*gdetc*vconld;
  t[284]=7695.*Bconrd*gdetc*vconld;
  t[285]=29241.*Bconru*gdetc*vconld;
  t[286]=t[281]+t[282];
  t[287]=t[283]+t[284]+t[285];
  t[288]=t[277]+t[278]+t[279];
  t[289]=t[280]+t[286]+t[287];
  t[290]=-76437.*Bconu*gdetc*vconld;
  t[291]=2025.*Bconl*gdetd*vconld;
  t[292]=90.*Bconld*gdetd*vconld;
  t[293]=-2430.*Bconlu*gdetd*vconld;
  t[294]=7695.*Bconr*gdetd*vconld;
  t[295]=342.*Bconrd*gdetd*vconld;
  t[296]=-9234.*Bconru*gdetd*vconld;
  t[297]=24138.*Bconu*gdetd*vconld;
  t[298]=-894.*Bconl*gdetl*vconld;
  t[299]=t[294]+t[295];
  t[300]=t[296]+t[297]+t[298];
  t[301]=t[290]+t[291]+t[292];
  t[302]=t[293]+t[299]+t[300];
  t[303]=90.*Bconld*gdetl*vconld;
  t[304]=342.*Bconlu*gdetl*vconld;
  t[305]=24138.*Bconr*gdetl*vconld;
  t[306]=-2430.*Bconrd*gdetl*vconld;
  t[307]=-9234.*Bconru*gdetl*vconld;
  t[308]=7695.*Bconu*gdetl*vconld;
  t[309]=90.*Bconl*gdetld*vconld;
  t[310]=-44096.*Bconld*gdetld*vconld;
  t[311]=-108.*Bconlu*gdetld*vconld;
  t[312]=t[307]+t[308];
  t[313]=t[309]+t[310]+t[311];
  t[314]=t[303]+t[304]+t[305];
  t[315]=t[306]+t[312]+t[313];
  t[316]=t[275]+t[276]+t[288]+t[289];
  t[317]=t[301]+t[302]+t[314]+t[315];
  t[318]=-2430.*Bconr*gdetld*vconld;
  t[319]=-108.*Bconrd*gdetld*vconld;
  t[320]=2916.*Bconru*gdetld*vconld;
  t[321]=-2430.*Bconu*gdetld*vconld;
  t[322]=342.*Bconl*gdetlu*vconld;
  t[323]=-108.*Bconld*gdetlu*vconld;
  t[324]=-108.*Bconlu*gdetlu*vconld;
  t[325]=-9234.*Bconr*gdetlu*vconld;
  t[326]=t[318]+t[319]+t[320]+t[321];
  t[327]=t[322]+t[323]+t[324]+t[325];
  t[328]=2916.*Bconrd*gdetlu*vconld;
  t[329]=2916.*Bconru*gdetlu*vconld;
  t[330]=-2430.*Bconu*gdetlu*vconld;
  t[331]=24138.*Bconl*gdetr*vconld;
  t[332]=-2430.*Bconld*gdetr*vconld;
  t[333]=-9234.*Bconlu*gdetr*vconld;
  t[334]=24138.*Bconr*gdetr*vconld;
  t[335]=-2430.*Bconrd*gdetr*vconld;
  t[336]=-9234.*Bconru*gdetr*vconld;
  t[337]=t[332]+t[333];
  t[338]=t[334]+t[335]+t[336];
  t[339]=t[328]+t[329]+t[330];
  t[340]=t[331]+t[337]+t[338];
  t[341]=29241.*Bconu*gdetr*vconld;
  t[342]=-2430.*Bconl*gdetrd*vconld;
  t[343]=-108.*Bconld*gdetrd*vconld;
  t[344]=2916.*Bconlu*gdetrd*vconld;
  t[345]=-2430.*Bconr*gdetrd*vconld;
  t[346]=-108.*Bconrd*gdetrd*vconld;
  t[347]=2916.*Bconru*gdetrd*vconld;
  t[348]=-9234.*Bconu*gdetrd*vconld;
  t[349]=-9234.*Bconl*gdetru*vconld;
  t[350]=t[345]+t[346];
  t[351]=t[347]+t[348]+t[349];
  t[352]=t[341]+t[342]+t[343];
  t[353]=t[344]+t[350]+t[351];
  t[354]=2916.*Bconld*gdetru*vconld;
  t[355]=2916.*Bconlu*gdetru*vconld;
  t[356]=-9234.*Bconr*gdetru*vconld;
  t[357]=2916.*Bconrd*gdetru*vconld;
  t[358]=2916.*Bconru*gdetru*vconld;
  t[359]=-9234.*Bconu*gdetru*vconld;
  t[360]=7695.*Bconl*gdetu*vconld;
  t[361]=-2430.*Bconld*gdetu*vconld;
  t[362]=-2430.*Bconlu*gdetu*vconld;
  t[363]=t[358]+t[359];
  t[364]=t[360]+t[361]+t[362];
  t[365]=t[354]+t[355]+t[356];
  t[366]=t[357]+t[363]+t[364];
  t[367]=t[326]+t[327]+t[339]+t[340];
  t[368]=t[352]+t[353]+t[365]+t[366];
  t[369]=29241.*Bconr*gdetu*vconld;
  t[370]=-9234.*Bconrd*gdetu*vconld;
  t[371]=-9234.*Bconru*gdetu*vconld;
  t[372]=24138.*Bconu*gdetu*vconld;
  t[373]=-10665.*Bconl*gdetc*vconlu;
  t[374]=7695.*Bconld*gdetc*vconlu;
  t[375]=2025.*Bconlu*gdetc*vconlu;
  t[376]=-40527.*Bconr*gdetc*vconlu;
  t[377]=29241.*Bconrd*gdetc*vconlu;
  t[378]=t[373]+t[374];
  t[379]=t[375]+t[376]+t[377];
  t[380]=t[369]+t[370]+t[371];
  t[381]=t[372]+t[378]+t[379];
  t[382]=7695.*Bconru*gdetc*vconlu;
  t[383]=-20115.*Bconu*gdetc*vconlu;
  t[384]=7695.*Bconl*gdetd*vconlu;
  t[385]=-2430.*Bconld*gdetd*vconlu;
  t[386]=-2430.*Bconlu*gdetd*vconlu;
  t[387]=29241.*Bconr*gdetd*vconlu;
  t[388]=-9234.*Bconrd*gdetd*vconlu;
  t[389]=-9234.*Bconru*gdetd*vconlu;
  t[390]=24138.*Bconu*gdetd*vconlu;
  t[391]=t[386]+t[387];
  t[392]=t[388]+t[389]+t[390];
  t[393]=t[382]+t[383]+t[384];
  t[394]=t[385]+t[391]+t[392];
  t[395]=-474.*Bconl*gdetl*vconlu;
  t[396]=342.*Bconld*gdetl*vconlu;
  t[397]=90.*Bconlu*gdetl*vconlu;
  t[398]=12798.*Bconr*gdetl*vconlu;
  t[399]=-9234.*Bconrd*gdetl*vconlu;
  t[400]=-2430.*Bconru*gdetl*vconlu;
  t[401]=2025.*Bconu*gdetl*vconlu;
  t[402]=342.*Bconl*gdetld*vconlu;
  t[403]=-108.*Bconld*gdetld*vconlu;
  t[404]=t[399]+t[400];
  t[405]=t[401]+t[402]+t[403];
  t[406]=t[395]+t[396]+t[397];
  t[407]=t[398]+t[404]+t[405];
  t[408]=-108.*Bconlu*gdetld*vconlu;
  t[409]=-9234.*Bconr*gdetld*vconlu;
  t[410]=2916.*Bconrd*gdetld*vconlu;
  t[411]=2916.*Bconru*gdetld*vconlu;
  t[412]=-2430.*Bconu*gdetld*vconlu;
  t[413]=90.*Bconl*gdetlu*vconlu;
  t[414]=-108.*Bconld*gdetlu*vconlu;
  t[415]=4.*Bconlu*gdetlu*vconlu;
  t[416]=-2430.*Bconr*gdetlu*vconlu;
  t[417]=t[412]+t[413];
  t[418]=t[414]+t[415]+t[416];
  t[419]=t[408]+t[409]+t[410];
  t[420]=t[411]+t[417]+t[418];
  t[421]=t[380]+t[381]+t[393]+t[394];
  t[422]=t[406]+t[407]+t[419]+t[420];
  t[423]=t[262]+t[263]+t[316]+t[317];
  t[424]=t[367]+t[368]+t[421]+t[422];
  t[425]=2916.*Bconrd*gdetlu*vconlu;
  t[426]=-108.*Bconru*gdetlu*vconlu;
  t[427]=90.*Bconu*gdetlu*vconlu;
  t[428]=12798.*Bconl*gdetr*vconlu;
  t[429]=-9234.*Bconld*gdetr*vconlu;
  t[430]=-2430.*Bconlu*gdetr*vconlu;
  t[431]=12798.*Bconr*gdetr*vconlu;
  t[432]=-9234.*Bconrd*gdetr*vconlu;
  t[433]=t[425]+t[426]+t[427]+t[428];
  t[434]=t[429]+t[430]+t[431]+t[432];
  t[435]=-2430.*Bconru*gdetr*vconlu;
  t[436]=7695.*Bconu*gdetr*vconlu;
  t[437]=-9234.*Bconl*gdetrd*vconlu;
  t[438]=2916.*Bconld*gdetrd*vconlu;
  t[439]=2916.*Bconlu*gdetrd*vconlu;
  t[440]=-9234.*Bconr*gdetrd*vconlu;
  t[441]=2916.*Bconrd*gdetrd*vconlu;
  t[442]=2916.*Bconru*gdetrd*vconlu;
  t[443]=-9234.*Bconu*gdetrd*vconlu;
  t[444]=t[439]+t[440];
  t[445]=t[441]+t[442]+t[443];
  t[446]=t[435]+t[436]+t[437];
  t[447]=t[438]+t[444]+t[445];
  t[448]=-2430.*Bconl*gdetru*vconlu;
  t[449]=2916.*Bconld*gdetru*vconlu;
  t[450]=-108.*Bconlu*gdetru*vconlu;
  t[451]=-2430.*Bconr*gdetru*vconlu;
  t[452]=2916.*Bconrd*gdetru*vconlu;
  t[453]=-108.*Bconru*gdetru*vconlu;
  t[454]=342.*Bconu*gdetru*vconlu;
  t[455]=2025.*Bconl*gdetu*vconlu;
  t[456]=-2430.*Bconld*gdetu*vconlu;
  t[457]=t[452]+t[453];
  t[458]=t[454]+t[455]+t[456];
  t[459]=t[448]+t[449]+t[450];
  t[460]=t[451]+t[457]+t[458];
  t[461]=90.*Bconlu*gdetu*vconlu;
  t[462]=7695.*Bconr*gdetu*vconlu;
  t[463]=-9234.*Bconrd*gdetu*vconlu;
  t[464]=342.*Bconru*gdetu*vconlu;
  t[465]=-894.*Bconu*gdetu*vconlu;
  t[466]=148086.*Bconl*gdetc*vconr;
  t[467]=-76437.*Bconld*gdetc*vconr;
  t[468]=-40527.*Bconlu*gdetc*vconr;
  t[469]=38970.*Bconr*gdetc*vconr;
  t[470]=t[465]+t[466];
  t[471]=t[467]+t[468]+t[469];
  t[472]=t[461]+t[462]+t[463];
  t[473]=t[464]+t[470]+t[471];
  t[474]=t[433]+t[434]+t[446]+t[447];
  t[475]=t[459]+t[460]+t[472]+t[473];
  t[476]=-20115.*Bconrd*gdetc*vconr;
  t[477]=-10665.*Bconru*gdetc*vconr;
  t[478]=56169.*Bconu*gdetc*vconr;
  t[479]=-76437.*Bconl*gdetd*vconr;
  t[480]=7695.*Bconld*gdetd*vconr;
  t[481]=29241.*Bconlu*gdetd*vconr;
  t[482]=-20115.*Bconr*gdetd*vconr;
  t[483]=2025.*Bconrd*gdetd*vconr;
  t[484]=7695.*Bconru*gdetd*vconr;
  t[485]=t[480]+t[481];
  t[486]=t[482]+t[483]+t[484];
  t[487]=t[476]+t[477]+t[478];
  t[488]=t[479]+t[485]+t[486];
  t[489]=-40527.*Bconu*gdetd*vconr;
  t[490]=-46764.*Bconl*gdetl*vconr;
  t[491]=24138.*Bconld*gdetl*vconr;
  t[492]=12798.*Bconlu*gdetl*vconr;
  t[493]=-46764.*Bconr*gdetl*vconr;
  t[494]=24138.*Bconrd*gdetl*vconr;
  t[495]=12798.*Bconru*gdetl*vconr;
  t[496]=-40527.*Bconu*gdetl*vconr;
  t[497]=24138.*Bconl*gdetld*vconr;
  t[498]=t[493]+t[494];
  t[499]=t[495]+t[496]+t[497];
  t[500]=t[489]+t[490]+t[491];
  t[501]=t[492]+t[498]+t[499];
  t[502]=-2430.*Bconld*gdetld*vconr;
  t[503]=-9234.*Bconlu*gdetld*vconr;
  t[504]=24138.*Bconr*gdetld*vconr;
  t[505]=-2430.*Bconrd*gdetld*vconr;
  t[506]=-9234.*Bconru*gdetld*vconr;
  t[507]=29241.*Bconu*gdetld*vconr;
  t[508]=12798.*Bconl*gdetlu*vconr;
  t[509]=-9234.*Bconld*gdetlu*vconr;
  t[510]=-2430.*Bconlu*gdetlu*vconr;
  t[511]=t[506]+t[507];
  t[512]=t[508]+t[509]+t[510];
  t[513]=t[502]+t[503]+t[504];
  t[514]=t[505]+t[511]+t[512];
  t[515]=12798.*Bconr*gdetlu*vconr;
  t[516]=-9234.*Bconrd*gdetlu*vconr;
  t[517]=-2430.*Bconru*gdetlu*vconr;
  t[518]=7695.*Bconu*gdetlu*vconr;
  t[519]=-46764.*Bconl*gdetr*vconr;
  t[520]=24138.*Bconld*gdetr*vconr;
  t[521]=12798.*Bconlu*gdetr*vconr;
  t[522]=1732.*Bconr*gdetr*vconr;
  t[523]=-894.*Bconrd*gdetr*vconr;
  t[524]=t[519]+t[520];
  t[525]=t[521]+t[522]+t[523];
  t[526]=t[515]+t[516]+t[517];
  t[527]=t[518]+t[524]+t[525];
  t[528]=t[487]+t[488]+t[500]+t[501];
  t[529]=t[513]+t[514]+t[526]+t[527];
  t[530]=-474.*Bconru*gdetr*vconr;
  t[531]=-10665.*Bconu*gdetr*vconr;
  t[532]=24138.*Bconl*gdetrd*vconr;
  t[533]=-2430.*Bconld*gdetrd*vconr;
  t[534]=-9234.*Bconlu*gdetrd*vconr;
  t[535]=-894.*Bconr*gdetrd*vconr;
  t[536]=90.*Bconrd*gdetrd*vconr;
  t[537]=342.*Bconru*gdetrd*vconr;
  t[538]=t[530]+t[531]+t[532]+t[533];
  t[539]=t[534]+t[535]+t[536]+t[537];
  t[540]=7695.*Bconu*gdetrd*vconr;
  t[541]=12798.*Bconl*gdetru*vconr;
  t[542]=-9234.*Bconld*gdetru*vconr;
  t[543]=-2430.*Bconlu*gdetru*vconr;
  t[544]=-474.*Bconr*gdetru*vconr;
  t[545]=342.*Bconrd*gdetru*vconr;
  t[546]=90.*Bconru*gdetru*vconr;
  t[547]=2025.*Bconu*gdetru*vconr;
  t[548]=-40527.*Bconl*gdetu*vconr;
  t[549]=t[544]+t[545];
  t[550]=t[546]+t[547]+t[548];
  t[551]=t[540]+t[541]+t[542];
  t[552]=t[543]+t[549]+t[550];
  t[553]=29241.*Bconld*gdetu*vconr;
  t[554]=7695.*Bconlu*gdetu*vconr;
  t[555]=-10665.*Bconr*gdetu*vconr;
  t[556]=7695.*Bconrd*gdetu*vconr;
  t[557]=2025.*Bconru*gdetu*vconr;
  t[558]=-10665.*Bconu*gdetu*vconr;
  t[559]=-76437.*Bconl*gdetc*vconrd;
  t[560]=7695.*Bconld*gdetc*vconrd;
  t[561]=29241.*Bconlu*gdetc*vconrd;
  t[562]=t[557]+t[558];
  t[563]=t[559]+t[560]+t[561];
  t[564]=t[553]+t[554]+t[555];
  t[565]=t[556]+t[562]+t[563];
  t[566]=-20115.*Bconr*gdetc*vconrd;
  t[567]=2025.*Bconrd*gdetc*vconrd;
  t[568]=7695.*Bconru*gdetc*vconrd;
  t[569]=-40527.*Bconu*gdetc*vconrd;
  t[570]=7695.*Bconl*gdetd*vconrd;
  t[571]=342.*Bconld*gdetd*vconrd;
  t[572]=-9234.*Bconlu*gdetd*vconrd;
  t[573]=2025.*Bconr*gdetd*vconrd;
  t[574]=90.*Bconrd*gdetd*vconrd;
  t[575]=t[570]+t[571];
  t[576]=t[572]+t[573]+t[574];
  t[577]=t[566]+t[567]+t[568];
  t[578]=t[569]+t[575]+t[576];
  t[579]=t[538]+t[539]+t[551]+t[552];
  t[580]=t[564]+t[565]+t[577]+t[578];
  t[581]=-2430.*Bconru*gdetd*vconrd;
  t[582]=12798.*Bconu*gdetd*vconrd;
  t[583]=24138.*Bconl*gdetl*vconrd;
  t[584]=-2430.*Bconld*gdetl*vconrd;
  t[585]=-9234.*Bconlu*gdetl*vconrd;
  t[586]=24138.*Bconr*gdetl*vconrd;
  t[587]=-2430.*Bconrd*gdetl*vconrd;
  t[588]=-9234.*Bconru*gdetl*vconrd;
  t[589]=29241.*Bconu*gdetl*vconrd;
  t[590]=t[585]+t[586];
  t[591]=t[587]+t[588]+t[589];
  t[592]=t[581]+t[582]+t[583];
  t[593]=t[584]+t[590]+t[591];
  t[594]=-2430.*Bconl*gdetld*vconrd;
  t[595]=-108.*Bconld*gdetld*vconrd;
  t[596]=2916.*Bconlu*gdetld*vconrd;
  t[597]=-2430.*Bconr*gdetld*vconrd;
  t[598]=-108.*Bconrd*gdetld*vconrd;
  t[599]=2916.*Bconru*gdetld*vconrd;
  t[600]=-9234.*Bconu*gdetld*vconrd;
  t[601]=-9234.*Bconl*gdetlu*vconrd;
  t[602]=2916.*Bconld*gdetlu*vconrd;
  t[603]=t[598]+t[599];
  t[604]=t[600]+t[601]+t[602];
  t[605]=t[594]+t[595]+t[596];
  t[606]=t[597]+t[603]+t[604];
  t[607]=2916.*Bconlu*gdetlu*vconrd;
  t[608]=-9234.*Bconr*gdetlu*vconrd;
  t[609]=2916.*Bconrd*gdetlu*vconrd;
  t[610]=2916.*Bconru*gdetlu*vconrd;
  t[611]=-9234.*Bconu*gdetlu*vconrd;
  t[612]=24138.*Bconl*gdetr*vconrd;
  t[613]=-2430.*Bconld*gdetr*vconrd;
  t[614]=-9234.*Bconlu*gdetr*vconrd;
  t[615]=-894.*Bconr*gdetr*vconrd;
  t[616]=t[611]+t[612];
  t[617]=t[613]+t[614]+t[615];
  t[618]=t[607]+t[608]+t[609];
  t[619]=t[610]+t[616]+t[617];
  t[620]=90.*Bconrd*gdetr*vconrd;
  t[621]=342.*Bconru*gdetr*vconrd;
  t[622]=7695.*Bconu*gdetr*vconrd;
  t[623]=-2430.*Bconl*gdetrd*vconrd;
  t[624]=-108.*Bconld*gdetrd*vconrd;
  t[625]=2916.*Bconlu*gdetrd*vconrd;
  t[626]=90.*Bconr*gdetrd*vconrd;
  t[627]=4.*Bconrd*gdetrd*vconrd;
  t[628]=-108.*Bconru*gdetrd*vconrd;
  t[629]=t[624]+t[625];
  t[630]=t[626]+t[627]+t[628];
  t[631]=t[620]+t[621]+t[622];
  t[632]=t[623]+t[629]+t[630];
  t[633]=t[592]+t[593]+t[605]+t[606];
  t[634]=t[618]+t[619]+t[631]+t[632];
  t[635]=t[474]+t[475]+t[528]+t[529];
  t[636]=t[579]+t[580]+t[633]+t[634];
  t[637]=-2430.*Bconu*gdetrd*vconrd;
  t[638]=-9234.*Bconl*gdetru*vconrd;
  t[639]=2916.*Bconld*gdetru*vconrd;
  t[640]=2916.*Bconlu*gdetru*vconrd;
  t[641]=342.*Bconr*gdetru*vconrd;
  t[642]=-108.*Bconrd*gdetru*vconrd;
  t[643]=-108.*Bconru*gdetru*vconrd;
  t[644]=-2430.*Bconu*gdetru*vconrd;
  t[645]=t[637]+t[638]+t[639]+t[640];
  t[646]=t[641]+t[642]+t[643]+t[644];
  t[647]=29241.*Bconl*gdetu*vconrd;
  t[648]=-9234.*Bconld*gdetu*vconrd;
  t[649]=-9234.*Bconlu*gdetu*vconrd;
  t[650]=7695.*Bconr*gdetu*vconrd;
  t[651]=-2430.*Bconrd*gdetu*vconrd;
  t[652]=-2430.*Bconru*gdetu*vconrd;
  t[653]=12798.*Bconu*gdetu*vconrd;
  t[654]=-40527.*Bconl*gdetc*vconru;
  t[655]=29241.*Bconld*gdetc*vconru;
  t[656]=t[651]+t[652];
  t[657]=t[653]+t[654]+t[655];
  t[658]=t[647]+t[648]+t[649];
  t[659]=t[650]+t[656]+t[657];
  t[660]=7695.*Bconlu*gdetc*vconru;
  t[661]=-10665.*Bconr*gdetc*vconru;
  t[662]=7695.*Bconrd*gdetc*vconru;
  t[663]=2025.*Bconru*gdetc*vconru;
  t[664]=-10665.*Bconu*gdetc*vconru;
  t[665]=29241.*Bconl*gdetd*vconru;
  t[666]=-9234.*Bconld*gdetd*vconru;
  t[667]=-9234.*Bconlu*gdetd*vconru;
  t[668]=7695.*Bconr*gdetd*vconru;
  t[669]=t[664]+t[665];
  t[670]=t[666]+t[667]+t[668];
  t[671]=t[660]+t[661]+t[662];
  t[672]=t[663]+t[669]+t[670];
  t[673]=-2430.*Bconrd*gdetd*vconru;
  t[674]=-2430.*Bconru*gdetd*vconru;
  t[675]=12798.*Bconu*gdetd*vconru;
  t[676]=12798.*Bconl*gdetl*vconru;
  t[677]=-9234.*Bconld*gdetl*vconru;
  t[678]=-2430.*Bconlu*gdetl*vconru;
  t[679]=12798.*Bconr*gdetl*vconru;
  t[680]=-9234.*Bconrd*gdetl*vconru;
  t[681]=-2430.*Bconru*gdetl*vconru;
  t[682]=t[677]+t[678];
  t[683]=t[679]+t[680]+t[681];
  t[684]=t[673]+t[674]+t[675];
  t[685]=t[676]+t[682]+t[683];
  t[686]=t[645]+t[646]+t[658]+t[659];
  t[687]=t[671]+t[672]+t[684]+t[685];
  t[688]=7695.*Bconu*gdetl*vconru;
  t[689]=-9234.*Bconl*gdetld*vconru;
  t[690]=2916.*Bconld*gdetld*vconru;
  t[691]=2916.*Bconlu*gdetld*vconru;
  t[692]=-9234.*Bconr*gdetld*vconru;
  t[693]=2916.*Bconrd*gdetld*vconru;
  t[694]=2916.*Bconru*gdetld*vconru;
  t[695]=-9234.*Bconu*gdetld*vconru;
  t[696]=-2430.*Bconl*gdetlu*vconru;
  t[697]=t[692]+t[693];
  t[698]=t[694]+t[695]+t[696];
  t[699]=t[688]+t[689]+t[690];
  t[700]=t[691]+t[697]+t[698];
  t[701]=2916.*Bconld*gdetlu*vconru;
  t[702]=-108.*Bconlu*gdetlu*vconru;
  t[703]=-2430.*Bconr*gdetlu*vconru;
  t[704]=2916.*Bconrd*gdetlu*vconru;
  t[705]=-108.*Bconru*gdetlu*vconru;
  t[706]=342.*Bconu*gdetlu*vconru;
  t[707]=12798.*Bconl*gdetr*vconru;
  t[708]=-9234.*Bconld*gdetr*vconru;
  t[709]=-2430.*Bconlu*gdetr*vconru;
  t[710]=t[705]+t[706];
  t[711]=t[707]+t[708]+t[709];
  t[712]=t[701]+t[702]+t[703];
  t[713]=t[704]+t[710]+t[711];
  t[714]=-474.*Bconr*gdetr*vconru;
  t[715]=342.*Bconrd*gdetr*vconru;
  t[716]=90.*Bconru*gdetr*vconru;
  t[717]=2025.*Bconu*gdetr*vconru;
  t[718]=-9234.*Bconl*gdetrd*vconru;
  t[719]=2916.*Bconld*gdetrd*vconru;
  t[720]=2916.*Bconlu*gdetrd*vconru;
  t[721]=342.*Bconr*gdetrd*vconru;
  t[722]=-108.*Bconrd*gdetrd*vconru;
  t[723]=t[718]+t[719];
  t[724]=t[720]+t[721]+t[722];
  t[725]=t[714]+t[715]+t[716];
  t[726]=t[717]+t[723]+t[724];
  t[727]=-108.*Bconru*gdetrd*vconru;
  t[728]=-2430.*Bconu*gdetrd*vconru;
  t[729]=-2430.*Bconl*gdetru*vconru;
  t[730]=2916.*Bconld*gdetru*vconru;
  t[731]=-108.*Bconlu*gdetru*vconru;
  t[732]=90.*Bconr*gdetru*vconru;
  t[733]=-108.*Bconrd*gdetru*vconru;
  t[734]=4.*Bconru*gdetru*vconru;
  t[735]=90.*Bconu*gdetru*vconru;
  t[736]=t[731]+t[732];
  t[737]=t[733]+t[734]+t[735];
  t[738]=t[727]+t[728]+t[729];
  t[739]=t[730]+t[736]+t[737];
  t[740]=t[699]+t[700]+t[712]+t[713];
  t[741]=t[725]+t[726]+t[738]+t[739];
  t[742]=7695.*Bconl*gdetu*vconru;
  t[743]=-9234.*Bconld*gdetu*vconru;
  t[744]=342.*Bconlu*gdetu*vconru;
  t[745]=2025.*Bconr*gdetu*vconru;
  t[746]=-2430.*Bconrd*gdetu*vconru;
  t[747]=90.*Bconru*gdetu*vconru;
  t[748]=-474.*Bconu*gdetu*vconru;
  t[749]=105939.*Bconl*gdetc*vconu;
  t[750]=-76437.*Bconld*gdetc*vconu;
  t[751]=t[746]+t[747];
  t[752]=t[748]+t[749]+t[750];
  t[753]=t[742]+t[743]+t[744];
  t[754]=t[745]+t[751]+t[752];
  t[755]=-20115.*Bconlu*gdetc*vconu;
  t[756]=56169.*Bconr*gdetc*vconu;
  t[757]=-40527.*Bconrd*gdetc*vconu;
  t[758]=-10665.*Bconru*gdetc*vconu;
  t[759]=38970.*Bconu*gdetc*vconu;
  t[760]=-76437.*Bconl*gdetd*vconu;
  t[761]=24138.*Bconld*gdetd*vconu;
  t[762]=24138.*Bconlu*gdetd*vconu;
  t[763]=-40527.*Bconr*gdetd*vconu;
  t[764]=t[759]+t[760];
  t[765]=t[761]+t[762]+t[763];
  t[766]=t[755]+t[756]+t[757];
  t[767]=t[758]+t[764]+t[765];
  t[768]=12798.*Bconrd*gdetd*vconu;
  t[769]=12798.*Bconru*gdetd*vconu;
  t[770]=-46764.*Bconu*gdetd*vconu;
  t[771]=-10665.*Bconl*gdetl*vconu;
  t[772]=7695.*Bconld*gdetl*vconu;
  t[773]=2025.*Bconlu*gdetl*vconu;
  t[774]=-40527.*Bconr*gdetl*vconu;
  t[775]=29241.*Bconrd*gdetl*vconu;
  t[776]=7695.*Bconru*gdetl*vconu;
  t[777]=t[772]+t[773];
  t[778]=t[774]+t[775]+t[776];
  t[779]=t[768]+t[769]+t[770];
  t[780]=t[771]+t[777]+t[778];
  t[781]=-20115.*Bconu*gdetl*vconu;
  t[782]=7695.*Bconl*gdetld*vconu;
  t[783]=-2430.*Bconld*gdetld*vconu;
  t[784]=-2430.*Bconlu*gdetld*vconu;
  t[785]=29241.*Bconr*gdetld*vconu;
  t[786]=-9234.*Bconrd*gdetld*vconu;
  t[787]=-9234.*Bconru*gdetld*vconu;
  t[788]=24138.*Bconu*gdetld*vconu;
  t[789]=2025.*Bconl*gdetlu*vconu;
  t[790]=t[785]+t[786];
  t[791]=t[787]+t[788]+t[789];
  t[792]=t[781]+t[782]+t[783];
  t[793]=t[784]+t[790]+t[791];
  t[794]=t[753]+t[754]+t[766]+t[767];
  t[795]=t[779]+t[780]+t[792]+t[793];
  t[796]=-2430.*Bconld*gdetlu*vconu;
  t[797]=90.*Bconlu*gdetlu*vconu;
  t[798]=7695.*Bconr*gdetlu*vconu;
  t[799]=-9234.*Bconrd*gdetlu*vconu;
  t[800]=342.*Bconru*gdetlu*vconu;
  t[801]=-894.*Bconu*gdetlu*vconu;
  t[802]=-40527.*Bconl*gdetr*vconu;
  t[803]=29241.*Bconld*gdetr*vconu;
  t[804]=7695.*Bconlu*gdetr*vconu;
  t[805]=t[800]+t[801];
  t[806]=t[802]+t[803]+t[804];
  t[807]=t[796]+t[797]+t[798];
  t[808]=t[799]+t[805]+t[806];
  t[809]=-10665.*Bconr*gdetr*vconu;
  t[810]=7695.*Bconrd*gdetr*vconu;
  t[811]=2025.*Bconru*gdetr*vconu;
  t[812]=-10665.*Bconu*gdetr*vconu;
  t[813]=29241.*Bconl*gdetrd*vconu;
  t[814]=-9234.*Bconld*gdetrd*vconu;
  t[815]=-9234.*Bconlu*gdetrd*vconu;
  t[816]=7695.*Bconr*gdetrd*vconu;
  t[817]=-2430.*Bconrd*gdetrd*vconu;
  t[818]=t[813]+t[814];
  t[819]=t[815]+t[816]+t[817];
  t[820]=t[809]+t[810]+t[811];
  t[821]=t[812]+t[818]+t[819];
  t[822]=-2430.*Bconru*gdetrd*vconu;
  t[823]=12798.*Bconu*gdetrd*vconu;
  t[824]=7695.*Bconl*gdetru*vconu;
  t[825]=-9234.*Bconld*gdetru*vconu;
  t[826]=342.*Bconlu*gdetru*vconu;
  t[827]=2025.*Bconr*gdetru*vconu;
  t[828]=-2430.*Bconrd*gdetru*vconu;
  t[829]=90.*Bconru*gdetru*vconu;
  t[830]=-474.*Bconu*gdetru*vconu;
  t[831]=t[826]+t[827];
  t[832]=t[828]+t[829]+t[830];
  t[833]=t[822]+t[823]+t[824];
  t[834]=t[825]+t[831]+t[832];
  t[835]=-20115.*Bconl*gdetu*vconu;
  t[836]=24138.*Bconld*gdetu*vconu;
  t[837]=-894.*Bconlu*gdetu*vconu;
  t[838]=-10665.*Bconr*gdetu*vconu;
  t[839]=12798.*Bconrd*gdetu*vconu;
  t[840]=-474.*Bconru*gdetu*vconu;
  t[841]=1732.*Bconu*gdetu*vconu;
  t[842]=-25479.*gdetlu*vconc;
  t[843]=35313.*gdetr*vconc;
  t[844]=-3555.*gdetrd*vconc;
  t[845]=-13509.*gdetru*vconc;
  t[846]=49362.*gdetu*vconc;
  t[847]=t[842]+t[843];
  t[848]=t[844]+t[845]+t[846];
  t[849]=8046.*gdetlu*vcond;
  t[850]=-3555.*gdetr*vcond;
  t[851]=-158.*gdetrd*vcond;
  t[852]=4266.*gdetru*vcond;
  t[853]=-15588.*gdetu*vcond;
  t[854]=2565.*gdetlu*vconl;
  t[855]=t[849]+t[850]+t[851];
  t[856]=t[852]+t[853]+t[854];
  t[857]=-25479.*gdetr*vconl;
  t[858]=2565.*gdetrd*vconl;
  t[859]=9747.*gdetru*vconl;
  t[860]=-25479.*gdetu*vconl;
  t[861]=-810.*gdetlu*vconld;
  t[862]=2565.*gdetr*vconld;
  t[863]=t[857]+t[858]+t[859];
  t[864]=t[860]+t[861]+t[862];
  t[865]=114.*gdetrd*vconld;
  t[866]=-3078.*gdetru*vconld;
  t[867]=8046.*gdetu*vconld;
  t[868]=-810.*gdetlu*vconlu;
  t[869]=9747.*gdetr*vconlu;
  t[870]=-3078.*gdetrd*vconlu;
  t[871]=t[865]+t[866]+t[867];
  t[872]=t[868]+t[869]+t[870];
  t[873]=t[847]+t[848]+t[855]+t[856];
  t[874]=t[863]+t[864]+t[871]+t[872];
  t[875]=-3078.*gdetru*vconlu;
  t[876]=8046.*gdetu*vconlu;
  t[877]=9747.*gdetlu*vconr;
  t[878]=-6705.*gdetr*vconr;
  t[879]=675.*gdetrd*vconr;
  t[880]=2565.*gdetru*vconr;
  t[881]=t[875]+t[876]+t[877];
  t[882]=t[878]+t[879]+t[880];
  t[883]=-13509.*gdetu*vconr;
  t[884]=-3078.*gdetlu*vconrd;
  t[885]=675.*gdetr*vconrd;
  t[886]=30.*gdetrd*vconrd;
  t[887]=-810.*gdetru*vconrd;
  t[888]=4266.*gdetu*vconrd;
  t[889]=t[883]+t[884]+t[885];
  t[890]=t[886]+t[887]+t[888];
  t[891]=-3078.*gdetlu*vconru;
  t[892]=2565.*gdetr*vconru;
  t[893]=-810.*gdetrd*vconru;
  t[894]=-810.*gdetru*vconru;
  t[895]=4266.*gdetu*vconru;
  t[896]=22201.*vconc;
  t[897]=-745.*vcond;
  t[898]=-745.*vconl;
  t[899]=75.*vconld;
  t[900]=285.*vconlu;
  t[901]=-2831.*vconr;
  t[902]=285.*vconrd;
  t[903]=1083.*vconru;
  t[904]=-2831.*vconu;
  t[905]=t[897]+t[898]+t[899]+t[900];
  t[906]=t[901]+t[902]+t[903]+t[904];
  t[907]=t[905]+t[906];
  t[908]=t[896];
  t[909]=3.*t[907];
  t[910]=gdetl*(t[908]+t[909]);
  t[911]=t[895];
  t[912]=3.*t[910];
  t[913]=t[891]+t[892]+t[893];
  t[914]=t[894]+t[911]+t[912];
  t[915]=8046.*gdetlu*vconu;
  t[916]=-13509.*gdetr*vconu;
  t[917]=4266.*gdetrd*vconu;
  t[918]=4266.*gdetru*vconu;
  t[919]=-15588.*gdetu*vconu;
  t[920]=-6705.*vconc;
  t[921]=-298.*vcond;
  t[922]=675.*vconl;
  t[923]=30.*vconld;
  t[924]=-810.*vconlu;
  t[925]=2565.*vconr;
  t[926]=114.*vconrd;
  t[927]=-3078.*vconru;
  t[928]=8046.*vconu;
  t[929]=t[924]+t[925];
  t[930]=t[926]+t[927]+t[928];
  t[931]=t[920]+t[921]+t[922];
  t[932]=t[923]+t[929]+t[930];
  t[933]=t[919];
  t[934]=gdetld*(t[931]+t[932]);
  t[935]=t[915]+t[916]+t[917];
  t[936]=t[918]+t[933]+t[934];
  t[937]=t[881]+t[882]+t[889]+t[890];
  t[938]=t[913]+t[914]+t[935]+t[936];
  t[939]=t[873]+t[874]+t[937]+t[938];
  t[940]=38970.*vconc;
  t[941]=1732.*vcond;
  t[942]=6705.*vconl;
  t[943]=298.*vconld;
  t[944]=-8046.*vconlu;
  t[945]=3555.*vconr;
  t[946]=158.*vconrd;
  t[947]=-4266.*vconru;
  t[948]=15588.*vconu;
  t[949]=t[942]+t[943]+t[944];
  t[950]=t[945]+t[946]+t[947]+t[948];
  t[951]=t[949]+t[950];
  t[952]=t[941];
  t[953]=-3.*t[951];
  t[954]=-3.*gdetc*vartemp8729;
  t[955]=gdetd*(t[940]+t[952]+t[953]);
  t[956]=t[954];
  t[957]=3.*t[939];
  t[958]=t[955]+t[956];
  t[959]=749956.*vconc;
  t[960]=129034.*vcond;
  t[961]=129034.*vconl;
  t[962]=-66603.*vconld;
  t[963]=-35313.*vconlu;
  t[964]=68414.*vconr;
  t[965]=-35313.*vconrd;
  t[966]=-18723.*vconru;
  t[967]=68414.*vconu;
  t[968]=t[960]+t[961]+t[962]+t[963];
  t[969]=t[964]+t[965]+t[966]+t[967];
  t[970]=t[968]+t[969];
  t[971]=t[959];
  t[972]=-3.*t[970];
  t[973]=-66603.*gdetld*vconc;
  t[974]=-35313.*gdetlu*vconc;
  t[975]=68414.*gdetr*vconc;
  t[976]=-35313.*gdetrd*vconc;
  t[977]=-18723.*gdetru*vconc;
  t[978]=68414.*gdetu*vconc;
  t[979]=6705.*gdetld*vcond;
  t[980]=t[973]+t[974]+t[975];
  t[981]=t[976]+t[977]+t[978]+t[979];
  t[982]=25479.*gdetlu*vcond;
  t[983]=-35313.*gdetr*vcond;
  t[984]=3555.*gdetrd*vcond;
  t[985]=13509.*gdetru*vcond;
  t[986]=-49362.*gdetu*vcond;
  t[987]=6705.*gdetld*vconl;
  t[988]=3555.*gdetlu*vconl;
  t[989]=t[982]+t[983]+t[984];
  t[990]=t[985]+t[986]+t[987]+t[988];
  t[991]=-49362.*gdetr*vconl;
  t[992]=25479.*gdetrd*vconl;
  t[993]=13509.*gdetru*vconl;
  t[994]=-35313.*gdetu*vconl;
  t[995]=-675.*gdetld*vconld;
  t[996]=-2565.*gdetlu*vconld;
  t[997]=25479.*gdetr*vconld;
  t[998]=t[991]+t[992]+t[993];
  t[999]=t[994]+t[995]+t[996]+t[997];
  t[1000]=-2565.*gdetrd*vconld;
  t[1001]=-9747.*gdetru*vconld;
  t[1002]=25479.*gdetu*vconld;
  t[1003]=-2565.*gdetld*vconlu;
  t[1004]=-675.*gdetlu*vconlu;
  t[1005]=13509.*gdetr*vconlu;
  t[1006]=-9747.*gdetrd*vconlu;
  t[1007]=t[1000]+t[1001]+t[1002];
  t[1008]=t[1003]+t[1004]+t[1005]+t[1006];
  t[1009]=t[980]+t[981]+t[989]+t[990];
  t[1010]=t[998]+t[999]+t[1007]+t[1008];
  t[1011]=-2565.*gdetru*vconlu;
  t[1012]=6705.*gdetu*vconlu;
  t[1013]=25479.*gdetld*vconr;
  t[1014]=13509.*gdetlu*vconr;
  t[1015]=-12990.*gdetr*vconr;
  t[1016]=6705.*gdetrd*vconr;
  t[1017]=3555.*gdetru*vconr;
  t[1018]=t[1011]+t[1012]+t[1013];
  t[1019]=t[1014]+t[1015]+t[1016]+t[1017];
  t[1020]=-18723.*gdetu*vconr;
  t[1021]=-2565.*gdetld*vconrd;
  t[1022]=-9747.*gdetlu*vconrd;
  t[1023]=6705.*gdetr*vconrd;
  t[1024]=-675.*gdetrd*vconrd;
  t[1025]=-2565.*gdetru*vconrd;
  t[1026]=13509.*gdetu*vconrd;
  t[1027]=t[1020]+t[1021]+t[1022];
  t[1028]=t[1023]+t[1024]+t[1025]+t[1026];
  t[1029]=-9747.*gdetld*vconru;
  t[1030]=-2565.*gdetlu*vconru;
  t[1031]=3555.*gdetr*vconru;
  t[1032]=-2565.*gdetrd*vconru;
  t[1033]=-675.*gdetru*vconru;
  t[1034]=3555.*gdetu*vconru;
  t[1035]=25479.*gdetld*vconu;
  t[1036]=t[1029]+t[1030]+t[1031];
  t[1037]=t[1032]+t[1033]+t[1034]+t[1035];
  t[1038]=6705.*gdetlu*vconu;
  t[1039]=-18723.*gdetr*vconu;
  t[1040]=13509.*gdetrd*vconu;
  t[1041]=3555.*gdetru*vconu;
  t[1042]=-12990.*gdetu*vconu;
  t[1043]=-66603.*vcond;
  t[1044]=-12990.*vconl;
  t[1045]=6705.*vconld;
  t[1046]=3555.*vconlu;
  t[1047]=-49362.*vconr;
  t[1048]=25479.*vconrd;
  t[1049]=13509.*vconru;
  t[1050]=-35313.*vconu+vartemp8718;
  t[1051]=t[1043]+t[1044]+t[1045]+t[1046];
  t[1052]=t[1047]+t[1048]+t[1049]+t[1050];
  t[1053]=gdetd*vartemp8729;
  t[1054]=gdetl*(t[1051]+t[1052]);
  t[1055]=t[1053];
  t[1056]=t[1038]+t[1039]+t[1040];
  t[1057]=t[1041]+t[1042]+t[1054]+t[1055];
  t[1058]=t[1018]+t[1019]+t[1027]+t[1028];
  t[1059]=t[1036]+t[1037]+t[1056]+t[1057];
  t[1060]=t[1009]+t[1010]+t[1058]+t[1059];
  t[1061]=gdetc*(t[971]+t[972]);
  t[1062]=-3.*t[1060];
  t[1063]=Bcond*(t[957]+t[958]);
  t[1064]=Bconc*(t[1061]+t[1062]);
  t[1065]=t[839]+t[840];
  t[1066]=t[841]+t[1063]+t[1064];
  t[1067]=t[835]+t[836]+t[837];
  t[1068]=t[838]+t[1065]+t[1066];
  t[1069]=t[807]+t[808]+t[820]+t[821];
  t[1070]=t[833]+t[834]+t[1067]+t[1068];
  t[1071]=t[686]+t[687]+t[740]+t[741];
  t[1072]=t[794]+t[795]+t[1069]+t[1070];
  t[1073]=t[211]+t[212]+t[423]+t[424];
  t[1074]=t[635]+t[636]+t[1071]+t[1072];
  t[1075]=t[1073]+t[1074];
  *Fld=0.000022675736961451247165532879818594*t[1075];









  FTYPE vartemp10574=129034.*vconc;
  FTYPE vartemp10575=-4330.*vcond;
  FTYPE vartemp10576=-11771.*vconl;
  FTYPE vartemp10577=1185.*vconld;
  FTYPE vartemp10578=4503.*vconlu;
  FTYPE vartemp10579=-22201.*vconr;
  FTYPE vartemp10580=2235.*vconrd;
  FTYPE vartemp10581=8493.*vconru;
  FTYPE vartemp10582=-16454.*vconu;
  FTYPE vartemp10583=vartemp10575+vartemp10576+vartemp10577+vartemp10578+vartemp10579+vartemp10580+vartemp10581+vartemp10582;
  FTYPE vartemp10584=3.*vartemp10583;
  FTYPE vartemp10585=vartemp10574+vartemp10584;
  t[1]=-205242.*Bconl*gdetc*vconc;
  t[2]=105939.*Bconld*gdetc*vconc;
  t[3]=56169.*Bconlu*gdetc*vconc;
  t[4]=-387102.*Bconr*gdetc*vconc;
  t[5]=199809.*Bconrd*gdetc*vconc;
  t[6]=105939.*Bconru*gdetc*vconc;
  t[7]=-205242.*Bconu*gdetc*vconc;
  t[8]=105939.*Bconl*gdetd*vconc;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=-10665.*Bconld*gdetd*vconc;
  t[12]=-40527.*Bconlu*gdetd*vconc;
  t[13]=199809.*Bconr*gdetd*vconc;
  t[14]=-20115.*Bconrd*gdetd*vconc;
  t[15]=-76437.*Bconru*gdetd*vconc;
  t[16]=148086.*Bconu*gdetd*vconc;
  t[17]=38970.*Bconl*gdetl*vconc;
  t[18]=-20115.*Bconld*gdetl*vconc;
  t[19]=-10665.*Bconlu*gdetl*vconc;
  t[20]=t[15]+t[16];
  t[21]=t[17]+t[18]+t[19];
  t[22]=t[11]+t[12]+t[13];
  t[23]=t[14]+t[20]+t[21];
  t[24]=148086.*Bconr*gdetl*vconc;
  t[25]=-76437.*Bconrd*gdetl*vconc;
  t[26]=-40527.*Bconru*gdetl*vconc;
  t[27]=56169.*Bconu*gdetl*vconc;
  t[28]=-20115.*Bconl*gdetld*vconc;
  t[29]=2025.*Bconld*gdetld*vconc;
  t[30]=7695.*Bconlu*gdetld*vconc;
  t[31]=-76437.*Bconr*gdetld*vconc;
  t[32]=7695.*Bconrd*gdetld*vconc;
  t[33]=t[28]+t[29];
  t[34]=t[30]+t[31]+t[32];
  t[35]=t[24]+t[25]+t[26];
  t[36]=t[27]+t[33]+t[34];
  t[37]=29241.*Bconru*gdetld*vconc;
  t[38]=-40527.*Bconu*gdetld*vconc;
  t[39]=-10665.*Bconl*gdetlu*vconc;
  t[40]=7695.*Bconld*gdetlu*vconc;
  t[41]=2025.*Bconlu*gdetlu*vconc;
  t[42]=-40527.*Bconr*gdetlu*vconc;
  t[43]=29241.*Bconrd*gdetlu*vconc;
  t[44]=7695.*Bconru*gdetlu*vconc;
  t[45]=-10665.*Bconu*gdetlu*vconc;
  t[46]=t[41]+t[42];
  t[47]=t[43]+t[44]+t[45];
  t[48]=t[37]+t[38]+t[39];
  t[49]=t[40]+t[46]+t[47];
  t[50]=t[9]+t[10]+t[22]+t[23];
  t[51]=t[35]+t[36]+t[48]+t[49];
  t[52]=148086.*Bconl*gdetr*vconc;
  t[53]=-76437.*Bconld*gdetr*vconc;
  t[54]=-40527.*Bconlu*gdetr*vconc;
  t[55]=38970.*Bconr*gdetr*vconc;
  t[56]=-20115.*Bconrd*gdetr*vconc;
  t[57]=-10665.*Bconru*gdetr*vconc;
  t[58]=105939.*Bconu*gdetr*vconc;
  t[59]=-76437.*Bconl*gdetrd*vconc;
  t[60]=7695.*Bconld*gdetrd*vconc;
  t[61]=t[56]+t[57];
  t[62]=t[58]+t[59]+t[60];
  t[63]=t[52]+t[53]+t[54];
  t[64]=t[55]+t[61]+t[62];
  t[65]=29241.*Bconlu*gdetrd*vconc;
  t[66]=-20115.*Bconr*gdetrd*vconc;
  t[67]=2025.*Bconrd*gdetrd*vconc;
  t[68]=7695.*Bconru*gdetrd*vconc;
  t[69]=-76437.*Bconu*gdetrd*vconc;
  t[70]=-40527.*Bconl*gdetru*vconc;
  t[71]=29241.*Bconld*gdetru*vconc;
  t[72]=7695.*Bconlu*gdetru*vconc;
  t[73]=-10665.*Bconr*gdetru*vconc;
  t[74]=t[69]+t[70];
  t[75]=t[71]+t[72]+t[73];
  t[76]=t[65]+t[66]+t[67];
  t[77]=t[68]+t[74]+t[75];
  t[78]=7695.*Bconrd*gdetru*vconc;
  t[79]=2025.*Bconru*gdetru*vconc;
  t[80]=-20115.*Bconu*gdetru*vconc;
  t[81]=56169.*Bconl*gdetu*vconc;
  t[82]=-40527.*Bconld*gdetu*vconc;
  t[83]=-10665.*Bconlu*gdetu*vconc;
  t[84]=105939.*Bconr*gdetu*vconc;
  t[85]=-76437.*Bconrd*gdetu*vconc;
  t[86]=-20115.*Bconru*gdetu*vconc;
  t[87]=t[82]+t[83];
  t[88]=t[84]+t[85]+t[86];
  t[89]=t[78]+t[79]+t[80];
  t[90]=t[81]+t[87]+t[88];
  t[91]=38970.*Bconu*gdetu*vconc;
  t[92]=105939.*Bconl*gdetc*vcond;
  t[93]=-10665.*Bconld*gdetc*vcond;
  t[94]=-40527.*Bconlu*gdetc*vcond;
  t[95]=199809.*Bconr*gdetc*vcond;
  t[96]=-20115.*Bconrd*gdetc*vcond;
  t[97]=-76437.*Bconru*gdetc*vcond;
  t[98]=148086.*Bconu*gdetc*vcond;
  t[99]=-10665.*Bconl*gdetd*vcond;
  t[100]=t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[91]+t[92]+t[93];
  t[103]=t[94]+t[100]+t[101];
  t[104]=t[63]+t[64]+t[76]+t[77];
  t[105]=t[89]+t[90]+t[102]+t[103];
  t[106]=-474.*Bconld*gdetd*vcond;
  t[107]=12798.*Bconlu*gdetd*vcond;
  t[108]=-20115.*Bconr*gdetd*vcond;
  t[109]=-894.*Bconrd*gdetd*vcond;
  t[110]=24138.*Bconru*gdetd*vcond;
  t[111]=-46764.*Bconu*gdetd*vcond;
  t[112]=-20115.*Bconl*gdetl*vcond;
  t[113]=2025.*Bconld*gdetl*vcond;
  t[114]=t[106]+t[107]+t[108]+t[109];
  t[115]=t[110]+t[111]+t[112]+t[113];
  t[116]=7695.*Bconlu*gdetl*vcond;
  t[117]=-76437.*Bconr*gdetl*vcond;
  t[118]=7695.*Bconrd*gdetl*vcond;
  t[119]=29241.*Bconru*gdetl*vcond;
  t[120]=-40527.*Bconu*gdetl*vcond;
  t[121]=2025.*Bconl*gdetld*vcond;
  t[122]=90.*Bconld*gdetld*vcond;
  t[123]=-2430.*Bconlu*gdetld*vcond;
  t[124]=7695.*Bconr*gdetld*vcond;
  t[125]=t[120]+t[121];
  t[126]=t[122]+t[123]+t[124];
  t[127]=t[116]+t[117]+t[118];
  t[128]=t[119]+t[125]+t[126];
  t[129]=342.*Bconrd*gdetld*vcond;
  t[130]=-9234.*Bconru*gdetld*vcond;
  t[131]=12798.*Bconu*gdetld*vcond;
  t[132]=7695.*Bconl*gdetlu*vcond;
  t[133]=-2430.*Bconld*gdetlu*vcond;
  t[134]=-2430.*Bconlu*gdetlu*vcond;
  t[135]=29241.*Bconr*gdetlu*vcond;
  t[136]=-9234.*Bconrd*gdetlu*vcond;
  t[137]=-9234.*Bconru*gdetlu*vcond;
  t[138]=t[133]+t[134];
  t[139]=t[135]+t[136]+t[137];
  t[140]=t[129]+t[130]+t[131];
  t[141]=t[132]+t[138]+t[139];
  t[142]=12798.*Bconu*gdetlu*vcond;
  t[143]=-76437.*Bconl*gdetr*vcond;
  t[144]=7695.*Bconld*gdetr*vcond;
  t[145]=29241.*Bconlu*gdetr*vcond;
  t[146]=-20115.*Bconr*gdetr*vcond;
  t[147]=2025.*Bconrd*gdetr*vcond;
  t[148]=7695.*Bconru*gdetr*vcond;
  t[149]=-76437.*Bconu*gdetr*vcond;
  t[150]=7695.*Bconl*gdetrd*vcond;
  t[151]=t[146]+t[147];
  t[152]=t[148]+t[149]+t[150];
  t[153]=t[142]+t[143]+t[144];
  t[154]=t[145]+t[151]+t[152];
  t[155]=t[114]+t[115]+t[127]+t[128];
  t[156]=t[140]+t[141]+t[153]+t[154];
  t[157]=342.*Bconld*gdetrd*vcond;
  t[158]=-9234.*Bconlu*gdetrd*vcond;
  t[159]=2025.*Bconr*gdetrd*vcond;
  t[160]=90.*Bconrd*gdetrd*vcond;
  t[161]=-2430.*Bconru*gdetrd*vcond;
  t[162]=24138.*Bconu*gdetrd*vcond;
  t[163]=29241.*Bconl*gdetru*vcond;
  t[164]=-9234.*Bconld*gdetru*vcond;
  t[165]=-9234.*Bconlu*gdetru*vcond;
  t[166]=t[161]+t[162];
  t[167]=t[163]+t[164]+t[165];
  t[168]=t[157]+t[158]+t[159];
  t[169]=t[160]+t[166]+t[167];
  t[170]=7695.*Bconr*gdetru*vcond;
  t[171]=-2430.*Bconrd*gdetru*vcond;
  t[172]=-2430.*Bconru*gdetru*vcond;
  t[173]=24138.*Bconu*gdetru*vcond;
  t[174]=-40527.*Bconl*gdetu*vcond;
  t[175]=12798.*Bconld*gdetu*vcond;
  t[176]=12798.*Bconlu*gdetu*vcond;
  t[177]=-76437.*Bconr*gdetu*vcond;
  t[178]=24138.*Bconrd*gdetu*vcond;
  t[179]=t[174]+t[175];
  t[180]=t[176]+t[177]+t[178];
  t[181]=t[170]+t[171]+t[172];
  t[182]=t[173]+t[179]+t[180];
  t[183]=24138.*Bconru*gdetu*vcond;
  t[184]=-46764.*Bconu*gdetu*vcond;
  t[185]=38970.*Bconl*gdetc*vconl;
  t[186]=-20115.*Bconld*gdetc*vconl;
  t[187]=-10665.*Bconlu*gdetc*vconl;
  t[188]=148086.*Bconr*gdetc*vconl;
  t[189]=-76437.*Bconrd*gdetc*vconl;
  t[190]=-40527.*Bconru*gdetc*vconl;
  t[191]=56169.*Bconu*gdetc*vconl;
  t[192]=t[187]+t[188];
  t[193]=t[189]+t[190]+t[191];
  t[194]=t[183]+t[184]+t[185];
  t[195]=t[186]+t[192]+t[193];
  t[196]=-20115.*Bconl*gdetd*vconl;
  t[197]=2025.*Bconld*gdetd*vconl;
  t[198]=7695.*Bconlu*gdetd*vconl;
  t[199]=-76437.*Bconr*gdetd*vconl;
  t[200]=7695.*Bconrd*gdetd*vconl;
  t[201]=29241.*Bconru*gdetd*vconl;
  t[202]=-40527.*Bconu*gdetd*vconl;
  t[203]=1732.*Bconl*gdetl*vconl;
  t[204]=-894.*Bconld*gdetl*vconl;
  t[205]=t[200]+t[201];
  t[206]=t[202]+t[203]+t[204];
  t[207]=t[196]+t[197]+t[198];
  t[208]=t[199]+t[205]+t[206];
  t[209]=t[168]+t[169]+t[181]+t[182];
  t[210]=t[194]+t[195]+t[207]+t[208];
  t[211]=t[50]+t[51]+t[104]+t[105];
  t[212]=t[155]+t[156]+t[209]+t[210];
  t[213]=-474.*Bconlu*gdetl*vconl;
  t[214]=-46764.*Bconr*gdetl*vconl;
  t[215]=24138.*Bconrd*gdetl*vconl;
  t[216]=12798.*Bconru*gdetl*vconl;
  t[217]=-10665.*Bconu*gdetl*vconl;
  t[218]=-894.*Bconl*gdetld*vconl;
  t[219]=90.*Bconld*gdetld*vconl;
  t[220]=342.*Bconlu*gdetld*vconl;
  t[221]=t[213]+t[214]+t[215]+t[216];
  t[222]=t[217]+t[218]+t[219]+t[220];
  t[223]=24138.*Bconr*gdetld*vconl;
  t[224]=-2430.*Bconrd*gdetld*vconl;
  t[225]=-9234.*Bconru*gdetld*vconl;
  t[226]=7695.*Bconu*gdetld*vconl;
  t[227]=-474.*Bconl*gdetlu*vconl;
  t[228]=342.*Bconld*gdetlu*vconl;
  t[229]=90.*Bconlu*gdetlu*vconl;
  t[230]=12798.*Bconr*gdetlu*vconl;
  t[231]=-9234.*Bconrd*gdetlu*vconl;
  t[232]=t[227]+t[228];
  t[233]=t[229]+t[230]+t[231];
  t[234]=t[223]+t[224]+t[225];
  t[235]=t[226]+t[232]+t[233];
  t[236]=-2430.*Bconru*gdetlu*vconl;
  t[237]=2025.*Bconu*gdetlu*vconl;
  t[238]=-46764.*Bconl*gdetr*vconl;
  t[239]=24138.*Bconld*gdetr*vconl;
  t[240]=12798.*Bconlu*gdetr*vconl;
  t[241]=-46764.*Bconr*gdetr*vconl;
  t[242]=24138.*Bconrd*gdetr*vconl;
  t[243]=12798.*Bconru*gdetr*vconl;
  t[244]=-40527.*Bconu*gdetr*vconl;
  t[245]=t[240]+t[241];
  t[246]=t[242]+t[243]+t[244];
  t[247]=t[236]+t[237]+t[238];
  t[248]=t[239]+t[245]+t[246];
  t[249]=24138.*Bconl*gdetrd*vconl;
  t[250]=-2430.*Bconld*gdetrd*vconl;
  t[251]=-9234.*Bconlu*gdetrd*vconl;
  t[252]=24138.*Bconr*gdetrd*vconl;
  t[253]=-2430.*Bconrd*gdetrd*vconl;
  t[254]=-9234.*Bconru*gdetrd*vconl;
  t[255]=29241.*Bconu*gdetrd*vconl;
  t[256]=12798.*Bconl*gdetru*vconl;
  t[257]=-9234.*Bconld*gdetru*vconl;
  t[258]=t[253]+t[254];
  t[259]=t[255]+t[256]+t[257];
  t[260]=t[249]+t[250]+t[251];
  t[261]=t[252]+t[258]+t[259];
  t[262]=t[221]+t[222]+t[234]+t[235];
  t[263]=t[247]+t[248]+t[260]+t[261];
  t[264]=-2430.*Bconlu*gdetru*vconl;
  t[265]=12798.*Bconr*gdetru*vconl;
  t[266]=-9234.*Bconrd*gdetru*vconl;
  t[267]=-2430.*Bconru*gdetru*vconl;
  t[268]=7695.*Bconu*gdetru*vconl;
  t[269]=-10665.*Bconl*gdetu*vconl;
  t[270]=7695.*Bconld*gdetu*vconl;
  t[271]=2025.*Bconlu*gdetu*vconl;
  t[272]=-40527.*Bconr*gdetu*vconl;
  t[273]=t[268]+t[269];
  t[274]=t[270]+t[271]+t[272];
  t[275]=t[264]+t[265]+t[266];
  t[276]=t[267]+t[273]+t[274];
  t[277]=29241.*Bconrd*gdetu*vconl;
  t[278]=7695.*Bconru*gdetu*vconl;
  t[279]=-10665.*Bconu*gdetu*vconl;
  t[280]=-20115.*Bconl*gdetc*vconld;
  t[281]=2025.*Bconld*gdetc*vconld;
  t[282]=7695.*Bconlu*gdetc*vconld;
  t[283]=-76437.*Bconr*gdetc*vconld;
  t[284]=7695.*Bconrd*gdetc*vconld;
  t[285]=29241.*Bconru*gdetc*vconld;
  t[286]=t[281]+t[282];
  t[287]=t[283]+t[284]+t[285];
  t[288]=t[277]+t[278]+t[279];
  t[289]=t[280]+t[286]+t[287];
  t[290]=-40527.*Bconu*gdetc*vconld;
  t[291]=2025.*Bconl*gdetd*vconld;
  t[292]=90.*Bconld*gdetd*vconld;
  t[293]=-2430.*Bconlu*gdetd*vconld;
  t[294]=7695.*Bconr*gdetd*vconld;
  t[295]=342.*Bconrd*gdetd*vconld;
  t[296]=-9234.*Bconru*gdetd*vconld;
  t[297]=12798.*Bconu*gdetd*vconld;
  t[298]=-894.*Bconl*gdetl*vconld;
  t[299]=t[294]+t[295];
  t[300]=t[296]+t[297]+t[298];
  t[301]=t[290]+t[291]+t[292];
  t[302]=t[293]+t[299]+t[300];
  t[303]=90.*Bconld*gdetl*vconld;
  t[304]=342.*Bconlu*gdetl*vconld;
  t[305]=24138.*Bconr*gdetl*vconld;
  t[306]=-2430.*Bconrd*gdetl*vconld;
  t[307]=-9234.*Bconru*gdetl*vconld;
  t[308]=7695.*Bconu*gdetl*vconld;
  t[309]=90.*Bconl*gdetld*vconld;
  t[310]=4.*Bconld*gdetld*vconld;
  t[311]=-108.*Bconlu*gdetld*vconld;
  t[312]=t[307]+t[308];
  t[313]=t[309]+t[310]+t[311];
  t[314]=t[303]+t[304]+t[305];
  t[315]=t[306]+t[312]+t[313];
  t[316]=t[275]+t[276]+t[288]+t[289];
  t[317]=t[301]+t[302]+t[314]+t[315];
  t[318]=-2430.*Bconr*gdetld*vconld;
  t[319]=-108.*Bconrd*gdetld*vconld;
  t[320]=2916.*Bconru*gdetld*vconld;
  t[321]=-2430.*Bconu*gdetld*vconld;
  t[322]=342.*Bconl*gdetlu*vconld;
  t[323]=-108.*Bconld*gdetlu*vconld;
  t[324]=-108.*Bconlu*gdetlu*vconld;
  t[325]=-9234.*Bconr*gdetlu*vconld;
  t[326]=t[318]+t[319]+t[320]+t[321];
  t[327]=t[322]+t[323]+t[324]+t[325];
  t[328]=2916.*Bconrd*gdetlu*vconld;
  t[329]=2916.*Bconru*gdetlu*vconld;
  t[330]=-2430.*Bconu*gdetlu*vconld;
  t[331]=24138.*Bconl*gdetr*vconld;
  t[332]=-2430.*Bconld*gdetr*vconld;
  t[333]=-9234.*Bconlu*gdetr*vconld;
  t[334]=24138.*Bconr*gdetr*vconld;
  t[335]=-2430.*Bconrd*gdetr*vconld;
  t[336]=-9234.*Bconru*gdetr*vconld;
  t[337]=t[332]+t[333];
  t[338]=t[334]+t[335]+t[336];
  t[339]=t[328]+t[329]+t[330];
  t[340]=t[331]+t[337]+t[338];
  t[341]=29241.*Bconu*gdetr*vconld;
  t[342]=-2430.*Bconl*gdetrd*vconld;
  t[343]=-108.*Bconld*gdetrd*vconld;
  t[344]=2916.*Bconlu*gdetrd*vconld;
  t[345]=-2430.*Bconr*gdetrd*vconld;
  t[346]=-108.*Bconrd*gdetrd*vconld;
  t[347]=2916.*Bconru*gdetrd*vconld;
  t[348]=-9234.*Bconu*gdetrd*vconld;
  t[349]=-9234.*Bconl*gdetru*vconld;
  t[350]=t[345]+t[346];
  t[351]=t[347]+t[348]+t[349];
  t[352]=t[341]+t[342]+t[343];
  t[353]=t[344]+t[350]+t[351];
  t[354]=2916.*Bconld*gdetru*vconld;
  t[355]=2916.*Bconlu*gdetru*vconld;
  t[356]=-9234.*Bconr*gdetru*vconld;
  t[357]=2916.*Bconrd*gdetru*vconld;
  t[358]=2916.*Bconru*gdetru*vconld;
  t[359]=-9234.*Bconu*gdetru*vconld;
  t[360]=7695.*Bconl*gdetu*vconld;
  t[361]=-2430.*Bconld*gdetu*vconld;
  t[362]=-2430.*Bconlu*gdetu*vconld;
  t[363]=t[358]+t[359];
  t[364]=t[360]+t[361]+t[362];
  t[365]=t[354]+t[355]+t[356];
  t[366]=t[357]+t[363]+t[364];
  t[367]=t[326]+t[327]+t[339]+t[340];
  t[368]=t[352]+t[353]+t[365]+t[366];
  t[369]=29241.*Bconr*gdetu*vconld;
  t[370]=-9234.*Bconrd*gdetu*vconld;
  t[371]=-9234.*Bconru*gdetu*vconld;
  t[372]=12798.*Bconu*gdetu*vconld;
  t[373]=-10665.*Bconl*gdetc*vconlu;
  t[374]=7695.*Bconld*gdetc*vconlu;
  t[375]=2025.*Bconlu*gdetc*vconlu;
  t[376]=-40527.*Bconr*gdetc*vconlu;
  t[377]=29241.*Bconrd*gdetc*vconlu;
  t[378]=t[373]+t[374];
  t[379]=t[375]+t[376]+t[377];
  t[380]=t[369]+t[370]+t[371];
  t[381]=t[372]+t[378]+t[379];
  t[382]=7695.*Bconru*gdetc*vconlu;
  t[383]=-10665.*Bconu*gdetc*vconlu;
  t[384]=7695.*Bconl*gdetd*vconlu;
  t[385]=-2430.*Bconld*gdetd*vconlu;
  t[386]=-2430.*Bconlu*gdetd*vconlu;
  t[387]=29241.*Bconr*gdetd*vconlu;
  t[388]=-9234.*Bconrd*gdetd*vconlu;
  t[389]=-9234.*Bconru*gdetd*vconlu;
  t[390]=12798.*Bconu*gdetd*vconlu;
  t[391]=t[386]+t[387];
  t[392]=t[388]+t[389]+t[390];
  t[393]=t[382]+t[383]+t[384];
  t[394]=t[385]+t[391]+t[392];
  t[395]=-474.*Bconl*gdetl*vconlu;
  t[396]=342.*Bconld*gdetl*vconlu;
  t[397]=90.*Bconlu*gdetl*vconlu;
  t[398]=12798.*Bconr*gdetl*vconlu;
  t[399]=-9234.*Bconrd*gdetl*vconlu;
  t[400]=-2430.*Bconru*gdetl*vconlu;
  t[401]=2025.*Bconu*gdetl*vconlu;
  t[402]=342.*Bconl*gdetld*vconlu;
  t[403]=-108.*Bconld*gdetld*vconlu;
  t[404]=t[399]+t[400];
  t[405]=t[401]+t[402]+t[403];
  t[406]=t[395]+t[396]+t[397];
  t[407]=t[398]+t[404]+t[405];
  t[408]=-108.*Bconlu*gdetld*vconlu;
  t[409]=-9234.*Bconr*gdetld*vconlu;
  t[410]=2916.*Bconrd*gdetld*vconlu;
  t[411]=2916.*Bconru*gdetld*vconlu;
  t[412]=-2430.*Bconu*gdetld*vconlu;
  t[413]=90.*Bconl*gdetlu*vconlu;
  t[414]=-108.*Bconld*gdetlu*vconlu;
  t[415]=4.*Bconlu*gdetlu*vconlu;
  t[416]=-2430.*Bconr*gdetlu*vconlu;
  t[417]=t[412]+t[413];
  t[418]=t[414]+t[415]+t[416];
  t[419]=t[408]+t[409]+t[410];
  t[420]=t[411]+t[417]+t[418];
  t[421]=t[380]+t[381]+t[393]+t[394];
  t[422]=t[406]+t[407]+t[419]+t[420];
  t[423]=t[262]+t[263]+t[316]+t[317];
  t[424]=t[367]+t[368]+t[421]+t[422];
  t[425]=2916.*Bconrd*gdetlu*vconlu;
  t[426]=-108.*Bconru*gdetlu*vconlu;
  t[427]=90.*Bconu*gdetlu*vconlu;
  t[428]=12798.*Bconl*gdetr*vconlu;
  t[429]=-9234.*Bconld*gdetr*vconlu;
  t[430]=-2430.*Bconlu*gdetr*vconlu;
  t[431]=12798.*Bconr*gdetr*vconlu;
  t[432]=-9234.*Bconrd*gdetr*vconlu;
  t[433]=t[425]+t[426]+t[427]+t[428];
  t[434]=t[429]+t[430]+t[431]+t[432];
  t[435]=-2430.*Bconru*gdetr*vconlu;
  t[436]=7695.*Bconu*gdetr*vconlu;
  t[437]=-9234.*Bconl*gdetrd*vconlu;
  t[438]=2916.*Bconld*gdetrd*vconlu;
  t[439]=2916.*Bconlu*gdetrd*vconlu;
  t[440]=-9234.*Bconr*gdetrd*vconlu;
  t[441]=2916.*Bconrd*gdetrd*vconlu;
  t[442]=2916.*Bconru*gdetrd*vconlu;
  t[443]=-9234.*Bconu*gdetrd*vconlu;
  t[444]=t[439]+t[440];
  t[445]=t[441]+t[442]+t[443];
  t[446]=t[435]+t[436]+t[437];
  t[447]=t[438]+t[444]+t[445];
  t[448]=-2430.*Bconl*gdetru*vconlu;
  t[449]=2916.*Bconld*gdetru*vconlu;
  t[450]=-108.*Bconlu*gdetru*vconlu;
  t[451]=-2430.*Bconr*gdetru*vconlu;
  t[452]=2916.*Bconrd*gdetru*vconlu;
  t[453]=-108.*Bconru*gdetru*vconlu;
  t[454]=342.*Bconu*gdetru*vconlu;
  t[455]=2025.*Bconl*gdetu*vconlu;
  t[456]=-2430.*Bconld*gdetu*vconlu;
  t[457]=t[452]+t[453];
  t[458]=t[454]+t[455]+t[456];
  t[459]=t[448]+t[449]+t[450];
  t[460]=t[451]+t[457]+t[458];
  t[461]=90.*Bconlu*gdetu*vconlu;
  t[462]=7695.*Bconr*gdetu*vconlu;
  t[463]=-9234.*Bconrd*gdetu*vconlu;
  t[464]=342.*Bconru*gdetu*vconlu;
  t[465]=-474.*Bconu*gdetu*vconlu;
  t[466]=148086.*Bconl*gdetc*vconr;
  t[467]=-76437.*Bconld*gdetc*vconr;
  t[468]=-40527.*Bconlu*gdetc*vconr;
  t[469]=38970.*Bconr*gdetc*vconr;
  t[470]=t[465]+t[466];
  t[471]=t[467]+t[468]+t[469];
  t[472]=t[461]+t[462]+t[463];
  t[473]=t[464]+t[470]+t[471];
  t[474]=t[433]+t[434]+t[446]+t[447];
  t[475]=t[459]+t[460]+t[472]+t[473];
  t[476]=-20115.*Bconrd*gdetc*vconr;
  t[477]=-10665.*Bconru*gdetc*vconr;
  t[478]=105939.*Bconu*gdetc*vconr;
  t[479]=-76437.*Bconl*gdetd*vconr;
  t[480]=7695.*Bconld*gdetd*vconr;
  t[481]=29241.*Bconlu*gdetd*vconr;
  t[482]=-20115.*Bconr*gdetd*vconr;
  t[483]=2025.*Bconrd*gdetd*vconr;
  t[484]=7695.*Bconru*gdetd*vconr;
  t[485]=t[480]+t[481];
  t[486]=t[482]+t[483]+t[484];
  t[487]=t[476]+t[477]+t[478];
  t[488]=t[479]+t[485]+t[486];
  t[489]=-76437.*Bconu*gdetd*vconr;
  t[490]=-46764.*Bconl*gdetl*vconr;
  t[491]=24138.*Bconld*gdetl*vconr;
  t[492]=12798.*Bconlu*gdetl*vconr;
  t[493]=-46764.*Bconr*gdetl*vconr;
  t[494]=24138.*Bconrd*gdetl*vconr;
  t[495]=12798.*Bconru*gdetl*vconr;
  t[496]=-40527.*Bconu*gdetl*vconr;
  t[497]=24138.*Bconl*gdetld*vconr;
  t[498]=t[493]+t[494];
  t[499]=t[495]+t[496]+t[497];
  t[500]=t[489]+t[490]+t[491];
  t[501]=t[492]+t[498]+t[499];
  t[502]=-2430.*Bconld*gdetld*vconr;
  t[503]=-9234.*Bconlu*gdetld*vconr;
  t[504]=24138.*Bconr*gdetld*vconr;
  t[505]=-2430.*Bconrd*gdetld*vconr;
  t[506]=-9234.*Bconru*gdetld*vconr;
  t[507]=29241.*Bconu*gdetld*vconr;
  t[508]=12798.*Bconl*gdetlu*vconr;
  t[509]=-9234.*Bconld*gdetlu*vconr;
  t[510]=-2430.*Bconlu*gdetlu*vconr;
  t[511]=t[506]+t[507];
  t[512]=t[508]+t[509]+t[510];
  t[513]=t[502]+t[503]+t[504];
  t[514]=t[505]+t[511]+t[512];
  t[515]=12798.*Bconr*gdetlu*vconr;
  t[516]=-9234.*Bconrd*gdetlu*vconr;
  t[517]=-2430.*Bconru*gdetlu*vconr;
  t[518]=7695.*Bconu*gdetlu*vconr;
  t[519]=-46764.*Bconl*gdetr*vconr;
  t[520]=24138.*Bconld*gdetr*vconr;
  t[521]=12798.*Bconlu*gdetr*vconr;
  t[522]=1732.*Bconr*gdetr*vconr;
  t[523]=-894.*Bconrd*gdetr*vconr;
  t[524]=t[519]+t[520];
  t[525]=t[521]+t[522]+t[523];
  t[526]=t[515]+t[516]+t[517];
  t[527]=t[518]+t[524]+t[525];
  t[528]=t[487]+t[488]+t[500]+t[501];
  t[529]=t[513]+t[514]+t[526]+t[527];
  t[530]=-474.*Bconru*gdetr*vconr;
  t[531]=-10665.*Bconu*gdetr*vconr;
  t[532]=24138.*Bconl*gdetrd*vconr;
  t[533]=-2430.*Bconld*gdetrd*vconr;
  t[534]=-9234.*Bconlu*gdetrd*vconr;
  t[535]=-894.*Bconr*gdetrd*vconr;
  t[536]=90.*Bconrd*gdetrd*vconr;
  t[537]=342.*Bconru*gdetrd*vconr;
  t[538]=t[530]+t[531]+t[532]+t[533];
  t[539]=t[534]+t[535]+t[536]+t[537];
  t[540]=7695.*Bconu*gdetrd*vconr;
  t[541]=12798.*Bconl*gdetru*vconr;
  t[542]=-9234.*Bconld*gdetru*vconr;
  t[543]=-2430.*Bconlu*gdetru*vconr;
  t[544]=-474.*Bconr*gdetru*vconr;
  t[545]=342.*Bconrd*gdetru*vconr;
  t[546]=90.*Bconru*gdetru*vconr;
  t[547]=2025.*Bconu*gdetru*vconr;
  t[548]=-40527.*Bconl*gdetu*vconr;
  t[549]=t[544]+t[545];
  t[550]=t[546]+t[547]+t[548];
  t[551]=t[540]+t[541]+t[542];
  t[552]=t[543]+t[549]+t[550];
  t[553]=29241.*Bconld*gdetu*vconr;
  t[554]=7695.*Bconlu*gdetu*vconr;
  t[555]=-10665.*Bconr*gdetu*vconr;
  t[556]=7695.*Bconrd*gdetu*vconr;
  t[557]=2025.*Bconru*gdetu*vconr;
  t[558]=-20115.*Bconu*gdetu*vconr;
  t[559]=-76437.*Bconl*gdetc*vconrd;
  t[560]=7695.*Bconld*gdetc*vconrd;
  t[561]=29241.*Bconlu*gdetc*vconrd;
  t[562]=t[557]+t[558];
  t[563]=t[559]+t[560]+t[561];
  t[564]=t[553]+t[554]+t[555];
  t[565]=t[556]+t[562]+t[563];
  t[566]=-20115.*Bconr*gdetc*vconrd;
  t[567]=2025.*Bconrd*gdetc*vconrd;
  t[568]=7695.*Bconru*gdetc*vconrd;
  t[569]=-76437.*Bconu*gdetc*vconrd;
  t[570]=7695.*Bconl*gdetd*vconrd;
  t[571]=342.*Bconld*gdetd*vconrd;
  t[572]=-9234.*Bconlu*gdetd*vconrd;
  t[573]=2025.*Bconr*gdetd*vconrd;
  t[574]=90.*Bconrd*gdetd*vconrd;
  t[575]=t[570]+t[571];
  t[576]=t[572]+t[573]+t[574];
  t[577]=t[566]+t[567]+t[568];
  t[578]=t[569]+t[575]+t[576];
  t[579]=t[538]+t[539]+t[551]+t[552];
  t[580]=t[564]+t[565]+t[577]+t[578];
  t[581]=-2430.*Bconru*gdetd*vconrd;
  t[582]=24138.*Bconu*gdetd*vconrd;
  t[583]=24138.*Bconl*gdetl*vconrd;
  t[584]=-2430.*Bconld*gdetl*vconrd;
  t[585]=-9234.*Bconlu*gdetl*vconrd;
  t[586]=24138.*Bconr*gdetl*vconrd;
  t[587]=-2430.*Bconrd*gdetl*vconrd;
  t[588]=-9234.*Bconru*gdetl*vconrd;
  t[589]=29241.*Bconu*gdetl*vconrd;
  t[590]=t[585]+t[586];
  t[591]=t[587]+t[588]+t[589];
  t[592]=t[581]+t[582]+t[583];
  t[593]=t[584]+t[590]+t[591];
  t[594]=-2430.*Bconl*gdetld*vconrd;
  t[595]=-108.*Bconld*gdetld*vconrd;
  t[596]=2916.*Bconlu*gdetld*vconrd;
  t[597]=-2430.*Bconr*gdetld*vconrd;
  t[598]=-108.*Bconrd*gdetld*vconrd;
  t[599]=2916.*Bconru*gdetld*vconrd;
  t[600]=-9234.*Bconu*gdetld*vconrd;
  t[601]=-9234.*Bconl*gdetlu*vconrd;
  t[602]=2916.*Bconld*gdetlu*vconrd;
  t[603]=t[598]+t[599];
  t[604]=t[600]+t[601]+t[602];
  t[605]=t[594]+t[595]+t[596];
  t[606]=t[597]+t[603]+t[604];
  t[607]=2916.*Bconlu*gdetlu*vconrd;
  t[608]=-9234.*Bconr*gdetlu*vconrd;
  t[609]=2916.*Bconrd*gdetlu*vconrd;
  t[610]=2916.*Bconru*gdetlu*vconrd;
  t[611]=-9234.*Bconu*gdetlu*vconrd;
  t[612]=24138.*Bconl*gdetr*vconrd;
  t[613]=-2430.*Bconld*gdetr*vconrd;
  t[614]=-9234.*Bconlu*gdetr*vconrd;
  t[615]=-894.*Bconr*gdetr*vconrd;
  t[616]=t[611]+t[612];
  t[617]=t[613]+t[614]+t[615];
  t[618]=t[607]+t[608]+t[609];
  t[619]=t[610]+t[616]+t[617];
  t[620]=90.*Bconrd*gdetr*vconrd;
  t[621]=342.*Bconru*gdetr*vconrd;
  t[622]=7695.*Bconu*gdetr*vconrd;
  t[623]=-2430.*Bconl*gdetrd*vconrd;
  t[624]=-108.*Bconld*gdetrd*vconrd;
  t[625]=2916.*Bconlu*gdetrd*vconrd;
  t[626]=90.*Bconr*gdetrd*vconrd;
  t[627]=-44096.*Bconrd*gdetrd*vconrd;
  t[628]=-108.*Bconru*gdetrd*vconrd;
  t[629]=t[624]+t[625];
  t[630]=t[626]+t[627]+t[628];
  t[631]=t[620]+t[621]+t[622];
  t[632]=t[623]+t[629]+t[630];
  t[633]=t[592]+t[593]+t[605]+t[606];
  t[634]=t[618]+t[619]+t[631]+t[632];
  t[635]=t[474]+t[475]+t[528]+t[529];
  t[636]=t[579]+t[580]+t[633]+t[634];
  t[637]=-2430.*Bconu*gdetrd*vconrd;
  t[638]=-9234.*Bconl*gdetru*vconrd;
  t[639]=2916.*Bconld*gdetru*vconrd;
  t[640]=2916.*Bconlu*gdetru*vconrd;
  t[641]=342.*Bconr*gdetru*vconrd;
  t[642]=-108.*Bconrd*gdetru*vconrd;
  t[643]=-108.*Bconru*gdetru*vconrd;
  t[644]=-2430.*Bconu*gdetru*vconrd;
  t[645]=t[637]+t[638]+t[639]+t[640];
  t[646]=t[641]+t[642]+t[643]+t[644];
  t[647]=29241.*Bconl*gdetu*vconrd;
  t[648]=-9234.*Bconld*gdetu*vconrd;
  t[649]=-9234.*Bconlu*gdetu*vconrd;
  t[650]=7695.*Bconr*gdetu*vconrd;
  t[651]=-2430.*Bconrd*gdetu*vconrd;
  t[652]=-2430.*Bconru*gdetu*vconrd;
  t[653]=24138.*Bconu*gdetu*vconrd;
  t[654]=-40527.*Bconl*gdetc*vconru;
  t[655]=29241.*Bconld*gdetc*vconru;
  t[656]=t[651]+t[652];
  t[657]=t[653]+t[654]+t[655];
  t[658]=t[647]+t[648]+t[649];
  t[659]=t[650]+t[656]+t[657];
  t[660]=7695.*Bconlu*gdetc*vconru;
  t[661]=-10665.*Bconr*gdetc*vconru;
  t[662]=7695.*Bconrd*gdetc*vconru;
  t[663]=2025.*Bconru*gdetc*vconru;
  t[664]=-20115.*Bconu*gdetc*vconru;
  t[665]=29241.*Bconl*gdetd*vconru;
  t[666]=-9234.*Bconld*gdetd*vconru;
  t[667]=-9234.*Bconlu*gdetd*vconru;
  t[668]=7695.*Bconr*gdetd*vconru;
  t[669]=t[664]+t[665];
  t[670]=t[666]+t[667]+t[668];
  t[671]=t[660]+t[661]+t[662];
  t[672]=t[663]+t[669]+t[670];
  t[673]=-2430.*Bconrd*gdetd*vconru;
  t[674]=-2430.*Bconru*gdetd*vconru;
  t[675]=24138.*Bconu*gdetd*vconru;
  t[676]=12798.*Bconl*gdetl*vconru;
  t[677]=-9234.*Bconld*gdetl*vconru;
  t[678]=-2430.*Bconlu*gdetl*vconru;
  t[679]=12798.*Bconr*gdetl*vconru;
  t[680]=-9234.*Bconrd*gdetl*vconru;
  t[681]=-2430.*Bconru*gdetl*vconru;
  t[682]=t[677]+t[678];
  t[683]=t[679]+t[680]+t[681];
  t[684]=t[673]+t[674]+t[675];
  t[685]=t[676]+t[682]+t[683];
  t[686]=t[645]+t[646]+t[658]+t[659];
  t[687]=t[671]+t[672]+t[684]+t[685];
  t[688]=7695.*Bconu*gdetl*vconru;
  t[689]=-9234.*Bconl*gdetld*vconru;
  t[690]=2916.*Bconld*gdetld*vconru;
  t[691]=2916.*Bconlu*gdetld*vconru;
  t[692]=-9234.*Bconr*gdetld*vconru;
  t[693]=2916.*Bconrd*gdetld*vconru;
  t[694]=2916.*Bconru*gdetld*vconru;
  t[695]=-9234.*Bconu*gdetld*vconru;
  t[696]=-2430.*Bconl*gdetlu*vconru;
  t[697]=t[692]+t[693];
  t[698]=t[694]+t[695]+t[696];
  t[699]=t[688]+t[689]+t[690];
  t[700]=t[691]+t[697]+t[698];
  t[701]=2916.*Bconld*gdetlu*vconru;
  t[702]=-108.*Bconlu*gdetlu*vconru;
  t[703]=-2430.*Bconr*gdetlu*vconru;
  t[704]=2916.*Bconrd*gdetlu*vconru;
  t[705]=-108.*Bconru*gdetlu*vconru;
  t[706]=342.*Bconu*gdetlu*vconru;
  t[707]=12798.*Bconl*gdetr*vconru;
  t[708]=-9234.*Bconld*gdetr*vconru;
  t[709]=-2430.*Bconlu*gdetr*vconru;
  t[710]=t[705]+t[706];
  t[711]=t[707]+t[708]+t[709];
  t[712]=t[701]+t[702]+t[703];
  t[713]=t[704]+t[710]+t[711];
  t[714]=-474.*Bconr*gdetr*vconru;
  t[715]=342.*Bconrd*gdetr*vconru;
  t[716]=90.*Bconru*gdetr*vconru;
  t[717]=2025.*Bconu*gdetr*vconru;
  t[718]=-9234.*Bconl*gdetrd*vconru;
  t[719]=2916.*Bconld*gdetrd*vconru;
  t[720]=2916.*Bconlu*gdetrd*vconru;
  t[721]=342.*Bconr*gdetrd*vconru;
  t[722]=-108.*Bconrd*gdetrd*vconru;
  t[723]=t[718]+t[719];
  t[724]=t[720]+t[721]+t[722];
  t[725]=t[714]+t[715]+t[716];
  t[726]=t[717]+t[723]+t[724];
  t[727]=-108.*Bconru*gdetrd*vconru;
  t[728]=-2430.*Bconu*gdetrd*vconru;
  t[729]=-2430.*Bconl*gdetru*vconru;
  t[730]=2916.*Bconld*gdetru*vconru;
  t[731]=-108.*Bconlu*gdetru*vconru;
  t[732]=90.*Bconr*gdetru*vconru;
  t[733]=-108.*Bconrd*gdetru*vconru;
  t[734]=4.*Bconru*gdetru*vconru;
  t[735]=90.*Bconu*gdetru*vconru;
  t[736]=t[731]+t[732];
  t[737]=t[733]+t[734]+t[735];
  t[738]=t[727]+t[728]+t[729];
  t[739]=t[730]+t[736]+t[737];
  t[740]=t[699]+t[700]+t[712]+t[713];
  t[741]=t[725]+t[726]+t[738]+t[739];
  t[742]=7695.*Bconl*gdetu*vconru;
  t[743]=-9234.*Bconld*gdetu*vconru;
  t[744]=342.*Bconlu*gdetu*vconru;
  t[745]=2025.*Bconr*gdetu*vconru;
  t[746]=-2430.*Bconrd*gdetu*vconru;
  t[747]=90.*Bconru*gdetu*vconru;
  t[748]=-894.*Bconu*gdetu*vconru;
  t[749]=56169.*Bconl*gdetc*vconu;
  t[750]=-40527.*Bconld*gdetc*vconu;
  t[751]=t[746]+t[747];
  t[752]=t[748]+t[749]+t[750];
  t[753]=t[742]+t[743]+t[744];
  t[754]=t[745]+t[751]+t[752];
  t[755]=-10665.*Bconlu*gdetc*vconu;
  t[756]=105939.*Bconr*gdetc*vconu;
  t[757]=-76437.*Bconrd*gdetc*vconu;
  t[758]=-20115.*Bconru*gdetc*vconu;
  t[759]=38970.*Bconu*gdetc*vconu;
  t[760]=-40527.*Bconl*gdetd*vconu;
  t[761]=12798.*Bconld*gdetd*vconu;
  t[762]=12798.*Bconlu*gdetd*vconu;
  t[763]=-76437.*Bconr*gdetd*vconu;
  t[764]=t[759]+t[760];
  t[765]=t[761]+t[762]+t[763];
  t[766]=t[755]+t[756]+t[757];
  t[767]=t[758]+t[764]+t[765];
  t[768]=24138.*Bconrd*gdetd*vconu;
  t[769]=24138.*Bconru*gdetd*vconu;
  t[770]=-46764.*Bconu*gdetd*vconu;
  t[771]=-10665.*Bconl*gdetl*vconu;
  t[772]=7695.*Bconld*gdetl*vconu;
  t[773]=2025.*Bconlu*gdetl*vconu;
  t[774]=-40527.*Bconr*gdetl*vconu;
  t[775]=29241.*Bconrd*gdetl*vconu;
  t[776]=7695.*Bconru*gdetl*vconu;
  t[777]=t[772]+t[773];
  t[778]=t[774]+t[775]+t[776];
  t[779]=t[768]+t[769]+t[770];
  t[780]=t[771]+t[777]+t[778];
  t[781]=-10665.*Bconu*gdetl*vconu;
  t[782]=7695.*Bconl*gdetld*vconu;
  t[783]=-2430.*Bconld*gdetld*vconu;
  t[784]=-2430.*Bconlu*gdetld*vconu;
  t[785]=29241.*Bconr*gdetld*vconu;
  t[786]=-9234.*Bconrd*gdetld*vconu;
  t[787]=-9234.*Bconru*gdetld*vconu;
  t[788]=12798.*Bconu*gdetld*vconu;
  t[789]=2025.*Bconl*gdetlu*vconu;
  t[790]=t[785]+t[786];
  t[791]=t[787]+t[788]+t[789];
  t[792]=t[781]+t[782]+t[783];
  t[793]=t[784]+t[790]+t[791];
  t[794]=t[753]+t[754]+t[766]+t[767];
  t[795]=t[779]+t[780]+t[792]+t[793];
  t[796]=-2430.*Bconld*gdetlu*vconu;
  t[797]=90.*Bconlu*gdetlu*vconu;
  t[798]=7695.*Bconr*gdetlu*vconu;
  t[799]=-9234.*Bconrd*gdetlu*vconu;
  t[800]=342.*Bconru*gdetlu*vconu;
  t[801]=-474.*Bconu*gdetlu*vconu;
  t[802]=-40527.*Bconl*gdetr*vconu;
  t[803]=29241.*Bconld*gdetr*vconu;
  t[804]=7695.*Bconlu*gdetr*vconu;
  t[805]=t[800]+t[801];
  t[806]=t[802]+t[803]+t[804];
  t[807]=t[796]+t[797]+t[798];
  t[808]=t[799]+t[805]+t[806];
  t[809]=-10665.*Bconr*gdetr*vconu;
  t[810]=7695.*Bconrd*gdetr*vconu;
  t[811]=2025.*Bconru*gdetr*vconu;
  t[812]=-20115.*Bconu*gdetr*vconu;
  t[813]=29241.*Bconl*gdetrd*vconu;
  t[814]=-9234.*Bconld*gdetrd*vconu;
  t[815]=-9234.*Bconlu*gdetrd*vconu;
  t[816]=7695.*Bconr*gdetrd*vconu;
  t[817]=-2430.*Bconrd*gdetrd*vconu;
  t[818]=t[813]+t[814];
  t[819]=t[815]+t[816]+t[817];
  t[820]=t[809]+t[810]+t[811];
  t[821]=t[812]+t[818]+t[819];
  t[822]=-2430.*Bconru*gdetrd*vconu;
  t[823]=24138.*Bconu*gdetrd*vconu;
  t[824]=7695.*Bconl*gdetru*vconu;
  t[825]=-9234.*Bconld*gdetru*vconu;
  t[826]=342.*Bconlu*gdetru*vconu;
  t[827]=2025.*Bconr*gdetru*vconu;
  t[828]=-2430.*Bconrd*gdetru*vconu;
  t[829]=90.*Bconru*gdetru*vconu;
  t[830]=-894.*Bconu*gdetru*vconu;
  t[831]=t[826]+t[827];
  t[832]=t[828]+t[829]+t[830];
  t[833]=t[822]+t[823]+t[824];
  t[834]=t[825]+t[831]+t[832];
  t[835]=-10665.*Bconl*gdetu*vconu;
  t[836]=12798.*Bconld*gdetu*vconu;
  t[837]=-474.*Bconlu*gdetu*vconu;
  t[838]=-20115.*Bconr*gdetu*vconu;
  t[839]=24138.*Bconrd*gdetu*vconu;
  t[840]=-894.*Bconru*gdetu*vconu;
  t[841]=1732.*Bconu*gdetu*vconu;
  t[842]=-13509.*gdetlu*vconc;
  t[843]=66603.*gdetr*vconc;
  t[844]=-6705.*gdetrd*vconc;
  t[845]=-25479.*gdetru*vconc;
  t[846]=49362.*gdetu*vconc;
  t[847]=t[842]+t[843];
  t[848]=t[844]+t[845]+t[846];
  t[849]=4266.*gdetlu*vcond;
  t[850]=-6705.*gdetr*vcond;
  t[851]=-298.*gdetrd*vcond;
  t[852]=8046.*gdetru*vcond;
  t[853]=-15588.*gdetu*vcond;
  t[854]=2565.*gdetlu*vconl;
  t[855]=t[849]+t[850]+t[851];
  t[856]=t[852]+t[853]+t[854];
  t[857]=-25479.*gdetr*vconl;
  t[858]=2565.*gdetrd*vconl;
  t[859]=9747.*gdetru*vconl;
  t[860]=-13509.*gdetu*vconl;
  t[861]=-810.*gdetlu*vconld;
  t[862]=2565.*gdetr*vconld;
  t[863]=t[857]+t[858]+t[859];
  t[864]=t[860]+t[861]+t[862];
  t[865]=114.*gdetrd*vconld;
  t[866]=-3078.*gdetru*vconld;
  t[867]=4266.*gdetu*vconld;
  t[868]=-810.*gdetlu*vconlu;
  t[869]=9747.*gdetr*vconlu;
  t[870]=-3078.*gdetrd*vconlu;
  t[871]=t[865]+t[866]+t[867];
  t[872]=t[868]+t[869]+t[870];
  t[873]=t[847]+t[848]+t[855]+t[856];
  t[874]=t[863]+t[864]+t[871]+t[872];
  t[875]=-3078.*gdetru*vconlu;
  t[876]=4266.*gdetu*vconlu;
  t[877]=9747.*gdetlu*vconr;
  t[878]=-6705.*gdetr*vconr;
  t[879]=675.*gdetrd*vconr;
  t[880]=2565.*gdetru*vconr;
  t[881]=t[875]+t[876]+t[877];
  t[882]=t[878]+t[879]+t[880];
  t[883]=-25479.*gdetu*vconr;
  t[884]=-3078.*gdetlu*vconrd;
  t[885]=675.*gdetr*vconrd;
  t[886]=30.*gdetrd*vconrd;
  t[887]=-810.*gdetru*vconrd;
  t[888]=8046.*gdetu*vconrd;
  t[889]=t[883]+t[884]+t[885];
  t[890]=t[886]+t[887]+t[888];
  t[891]=-3078.*gdetlu*vconru;
  t[892]=2565.*gdetr*vconru;
  t[893]=-810.*gdetrd*vconru;
  t[894]=-810.*gdetru*vconru;
  t[895]=8046.*gdetu*vconru;
  t[896]=11771.*vconc;
  t[897]=-395.*vcond;
  t[898]=-745.*vconl;
  t[899]=75.*vconld;
  t[900]=285.*vconlu;
  t[901]=-2831.*vconr;
  t[902]=285.*vconrd;
  t[903]=1083.*vconru;
  t[904]=-1501.*vconu;
  t[905]=t[897]+t[898]+t[899]+t[900];
  t[906]=t[901]+t[902]+t[903]+t[904];
  t[907]=t[905]+t[906];
  t[908]=t[896];
  t[909]=3.*t[907];
  t[910]=gdetl*(t[908]+t[909]);
  t[911]=t[895];
  t[912]=3.*t[910];
  t[913]=t[891]+t[892]+t[893];
  t[914]=t[894]+t[911]+t[912];
  t[915]=4266.*gdetlu*vconu;
  t[916]=-25479.*gdetr*vconu;
  t[917]=8046.*gdetrd*vconu;
  t[918]=8046.*gdetru*vconu;
  t[919]=-15588.*gdetu*vconu;
  t[920]=-3555.*vconc;
  t[921]=-158.*vcond;
  t[922]=675.*vconl;
  t[923]=30.*vconld;
  t[924]=-810.*vconlu;
  t[925]=2565.*vconr;
  t[926]=114.*vconrd;
  t[927]=-3078.*vconru;
  t[928]=4266.*vconu;
  t[929]=t[924]+t[925];
  t[930]=t[926]+t[927]+t[928];
  t[931]=t[920]+t[921]+t[922];
  t[932]=t[923]+t[929]+t[930];
  t[933]=t[919];
  t[934]=gdetld*(t[931]+t[932]);
  t[935]=t[915]+t[916]+t[917];
  t[936]=t[918]+t[933]+t[934];
  t[937]=t[881]+t[882]+t[889]+t[890];
  t[938]=t[913]+t[914]+t[935]+t[936];
  t[939]=t[873]+t[874]+t[937]+t[938];
  t[940]=38970.*vconc;
  t[941]=1732.*vcond;
  t[942]=3555.*vconl;
  t[943]=158.*vconld;
  t[944]=-4266.*vconlu;
  t[945]=6705.*vconr;
  t[946]=298.*vconrd;
  t[947]=-8046.*vconru;
  t[948]=15588.*vconu;
  t[949]=t[942]+t[943]+t[944];
  t[950]=t[945]+t[946]+t[947]+t[948];
  t[951]=t[949]+t[950];
  t[952]=t[941];
  t[953]=-3.*t[951];
  t[954]=-3.*gdetc*vartemp10585;
  t[955]=gdetd*(t[940]+t[952]+t[953]);
  t[956]=t[954];
  t[957]=3.*t[939];
  t[958]=t[955]+t[956];
  t[959]=749956.*vconc;
  t[960]=129034.*vcond;
  t[961]=68414.*vconl;
  t[962]=-35313.*vconld;
  t[963]=-18723.*vconlu;
  t[964]=129034.*vconr;
  t[965]=-66603.*vconrd;
  t[966]=-35313.*vconru;
  t[967]=68414.*vconu;
  t[968]=t[960]+t[961]+t[962]+t[963];
  t[969]=t[964]+t[965]+t[966]+t[967];
  t[970]=t[968]+t[969];
  t[971]=t[959];
  t[972]=-3.*t[970];
  t[973]=-35313.*gdetld*vconc;
  t[974]=-18723.*gdetlu*vconc;
  t[975]=129034.*gdetr*vconc;
  t[976]=-66603.*gdetrd*vconc;
  t[977]=-35313.*gdetru*vconc;
  t[978]=68414.*gdetu*vconc;
  t[979]=3555.*gdetld*vcond;
  t[980]=t[973]+t[974]+t[975];
  t[981]=t[976]+t[977]+t[978]+t[979];
  t[982]=13509.*gdetlu*vcond;
  t[983]=-66603.*gdetr*vcond;
  t[984]=6705.*gdetrd*vcond;
  t[985]=25479.*gdetru*vcond;
  t[986]=-49362.*gdetu*vcond;
  t[987]=6705.*gdetld*vconl;
  t[988]=3555.*gdetlu*vconl;
  t[989]=t[982]+t[983]+t[984];
  t[990]=t[985]+t[986]+t[987]+t[988];
  t[991]=-49362.*gdetr*vconl;
  t[992]=25479.*gdetrd*vconl;
  t[993]=13509.*gdetru*vconl;
  t[994]=-18723.*gdetu*vconl;
  t[995]=-675.*gdetld*vconld;
  t[996]=-2565.*gdetlu*vconld;
  t[997]=25479.*gdetr*vconld;
  t[998]=t[991]+t[992]+t[993];
  t[999]=t[994]+t[995]+t[996]+t[997];
  t[1000]=-2565.*gdetrd*vconld;
  t[1001]=-9747.*gdetru*vconld;
  t[1002]=13509.*gdetu*vconld;
  t[1003]=-2565.*gdetld*vconlu;
  t[1004]=-675.*gdetlu*vconlu;
  t[1005]=13509.*gdetr*vconlu;
  t[1006]=-9747.*gdetrd*vconlu;
  t[1007]=t[1000]+t[1001]+t[1002];
  t[1008]=t[1003]+t[1004]+t[1005]+t[1006];
  t[1009]=t[980]+t[981]+t[989]+t[990];
  t[1010]=t[998]+t[999]+t[1007]+t[1008];
  t[1011]=-2565.*gdetru*vconlu;
  t[1012]=3555.*gdetu*vconlu;
  t[1013]=25479.*gdetld*vconr;
  t[1014]=13509.*gdetlu*vconr;
  t[1015]=-12990.*gdetr*vconr;
  t[1016]=6705.*gdetrd*vconr;
  t[1017]=3555.*gdetru*vconr;
  t[1018]=t[1011]+t[1012]+t[1013];
  t[1019]=t[1014]+t[1015]+t[1016]+t[1017];
  t[1020]=-35313.*gdetu*vconr;
  t[1021]=-2565.*gdetld*vconrd;
  t[1022]=-9747.*gdetlu*vconrd;
  t[1023]=6705.*gdetr*vconrd;
  t[1024]=-675.*gdetrd*vconrd;
  t[1025]=-2565.*gdetru*vconrd;
  t[1026]=25479.*gdetu*vconrd;
  t[1027]=t[1020]+t[1021]+t[1022];
  t[1028]=t[1023]+t[1024]+t[1025]+t[1026];
  t[1029]=-9747.*gdetld*vconru;
  t[1030]=-2565.*gdetlu*vconru;
  t[1031]=3555.*gdetr*vconru;
  t[1032]=-2565.*gdetrd*vconru;
  t[1033]=-675.*gdetru*vconru;
  t[1034]=6705.*gdetu*vconru;
  t[1035]=68414.*vconc;
  t[1036]=-35313.*vcond;
  t[1037]=-12990.*vconl;
  t[1038]=6705.*vconld;
  t[1039]=3555.*vconlu;
  t[1040]=-49362.*vconr;
  t[1041]=25479.*vconrd;
  t[1042]=13509.*vconru;
  t[1043]=-18723.*vconu;
  t[1044]=t[1039]+t[1040];
  t[1045]=t[1041]+t[1042]+t[1043];
  t[1046]=t[1035]+t[1036]+t[1037];
  t[1047]=t[1038]+t[1044]+t[1045];
  t[1048]=t[1034];
  t[1049]=gdetl*(t[1046]+t[1047]);
  t[1050]=t[1029]+t[1030]+t[1031];
  t[1051]=t[1032]+t[1033]+t[1048]+t[1049];
  t[1052]=13509.*gdetld*vconu;
  t[1053]=3555.*gdetlu*vconu;
  t[1054]=-35313.*gdetr*vconu;
  t[1055]=25479.*gdetrd*vconu;
  t[1056]=6705.*gdetru*vconu;
  t[1057]=-12990.*gdetu*vconu+gdetd*vartemp10585;
  t[1058]=t[1052]+t[1053]+t[1054];
  t[1059]=t[1055]+t[1056]+t[1057];
  t[1060]=t[1018]+t[1019]+t[1027]+t[1028];
  t[1061]=t[1050]+t[1051]+t[1058]+t[1059];
  t[1062]=t[1009]+t[1010]+t[1060]+t[1061];
  t[1063]=gdetc*(t[971]+t[972]);
  t[1064]=-3.*t[1062];
  t[1065]=Bcond*(t[957]+t[958]);
  t[1066]=Bconc*(t[1063]+t[1064]);
  t[1067]=t[839]+t[840];
  t[1068]=t[841]+t[1065]+t[1066];
  t[1069]=t[835]+t[836]+t[837];
  t[1070]=t[838]+t[1067]+t[1068];
  t[1071]=t[807]+t[808]+t[820]+t[821];
  t[1072]=t[833]+t[834]+t[1069]+t[1070];
  t[1073]=t[686]+t[687]+t[740]+t[741];
  t[1074]=t[794]+t[795]+t[1071]+t[1072];
  t[1075]=t[211]+t[212]+t[423]+t[424];
  t[1076]=t[635]+t[636]+t[1073]+t[1074];
  t[1077]=t[1075]+t[1076];
  *Frd=0.000022675736961451247165532879818594*t[1077];



  FTYPE vartemp12517=68414.*vconc;
  FTYPE vartemp12518=4330.*vcond;
  FTYPE vartemp12519=11771.*vconl;
  FTYPE vartemp12520=-2235.*vconld;
  FTYPE vartemp12521=-8493.*vconlu;
  FTYPE vartemp12522=6241.*vconr;
  FTYPE vartemp12523=-1185.*vconrd;
  FTYPE vartemp12524=-4503.*vconru;
  FTYPE vartemp12525=16454.*vconu;
  FTYPE vartemp12526=vartemp12518+vartemp12519+vartemp12520+vartemp12521+vartemp12522+vartemp12523+vartemp12524+vartemp12525;
  FTYPE vartemp12527=-3.*vartemp12526;
  FTYPE vartemp12528=vartemp12517+vartemp12527;
  t[1]=-387102.*Bconl*gdetc*vconc;
  t[2]=105939.*Bconld*gdetc*vconc;
  t[3]=199809.*Bconlu*gdetc*vconc;
  t[4]=-205242.*Bconr*gdetc*vconc;
  t[5]=56169.*Bconrd*gdetc*vconc;
  t[6]=105939.*Bconru*gdetc*vconc;
  t[7]=-387102.*Bconu*gdetc*vconc;
  t[8]=105939.*Bconl*gdetd*vconc;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=-20115.*Bconld*gdetd*vconc;
  t[12]=-76437.*Bconlu*gdetd*vconc;
  t[13]=56169.*Bconr*gdetd*vconc;
  t[14]=-10665.*Bconrd*gdetd*vconc;
  t[15]=-40527.*Bconru*gdetd*vconc;
  t[16]=148086.*Bconu*gdetd*vconc;
  t[17]=38970.*Bconl*gdetl*vconc;
  t[18]=-10665.*Bconld*gdetl*vconc;
  t[19]=-20115.*Bconlu*gdetl*vconc;
  t[20]=t[15]+t[16];
  t[21]=t[17]+t[18]+t[19];
  t[22]=t[11]+t[12]+t[13];
  t[23]=t[14]+t[20]+t[21];
  t[24]=148086.*Bconr*gdetl*vconc;
  t[25]=-40527.*Bconrd*gdetl*vconc;
  t[26]=-76437.*Bconru*gdetl*vconc;
  t[27]=199809.*Bconu*gdetl*vconc;
  t[28]=-10665.*Bconl*gdetld*vconc;
  t[29]=2025.*Bconld*gdetld*vconc;
  t[30]=7695.*Bconlu*gdetld*vconc;
  t[31]=-40527.*Bconr*gdetld*vconc;
  t[32]=7695.*Bconrd*gdetld*vconc;
  t[33]=t[28]+t[29];
  t[34]=t[30]+t[31]+t[32];
  t[35]=t[24]+t[25]+t[26];
  t[36]=t[27]+t[33]+t[34];
  t[37]=29241.*Bconru*gdetld*vconc;
  t[38]=-76437.*Bconu*gdetld*vconc;
  t[39]=-20115.*Bconl*gdetlu*vconc;
  t[40]=7695.*Bconld*gdetlu*vconc;
  t[41]=2025.*Bconlu*gdetlu*vconc;
  t[42]=-76437.*Bconr*gdetlu*vconc;
  t[43]=29241.*Bconrd*gdetlu*vconc;
  t[44]=7695.*Bconru*gdetlu*vconc;
  t[45]=-20115.*Bconu*gdetlu*vconc;
  t[46]=t[41]+t[42];
  t[47]=t[43]+t[44]+t[45];
  t[48]=t[37]+t[38]+t[39];
  t[49]=t[40]+t[46]+t[47];
  t[50]=t[9]+t[10]+t[22]+t[23];
  t[51]=t[35]+t[36]+t[48]+t[49];
  t[52]=148086.*Bconl*gdetr*vconc;
  t[53]=-40527.*Bconld*gdetr*vconc;
  t[54]=-76437.*Bconlu*gdetr*vconc;
  t[55]=38970.*Bconr*gdetr*vconc;
  t[56]=-10665.*Bconrd*gdetr*vconc;
  t[57]=-20115.*Bconru*gdetr*vconc;
  t[58]=105939.*Bconu*gdetr*vconc;
  t[59]=-40527.*Bconl*gdetrd*vconc;
  t[60]=7695.*Bconld*gdetrd*vconc;
  t[61]=t[56]+t[57];
  t[62]=t[58]+t[59]+t[60];
  t[63]=t[52]+t[53]+t[54];
  t[64]=t[55]+t[61]+t[62];
  t[65]=29241.*Bconlu*gdetrd*vconc;
  t[66]=-10665.*Bconr*gdetrd*vconc;
  t[67]=2025.*Bconrd*gdetrd*vconc;
  t[68]=7695.*Bconru*gdetrd*vconc;
  t[69]=-40527.*Bconu*gdetrd*vconc;
  t[70]=-76437.*Bconl*gdetru*vconc;
  t[71]=29241.*Bconld*gdetru*vconc;
  t[72]=7695.*Bconlu*gdetru*vconc;
  t[73]=-20115.*Bconr*gdetru*vconc;
  t[74]=t[69]+t[70];
  t[75]=t[71]+t[72]+t[73];
  t[76]=t[65]+t[66]+t[67];
  t[77]=t[68]+t[74]+t[75];
  t[78]=7695.*Bconrd*gdetru*vconc;
  t[79]=2025.*Bconru*gdetru*vconc;
  t[80]=-10665.*Bconu*gdetru*vconc;
  t[81]=199809.*Bconl*gdetu*vconc;
  t[82]=-76437.*Bconld*gdetu*vconc;
  t[83]=-20115.*Bconlu*gdetu*vconc;
  t[84]=105939.*Bconr*gdetu*vconc;
  t[85]=-40527.*Bconrd*gdetu*vconc;
  t[86]=-10665.*Bconru*gdetu*vconc;
  t[87]=t[82]+t[83];
  t[88]=t[84]+t[85]+t[86];
  t[89]=t[78]+t[79]+t[80];
  t[90]=t[81]+t[87]+t[88];
  t[91]=38970.*Bconu*gdetu*vconc;
  t[92]=105939.*Bconl*gdetc*vcond;
  t[93]=-20115.*Bconld*gdetc*vcond;
  t[94]=-76437.*Bconlu*gdetc*vcond;
  t[95]=56169.*Bconr*gdetc*vcond;
  t[96]=-10665.*Bconrd*gdetc*vcond;
  t[97]=-40527.*Bconru*gdetc*vcond;
  t[98]=148086.*Bconu*gdetc*vcond;
  t[99]=-20115.*Bconl*gdetd*vcond;
  t[100]=t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[91]+t[92]+t[93];
  t[103]=t[94]+t[100]+t[101];
  t[104]=t[63]+t[64]+t[76]+t[77];
  t[105]=t[89]+t[90]+t[102]+t[103];
  t[106]=-894.*Bconld*gdetd*vcond;
  t[107]=24138.*Bconlu*gdetd*vcond;
  t[108]=-10665.*Bconr*gdetd*vcond;
  t[109]=-474.*Bconrd*gdetd*vcond;
  t[110]=12798.*Bconru*gdetd*vcond;
  t[111]=-46764.*Bconu*gdetd*vcond;
  t[112]=-10665.*Bconl*gdetl*vcond;
  t[113]=2025.*Bconld*gdetl*vcond;
  t[114]=t[106]+t[107]+t[108]+t[109];
  t[115]=t[110]+t[111]+t[112]+t[113];
  t[116]=7695.*Bconlu*gdetl*vcond;
  t[117]=-40527.*Bconr*gdetl*vcond;
  t[118]=7695.*Bconrd*gdetl*vcond;
  t[119]=29241.*Bconru*gdetl*vcond;
  t[120]=-76437.*Bconu*gdetl*vcond;
  t[121]=2025.*Bconl*gdetld*vcond;
  t[122]=90.*Bconld*gdetld*vcond;
  t[123]=-2430.*Bconlu*gdetld*vcond;
  t[124]=7695.*Bconr*gdetld*vcond;
  t[125]=t[120]+t[121];
  t[126]=t[122]+t[123]+t[124];
  t[127]=t[116]+t[117]+t[118];
  t[128]=t[119]+t[125]+t[126];
  t[129]=342.*Bconrd*gdetld*vcond;
  t[130]=-9234.*Bconru*gdetld*vcond;
  t[131]=24138.*Bconu*gdetld*vcond;
  t[132]=7695.*Bconl*gdetlu*vcond;
  t[133]=-2430.*Bconld*gdetlu*vcond;
  t[134]=-2430.*Bconlu*gdetlu*vcond;
  t[135]=29241.*Bconr*gdetlu*vcond;
  t[136]=-9234.*Bconrd*gdetlu*vcond;
  t[137]=-9234.*Bconru*gdetlu*vcond;
  t[138]=t[133]+t[134];
  t[139]=t[135]+t[136]+t[137];
  t[140]=t[129]+t[130]+t[131];
  t[141]=t[132]+t[138]+t[139];
  t[142]=24138.*Bconu*gdetlu*vcond;
  t[143]=-40527.*Bconl*gdetr*vcond;
  t[144]=7695.*Bconld*gdetr*vcond;
  t[145]=29241.*Bconlu*gdetr*vcond;
  t[146]=-10665.*Bconr*gdetr*vcond;
  t[147]=2025.*Bconrd*gdetr*vcond;
  t[148]=7695.*Bconru*gdetr*vcond;
  t[149]=-40527.*Bconu*gdetr*vcond;
  t[150]=7695.*Bconl*gdetrd*vcond;
  t[151]=t[146]+t[147];
  t[152]=t[148]+t[149]+t[150];
  t[153]=t[142]+t[143]+t[144];
  t[154]=t[145]+t[151]+t[152];
  t[155]=t[114]+t[115]+t[127]+t[128];
  t[156]=t[140]+t[141]+t[153]+t[154];
  t[157]=342.*Bconld*gdetrd*vcond;
  t[158]=-9234.*Bconlu*gdetrd*vcond;
  t[159]=2025.*Bconr*gdetrd*vcond;
  t[160]=90.*Bconrd*gdetrd*vcond;
  t[161]=-2430.*Bconru*gdetrd*vcond;
  t[162]=12798.*Bconu*gdetrd*vcond;
  t[163]=29241.*Bconl*gdetru*vcond;
  t[164]=-9234.*Bconld*gdetru*vcond;
  t[165]=-9234.*Bconlu*gdetru*vcond;
  t[166]=t[161]+t[162];
  t[167]=t[163]+t[164]+t[165];
  t[168]=t[157]+t[158]+t[159];
  t[169]=t[160]+t[166]+t[167];
  t[170]=7695.*Bconr*gdetru*vcond;
  t[171]=-2430.*Bconrd*gdetru*vcond;
  t[172]=-2430.*Bconru*gdetru*vcond;
  t[173]=12798.*Bconu*gdetru*vcond;
  t[174]=-76437.*Bconl*gdetu*vcond;
  t[175]=24138.*Bconld*gdetu*vcond;
  t[176]=24138.*Bconlu*gdetu*vcond;
  t[177]=-40527.*Bconr*gdetu*vcond;
  t[178]=12798.*Bconrd*gdetu*vcond;
  t[179]=t[174]+t[175];
  t[180]=t[176]+t[177]+t[178];
  t[181]=t[170]+t[171]+t[172];
  t[182]=t[173]+t[179]+t[180];
  t[183]=12798.*Bconru*gdetu*vcond;
  t[184]=-46764.*Bconu*gdetu*vcond;
  t[185]=38970.*Bconl*gdetc*vconl;
  t[186]=-10665.*Bconld*gdetc*vconl;
  t[187]=-20115.*Bconlu*gdetc*vconl;
  t[188]=148086.*Bconr*gdetc*vconl;
  t[189]=-40527.*Bconrd*gdetc*vconl;
  t[190]=-76437.*Bconru*gdetc*vconl;
  t[191]=199809.*Bconu*gdetc*vconl;
  t[192]=t[187]+t[188];
  t[193]=t[189]+t[190]+t[191];
  t[194]=t[183]+t[184]+t[185];
  t[195]=t[186]+t[192]+t[193];
  t[196]=-10665.*Bconl*gdetd*vconl;
  t[197]=2025.*Bconld*gdetd*vconl;
  t[198]=7695.*Bconlu*gdetd*vconl;
  t[199]=-40527.*Bconr*gdetd*vconl;
  t[200]=7695.*Bconrd*gdetd*vconl;
  t[201]=29241.*Bconru*gdetd*vconl;
  t[202]=-76437.*Bconu*gdetd*vconl;
  t[203]=1732.*Bconl*gdetl*vconl;
  t[204]=-474.*Bconld*gdetl*vconl;
  t[205]=t[200]+t[201];
  t[206]=t[202]+t[203]+t[204];
  t[207]=t[196]+t[197]+t[198];
  t[208]=t[199]+t[205]+t[206];
  t[209]=t[168]+t[169]+t[181]+t[182];
  t[210]=t[194]+t[195]+t[207]+t[208];
  t[211]=t[50]+t[51]+t[104]+t[105];
  t[212]=t[155]+t[156]+t[209]+t[210];
  t[213]=-894.*Bconlu*gdetl*vconl;
  t[214]=-46764.*Bconr*gdetl*vconl;
  t[215]=12798.*Bconrd*gdetl*vconl;
  t[216]=24138.*Bconru*gdetl*vconl;
  t[217]=-20115.*Bconu*gdetl*vconl;
  t[218]=-474.*Bconl*gdetld*vconl;
  t[219]=90.*Bconld*gdetld*vconl;
  t[220]=342.*Bconlu*gdetld*vconl;
  t[221]=t[213]+t[214]+t[215]+t[216];
  t[222]=t[217]+t[218]+t[219]+t[220];
  t[223]=12798.*Bconr*gdetld*vconl;
  t[224]=-2430.*Bconrd*gdetld*vconl;
  t[225]=-9234.*Bconru*gdetld*vconl;
  t[226]=7695.*Bconu*gdetld*vconl;
  t[227]=-894.*Bconl*gdetlu*vconl;
  t[228]=342.*Bconld*gdetlu*vconl;
  t[229]=90.*Bconlu*gdetlu*vconl;
  t[230]=24138.*Bconr*gdetlu*vconl;
  t[231]=-9234.*Bconrd*gdetlu*vconl;
  t[232]=t[227]+t[228];
  t[233]=t[229]+t[230]+t[231];
  t[234]=t[223]+t[224]+t[225];
  t[235]=t[226]+t[232]+t[233];
  t[236]=-2430.*Bconru*gdetlu*vconl;
  t[237]=2025.*Bconu*gdetlu*vconl;
  t[238]=-46764.*Bconl*gdetr*vconl;
  t[239]=12798.*Bconld*gdetr*vconl;
  t[240]=24138.*Bconlu*gdetr*vconl;
  t[241]=-46764.*Bconr*gdetr*vconl;
  t[242]=12798.*Bconrd*gdetr*vconl;
  t[243]=24138.*Bconru*gdetr*vconl;
  t[244]=-76437.*Bconu*gdetr*vconl;
  t[245]=t[240]+t[241];
  t[246]=t[242]+t[243]+t[244];
  t[247]=t[236]+t[237]+t[238];
  t[248]=t[239]+t[245]+t[246];
  t[249]=12798.*Bconl*gdetrd*vconl;
  t[250]=-2430.*Bconld*gdetrd*vconl;
  t[251]=-9234.*Bconlu*gdetrd*vconl;
  t[252]=12798.*Bconr*gdetrd*vconl;
  t[253]=-2430.*Bconrd*gdetrd*vconl;
  t[254]=-9234.*Bconru*gdetrd*vconl;
  t[255]=29241.*Bconu*gdetrd*vconl;
  t[256]=24138.*Bconl*gdetru*vconl;
  t[257]=-9234.*Bconld*gdetru*vconl;
  t[258]=t[253]+t[254];
  t[259]=t[255]+t[256]+t[257];
  t[260]=t[249]+t[250]+t[251];
  t[261]=t[252]+t[258]+t[259];
  t[262]=t[221]+t[222]+t[234]+t[235];
  t[263]=t[247]+t[248]+t[260]+t[261];
  t[264]=-2430.*Bconlu*gdetru*vconl;
  t[265]=24138.*Bconr*gdetru*vconl;
  t[266]=-9234.*Bconrd*gdetru*vconl;
  t[267]=-2430.*Bconru*gdetru*vconl;
  t[268]=7695.*Bconu*gdetru*vconl;
  t[269]=-20115.*Bconl*gdetu*vconl;
  t[270]=7695.*Bconld*gdetu*vconl;
  t[271]=2025.*Bconlu*gdetu*vconl;
  t[272]=-76437.*Bconr*gdetu*vconl;
  t[273]=t[268]+t[269];
  t[274]=t[270]+t[271]+t[272];
  t[275]=t[264]+t[265]+t[266];
  t[276]=t[267]+t[273]+t[274];
  t[277]=29241.*Bconrd*gdetu*vconl;
  t[278]=7695.*Bconru*gdetu*vconl;
  t[279]=-20115.*Bconu*gdetu*vconl;
  t[280]=-10665.*Bconl*gdetc*vconld;
  t[281]=2025.*Bconld*gdetc*vconld;
  t[282]=7695.*Bconlu*gdetc*vconld;
  t[283]=-40527.*Bconr*gdetc*vconld;
  t[284]=7695.*Bconrd*gdetc*vconld;
  t[285]=29241.*Bconru*gdetc*vconld;
  t[286]=t[281]+t[282];
  t[287]=t[283]+t[284]+t[285];
  t[288]=t[277]+t[278]+t[279];
  t[289]=t[280]+t[286]+t[287];
  t[290]=-76437.*Bconu*gdetc*vconld;
  t[291]=2025.*Bconl*gdetd*vconld;
  t[292]=90.*Bconld*gdetd*vconld;
  t[293]=-2430.*Bconlu*gdetd*vconld;
  t[294]=7695.*Bconr*gdetd*vconld;
  t[295]=342.*Bconrd*gdetd*vconld;
  t[296]=-9234.*Bconru*gdetd*vconld;
  t[297]=24138.*Bconu*gdetd*vconld;
  t[298]=-474.*Bconl*gdetl*vconld;
  t[299]=t[294]+t[295];
  t[300]=t[296]+t[297]+t[298];
  t[301]=t[290]+t[291]+t[292];
  t[302]=t[293]+t[299]+t[300];
  t[303]=90.*Bconld*gdetl*vconld;
  t[304]=342.*Bconlu*gdetl*vconld;
  t[305]=12798.*Bconr*gdetl*vconld;
  t[306]=-2430.*Bconrd*gdetl*vconld;
  t[307]=-9234.*Bconru*gdetl*vconld;
  t[308]=7695.*Bconu*gdetl*vconld;
  t[309]=90.*Bconl*gdetld*vconld;
  t[310]=4.*Bconld*gdetld*vconld;
  t[311]=-108.*Bconlu*gdetld*vconld;
  t[312]=t[307]+t[308];
  t[313]=t[309]+t[310]+t[311];
  t[314]=t[303]+t[304]+t[305];
  t[315]=t[306]+t[312]+t[313];
  t[316]=t[275]+t[276]+t[288]+t[289];
  t[317]=t[301]+t[302]+t[314]+t[315];
  t[318]=-2430.*Bconr*gdetld*vconld;
  t[319]=-108.*Bconrd*gdetld*vconld;
  t[320]=2916.*Bconru*gdetld*vconld;
  t[321]=-2430.*Bconu*gdetld*vconld;
  t[322]=342.*Bconl*gdetlu*vconld;
  t[323]=-108.*Bconld*gdetlu*vconld;
  t[324]=-108.*Bconlu*gdetlu*vconld;
  t[325]=-9234.*Bconr*gdetlu*vconld;
  t[326]=t[318]+t[319]+t[320]+t[321];
  t[327]=t[322]+t[323]+t[324]+t[325];
  t[328]=2916.*Bconrd*gdetlu*vconld;
  t[329]=2916.*Bconru*gdetlu*vconld;
  t[330]=-2430.*Bconu*gdetlu*vconld;
  t[331]=12798.*Bconl*gdetr*vconld;
  t[332]=-2430.*Bconld*gdetr*vconld;
  t[333]=-9234.*Bconlu*gdetr*vconld;
  t[334]=12798.*Bconr*gdetr*vconld;
  t[335]=-2430.*Bconrd*gdetr*vconld;
  t[336]=-9234.*Bconru*gdetr*vconld;
  t[337]=t[332]+t[333];
  t[338]=t[334]+t[335]+t[336];
  t[339]=t[328]+t[329]+t[330];
  t[340]=t[331]+t[337]+t[338];
  t[341]=29241.*Bconu*gdetr*vconld;
  t[342]=-2430.*Bconl*gdetrd*vconld;
  t[343]=-108.*Bconld*gdetrd*vconld;
  t[344]=2916.*Bconlu*gdetrd*vconld;
  t[345]=-2430.*Bconr*gdetrd*vconld;
  t[346]=-108.*Bconrd*gdetrd*vconld;
  t[347]=2916.*Bconru*gdetrd*vconld;
  t[348]=-9234.*Bconu*gdetrd*vconld;
  t[349]=-9234.*Bconl*gdetru*vconld;
  t[350]=t[345]+t[346];
  t[351]=t[347]+t[348]+t[349];
  t[352]=t[341]+t[342]+t[343];
  t[353]=t[344]+t[350]+t[351];
  t[354]=2916.*Bconld*gdetru*vconld;
  t[355]=2916.*Bconlu*gdetru*vconld;
  t[356]=-9234.*Bconr*gdetru*vconld;
  t[357]=2916.*Bconrd*gdetru*vconld;
  t[358]=2916.*Bconru*gdetru*vconld;
  t[359]=-9234.*Bconu*gdetru*vconld;
  t[360]=7695.*Bconl*gdetu*vconld;
  t[361]=-2430.*Bconld*gdetu*vconld;
  t[362]=-2430.*Bconlu*gdetu*vconld;
  t[363]=t[358]+t[359];
  t[364]=t[360]+t[361]+t[362];
  t[365]=t[354]+t[355]+t[356];
  t[366]=t[357]+t[363]+t[364];
  t[367]=t[326]+t[327]+t[339]+t[340];
  t[368]=t[352]+t[353]+t[365]+t[366];
  t[369]=29241.*Bconr*gdetu*vconld;
  t[370]=-9234.*Bconrd*gdetu*vconld;
  t[371]=-9234.*Bconru*gdetu*vconld;
  t[372]=24138.*Bconu*gdetu*vconld;
  t[373]=-20115.*Bconl*gdetc*vconlu;
  t[374]=7695.*Bconld*gdetc*vconlu;
  t[375]=2025.*Bconlu*gdetc*vconlu;
  t[376]=-76437.*Bconr*gdetc*vconlu;
  t[377]=29241.*Bconrd*gdetc*vconlu;
  t[378]=t[373]+t[374];
  t[379]=t[375]+t[376]+t[377];
  t[380]=t[369]+t[370]+t[371];
  t[381]=t[372]+t[378]+t[379];
  t[382]=7695.*Bconru*gdetc*vconlu;
  t[383]=-20115.*Bconu*gdetc*vconlu;
  t[384]=7695.*Bconl*gdetd*vconlu;
  t[385]=-2430.*Bconld*gdetd*vconlu;
  t[386]=-2430.*Bconlu*gdetd*vconlu;
  t[387]=29241.*Bconr*gdetd*vconlu;
  t[388]=-9234.*Bconrd*gdetd*vconlu;
  t[389]=-9234.*Bconru*gdetd*vconlu;
  t[390]=24138.*Bconu*gdetd*vconlu;
  t[391]=t[386]+t[387];
  t[392]=t[388]+t[389]+t[390];
  t[393]=t[382]+t[383]+t[384];
  t[394]=t[385]+t[391]+t[392];
  t[395]=-894.*Bconl*gdetl*vconlu;
  t[396]=342.*Bconld*gdetl*vconlu;
  t[397]=90.*Bconlu*gdetl*vconlu;
  t[398]=24138.*Bconr*gdetl*vconlu;
  t[399]=-9234.*Bconrd*gdetl*vconlu;
  t[400]=-2430.*Bconru*gdetl*vconlu;
  t[401]=2025.*Bconu*gdetl*vconlu;
  t[402]=342.*Bconl*gdetld*vconlu;
  t[403]=-108.*Bconld*gdetld*vconlu;
  t[404]=t[399]+t[400];
  t[405]=t[401]+t[402]+t[403];
  t[406]=t[395]+t[396]+t[397];
  t[407]=t[398]+t[404]+t[405];
  t[408]=-108.*Bconlu*gdetld*vconlu;
  t[409]=-9234.*Bconr*gdetld*vconlu;
  t[410]=2916.*Bconrd*gdetld*vconlu;
  t[411]=2916.*Bconru*gdetld*vconlu;
  t[412]=-2430.*Bconu*gdetld*vconlu;
  t[413]=90.*Bconl*gdetlu*vconlu;
  t[414]=-108.*Bconld*gdetlu*vconlu;
  t[415]=-44096.*Bconlu*gdetlu*vconlu;
  t[416]=-2430.*Bconr*gdetlu*vconlu;
  t[417]=t[412]+t[413];
  t[418]=t[414]+t[415]+t[416];
  t[419]=t[408]+t[409]+t[410];
  t[420]=t[411]+t[417]+t[418];
  t[421]=t[380]+t[381]+t[393]+t[394];
  t[422]=t[406]+t[407]+t[419]+t[420];
  t[423]=t[262]+t[263]+t[316]+t[317];
  t[424]=t[367]+t[368]+t[421]+t[422];
  t[425]=2916.*Bconrd*gdetlu*vconlu;
  t[426]=-108.*Bconru*gdetlu*vconlu;
  t[427]=90.*Bconu*gdetlu*vconlu;
  t[428]=24138.*Bconl*gdetr*vconlu;
  t[429]=-9234.*Bconld*gdetr*vconlu;
  t[430]=-2430.*Bconlu*gdetr*vconlu;
  t[431]=24138.*Bconr*gdetr*vconlu;
  t[432]=-9234.*Bconrd*gdetr*vconlu;
  t[433]=t[425]+t[426]+t[427]+t[428];
  t[434]=t[429]+t[430]+t[431]+t[432];
  t[435]=-2430.*Bconru*gdetr*vconlu;
  t[436]=7695.*Bconu*gdetr*vconlu;
  t[437]=-9234.*Bconl*gdetrd*vconlu;
  t[438]=2916.*Bconld*gdetrd*vconlu;
  t[439]=2916.*Bconlu*gdetrd*vconlu;
  t[440]=-9234.*Bconr*gdetrd*vconlu;
  t[441]=2916.*Bconrd*gdetrd*vconlu;
  t[442]=2916.*Bconru*gdetrd*vconlu;
  t[443]=-9234.*Bconu*gdetrd*vconlu;
  t[444]=t[439]+t[440];
  t[445]=t[441]+t[442]+t[443];
  t[446]=t[435]+t[436]+t[437];
  t[447]=t[438]+t[444]+t[445];
  t[448]=-2430.*Bconl*gdetru*vconlu;
  t[449]=2916.*Bconld*gdetru*vconlu;
  t[450]=-108.*Bconlu*gdetru*vconlu;
  t[451]=-2430.*Bconr*gdetru*vconlu;
  t[452]=2916.*Bconrd*gdetru*vconlu;
  t[453]=-108.*Bconru*gdetru*vconlu;
  t[454]=342.*Bconu*gdetru*vconlu;
  t[455]=2025.*Bconl*gdetu*vconlu;
  t[456]=-2430.*Bconld*gdetu*vconlu;
  t[457]=t[452]+t[453];
  t[458]=t[454]+t[455]+t[456];
  t[459]=t[448]+t[449]+t[450];
  t[460]=t[451]+t[457]+t[458];
  t[461]=90.*Bconlu*gdetu*vconlu;
  t[462]=7695.*Bconr*gdetu*vconlu;
  t[463]=-9234.*Bconrd*gdetu*vconlu;
  t[464]=342.*Bconru*gdetu*vconlu;
  t[465]=-894.*Bconu*gdetu*vconlu;
  t[466]=148086.*Bconl*gdetc*vconr;
  t[467]=-40527.*Bconld*gdetc*vconr;
  t[468]=-76437.*Bconlu*gdetc*vconr;
  t[469]=38970.*Bconr*gdetc*vconr;
  t[470]=t[465]+t[466];
  t[471]=t[467]+t[468]+t[469];
  t[472]=t[461]+t[462]+t[463];
  t[473]=t[464]+t[470]+t[471];
  t[474]=t[433]+t[434]+t[446]+t[447];
  t[475]=t[459]+t[460]+t[472]+t[473];
  t[476]=-10665.*Bconrd*gdetc*vconr;
  t[477]=-20115.*Bconru*gdetc*vconr;
  t[478]=105939.*Bconu*gdetc*vconr;
  t[479]=-40527.*Bconl*gdetd*vconr;
  t[480]=7695.*Bconld*gdetd*vconr;
  t[481]=29241.*Bconlu*gdetd*vconr;
  t[482]=-10665.*Bconr*gdetd*vconr;
  t[483]=2025.*Bconrd*gdetd*vconr;
  t[484]=7695.*Bconru*gdetd*vconr;
  t[485]=t[480]+t[481];
  t[486]=t[482]+t[483]+t[484];
  t[487]=t[476]+t[477]+t[478];
  t[488]=t[479]+t[485]+t[486];
  t[489]=-40527.*Bconu*gdetd*vconr;
  t[490]=-46764.*Bconl*gdetl*vconr;
  t[491]=12798.*Bconld*gdetl*vconr;
  t[492]=24138.*Bconlu*gdetl*vconr;
  t[493]=-46764.*Bconr*gdetl*vconr;
  t[494]=12798.*Bconrd*gdetl*vconr;
  t[495]=24138.*Bconru*gdetl*vconr;
  t[496]=-76437.*Bconu*gdetl*vconr;
  t[497]=12798.*Bconl*gdetld*vconr;
  t[498]=t[493]+t[494];
  t[499]=t[495]+t[496]+t[497];
  t[500]=t[489]+t[490]+t[491];
  t[501]=t[492]+t[498]+t[499];
  t[502]=-2430.*Bconld*gdetld*vconr;
  t[503]=-9234.*Bconlu*gdetld*vconr;
  t[504]=12798.*Bconr*gdetld*vconr;
  t[505]=-2430.*Bconrd*gdetld*vconr;
  t[506]=-9234.*Bconru*gdetld*vconr;
  t[507]=29241.*Bconu*gdetld*vconr;
  t[508]=24138.*Bconl*gdetlu*vconr;
  t[509]=-9234.*Bconld*gdetlu*vconr;
  t[510]=-2430.*Bconlu*gdetlu*vconr;
  t[511]=t[506]+t[507];
  t[512]=t[508]+t[509]+t[510];
  t[513]=t[502]+t[503]+t[504];
  t[514]=t[505]+t[511]+t[512];
  t[515]=24138.*Bconr*gdetlu*vconr;
  t[516]=-9234.*Bconrd*gdetlu*vconr;
  t[517]=-2430.*Bconru*gdetlu*vconr;
  t[518]=7695.*Bconu*gdetlu*vconr;
  t[519]=-46764.*Bconl*gdetr*vconr;
  t[520]=12798.*Bconld*gdetr*vconr;
  t[521]=24138.*Bconlu*gdetr*vconr;
  t[522]=1732.*Bconr*gdetr*vconr;
  t[523]=-474.*Bconrd*gdetr*vconr;
  t[524]=t[519]+t[520];
  t[525]=t[521]+t[522]+t[523];
  t[526]=t[515]+t[516]+t[517];
  t[527]=t[518]+t[524]+t[525];
  t[528]=t[487]+t[488]+t[500]+t[501];
  t[529]=t[513]+t[514]+t[526]+t[527];
  t[530]=-894.*Bconru*gdetr*vconr;
  t[531]=-20115.*Bconu*gdetr*vconr;
  t[532]=12798.*Bconl*gdetrd*vconr;
  t[533]=-2430.*Bconld*gdetrd*vconr;
  t[534]=-9234.*Bconlu*gdetrd*vconr;
  t[535]=-474.*Bconr*gdetrd*vconr;
  t[536]=90.*Bconrd*gdetrd*vconr;
  t[537]=342.*Bconru*gdetrd*vconr;
  t[538]=t[530]+t[531]+t[532]+t[533];
  t[539]=t[534]+t[535]+t[536]+t[537];
  t[540]=7695.*Bconu*gdetrd*vconr;
  t[541]=24138.*Bconl*gdetru*vconr;
  t[542]=-9234.*Bconld*gdetru*vconr;
  t[543]=-2430.*Bconlu*gdetru*vconr;
  t[544]=-894.*Bconr*gdetru*vconr;
  t[545]=342.*Bconrd*gdetru*vconr;
  t[546]=90.*Bconru*gdetru*vconr;
  t[547]=2025.*Bconu*gdetru*vconr;
  t[548]=-76437.*Bconl*gdetu*vconr;
  t[549]=t[544]+t[545];
  t[550]=t[546]+t[547]+t[548];
  t[551]=t[540]+t[541]+t[542];
  t[552]=t[543]+t[549]+t[550];
  t[553]=29241.*Bconld*gdetu*vconr;
  t[554]=7695.*Bconlu*gdetu*vconr;
  t[555]=-20115.*Bconr*gdetu*vconr;
  t[556]=7695.*Bconrd*gdetu*vconr;
  t[557]=2025.*Bconru*gdetu*vconr;
  t[558]=-10665.*Bconu*gdetu*vconr;
  t[559]=-40527.*Bconl*gdetc*vconrd;
  t[560]=7695.*Bconld*gdetc*vconrd;
  t[561]=29241.*Bconlu*gdetc*vconrd;
  t[562]=t[557]+t[558];
  t[563]=t[559]+t[560]+t[561];
  t[564]=t[553]+t[554]+t[555];
  t[565]=t[556]+t[562]+t[563];
  t[566]=-10665.*Bconr*gdetc*vconrd;
  t[567]=2025.*Bconrd*gdetc*vconrd;
  t[568]=7695.*Bconru*gdetc*vconrd;
  t[569]=-40527.*Bconu*gdetc*vconrd;
  t[570]=7695.*Bconl*gdetd*vconrd;
  t[571]=342.*Bconld*gdetd*vconrd;
  t[572]=-9234.*Bconlu*gdetd*vconrd;
  t[573]=2025.*Bconr*gdetd*vconrd;
  t[574]=90.*Bconrd*gdetd*vconrd;
  t[575]=t[570]+t[571];
  t[576]=t[572]+t[573]+t[574];
  t[577]=t[566]+t[567]+t[568];
  t[578]=t[569]+t[575]+t[576];
  t[579]=t[538]+t[539]+t[551]+t[552];
  t[580]=t[564]+t[565]+t[577]+t[578];
  t[581]=-2430.*Bconru*gdetd*vconrd;
  t[582]=12798.*Bconu*gdetd*vconrd;
  t[583]=12798.*Bconl*gdetl*vconrd;
  t[584]=-2430.*Bconld*gdetl*vconrd;
  t[585]=-9234.*Bconlu*gdetl*vconrd;
  t[586]=12798.*Bconr*gdetl*vconrd;
  t[587]=-2430.*Bconrd*gdetl*vconrd;
  t[588]=-9234.*Bconru*gdetl*vconrd;
  t[589]=29241.*Bconu*gdetl*vconrd;
  t[590]=t[585]+t[586];
  t[591]=t[587]+t[588]+t[589];
  t[592]=t[581]+t[582]+t[583];
  t[593]=t[584]+t[590]+t[591];
  t[594]=-2430.*Bconl*gdetld*vconrd;
  t[595]=-108.*Bconld*gdetld*vconrd;
  t[596]=2916.*Bconlu*gdetld*vconrd;
  t[597]=-2430.*Bconr*gdetld*vconrd;
  t[598]=-108.*Bconrd*gdetld*vconrd;
  t[599]=2916.*Bconru*gdetld*vconrd;
  t[600]=-9234.*Bconu*gdetld*vconrd;
  t[601]=-9234.*Bconl*gdetlu*vconrd;
  t[602]=2916.*Bconld*gdetlu*vconrd;
  t[603]=t[598]+t[599];
  t[604]=t[600]+t[601]+t[602];
  t[605]=t[594]+t[595]+t[596];
  t[606]=t[597]+t[603]+t[604];
  t[607]=2916.*Bconlu*gdetlu*vconrd;
  t[608]=-9234.*Bconr*gdetlu*vconrd;
  t[609]=2916.*Bconrd*gdetlu*vconrd;
  t[610]=2916.*Bconru*gdetlu*vconrd;
  t[611]=-9234.*Bconu*gdetlu*vconrd;
  t[612]=12798.*Bconl*gdetr*vconrd;
  t[613]=-2430.*Bconld*gdetr*vconrd;
  t[614]=-9234.*Bconlu*gdetr*vconrd;
  t[615]=-474.*Bconr*gdetr*vconrd;
  t[616]=t[611]+t[612];
  t[617]=t[613]+t[614]+t[615];
  t[618]=t[607]+t[608]+t[609];
  t[619]=t[610]+t[616]+t[617];
  t[620]=90.*Bconrd*gdetr*vconrd;
  t[621]=342.*Bconru*gdetr*vconrd;
  t[622]=7695.*Bconu*gdetr*vconrd;
  t[623]=-2430.*Bconl*gdetrd*vconrd;
  t[624]=-108.*Bconld*gdetrd*vconrd;
  t[625]=2916.*Bconlu*gdetrd*vconrd;
  t[626]=90.*Bconr*gdetrd*vconrd;
  t[627]=4.*Bconrd*gdetrd*vconrd;
  t[628]=-108.*Bconru*gdetrd*vconrd;
  t[629]=t[624]+t[625];
  t[630]=t[626]+t[627]+t[628];
  t[631]=t[620]+t[621]+t[622];
  t[632]=t[623]+t[629]+t[630];
  t[633]=t[592]+t[593]+t[605]+t[606];
  t[634]=t[618]+t[619]+t[631]+t[632];
  t[635]=t[474]+t[475]+t[528]+t[529];
  t[636]=t[579]+t[580]+t[633]+t[634];
  t[637]=-2430.*Bconu*gdetrd*vconrd;
  t[638]=-9234.*Bconl*gdetru*vconrd;
  t[639]=2916.*Bconld*gdetru*vconrd;
  t[640]=2916.*Bconlu*gdetru*vconrd;
  t[641]=342.*Bconr*gdetru*vconrd;
  t[642]=-108.*Bconrd*gdetru*vconrd;
  t[643]=-108.*Bconru*gdetru*vconrd;
  t[644]=-2430.*Bconu*gdetru*vconrd;
  t[645]=t[637]+t[638]+t[639]+t[640];
  t[646]=t[641]+t[642]+t[643]+t[644];
  t[647]=29241.*Bconl*gdetu*vconrd;
  t[648]=-9234.*Bconld*gdetu*vconrd;
  t[649]=-9234.*Bconlu*gdetu*vconrd;
  t[650]=7695.*Bconr*gdetu*vconrd;
  t[651]=-2430.*Bconrd*gdetu*vconrd;
  t[652]=-2430.*Bconru*gdetu*vconrd;
  t[653]=12798.*Bconu*gdetu*vconrd;
  t[654]=-76437.*Bconl*gdetc*vconru;
  t[655]=29241.*Bconld*gdetc*vconru;
  t[656]=t[651]+t[652];
  t[657]=t[653]+t[654]+t[655];
  t[658]=t[647]+t[648]+t[649];
  t[659]=t[650]+t[656]+t[657];
  t[660]=7695.*Bconlu*gdetc*vconru;
  t[661]=-20115.*Bconr*gdetc*vconru;
  t[662]=7695.*Bconrd*gdetc*vconru;
  t[663]=2025.*Bconru*gdetc*vconru;
  t[664]=-10665.*Bconu*gdetc*vconru;
  t[665]=29241.*Bconl*gdetd*vconru;
  t[666]=-9234.*Bconld*gdetd*vconru;
  t[667]=-9234.*Bconlu*gdetd*vconru;
  t[668]=7695.*Bconr*gdetd*vconru;
  t[669]=t[664]+t[665];
  t[670]=t[666]+t[667]+t[668];
  t[671]=t[660]+t[661]+t[662];
  t[672]=t[663]+t[669]+t[670];
  t[673]=-2430.*Bconrd*gdetd*vconru;
  t[674]=-2430.*Bconru*gdetd*vconru;
  t[675]=12798.*Bconu*gdetd*vconru;
  t[676]=24138.*Bconl*gdetl*vconru;
  t[677]=-9234.*Bconld*gdetl*vconru;
  t[678]=-2430.*Bconlu*gdetl*vconru;
  t[679]=24138.*Bconr*gdetl*vconru;
  t[680]=-9234.*Bconrd*gdetl*vconru;
  t[681]=-2430.*Bconru*gdetl*vconru;
  t[682]=t[677]+t[678];
  t[683]=t[679]+t[680]+t[681];
  t[684]=t[673]+t[674]+t[675];
  t[685]=t[676]+t[682]+t[683];
  t[686]=t[645]+t[646]+t[658]+t[659];
  t[687]=t[671]+t[672]+t[684]+t[685];
  t[688]=7695.*Bconu*gdetl*vconru;
  t[689]=-9234.*Bconl*gdetld*vconru;
  t[690]=2916.*Bconld*gdetld*vconru;
  t[691]=2916.*Bconlu*gdetld*vconru;
  t[692]=-9234.*Bconr*gdetld*vconru;
  t[693]=2916.*Bconrd*gdetld*vconru;
  t[694]=2916.*Bconru*gdetld*vconru;
  t[695]=-9234.*Bconu*gdetld*vconru;
  t[696]=-2430.*Bconl*gdetlu*vconru;
  t[697]=t[692]+t[693];
  t[698]=t[694]+t[695]+t[696];
  t[699]=t[688]+t[689]+t[690];
  t[700]=t[691]+t[697]+t[698];
  t[701]=2916.*Bconld*gdetlu*vconru;
  t[702]=-108.*Bconlu*gdetlu*vconru;
  t[703]=-2430.*Bconr*gdetlu*vconru;
  t[704]=2916.*Bconrd*gdetlu*vconru;
  t[705]=-108.*Bconru*gdetlu*vconru;
  t[706]=342.*Bconu*gdetlu*vconru;
  t[707]=24138.*Bconl*gdetr*vconru;
  t[708]=-9234.*Bconld*gdetr*vconru;
  t[709]=-2430.*Bconlu*gdetr*vconru;
  t[710]=t[705]+t[706];
  t[711]=t[707]+t[708]+t[709];
  t[712]=t[701]+t[702]+t[703];
  t[713]=t[704]+t[710]+t[711];
  t[714]=-894.*Bconr*gdetr*vconru;
  t[715]=342.*Bconrd*gdetr*vconru;
  t[716]=90.*Bconru*gdetr*vconru;
  t[717]=2025.*Bconu*gdetr*vconru;
  t[718]=-9234.*Bconl*gdetrd*vconru;
  t[719]=2916.*Bconld*gdetrd*vconru;
  t[720]=2916.*Bconlu*gdetrd*vconru;
  t[721]=342.*Bconr*gdetrd*vconru;
  t[722]=-108.*Bconrd*gdetrd*vconru;
  t[723]=t[718]+t[719];
  t[724]=t[720]+t[721]+t[722];
  t[725]=t[714]+t[715]+t[716];
  t[726]=t[717]+t[723]+t[724];
  t[727]=-108.*Bconru*gdetrd*vconru;
  t[728]=-2430.*Bconu*gdetrd*vconru;
  t[729]=-2430.*Bconl*gdetru*vconru;
  t[730]=2916.*Bconld*gdetru*vconru;
  t[731]=-108.*Bconlu*gdetru*vconru;
  t[732]=90.*Bconr*gdetru*vconru;
  t[733]=-108.*Bconrd*gdetru*vconru;
  t[734]=4.*Bconru*gdetru*vconru;
  t[735]=90.*Bconu*gdetru*vconru;
  t[736]=t[731]+t[732];
  t[737]=t[733]+t[734]+t[735];
  t[738]=t[727]+t[728]+t[729];
  t[739]=t[730]+t[736]+t[737];
  t[740]=t[699]+t[700]+t[712]+t[713];
  t[741]=t[725]+t[726]+t[738]+t[739];
  t[742]=7695.*Bconl*gdetu*vconru;
  t[743]=-9234.*Bconld*gdetu*vconru;
  t[744]=342.*Bconlu*gdetu*vconru;
  t[745]=2025.*Bconr*gdetu*vconru;
  t[746]=-2430.*Bconrd*gdetu*vconru;
  t[747]=90.*Bconru*gdetu*vconru;
  t[748]=-474.*Bconu*gdetu*vconru;
  t[749]=199809.*Bconl*gdetc*vconu;
  t[750]=-76437.*Bconld*gdetc*vconu;
  t[751]=t[746]+t[747];
  t[752]=t[748]+t[749]+t[750];
  t[753]=t[742]+t[743]+t[744];
  t[754]=t[745]+t[751]+t[752];
  t[755]=-20115.*Bconlu*gdetc*vconu;
  t[756]=105939.*Bconr*gdetc*vconu;
  t[757]=-40527.*Bconrd*gdetc*vconu;
  t[758]=-10665.*Bconru*gdetc*vconu;
  t[759]=38970.*Bconu*gdetc*vconu;
  t[760]=-76437.*Bconl*gdetd*vconu;
  t[761]=24138.*Bconld*gdetd*vconu;
  t[762]=24138.*Bconlu*gdetd*vconu;
  t[763]=-40527.*Bconr*gdetd*vconu;
  t[764]=t[759]+t[760];
  t[765]=t[761]+t[762]+t[763];
  t[766]=t[755]+t[756]+t[757];
  t[767]=t[758]+t[764]+t[765];
  t[768]=12798.*Bconrd*gdetd*vconu;
  t[769]=12798.*Bconru*gdetd*vconu;
  t[770]=-46764.*Bconu*gdetd*vconu;
  t[771]=-20115.*Bconl*gdetl*vconu;
  t[772]=7695.*Bconld*gdetl*vconu;
  t[773]=2025.*Bconlu*gdetl*vconu;
  t[774]=-76437.*Bconr*gdetl*vconu;
  t[775]=29241.*Bconrd*gdetl*vconu;
  t[776]=7695.*Bconru*gdetl*vconu;
  t[777]=t[772]+t[773];
  t[778]=t[774]+t[775]+t[776];
  t[779]=t[768]+t[769]+t[770];
  t[780]=t[771]+t[777]+t[778];
  t[781]=-20115.*Bconu*gdetl*vconu;
  t[782]=7695.*Bconl*gdetld*vconu;
  t[783]=-2430.*Bconld*gdetld*vconu;
  t[784]=-2430.*Bconlu*gdetld*vconu;
  t[785]=29241.*Bconr*gdetld*vconu;
  t[786]=-9234.*Bconrd*gdetld*vconu;
  t[787]=-9234.*Bconru*gdetld*vconu;
  t[788]=24138.*Bconu*gdetld*vconu;
  t[789]=2025.*Bconl*gdetlu*vconu;
  t[790]=t[785]+t[786];
  t[791]=t[787]+t[788]+t[789];
  t[792]=t[781]+t[782]+t[783];
  t[793]=t[784]+t[790]+t[791];
  t[794]=t[753]+t[754]+t[766]+t[767];
  t[795]=t[779]+t[780]+t[792]+t[793];
  t[796]=-2430.*Bconld*gdetlu*vconu;
  t[797]=90.*Bconlu*gdetlu*vconu;
  t[798]=7695.*Bconr*gdetlu*vconu;
  t[799]=-9234.*Bconrd*gdetlu*vconu;
  t[800]=342.*Bconru*gdetlu*vconu;
  t[801]=-894.*Bconu*gdetlu*vconu;
  t[802]=-76437.*Bconl*gdetr*vconu;
  t[803]=29241.*Bconld*gdetr*vconu;
  t[804]=7695.*Bconlu*gdetr*vconu;
  t[805]=t[800]+t[801];
  t[806]=t[802]+t[803]+t[804];
  t[807]=t[796]+t[797]+t[798];
  t[808]=t[799]+t[805]+t[806];
  t[809]=-20115.*Bconr*gdetr*vconu;
  t[810]=7695.*Bconrd*gdetr*vconu;
  t[811]=2025.*Bconru*gdetr*vconu;
  t[812]=-10665.*Bconu*gdetr*vconu;
  t[813]=29241.*Bconl*gdetrd*vconu;
  t[814]=-9234.*Bconld*gdetrd*vconu;
  t[815]=-9234.*Bconlu*gdetrd*vconu;
  t[816]=7695.*Bconr*gdetrd*vconu;
  t[817]=-2430.*Bconrd*gdetrd*vconu;
  t[818]=t[813]+t[814];
  t[819]=t[815]+t[816]+t[817];
  t[820]=t[809]+t[810]+t[811];
  t[821]=t[812]+t[818]+t[819];
  t[822]=-2430.*Bconru*gdetrd*vconu;
  t[823]=12798.*Bconu*gdetrd*vconu;
  t[824]=7695.*Bconl*gdetru*vconu;
  t[825]=-9234.*Bconld*gdetru*vconu;
  t[826]=342.*Bconlu*gdetru*vconu;
  t[827]=2025.*Bconr*gdetru*vconu;
  t[828]=-2430.*Bconrd*gdetru*vconu;
  t[829]=90.*Bconru*gdetru*vconu;
  t[830]=-474.*Bconu*gdetru*vconu;
  t[831]=t[826]+t[827];
  t[832]=t[828]+t[829]+t[830];
  t[833]=t[822]+t[823]+t[824];
  t[834]=t[825]+t[831]+t[832];
  t[835]=-20115.*Bconl*gdetu*vconu;
  t[836]=24138.*Bconld*gdetu*vconu;
  t[837]=-894.*Bconlu*gdetu*vconu;
  t[838]=-10665.*Bconr*gdetu*vconu;
  t[839]=12798.*Bconrd*gdetu*vconu;
  t[840]=-474.*Bconru*gdetu*vconu;
  t[841]=1732.*Bconu*gdetu*vconu;
  t[842]=-25479.*gdetlu*vconc;
  t[843]=18723.*gdetr*vconc;
  t[844]=-3555.*gdetrd*vconc;
  t[845]=-13509.*gdetru*vconc;
  t[846]=49362.*gdetu*vconc;
  t[847]=t[842]+t[843];
  t[848]=t[844]+t[845]+t[846];
  t[849]=8046.*gdetlu*vcond;
  t[850]=-3555.*gdetr*vcond;
  t[851]=-158.*gdetrd*vcond;
  t[852]=4266.*gdetru*vcond;
  t[853]=-15588.*gdetu*vcond;
  t[854]=2565.*gdetlu*vconl;
  t[855]=t[849]+t[850]+t[851];
  t[856]=t[852]+t[853]+t[854];
  t[857]=-13509.*gdetr*vconl;
  t[858]=2565.*gdetrd*vconl;
  t[859]=9747.*gdetru*vconl;
  t[860]=-25479.*gdetu*vconl;
  t[861]=-810.*gdetlu*vconld;
  t[862]=2565.*gdetr*vconld;
  t[863]=t[857]+t[858]+t[859];
  t[864]=t[860]+t[861]+t[862];
  t[865]=114.*gdetrd*vconld;
  t[866]=-3078.*gdetru*vconld;
  t[867]=8046.*gdetu*vconld;
  t[868]=-810.*gdetlu*vconlu;
  t[869]=9747.*gdetr*vconlu;
  t[870]=-3078.*gdetrd*vconlu;
  t[871]=t[865]+t[866]+t[867];
  t[872]=t[868]+t[869]+t[870];
  t[873]=t[847]+t[848]+t[855]+t[856];
  t[874]=t[863]+t[864]+t[871]+t[872];
  t[875]=-3078.*gdetru*vconlu;
  t[876]=8046.*gdetu*vconlu;
  t[877]=9747.*gdetlu*vconr;
  t[878]=-3555.*gdetr*vconr;
  t[879]=675.*gdetrd*vconr;
  t[880]=2565.*gdetru*vconr;
  t[881]=t[875]+t[876]+t[877];
  t[882]=t[878]+t[879]+t[880];
  t[883]=-13509.*gdetu*vconr;
  t[884]=-3078.*gdetlu*vconrd;
  t[885]=675.*gdetr*vconrd;
  t[886]=30.*gdetrd*vconrd;
  t[887]=-810.*gdetru*vconrd;
  t[888]=4266.*gdetu*vconrd;
  t[889]=t[883]+t[884]+t[885];
  t[890]=t[886]+t[887]+t[888];
  t[891]=-3078.*gdetlu*vconru;
  t[892]=2565.*gdetr*vconru;
  t[893]=-810.*gdetrd*vconru;
  t[894]=-810.*gdetru*vconru;
  t[895]=4266.*gdetu*vconru;
  t[896]=11771.*vconc;
  t[897]=-745.*vcond;
  t[898]=-395.*vconl;
  t[899]=75.*vconld;
  t[900]=285.*vconlu;
  t[901]=-1501.*vconr;
  t[902]=285.*vconrd;
  t[903]=1083.*vconru;
  t[904]=-2831.*vconu;
  t[905]=t[897]+t[898]+t[899]+t[900];
  t[906]=t[901]+t[902]+t[903]+t[904];
  t[907]=t[905]+t[906];
  t[908]=t[896];
  t[909]=3.*t[907];
  t[910]=gdetl*(t[908]+t[909]);
  t[911]=t[895];
  t[912]=3.*t[910];
  t[913]=t[891]+t[892]+t[893];
  t[914]=t[894]+t[911]+t[912];
  t[915]=8046.*gdetlu*vconu;
  t[916]=-13509.*gdetr*vconu;
  t[917]=4266.*gdetrd*vconu;
  t[918]=4266.*gdetru*vconu;
  t[919]=-15588.*gdetu*vconu;
  t[920]=-6705.*vconc;
  t[921]=-298.*vcond;
  t[922]=675.*vconl;
  t[923]=30.*vconld;
  t[924]=-810.*vconlu;
  t[925]=2565.*vconr;
  t[926]=114.*vconrd;
  t[927]=-3078.*vconru;
  t[928]=8046.*vconu;
  t[929]=t[924]+t[925];
  t[930]=t[926]+t[927]+t[928];
  t[931]=t[920]+t[921]+t[922];
  t[932]=t[923]+t[929]+t[930];
  t[933]=t[919];
  t[934]=gdetld*(t[931]+t[932]);
  t[935]=t[915]+t[916]+t[917];
  t[936]=t[918]+t[933]+t[934];
  t[937]=t[881]+t[882]+t[889]+t[890];
  t[938]=t[913]+t[914]+t[935]+t[936];
  t[939]=t[873]+t[874]+t[937]+t[938];
  t[940]=38970.*vconc;
  t[941]=1732.*vcond;
  t[942]=6705.*vconl;
  t[943]=298.*vconld;
  t[944]=-8046.*vconlu;
  t[945]=3555.*vconr;
  t[946]=158.*vconrd;
  t[947]=-4266.*vconru;
  t[948]=15588.*vconu;
  t[949]=t[942]+t[943]+t[944];
  t[950]=t[945]+t[946]+t[947]+t[948];
  t[951]=t[949]+t[950];
  t[952]=t[941];
  t[953]=-3.*t[951];
  t[954]=-3.*gdetc*vartemp12528;
  t[955]=gdetd*(t[940]+t[952]+t[953]);
  t[956]=t[954];
  t[957]=3.*t[939];
  t[958]=t[955]+t[956];
  t[959]=749956.*vconc;
  t[960]=68414.*vcond;
  t[961]=129034.*vconl;
  t[962]=-35313.*vconld;
  t[963]=-66603.*vconlu;
  t[964]=68414.*vconr;
  t[965]=-18723.*vconrd;
  t[966]=-35313.*vconru;
  t[967]=129034.*vconu;
  t[968]=t[960]+t[961]+t[962]+t[963];
  t[969]=t[964]+t[965]+t[966]+t[967];
  t[970]=t[968]+t[969];
  t[971]=t[959];
  t[972]=-3.*t[970];
  t[973]=-35313.*gdetld*vconc;
  t[974]=-66603.*gdetlu*vconc;
  t[975]=68414.*gdetr*vconc;
  t[976]=-18723.*gdetrd*vconc;
  t[977]=-35313.*gdetru*vconc;
  t[978]=129034.*gdetu*vconc;
  t[979]=6705.*gdetld*vcond;
  t[980]=t[973]+t[974]+t[975];
  t[981]=t[976]+t[977]+t[978]+t[979];
  t[982]=25479.*gdetlu*vcond;
  t[983]=-18723.*gdetr*vcond;
  t[984]=3555.*gdetrd*vcond;
  t[985]=13509.*gdetru*vcond;
  t[986]=-49362.*gdetu*vcond;
  t[987]=3555.*gdetld*vconl;
  t[988]=6705.*gdetlu*vconl;
  t[989]=t[982]+t[983]+t[984];
  t[990]=t[985]+t[986]+t[987]+t[988];
  t[991]=-49362.*gdetr*vconl;
  t[992]=13509.*gdetrd*vconl;
  t[993]=25479.*gdetru*vconl;
  t[994]=-66603.*gdetu*vconl;
  t[995]=-675.*gdetld*vconld;
  t[996]=-2565.*gdetlu*vconld;
  t[997]=13509.*gdetr*vconld;
  t[998]=t[991]+t[992]+t[993];
  t[999]=t[994]+t[995]+t[996]+t[997];
  t[1000]=-2565.*gdetrd*vconld;
  t[1001]=-9747.*gdetru*vconld;
  t[1002]=25479.*gdetu*vconld;
  t[1003]=-2565.*gdetld*vconlu;
  t[1004]=-675.*gdetlu*vconlu;
  t[1005]=25479.*gdetr*vconlu;
  t[1006]=-9747.*gdetrd*vconlu;
  t[1007]=t[1000]+t[1001]+t[1002];
  t[1008]=t[1003]+t[1004]+t[1005]+t[1006];
  t[1009]=t[980]+t[981]+t[989]+t[990];
  t[1010]=t[998]+t[999]+t[1007]+t[1008];
  t[1011]=-2565.*gdetru*vconlu;
  t[1012]=6705.*gdetu*vconlu;
  t[1013]=13509.*gdetld*vconr;
  t[1014]=25479.*gdetlu*vconr;
  t[1015]=-12990.*gdetr*vconr;
  t[1016]=3555.*gdetrd*vconr;
  t[1017]=6705.*gdetru*vconr;
  t[1018]=t[1011]+t[1012]+t[1013];
  t[1019]=t[1014]+t[1015]+t[1016]+t[1017];
  t[1020]=-35313.*gdetu*vconr;
  t[1021]=-2565.*gdetld*vconrd;
  t[1022]=-9747.*gdetlu*vconrd;
  t[1023]=3555.*gdetr*vconrd;
  t[1024]=-675.*gdetrd*vconrd;
  t[1025]=-2565.*gdetru*vconrd;
  t[1026]=13509.*gdetu*vconrd;
  t[1027]=t[1020]+t[1021]+t[1022];
  t[1028]=t[1023]+t[1024]+t[1025]+t[1026];
  t[1029]=-9747.*gdetld*vconru;
  t[1030]=-2565.*gdetlu*vconru;
  t[1031]=6705.*gdetr*vconru;
  t[1032]=-2565.*gdetrd*vconru;
  t[1033]=-675.*gdetru*vconru;
  t[1034]=3555.*gdetu*vconru;
  t[1035]=129034.*vconc;
  t[1036]=-35313.*vcond;
  t[1037]=-12990.*vconl;
  t[1038]=3555.*vconld;
  t[1039]=6705.*vconlu;
  t[1040]=-49362.*vconr;
  t[1041]=13509.*vconrd;
  t[1042]=25479.*vconru;
  t[1043]=-66603.*vconu;
  t[1044]=t[1039]+t[1040];
  t[1045]=t[1041]+t[1042]+t[1043];
  t[1046]=t[1035]+t[1036]+t[1037];
  t[1047]=t[1038]+t[1044]+t[1045];
  t[1048]=t[1034];
  t[1049]=gdetl*(t[1046]+t[1047]);
  t[1050]=t[1029]+t[1030]+t[1031];
  t[1051]=t[1032]+t[1033]+t[1048]+t[1049];
  t[1052]=25479.*gdetld*vconu;
  t[1053]=6705.*gdetlu*vconu;
  t[1054]=-35313.*gdetr*vconu;
  t[1055]=13509.*gdetrd*vconu;
  t[1056]=3555.*gdetru*vconu;
  t[1057]=-12990.*gdetu*vconu+gdetd*vartemp12528;
  t[1058]=t[1052]+t[1053]+t[1054];
  t[1059]=t[1055]+t[1056]+t[1057];
  t[1060]=t[1018]+t[1019]+t[1027]+t[1028];
  t[1061]=t[1050]+t[1051]+t[1058]+t[1059];
  t[1062]=t[1009]+t[1010]+t[1060]+t[1061];
  t[1063]=gdetc*(t[971]+t[972]);
  t[1064]=-3.*t[1062];
  t[1065]=Bcond*(t[957]+t[958]);
  t[1066]=Bconc*(t[1063]+t[1064]);
  t[1067]=t[839]+t[840];
  t[1068]=t[841]+t[1065]+t[1066];
  t[1069]=t[835]+t[836]+t[837];
  t[1070]=t[838]+t[1067]+t[1068];
  t[1071]=t[807]+t[808]+t[820]+t[821];
  t[1072]=t[833]+t[834]+t[1069]+t[1070];
  t[1073]=t[686]+t[687]+t[740]+t[741];
  t[1074]=t[794]+t[795]+t[1071]+t[1072];
  t[1075]=t[211]+t[212]+t[423]+t[424];
  t[1076]=t[635]+t[636]+t[1073]+t[1074];
  t[1077]=t[1075]+t[1076];
  *Flu=0.000022675736961451247165532879818594*t[1077];









  FTYPE vartemp14376=68414.*vconc;
  FTYPE vartemp14377=4330.*vcond;
  FTYPE vartemp14378=6241.*vconl;
  FTYPE vartemp14379=-1185.*vconld;
  FTYPE vartemp14380=-4503.*vconlu;
  FTYPE vartemp14381=11771.*vconr;
  FTYPE vartemp14382=-2235.*vconrd;
  FTYPE vartemp14383=-8493.*vconru;
  FTYPE vartemp14384=16454.*vconu;
  FTYPE vartemp14385=vartemp14377+vartemp14378+vartemp14379+vartemp14380+vartemp14381+vartemp14382+vartemp14383+vartemp14384;
  FTYPE vartemp14386=-3.*vartemp14385;
  FTYPE vartemp14387=vartemp14376+vartemp14386;
  t[1]=-205242.*Bconl*gdetc*vconc;
  t[2]=56169.*Bconld*gdetc*vconc;
  t[3]=105939.*Bconlu*gdetc*vconc;
  t[4]=-387102.*Bconr*gdetc*vconc;
  t[5]=105939.*Bconrd*gdetc*vconc;
  t[6]=199809.*Bconru*gdetc*vconc;
  t[7]=-387102.*Bconu*gdetc*vconc;
  t[8]=56169.*Bconl*gdetd*vconc;
  t[9]=t[1]+t[2]+t[3]+t[4];
  t[10]=t[5]+t[6]+t[7]+t[8];
  t[11]=-10665.*Bconld*gdetd*vconc;
  t[12]=-40527.*Bconlu*gdetd*vconc;
  t[13]=105939.*Bconr*gdetd*vconc;
  t[14]=-20115.*Bconrd*gdetd*vconc;
  t[15]=-76437.*Bconru*gdetd*vconc;
  t[16]=148086.*Bconu*gdetd*vconc;
  t[17]=38970.*Bconl*gdetl*vconc;
  t[18]=-10665.*Bconld*gdetl*vconc;
  t[19]=-20115.*Bconlu*gdetl*vconc;
  t[20]=t[15]+t[16];
  t[21]=t[17]+t[18]+t[19];
  t[22]=t[11]+t[12]+t[13];
  t[23]=t[14]+t[20]+t[21];
  t[24]=148086.*Bconr*gdetl*vconc;
  t[25]=-40527.*Bconrd*gdetl*vconc;
  t[26]=-76437.*Bconru*gdetl*vconc;
  t[27]=105939.*Bconu*gdetl*vconc;
  t[28]=-10665.*Bconl*gdetld*vconc;
  t[29]=2025.*Bconld*gdetld*vconc;
  t[30]=7695.*Bconlu*gdetld*vconc;
  t[31]=-40527.*Bconr*gdetld*vconc;
  t[32]=7695.*Bconrd*gdetld*vconc;
  t[33]=t[28]+t[29];
  t[34]=t[30]+t[31]+t[32];
  t[35]=t[24]+t[25]+t[26];
  t[36]=t[27]+t[33]+t[34];
  t[37]=29241.*Bconru*gdetld*vconc;
  t[38]=-40527.*Bconu*gdetld*vconc;
  t[39]=-20115.*Bconl*gdetlu*vconc;
  t[40]=7695.*Bconld*gdetlu*vconc;
  t[41]=2025.*Bconlu*gdetlu*vconc;
  t[42]=-76437.*Bconr*gdetlu*vconc;
  t[43]=29241.*Bconrd*gdetlu*vconc;
  t[44]=7695.*Bconru*gdetlu*vconc;
  t[45]=-10665.*Bconu*gdetlu*vconc;
  t[46]=t[41]+t[42];
  t[47]=t[43]+t[44]+t[45];
  t[48]=t[37]+t[38]+t[39];
  t[49]=t[40]+t[46]+t[47];
  t[50]=t[9]+t[10]+t[22]+t[23];
  t[51]=t[35]+t[36]+t[48]+t[49];
  t[52]=148086.*Bconl*gdetr*vconc;
  t[53]=-40527.*Bconld*gdetr*vconc;
  t[54]=-76437.*Bconlu*gdetr*vconc;
  t[55]=38970.*Bconr*gdetr*vconc;
  t[56]=-10665.*Bconrd*gdetr*vconc;
  t[57]=-20115.*Bconru*gdetr*vconc;
  t[58]=199809.*Bconu*gdetr*vconc;
  t[59]=-40527.*Bconl*gdetrd*vconc;
  t[60]=7695.*Bconld*gdetrd*vconc;
  t[61]=t[56]+t[57];
  t[62]=t[58]+t[59]+t[60];
  t[63]=t[52]+t[53]+t[54];
  t[64]=t[55]+t[61]+t[62];
  t[65]=29241.*Bconlu*gdetrd*vconc;
  t[66]=-10665.*Bconr*gdetrd*vconc;
  t[67]=2025.*Bconrd*gdetrd*vconc;
  t[68]=7695.*Bconru*gdetrd*vconc;
  t[69]=-76437.*Bconu*gdetrd*vconc;
  t[70]=-76437.*Bconl*gdetru*vconc;
  t[71]=29241.*Bconld*gdetru*vconc;
  t[72]=7695.*Bconlu*gdetru*vconc;
  t[73]=-20115.*Bconr*gdetru*vconc;
  t[74]=t[69]+t[70];
  t[75]=t[71]+t[72]+t[73];
  t[76]=t[65]+t[66]+t[67];
  t[77]=t[68]+t[74]+t[75];
  t[78]=7695.*Bconrd*gdetru*vconc;
  t[79]=2025.*Bconru*gdetru*vconc;
  t[80]=-20115.*Bconu*gdetru*vconc;
  t[81]=105939.*Bconl*gdetu*vconc;
  t[82]=-40527.*Bconld*gdetu*vconc;
  t[83]=-10665.*Bconlu*gdetu*vconc;
  t[84]=199809.*Bconr*gdetu*vconc;
  t[85]=-76437.*Bconrd*gdetu*vconc;
  t[86]=-20115.*Bconru*gdetu*vconc;
  t[87]=t[82]+t[83];
  t[88]=t[84]+t[85]+t[86];
  t[89]=t[78]+t[79]+t[80];
  t[90]=t[81]+t[87]+t[88];
  t[91]=38970.*Bconu*gdetu*vconc;
  t[92]=56169.*Bconl*gdetc*vcond;
  t[93]=-10665.*Bconld*gdetc*vcond;
  t[94]=-40527.*Bconlu*gdetc*vcond;
  t[95]=105939.*Bconr*gdetc*vcond;
  t[96]=-20115.*Bconrd*gdetc*vcond;
  t[97]=-76437.*Bconru*gdetc*vcond;
  t[98]=148086.*Bconu*gdetc*vcond;
  t[99]=-10665.*Bconl*gdetd*vcond;
  t[100]=t[95]+t[96];
  t[101]=t[97]+t[98]+t[99];
  t[102]=t[91]+t[92]+t[93];
  t[103]=t[94]+t[100]+t[101];
  t[104]=t[63]+t[64]+t[76]+t[77];
  t[105]=t[89]+t[90]+t[102]+t[103];
  t[106]=-474.*Bconld*gdetd*vcond;
  t[107]=12798.*Bconlu*gdetd*vcond;
  t[108]=-20115.*Bconr*gdetd*vcond;
  t[109]=-894.*Bconrd*gdetd*vcond;
  t[110]=24138.*Bconru*gdetd*vcond;
  t[111]=-46764.*Bconu*gdetd*vcond;
  t[112]=-10665.*Bconl*gdetl*vcond;
  t[113]=2025.*Bconld*gdetl*vcond;
  t[114]=t[106]+t[107]+t[108]+t[109];
  t[115]=t[110]+t[111]+t[112]+t[113];
  t[116]=7695.*Bconlu*gdetl*vcond;
  t[117]=-40527.*Bconr*gdetl*vcond;
  t[118]=7695.*Bconrd*gdetl*vcond;
  t[119]=29241.*Bconru*gdetl*vcond;
  t[120]=-40527.*Bconu*gdetl*vcond;
  t[121]=2025.*Bconl*gdetld*vcond;
  t[122]=90.*Bconld*gdetld*vcond;
  t[123]=-2430.*Bconlu*gdetld*vcond;
  t[124]=7695.*Bconr*gdetld*vcond;
  t[125]=t[120]+t[121];
  t[126]=t[122]+t[123]+t[124];
  t[127]=t[116]+t[117]+t[118];
  t[128]=t[119]+t[125]+t[126];
  t[129]=342.*Bconrd*gdetld*vcond;
  t[130]=-9234.*Bconru*gdetld*vcond;
  t[131]=12798.*Bconu*gdetld*vcond;
  t[132]=7695.*Bconl*gdetlu*vcond;
  t[133]=-2430.*Bconld*gdetlu*vcond;
  t[134]=-2430.*Bconlu*gdetlu*vcond;
  t[135]=29241.*Bconr*gdetlu*vcond;
  t[136]=-9234.*Bconrd*gdetlu*vcond;
  t[137]=-9234.*Bconru*gdetlu*vcond;
  t[138]=t[133]+t[134];
  t[139]=t[135]+t[136]+t[137];
  t[140]=t[129]+t[130]+t[131];
  t[141]=t[132]+t[138]+t[139];
  t[142]=12798.*Bconu*gdetlu*vcond;
  t[143]=-40527.*Bconl*gdetr*vcond;
  t[144]=7695.*Bconld*gdetr*vcond;
  t[145]=29241.*Bconlu*gdetr*vcond;
  t[146]=-10665.*Bconr*gdetr*vcond;
  t[147]=2025.*Bconrd*gdetr*vcond;
  t[148]=7695.*Bconru*gdetr*vcond;
  t[149]=-76437.*Bconu*gdetr*vcond;
  t[150]=7695.*Bconl*gdetrd*vcond;
  t[151]=t[146]+t[147];
  t[152]=t[148]+t[149]+t[150];
  t[153]=t[142]+t[143]+t[144];
  t[154]=t[145]+t[151]+t[152];
  t[155]=t[114]+t[115]+t[127]+t[128];
  t[156]=t[140]+t[141]+t[153]+t[154];
  t[157]=342.*Bconld*gdetrd*vcond;
  t[158]=-9234.*Bconlu*gdetrd*vcond;
  t[159]=2025.*Bconr*gdetrd*vcond;
  t[160]=90.*Bconrd*gdetrd*vcond;
  t[161]=-2430.*Bconru*gdetrd*vcond;
  t[162]=24138.*Bconu*gdetrd*vcond;
  t[163]=29241.*Bconl*gdetru*vcond;
  t[164]=-9234.*Bconld*gdetru*vcond;
  t[165]=-9234.*Bconlu*gdetru*vcond;
  t[166]=t[161]+t[162];
  t[167]=t[163]+t[164]+t[165];
  t[168]=t[157]+t[158]+t[159];
  t[169]=t[160]+t[166]+t[167];
  t[170]=7695.*Bconr*gdetru*vcond;
  t[171]=-2430.*Bconrd*gdetru*vcond;
  t[172]=-2430.*Bconru*gdetru*vcond;
  t[173]=24138.*Bconu*gdetru*vcond;
  t[174]=-40527.*Bconl*gdetu*vcond;
  t[175]=12798.*Bconld*gdetu*vcond;
  t[176]=12798.*Bconlu*gdetu*vcond;
  t[177]=-76437.*Bconr*gdetu*vcond;
  t[178]=24138.*Bconrd*gdetu*vcond;
  t[179]=t[174]+t[175];
  t[180]=t[176]+t[177]+t[178];
  t[181]=t[170]+t[171]+t[172];
  t[182]=t[173]+t[179]+t[180];
  t[183]=24138.*Bconru*gdetu*vcond;
  t[184]=-46764.*Bconu*gdetu*vcond;
  t[185]=38970.*Bconl*gdetc*vconl;
  t[186]=-10665.*Bconld*gdetc*vconl;
  t[187]=-20115.*Bconlu*gdetc*vconl;
  t[188]=148086.*Bconr*gdetc*vconl;
  t[189]=-40527.*Bconrd*gdetc*vconl;
  t[190]=-76437.*Bconru*gdetc*vconl;
  t[191]=105939.*Bconu*gdetc*vconl;
  t[192]=t[187]+t[188];
  t[193]=t[189]+t[190]+t[191];
  t[194]=t[183]+t[184]+t[185];
  t[195]=t[186]+t[192]+t[193];
  t[196]=-10665.*Bconl*gdetd*vconl;
  t[197]=2025.*Bconld*gdetd*vconl;
  t[198]=7695.*Bconlu*gdetd*vconl;
  t[199]=-40527.*Bconr*gdetd*vconl;
  t[200]=7695.*Bconrd*gdetd*vconl;
  t[201]=29241.*Bconru*gdetd*vconl;
  t[202]=-40527.*Bconu*gdetd*vconl;
  t[203]=1732.*Bconl*gdetl*vconl;
  t[204]=-474.*Bconld*gdetl*vconl;
  t[205]=t[200]+t[201];
  t[206]=t[202]+t[203]+t[204];
  t[207]=t[196]+t[197]+t[198];
  t[208]=t[199]+t[205]+t[206];
  t[209]=t[168]+t[169]+t[181]+t[182];
  t[210]=t[194]+t[195]+t[207]+t[208];
  t[211]=t[50]+t[51]+t[104]+t[105];
  t[212]=t[155]+t[156]+t[209]+t[210];
  t[213]=-894.*Bconlu*gdetl*vconl;
  t[214]=-46764.*Bconr*gdetl*vconl;
  t[215]=12798.*Bconrd*gdetl*vconl;
  t[216]=24138.*Bconru*gdetl*vconl;
  t[217]=-20115.*Bconu*gdetl*vconl;
  t[218]=-474.*Bconl*gdetld*vconl;
  t[219]=90.*Bconld*gdetld*vconl;
  t[220]=342.*Bconlu*gdetld*vconl;
  t[221]=t[213]+t[214]+t[215]+t[216];
  t[222]=t[217]+t[218]+t[219]+t[220];
  t[223]=12798.*Bconr*gdetld*vconl;
  t[224]=-2430.*Bconrd*gdetld*vconl;
  t[225]=-9234.*Bconru*gdetld*vconl;
  t[226]=7695.*Bconu*gdetld*vconl;
  t[227]=-894.*Bconl*gdetlu*vconl;
  t[228]=342.*Bconld*gdetlu*vconl;
  t[229]=90.*Bconlu*gdetlu*vconl;
  t[230]=24138.*Bconr*gdetlu*vconl;
  t[231]=-9234.*Bconrd*gdetlu*vconl;
  t[232]=t[227]+t[228];
  t[233]=t[229]+t[230]+t[231];
  t[234]=t[223]+t[224]+t[225];
  t[235]=t[226]+t[232]+t[233];
  t[236]=-2430.*Bconru*gdetlu*vconl;
  t[237]=2025.*Bconu*gdetlu*vconl;
  t[238]=-46764.*Bconl*gdetr*vconl;
  t[239]=12798.*Bconld*gdetr*vconl;
  t[240]=24138.*Bconlu*gdetr*vconl;
  t[241]=-46764.*Bconr*gdetr*vconl;
  t[242]=12798.*Bconrd*gdetr*vconl;
  t[243]=24138.*Bconru*gdetr*vconl;
  t[244]=-76437.*Bconu*gdetr*vconl;
  t[245]=t[240]+t[241];
  t[246]=t[242]+t[243]+t[244];
  t[247]=t[236]+t[237]+t[238];
  t[248]=t[239]+t[245]+t[246];
  t[249]=12798.*Bconl*gdetrd*vconl;
  t[250]=-2430.*Bconld*gdetrd*vconl;
  t[251]=-9234.*Bconlu*gdetrd*vconl;
  t[252]=12798.*Bconr*gdetrd*vconl;
  t[253]=-2430.*Bconrd*gdetrd*vconl;
  t[254]=-9234.*Bconru*gdetrd*vconl;
  t[255]=29241.*Bconu*gdetrd*vconl;
  t[256]=24138.*Bconl*gdetru*vconl;
  t[257]=-9234.*Bconld*gdetru*vconl;
  t[258]=t[253]+t[254];
  t[259]=t[255]+t[256]+t[257];
  t[260]=t[249]+t[250]+t[251];
  t[261]=t[252]+t[258]+t[259];
  t[262]=t[221]+t[222]+t[234]+t[235];
  t[263]=t[247]+t[248]+t[260]+t[261];
  t[264]=-2430.*Bconlu*gdetru*vconl;
  t[265]=24138.*Bconr*gdetru*vconl;
  t[266]=-9234.*Bconrd*gdetru*vconl;
  t[267]=-2430.*Bconru*gdetru*vconl;
  t[268]=7695.*Bconu*gdetru*vconl;
  t[269]=-20115.*Bconl*gdetu*vconl;
  t[270]=7695.*Bconld*gdetu*vconl;
  t[271]=2025.*Bconlu*gdetu*vconl;
  t[272]=-76437.*Bconr*gdetu*vconl;
  t[273]=t[268]+t[269];
  t[274]=t[270]+t[271]+t[272];
  t[275]=t[264]+t[265]+t[266];
  t[276]=t[267]+t[273]+t[274];
  t[277]=29241.*Bconrd*gdetu*vconl;
  t[278]=7695.*Bconru*gdetu*vconl;
  t[279]=-10665.*Bconu*gdetu*vconl;
  t[280]=-10665.*Bconl*gdetc*vconld;
  t[281]=2025.*Bconld*gdetc*vconld;
  t[282]=7695.*Bconlu*gdetc*vconld;
  t[283]=-40527.*Bconr*gdetc*vconld;
  t[284]=7695.*Bconrd*gdetc*vconld;
  t[285]=29241.*Bconru*gdetc*vconld;
  t[286]=t[281]+t[282];
  t[287]=t[283]+t[284]+t[285];
  t[288]=t[277]+t[278]+t[279];
  t[289]=t[280]+t[286]+t[287];
  t[290]=-40527.*Bconu*gdetc*vconld;
  t[291]=2025.*Bconl*gdetd*vconld;
  t[292]=90.*Bconld*gdetd*vconld;
  t[293]=-2430.*Bconlu*gdetd*vconld;
  t[294]=7695.*Bconr*gdetd*vconld;
  t[295]=342.*Bconrd*gdetd*vconld;
  t[296]=-9234.*Bconru*gdetd*vconld;
  t[297]=12798.*Bconu*gdetd*vconld;
  t[298]=-474.*Bconl*gdetl*vconld;
  t[299]=t[294]+t[295];
  t[300]=t[296]+t[297]+t[298];
  t[301]=t[290]+t[291]+t[292];
  t[302]=t[293]+t[299]+t[300];
  t[303]=90.*Bconld*gdetl*vconld;
  t[304]=342.*Bconlu*gdetl*vconld;
  t[305]=12798.*Bconr*gdetl*vconld;
  t[306]=-2430.*Bconrd*gdetl*vconld;
  t[307]=-9234.*Bconru*gdetl*vconld;
  t[308]=7695.*Bconu*gdetl*vconld;
  t[309]=90.*Bconl*gdetld*vconld;
  t[310]=4.*Bconld*gdetld*vconld;
  t[311]=-108.*Bconlu*gdetld*vconld;
  t[312]=t[307]+t[308];
  t[313]=t[309]+t[310]+t[311];
  t[314]=t[303]+t[304]+t[305];
  t[315]=t[306]+t[312]+t[313];
  t[316]=t[275]+t[276]+t[288]+t[289];
  t[317]=t[301]+t[302]+t[314]+t[315];
  t[318]=-2430.*Bconr*gdetld*vconld;
  t[319]=-108.*Bconrd*gdetld*vconld;
  t[320]=2916.*Bconru*gdetld*vconld;
  t[321]=-2430.*Bconu*gdetld*vconld;
  t[322]=342.*Bconl*gdetlu*vconld;
  t[323]=-108.*Bconld*gdetlu*vconld;
  t[324]=-108.*Bconlu*gdetlu*vconld;
  t[325]=-9234.*Bconr*gdetlu*vconld;
  t[326]=t[318]+t[319]+t[320]+t[321];
  t[327]=t[322]+t[323]+t[324]+t[325];
  t[328]=2916.*Bconrd*gdetlu*vconld;
  t[329]=2916.*Bconru*gdetlu*vconld;
  t[330]=-2430.*Bconu*gdetlu*vconld;
  t[331]=12798.*Bconl*gdetr*vconld;
  t[332]=-2430.*Bconld*gdetr*vconld;
  t[333]=-9234.*Bconlu*gdetr*vconld;
  t[334]=12798.*Bconr*gdetr*vconld;
  t[335]=-2430.*Bconrd*gdetr*vconld;
  t[336]=-9234.*Bconru*gdetr*vconld;
  t[337]=t[332]+t[333];
  t[338]=t[334]+t[335]+t[336];
  t[339]=t[328]+t[329]+t[330];
  t[340]=t[331]+t[337]+t[338];
  t[341]=29241.*Bconu*gdetr*vconld;
  t[342]=-2430.*Bconl*gdetrd*vconld;
  t[343]=-108.*Bconld*gdetrd*vconld;
  t[344]=2916.*Bconlu*gdetrd*vconld;
  t[345]=-2430.*Bconr*gdetrd*vconld;
  t[346]=-108.*Bconrd*gdetrd*vconld;
  t[347]=2916.*Bconru*gdetrd*vconld;
  t[348]=-9234.*Bconu*gdetrd*vconld;
  t[349]=-9234.*Bconl*gdetru*vconld;
  t[350]=t[345]+t[346];
  t[351]=t[347]+t[348]+t[349];
  t[352]=t[341]+t[342]+t[343];
  t[353]=t[344]+t[350]+t[351];
  t[354]=2916.*Bconld*gdetru*vconld;
  t[355]=2916.*Bconlu*gdetru*vconld;
  t[356]=-9234.*Bconr*gdetru*vconld;
  t[357]=2916.*Bconrd*gdetru*vconld;
  t[358]=2916.*Bconru*gdetru*vconld;
  t[359]=-9234.*Bconu*gdetru*vconld;
  t[360]=7695.*Bconl*gdetu*vconld;
  t[361]=-2430.*Bconld*gdetu*vconld;
  t[362]=-2430.*Bconlu*gdetu*vconld;
  t[363]=t[358]+t[359];
  t[364]=t[360]+t[361]+t[362];
  t[365]=t[354]+t[355]+t[356];
  t[366]=t[357]+t[363]+t[364];
  t[367]=t[326]+t[327]+t[339]+t[340];
  t[368]=t[352]+t[353]+t[365]+t[366];
  t[369]=29241.*Bconr*gdetu*vconld;
  t[370]=-9234.*Bconrd*gdetu*vconld;
  t[371]=-9234.*Bconru*gdetu*vconld;
  t[372]=12798.*Bconu*gdetu*vconld;
  t[373]=-20115.*Bconl*gdetc*vconlu;
  t[374]=7695.*Bconld*gdetc*vconlu;
  t[375]=2025.*Bconlu*gdetc*vconlu;
  t[376]=-76437.*Bconr*gdetc*vconlu;
  t[377]=29241.*Bconrd*gdetc*vconlu;
  t[378]=t[373]+t[374];
  t[379]=t[375]+t[376]+t[377];
  t[380]=t[369]+t[370]+t[371];
  t[381]=t[372]+t[378]+t[379];
  t[382]=7695.*Bconru*gdetc*vconlu;
  t[383]=-10665.*Bconu*gdetc*vconlu;
  t[384]=7695.*Bconl*gdetd*vconlu;
  t[385]=-2430.*Bconld*gdetd*vconlu;
  t[386]=-2430.*Bconlu*gdetd*vconlu;
  t[387]=29241.*Bconr*gdetd*vconlu;
  t[388]=-9234.*Bconrd*gdetd*vconlu;
  t[389]=-9234.*Bconru*gdetd*vconlu;
  t[390]=12798.*Bconu*gdetd*vconlu;
  t[391]=t[386]+t[387];
  t[392]=t[388]+t[389]+t[390];
  t[393]=t[382]+t[383]+t[384];
  t[394]=t[385]+t[391]+t[392];
  t[395]=-894.*Bconl*gdetl*vconlu;
  t[396]=342.*Bconld*gdetl*vconlu;
  t[397]=90.*Bconlu*gdetl*vconlu;
  t[398]=24138.*Bconr*gdetl*vconlu;
  t[399]=-9234.*Bconrd*gdetl*vconlu;
  t[400]=-2430.*Bconru*gdetl*vconlu;
  t[401]=2025.*Bconu*gdetl*vconlu;
  t[402]=342.*Bconl*gdetld*vconlu;
  t[403]=-108.*Bconld*gdetld*vconlu;
  t[404]=t[399]+t[400];
  t[405]=t[401]+t[402]+t[403];
  t[406]=t[395]+t[396]+t[397];
  t[407]=t[398]+t[404]+t[405];
  t[408]=-108.*Bconlu*gdetld*vconlu;
  t[409]=-9234.*Bconr*gdetld*vconlu;
  t[410]=2916.*Bconrd*gdetld*vconlu;
  t[411]=2916.*Bconru*gdetld*vconlu;
  t[412]=-2430.*Bconu*gdetld*vconlu;
  t[413]=90.*Bconl*gdetlu*vconlu;
  t[414]=-108.*Bconld*gdetlu*vconlu;
  t[415]=4.*Bconlu*gdetlu*vconlu;
  t[416]=-2430.*Bconr*gdetlu*vconlu;
  t[417]=t[412]+t[413];
  t[418]=t[414]+t[415]+t[416];
  t[419]=t[408]+t[409]+t[410];
  t[420]=t[411]+t[417]+t[418];
  t[421]=t[380]+t[381]+t[393]+t[394];
  t[422]=t[406]+t[407]+t[419]+t[420];
  t[423]=t[262]+t[263]+t[316]+t[317];
  t[424]=t[367]+t[368]+t[421]+t[422];
  t[425]=2916.*Bconrd*gdetlu*vconlu;
  t[426]=-108.*Bconru*gdetlu*vconlu;
  t[427]=90.*Bconu*gdetlu*vconlu;
  t[428]=24138.*Bconl*gdetr*vconlu;
  t[429]=-9234.*Bconld*gdetr*vconlu;
  t[430]=-2430.*Bconlu*gdetr*vconlu;
  t[431]=24138.*Bconr*gdetr*vconlu;
  t[432]=-9234.*Bconrd*gdetr*vconlu;
  t[433]=t[425]+t[426]+t[427]+t[428];
  t[434]=t[429]+t[430]+t[431]+t[432];
  t[435]=-2430.*Bconru*gdetr*vconlu;
  t[436]=7695.*Bconu*gdetr*vconlu;
  t[437]=-9234.*Bconl*gdetrd*vconlu;
  t[438]=2916.*Bconld*gdetrd*vconlu;
  t[439]=2916.*Bconlu*gdetrd*vconlu;
  t[440]=-9234.*Bconr*gdetrd*vconlu;
  t[441]=2916.*Bconrd*gdetrd*vconlu;
  t[442]=2916.*Bconru*gdetrd*vconlu;
  t[443]=-9234.*Bconu*gdetrd*vconlu;
  t[444]=t[439]+t[440];
  t[445]=t[441]+t[442]+t[443];
  t[446]=t[435]+t[436]+t[437];
  t[447]=t[438]+t[444]+t[445];
  t[448]=-2430.*Bconl*gdetru*vconlu;
  t[449]=2916.*Bconld*gdetru*vconlu;
  t[450]=-108.*Bconlu*gdetru*vconlu;
  t[451]=-2430.*Bconr*gdetru*vconlu;
  t[452]=2916.*Bconrd*gdetru*vconlu;
  t[453]=-108.*Bconru*gdetru*vconlu;
  t[454]=342.*Bconu*gdetru*vconlu;
  t[455]=2025.*Bconl*gdetu*vconlu;
  t[456]=-2430.*Bconld*gdetu*vconlu;
  t[457]=t[452]+t[453];
  t[458]=t[454]+t[455]+t[456];
  t[459]=t[448]+t[449]+t[450];
  t[460]=t[451]+t[457]+t[458];
  t[461]=90.*Bconlu*gdetu*vconlu;
  t[462]=7695.*Bconr*gdetu*vconlu;
  t[463]=-9234.*Bconrd*gdetu*vconlu;
  t[464]=342.*Bconru*gdetu*vconlu;
  t[465]=-474.*Bconu*gdetu*vconlu;
  t[466]=148086.*Bconl*gdetc*vconr;
  t[467]=-40527.*Bconld*gdetc*vconr;
  t[468]=-76437.*Bconlu*gdetc*vconr;
  t[469]=38970.*Bconr*gdetc*vconr;
  t[470]=t[465]+t[466];
  t[471]=t[467]+t[468]+t[469];
  t[472]=t[461]+t[462]+t[463];
  t[473]=t[464]+t[470]+t[471];
  t[474]=t[433]+t[434]+t[446]+t[447];
  t[475]=t[459]+t[460]+t[472]+t[473];
  t[476]=-10665.*Bconrd*gdetc*vconr;
  t[477]=-20115.*Bconru*gdetc*vconr;
  t[478]=199809.*Bconu*gdetc*vconr;
  t[479]=-40527.*Bconl*gdetd*vconr;
  t[480]=7695.*Bconld*gdetd*vconr;
  t[481]=29241.*Bconlu*gdetd*vconr;
  t[482]=-10665.*Bconr*gdetd*vconr;
  t[483]=2025.*Bconrd*gdetd*vconr;
  t[484]=7695.*Bconru*gdetd*vconr;
  t[485]=t[480]+t[481];
  t[486]=t[482]+t[483]+t[484];
  t[487]=t[476]+t[477]+t[478];
  t[488]=t[479]+t[485]+t[486];
  t[489]=-76437.*Bconu*gdetd*vconr;
  t[490]=-46764.*Bconl*gdetl*vconr;
  t[491]=12798.*Bconld*gdetl*vconr;
  t[492]=24138.*Bconlu*gdetl*vconr;
  t[493]=-46764.*Bconr*gdetl*vconr;
  t[494]=12798.*Bconrd*gdetl*vconr;
  t[495]=24138.*Bconru*gdetl*vconr;
  t[496]=-76437.*Bconu*gdetl*vconr;
  t[497]=12798.*Bconl*gdetld*vconr;
  t[498]=t[493]+t[494];
  t[499]=t[495]+t[496]+t[497];
  t[500]=t[489]+t[490]+t[491];
  t[501]=t[492]+t[498]+t[499];
  t[502]=-2430.*Bconld*gdetld*vconr;
  t[503]=-9234.*Bconlu*gdetld*vconr;
  t[504]=12798.*Bconr*gdetld*vconr;
  t[505]=-2430.*Bconrd*gdetld*vconr;
  t[506]=-9234.*Bconru*gdetld*vconr;
  t[507]=29241.*Bconu*gdetld*vconr;
  t[508]=24138.*Bconl*gdetlu*vconr;
  t[509]=-9234.*Bconld*gdetlu*vconr;
  t[510]=-2430.*Bconlu*gdetlu*vconr;
  t[511]=t[506]+t[507];
  t[512]=t[508]+t[509]+t[510];
  t[513]=t[502]+t[503]+t[504];
  t[514]=t[505]+t[511]+t[512];
  t[515]=24138.*Bconr*gdetlu*vconr;
  t[516]=-9234.*Bconrd*gdetlu*vconr;
  t[517]=-2430.*Bconru*gdetlu*vconr;
  t[518]=7695.*Bconu*gdetlu*vconr;
  t[519]=-46764.*Bconl*gdetr*vconr;
  t[520]=12798.*Bconld*gdetr*vconr;
  t[521]=24138.*Bconlu*gdetr*vconr;
  t[522]=1732.*Bconr*gdetr*vconr;
  t[523]=-474.*Bconrd*gdetr*vconr;
  t[524]=t[519]+t[520];
  t[525]=t[521]+t[522]+t[523];
  t[526]=t[515]+t[516]+t[517];
  t[527]=t[518]+t[524]+t[525];
  t[528]=t[487]+t[488]+t[500]+t[501];
  t[529]=t[513]+t[514]+t[526]+t[527];
  t[530]=-894.*Bconru*gdetr*vconr;
  t[531]=-20115.*Bconu*gdetr*vconr;
  t[532]=12798.*Bconl*gdetrd*vconr;
  t[533]=-2430.*Bconld*gdetrd*vconr;
  t[534]=-9234.*Bconlu*gdetrd*vconr;
  t[535]=-474.*Bconr*gdetrd*vconr;
  t[536]=90.*Bconrd*gdetrd*vconr;
  t[537]=342.*Bconru*gdetrd*vconr;
  t[538]=t[530]+t[531]+t[532]+t[533];
  t[539]=t[534]+t[535]+t[536]+t[537];
  t[540]=7695.*Bconu*gdetrd*vconr;
  t[541]=24138.*Bconl*gdetru*vconr;
  t[542]=-9234.*Bconld*gdetru*vconr;
  t[543]=-2430.*Bconlu*gdetru*vconr;
  t[544]=-894.*Bconr*gdetru*vconr;
  t[545]=342.*Bconrd*gdetru*vconr;
  t[546]=90.*Bconru*gdetru*vconr;
  t[547]=2025.*Bconu*gdetru*vconr;
  t[548]=-76437.*Bconl*gdetu*vconr;
  t[549]=t[544]+t[545];
  t[550]=t[546]+t[547]+t[548];
  t[551]=t[540]+t[541]+t[542];
  t[552]=t[543]+t[549]+t[550];
  t[553]=29241.*Bconld*gdetu*vconr;
  t[554]=7695.*Bconlu*gdetu*vconr;
  t[555]=-20115.*Bconr*gdetu*vconr;
  t[556]=7695.*Bconrd*gdetu*vconr;
  t[557]=2025.*Bconru*gdetu*vconr;
  t[558]=-20115.*Bconu*gdetu*vconr;
  t[559]=-40527.*Bconl*gdetc*vconrd;
  t[560]=7695.*Bconld*gdetc*vconrd;
  t[561]=29241.*Bconlu*gdetc*vconrd;
  t[562]=t[557]+t[558];
  t[563]=t[559]+t[560]+t[561];
  t[564]=t[553]+t[554]+t[555];
  t[565]=t[556]+t[562]+t[563];
  t[566]=-10665.*Bconr*gdetc*vconrd;
  t[567]=2025.*Bconrd*gdetc*vconrd;
  t[568]=7695.*Bconru*gdetc*vconrd;
  t[569]=-76437.*Bconu*gdetc*vconrd;
  t[570]=7695.*Bconl*gdetd*vconrd;
  t[571]=342.*Bconld*gdetd*vconrd;
  t[572]=-9234.*Bconlu*gdetd*vconrd;
  t[573]=2025.*Bconr*gdetd*vconrd;
  t[574]=90.*Bconrd*gdetd*vconrd;
  t[575]=t[570]+t[571];
  t[576]=t[572]+t[573]+t[574];
  t[577]=t[566]+t[567]+t[568];
  t[578]=t[569]+t[575]+t[576];
  t[579]=t[538]+t[539]+t[551]+t[552];
  t[580]=t[564]+t[565]+t[577]+t[578];
  t[581]=-2430.*Bconru*gdetd*vconrd;
  t[582]=24138.*Bconu*gdetd*vconrd;
  t[583]=12798.*Bconl*gdetl*vconrd;
  t[584]=-2430.*Bconld*gdetl*vconrd;
  t[585]=-9234.*Bconlu*gdetl*vconrd;
  t[586]=12798.*Bconr*gdetl*vconrd;
  t[587]=-2430.*Bconrd*gdetl*vconrd;
  t[588]=-9234.*Bconru*gdetl*vconrd;
  t[589]=29241.*Bconu*gdetl*vconrd;
  t[590]=t[585]+t[586];
  t[591]=t[587]+t[588]+t[589];
  t[592]=t[581]+t[582]+t[583];
  t[593]=t[584]+t[590]+t[591];
  t[594]=-2430.*Bconl*gdetld*vconrd;
  t[595]=-108.*Bconld*gdetld*vconrd;
  t[596]=2916.*Bconlu*gdetld*vconrd;
  t[597]=-2430.*Bconr*gdetld*vconrd;
  t[598]=-108.*Bconrd*gdetld*vconrd;
  t[599]=2916.*Bconru*gdetld*vconrd;
  t[600]=-9234.*Bconu*gdetld*vconrd;
  t[601]=-9234.*Bconl*gdetlu*vconrd;
  t[602]=2916.*Bconld*gdetlu*vconrd;
  t[603]=t[598]+t[599];
  t[604]=t[600]+t[601]+t[602];
  t[605]=t[594]+t[595]+t[596];
  t[606]=t[597]+t[603]+t[604];
  t[607]=2916.*Bconlu*gdetlu*vconrd;
  t[608]=-9234.*Bconr*gdetlu*vconrd;
  t[609]=2916.*Bconrd*gdetlu*vconrd;
  t[610]=2916.*Bconru*gdetlu*vconrd;
  t[611]=-9234.*Bconu*gdetlu*vconrd;
  t[612]=12798.*Bconl*gdetr*vconrd;
  t[613]=-2430.*Bconld*gdetr*vconrd;
  t[614]=-9234.*Bconlu*gdetr*vconrd;
  t[615]=-474.*Bconr*gdetr*vconrd;
  t[616]=t[611]+t[612];
  t[617]=t[613]+t[614]+t[615];
  t[618]=t[607]+t[608]+t[609];
  t[619]=t[610]+t[616]+t[617];
  t[620]=90.*Bconrd*gdetr*vconrd;
  t[621]=342.*Bconru*gdetr*vconrd;
  t[622]=7695.*Bconu*gdetr*vconrd;
  t[623]=-2430.*Bconl*gdetrd*vconrd;
  t[624]=-108.*Bconld*gdetrd*vconrd;
  t[625]=2916.*Bconlu*gdetrd*vconrd;
  t[626]=90.*Bconr*gdetrd*vconrd;
  t[627]=4.*Bconrd*gdetrd*vconrd;
  t[628]=-108.*Bconru*gdetrd*vconrd;
  t[629]=t[624]+t[625];
  t[630]=t[626]+t[627]+t[628];
  t[631]=t[620]+t[621]+t[622];
  t[632]=t[623]+t[629]+t[630];
  t[633]=t[592]+t[593]+t[605]+t[606];
  t[634]=t[618]+t[619]+t[631]+t[632];
  t[635]=t[474]+t[475]+t[528]+t[529];
  t[636]=t[579]+t[580]+t[633]+t[634];
  t[637]=-2430.*Bconu*gdetrd*vconrd;
  t[638]=-9234.*Bconl*gdetru*vconrd;
  t[639]=2916.*Bconld*gdetru*vconrd;
  t[640]=2916.*Bconlu*gdetru*vconrd;
  t[641]=342.*Bconr*gdetru*vconrd;
  t[642]=-108.*Bconrd*gdetru*vconrd;
  t[643]=-108.*Bconru*gdetru*vconrd;
  t[644]=-2430.*Bconu*gdetru*vconrd;
  t[645]=t[637]+t[638]+t[639]+t[640];
  t[646]=t[641]+t[642]+t[643]+t[644];
  t[647]=29241.*Bconl*gdetu*vconrd;
  t[648]=-9234.*Bconld*gdetu*vconrd;
  t[649]=-9234.*Bconlu*gdetu*vconrd;
  t[650]=7695.*Bconr*gdetu*vconrd;
  t[651]=-2430.*Bconrd*gdetu*vconrd;
  t[652]=-2430.*Bconru*gdetu*vconrd;
  t[653]=24138.*Bconu*gdetu*vconrd;
  t[654]=-76437.*Bconl*gdetc*vconru;
  t[655]=29241.*Bconld*gdetc*vconru;
  t[656]=t[651]+t[652];
  t[657]=t[653]+t[654]+t[655];
  t[658]=t[647]+t[648]+t[649];
  t[659]=t[650]+t[656]+t[657];
  t[660]=7695.*Bconlu*gdetc*vconru;
  t[661]=-20115.*Bconr*gdetc*vconru;
  t[662]=7695.*Bconrd*gdetc*vconru;
  t[663]=2025.*Bconru*gdetc*vconru;
  t[664]=-20115.*Bconu*gdetc*vconru;
  t[665]=29241.*Bconl*gdetd*vconru;
  t[666]=-9234.*Bconld*gdetd*vconru;
  t[667]=-9234.*Bconlu*gdetd*vconru;
  t[668]=7695.*Bconr*gdetd*vconru;
  t[669]=t[664]+t[665];
  t[670]=t[666]+t[667]+t[668];
  t[671]=t[660]+t[661]+t[662];
  t[672]=t[663]+t[669]+t[670];
  t[673]=-2430.*Bconrd*gdetd*vconru;
  t[674]=-2430.*Bconru*gdetd*vconru;
  t[675]=24138.*Bconu*gdetd*vconru;
  t[676]=24138.*Bconl*gdetl*vconru;
  t[677]=-9234.*Bconld*gdetl*vconru;
  t[678]=-2430.*Bconlu*gdetl*vconru;
  t[679]=24138.*Bconr*gdetl*vconru;
  t[680]=-9234.*Bconrd*gdetl*vconru;
  t[681]=-2430.*Bconru*gdetl*vconru;
  t[682]=t[677]+t[678];
  t[683]=t[679]+t[680]+t[681];
  t[684]=t[673]+t[674]+t[675];
  t[685]=t[676]+t[682]+t[683];
  t[686]=t[645]+t[646]+t[658]+t[659];
  t[687]=t[671]+t[672]+t[684]+t[685];
  t[688]=7695.*Bconu*gdetl*vconru;
  t[689]=-9234.*Bconl*gdetld*vconru;
  t[690]=2916.*Bconld*gdetld*vconru;
  t[691]=2916.*Bconlu*gdetld*vconru;
  t[692]=-9234.*Bconr*gdetld*vconru;
  t[693]=2916.*Bconrd*gdetld*vconru;
  t[694]=2916.*Bconru*gdetld*vconru;
  t[695]=-9234.*Bconu*gdetld*vconru;
  t[696]=-2430.*Bconl*gdetlu*vconru;
  t[697]=t[692]+t[693];
  t[698]=t[694]+t[695]+t[696];
  t[699]=t[688]+t[689]+t[690];
  t[700]=t[691]+t[697]+t[698];
  t[701]=2916.*Bconld*gdetlu*vconru;
  t[702]=-108.*Bconlu*gdetlu*vconru;
  t[703]=-2430.*Bconr*gdetlu*vconru;
  t[704]=2916.*Bconrd*gdetlu*vconru;
  t[705]=-108.*Bconru*gdetlu*vconru;
  t[706]=342.*Bconu*gdetlu*vconru;
  t[707]=24138.*Bconl*gdetr*vconru;
  t[708]=-9234.*Bconld*gdetr*vconru;
  t[709]=-2430.*Bconlu*gdetr*vconru;
  t[710]=t[705]+t[706];
  t[711]=t[707]+t[708]+t[709];
  t[712]=t[701]+t[702]+t[703];
  t[713]=t[704]+t[710]+t[711];
  t[714]=-894.*Bconr*gdetr*vconru;
  t[715]=342.*Bconrd*gdetr*vconru;
  t[716]=90.*Bconru*gdetr*vconru;
  t[717]=2025.*Bconu*gdetr*vconru;
  t[718]=-9234.*Bconl*gdetrd*vconru;
  t[719]=2916.*Bconld*gdetrd*vconru;
  t[720]=2916.*Bconlu*gdetrd*vconru;
  t[721]=342.*Bconr*gdetrd*vconru;
  t[722]=-108.*Bconrd*gdetrd*vconru;
  t[723]=t[718]+t[719];
  t[724]=t[720]+t[721]+t[722];
  t[725]=t[714]+t[715]+t[716];
  t[726]=t[717]+t[723]+t[724];
  t[727]=-108.*Bconru*gdetrd*vconru;
  t[728]=-2430.*Bconu*gdetrd*vconru;
  t[729]=-2430.*Bconl*gdetru*vconru;
  t[730]=2916.*Bconld*gdetru*vconru;
  t[731]=-108.*Bconlu*gdetru*vconru;
  t[732]=90.*Bconr*gdetru*vconru;
  t[733]=-108.*Bconrd*gdetru*vconru;
  t[734]=-44096.*Bconru*gdetru*vconru;
  t[735]=90.*Bconu*gdetru*vconru;
  t[736]=t[731]+t[732];
  t[737]=t[733]+t[734]+t[735];
  t[738]=t[727]+t[728]+t[729];
  t[739]=t[730]+t[736]+t[737];
  t[740]=t[699]+t[700]+t[712]+t[713];
  t[741]=t[725]+t[726]+t[738]+t[739];
  t[742]=7695.*Bconl*gdetu*vconru;
  t[743]=-9234.*Bconld*gdetu*vconru;
  t[744]=342.*Bconlu*gdetu*vconru;
  t[745]=2025.*Bconr*gdetu*vconru;
  t[746]=-2430.*Bconrd*gdetu*vconru;
  t[747]=90.*Bconru*gdetu*vconru;
  t[748]=-894.*Bconu*gdetu*vconru;
  t[749]=105939.*Bconl*gdetc*vconu;
  t[750]=-40527.*Bconld*gdetc*vconu;
  t[751]=t[746]+t[747];
  t[752]=t[748]+t[749]+t[750];
  t[753]=t[742]+t[743]+t[744];
  t[754]=t[745]+t[751]+t[752];
  t[755]=-10665.*Bconlu*gdetc*vconu;
  t[756]=199809.*Bconr*gdetc*vconu;
  t[757]=-76437.*Bconrd*gdetc*vconu;
  t[758]=-20115.*Bconru*gdetc*vconu;
  t[759]=38970.*Bconu*gdetc*vconu;
  t[760]=-40527.*Bconl*gdetd*vconu;
  t[761]=12798.*Bconld*gdetd*vconu;
  t[762]=12798.*Bconlu*gdetd*vconu;
  t[763]=-76437.*Bconr*gdetd*vconu;
  t[764]=t[759]+t[760];
  t[765]=t[761]+t[762]+t[763];
  t[766]=t[755]+t[756]+t[757];
  t[767]=t[758]+t[764]+t[765];
  t[768]=24138.*Bconrd*gdetd*vconu;
  t[769]=24138.*Bconru*gdetd*vconu;
  t[770]=-46764.*Bconu*gdetd*vconu;
  t[771]=-20115.*Bconl*gdetl*vconu;
  t[772]=7695.*Bconld*gdetl*vconu;
  t[773]=2025.*Bconlu*gdetl*vconu;
  t[774]=-76437.*Bconr*gdetl*vconu;
  t[775]=29241.*Bconrd*gdetl*vconu;
  t[776]=7695.*Bconru*gdetl*vconu;
  t[777]=t[772]+t[773];
  t[778]=t[774]+t[775]+t[776];
  t[779]=t[768]+t[769]+t[770];
  t[780]=t[771]+t[777]+t[778];
  t[781]=-10665.*Bconu*gdetl*vconu;
  t[782]=7695.*Bconl*gdetld*vconu;
  t[783]=-2430.*Bconld*gdetld*vconu;
  t[784]=-2430.*Bconlu*gdetld*vconu;
  t[785]=29241.*Bconr*gdetld*vconu;
  t[786]=-9234.*Bconrd*gdetld*vconu;
  t[787]=-9234.*Bconru*gdetld*vconu;
  t[788]=12798.*Bconu*gdetld*vconu;
  t[789]=2025.*Bconl*gdetlu*vconu;
  t[790]=t[785]+t[786];
  t[791]=t[787]+t[788]+t[789];
  t[792]=t[781]+t[782]+t[783];
  t[793]=t[784]+t[790]+t[791];
  t[794]=t[753]+t[754]+t[766]+t[767];
  t[795]=t[779]+t[780]+t[792]+t[793];
  t[796]=-2430.*Bconld*gdetlu*vconu;
  t[797]=90.*Bconlu*gdetlu*vconu;
  t[798]=7695.*Bconr*gdetlu*vconu;
  t[799]=-9234.*Bconrd*gdetlu*vconu;
  t[800]=342.*Bconru*gdetlu*vconu;
  t[801]=-474.*Bconu*gdetlu*vconu;
  t[802]=-76437.*Bconl*gdetr*vconu;
  t[803]=29241.*Bconld*gdetr*vconu;
  t[804]=7695.*Bconlu*gdetr*vconu;
  t[805]=t[800]+t[801];
  t[806]=t[802]+t[803]+t[804];
  t[807]=t[796]+t[797]+t[798];
  t[808]=t[799]+t[805]+t[806];
  t[809]=-20115.*Bconr*gdetr*vconu;
  t[810]=7695.*Bconrd*gdetr*vconu;
  t[811]=2025.*Bconru*gdetr*vconu;
  t[812]=-20115.*Bconu*gdetr*vconu;
  t[813]=29241.*Bconl*gdetrd*vconu;
  t[814]=-9234.*Bconld*gdetrd*vconu;
  t[815]=-9234.*Bconlu*gdetrd*vconu;
  t[816]=7695.*Bconr*gdetrd*vconu;
  t[817]=-2430.*Bconrd*gdetrd*vconu;
  t[818]=t[813]+t[814];
  t[819]=t[815]+t[816]+t[817];
  t[820]=t[809]+t[810]+t[811];
  t[821]=t[812]+t[818]+t[819];
  t[822]=-2430.*Bconru*gdetrd*vconu;
  t[823]=24138.*Bconu*gdetrd*vconu;
  t[824]=7695.*Bconl*gdetru*vconu;
  t[825]=-9234.*Bconld*gdetru*vconu;
  t[826]=342.*Bconlu*gdetru*vconu;
  t[827]=2025.*Bconr*gdetru*vconu;
  t[828]=-2430.*Bconrd*gdetru*vconu;
  t[829]=90.*Bconru*gdetru*vconu;
  t[830]=-894.*Bconu*gdetru*vconu;
  t[831]=t[826]+t[827];
  t[832]=t[828]+t[829]+t[830];
  t[833]=t[822]+t[823]+t[824];
  t[834]=t[825]+t[831]+t[832];
  t[835]=-10665.*Bconl*gdetu*vconu;
  t[836]=12798.*Bconld*gdetu*vconu;
  t[837]=-474.*Bconlu*gdetu*vconu;
  t[838]=-20115.*Bconr*gdetu*vconu;
  t[839]=24138.*Bconrd*gdetu*vconu;
  t[840]=-894.*Bconru*gdetu*vconu;
  t[841]=1732.*Bconu*gdetu*vconu;
  t[842]=-13509.*gdetlu*vconc;
  t[843]=35313.*gdetr*vconc;
  t[844]=-6705.*gdetrd*vconc;
  t[845]=-25479.*gdetru*vconc;
  t[846]=49362.*gdetu*vconc;
  t[847]=t[842]+t[843];
  t[848]=t[844]+t[845]+t[846];
  t[849]=4266.*gdetlu*vcond;
  t[850]=-6705.*gdetr*vcond;
  t[851]=-298.*gdetrd*vcond;
  t[852]=8046.*gdetru*vcond;
  t[853]=-15588.*gdetu*vcond;
  t[854]=2565.*gdetlu*vconl;
  t[855]=t[849]+t[850]+t[851];
  t[856]=t[852]+t[853]+t[854];
  t[857]=-13509.*gdetr*vconl;
  t[858]=2565.*gdetrd*vconl;
  t[859]=9747.*gdetru*vconl;
  t[860]=-13509.*gdetu*vconl;
  t[861]=-810.*gdetlu*vconld;
  t[862]=2565.*gdetr*vconld;
  t[863]=t[857]+t[858]+t[859];
  t[864]=t[860]+t[861]+t[862];
  t[865]=114.*gdetrd*vconld;
  t[866]=-3078.*gdetru*vconld;
  t[867]=4266.*gdetu*vconld;
  t[868]=-810.*gdetlu*vconlu;
  t[869]=9747.*gdetr*vconlu;
  t[870]=-3078.*gdetrd*vconlu;
  t[871]=t[865]+t[866]+t[867];
  t[872]=t[868]+t[869]+t[870];
  t[873]=t[847]+t[848]+t[855]+t[856];
  t[874]=t[863]+t[864]+t[871]+t[872];
  t[875]=-3078.*gdetru*vconlu;
  t[876]=4266.*gdetu*vconlu;
  t[877]=9747.*gdetlu*vconr;
  t[878]=-3555.*gdetr*vconr;
  t[879]=675.*gdetrd*vconr;
  t[880]=2565.*gdetru*vconr;
  t[881]=t[875]+t[876]+t[877];
  t[882]=t[878]+t[879]+t[880];
  t[883]=-25479.*gdetu*vconr;
  t[884]=-3078.*gdetlu*vconrd;
  t[885]=675.*gdetr*vconrd;
  t[886]=30.*gdetrd*vconrd;
  t[887]=-810.*gdetru*vconrd;
  t[888]=8046.*gdetu*vconrd;
  t[889]=t[883]+t[884]+t[885];
  t[890]=t[886]+t[887]+t[888];
  t[891]=-3078.*gdetlu*vconru;
  t[892]=2565.*gdetr*vconru;
  t[893]=-810.*gdetrd*vconru;
  t[894]=-810.*gdetru*vconru;
  t[895]=8046.*gdetu*vconru;
  t[896]=6241.*vconc;
  t[897]=-395.*vcond;
  t[898]=-395.*vconl;
  t[899]=75.*vconld;
  t[900]=285.*vconlu;
  t[901]=-1501.*vconr;
  t[902]=285.*vconrd;
  t[903]=1083.*vconru;
  t[904]=-1501.*vconu;
  t[905]=t[897]+t[898]+t[899]+t[900];
  t[906]=t[901]+t[902]+t[903]+t[904];
  t[907]=t[905]+t[906];
  t[908]=t[896];
  t[909]=3.*t[907];
  t[910]=gdetl*(t[908]+t[909]);
  t[911]=t[895];
  t[912]=3.*t[910];
  t[913]=t[891]+t[892]+t[893];
  t[914]=t[894]+t[911]+t[912];
  t[915]=4266.*gdetlu*vconu;
  t[916]=-25479.*gdetr*vconu;
  t[917]=8046.*gdetrd*vconu;
  t[918]=8046.*gdetru*vconu;
  t[919]=-15588.*gdetu*vconu;
  t[920]=-3555.*vconc;
  t[921]=-158.*vcond;
  t[922]=675.*vconl;
  t[923]=30.*vconld;
  t[924]=-810.*vconlu;
  t[925]=2565.*vconr;
  t[926]=114.*vconrd;
  t[927]=-3078.*vconru;
  t[928]=4266.*vconu;
  t[929]=t[924]+t[925];
  t[930]=t[926]+t[927]+t[928];
  t[931]=t[920]+t[921]+t[922];
  t[932]=t[923]+t[929]+t[930];
  t[933]=t[919];
  t[934]=gdetld*(t[931]+t[932]);
  t[935]=t[915]+t[916]+t[917];
  t[936]=t[918]+t[933]+t[934];
  t[937]=t[881]+t[882]+t[889]+t[890];
  t[938]=t[913]+t[914]+t[935]+t[936];
  t[939]=t[873]+t[874]+t[937]+t[938];
  t[940]=38970.*vconc;
  t[941]=1732.*vcond;
  t[942]=3555.*vconl;
  t[943]=158.*vconld;
  t[944]=-4266.*vconlu;
  t[945]=6705.*vconr;
  t[946]=298.*vconrd;
  t[947]=-8046.*vconru;
  t[948]=15588.*vconu;
  t[949]=t[942]+t[943]+t[944];
  t[950]=t[945]+t[946]+t[947]+t[948];
  t[951]=t[949]+t[950];
  t[952]=t[941];
  t[953]=-3.*t[951];
  t[954]=-3.*gdetc*vartemp14387;
  t[955]=gdetd*(t[940]+t[952]+t[953]);
  t[956]=t[954];
  t[957]=3.*t[939];
  t[958]=t[955]+t[956];
  t[959]=749956.*vconc;
  t[960]=68414.*vcond;
  t[961]=68414.*vconl;
  t[962]=-18723.*vconld;
  t[963]=-35313.*vconlu;
  t[964]=129034.*vconr;
  t[965]=-35313.*vconrd;
  t[966]=-66603.*vconru;
  t[967]=129034.*vconu;
  t[968]=t[960]+t[961]+t[962]+t[963];
  t[969]=t[964]+t[965]+t[966]+t[967];
  t[970]=t[968]+t[969];
  t[971]=t[959];
  t[972]=-3.*t[970];
  t[973]=-18723.*gdetld*vconc;
  t[974]=-35313.*gdetlu*vconc;
  t[975]=129034.*gdetr*vconc;
  t[976]=-35313.*gdetrd*vconc;
  t[977]=-66603.*gdetru*vconc;
  t[978]=129034.*gdetu*vconc;
  t[979]=3555.*gdetld*vcond;
  t[980]=t[973]+t[974]+t[975];
  t[981]=t[976]+t[977]+t[978]+t[979];
  t[982]=13509.*gdetlu*vcond;
  t[983]=-35313.*gdetr*vcond;
  t[984]=6705.*gdetrd*vcond;
  t[985]=25479.*gdetru*vcond;
  t[986]=-49362.*gdetu*vcond;
  t[987]=3555.*gdetld*vconl;
  t[988]=6705.*gdetlu*vconl;
  t[989]=t[982]+t[983]+t[984];
  t[990]=t[985]+t[986]+t[987]+t[988];
  t[991]=-49362.*gdetr*vconl;
  t[992]=13509.*gdetrd*vconl;
  t[993]=25479.*gdetru*vconl;
  t[994]=-35313.*gdetu*vconl;
  t[995]=-675.*gdetld*vconld;
  t[996]=-2565.*gdetlu*vconld;
  t[997]=13509.*gdetr*vconld;
  t[998]=t[991]+t[992]+t[993];
  t[999]=t[994]+t[995]+t[996]+t[997];
  t[1000]=-2565.*gdetrd*vconld;
  t[1001]=-9747.*gdetru*vconld;
  t[1002]=13509.*gdetu*vconld;
  t[1003]=-2565.*gdetld*vconlu;
  t[1004]=-675.*gdetlu*vconlu;
  t[1005]=25479.*gdetr*vconlu;
  t[1006]=-9747.*gdetrd*vconlu;
  t[1007]=t[1000]+t[1001]+t[1002];
  t[1008]=t[1003]+t[1004]+t[1005]+t[1006];
  t[1009]=t[980]+t[981]+t[989]+t[990];
  t[1010]=t[998]+t[999]+t[1007]+t[1008];
  t[1011]=-2565.*gdetru*vconlu;
  t[1012]=3555.*gdetu*vconlu;
  t[1013]=13509.*gdetld*vconr;
  t[1014]=25479.*gdetlu*vconr;
  t[1015]=-12990.*gdetr*vconr;
  t[1016]=3555.*gdetrd*vconr;
  t[1017]=6705.*gdetru*vconr;
  t[1018]=t[1011]+t[1012]+t[1013];
  t[1019]=t[1014]+t[1015]+t[1016]+t[1017];
  t[1020]=-66603.*gdetu*vconr;
  t[1021]=-2565.*gdetld*vconrd;
  t[1022]=-9747.*gdetlu*vconrd;
  t[1023]=3555.*gdetr*vconrd;
  t[1024]=-675.*gdetrd*vconrd;
  t[1025]=-2565.*gdetru*vconrd;
  t[1026]=25479.*gdetu*vconrd;
  t[1027]=t[1020]+t[1021]+t[1022];
  t[1028]=t[1023]+t[1024]+t[1025]+t[1026];
  t[1029]=-9747.*gdetld*vconru;
  t[1030]=-2565.*gdetlu*vconru;
  t[1031]=6705.*gdetr*vconru;
  t[1032]=-2565.*gdetrd*vconru;
  t[1033]=-675.*gdetru*vconru;
  t[1034]=6705.*gdetu*vconru;
  t[1035]=13509.*gdetld*vconu;
  t[1036]=t[1029]+t[1030]+t[1031];
  t[1037]=t[1032]+t[1033]+t[1034]+t[1035];
  t[1038]=3555.*gdetlu*vconu;
  t[1039]=-66603.*gdetr*vconu;
  t[1040]=25479.*gdetrd*vconu;
  t[1041]=6705.*gdetru*vconu;
  t[1042]=-12990.*gdetu*vconu;
  t[1043]=-18723.*vcond;
  t[1044]=-12990.*vconl;
  t[1045]=3555.*vconld;
  t[1046]=6705.*vconlu;
  t[1047]=-49362.*vconr;
  t[1048]=13509.*vconrd;
  t[1049]=25479.*vconru;
  t[1050]=-35313.*vconu+vartemp14376;
  t[1051]=t[1043]+t[1044]+t[1045]+t[1046];
  t[1052]=t[1047]+t[1048]+t[1049]+t[1050];
  t[1053]=gdetd*vartemp14387;
  t[1054]=gdetl*(t[1051]+t[1052]);
  t[1055]=t[1053];
  t[1056]=t[1038]+t[1039]+t[1040];
  t[1057]=t[1041]+t[1042]+t[1054]+t[1055];
  t[1058]=t[1018]+t[1019]+t[1027]+t[1028];
  t[1059]=t[1036]+t[1037]+t[1056]+t[1057];
  t[1060]=t[1009]+t[1010]+t[1058]+t[1059];
  t[1061]=gdetc*(t[971]+t[972]);
  t[1062]=-3.*t[1060];
  t[1063]=Bcond*(t[957]+t[958]);
  t[1064]=Bconc*(t[1061]+t[1062]);
  t[1065]=t[839]+t[840];
  t[1066]=t[841]+t[1063]+t[1064];
  t[1067]=t[835]+t[836]+t[837];
  t[1068]=t[838]+t[1065]+t[1066];
  t[1069]=t[807]+t[808]+t[820]+t[821];
  t[1070]=t[833]+t[834]+t[1067]+t[1068];
  t[1071]=t[686]+t[687]+t[740]+t[741];
  t[1072]=t[794]+t[795]+t[1069]+t[1070];
  t[1073]=t[211]+t[212]+t[423]+t[424];
  t[1074]=t[635]+t[636]+t[1071]+t[1072];
  t[1075]=t[1073]+t[1074];
  *Fru=0.000022675736961451247165532879818594*t[1075];






  return(0);


}




