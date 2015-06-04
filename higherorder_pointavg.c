#include "decs.h"


// Notes:
//
// 1) For FV or FLUXRECON(non-point for field) methods, c2a is not direct inverse of a2c.  This leads to many issues, but for example the point field will not preserve the correct divb=0 when computing the flux differences.  Del Zanna et al. (2003) apply an iterative approach to inverting the stencil so the point value still gives a2c(c2a(f))=f and a2c(f+g)=a2cf+a2cg that are required to preserve divb=0.

// 2) Point field version of FLUXRECON method is bad since forces fixed a2c stencil.



// below used by flux.c before any flux calls
void preinterp_flux_point2avg(void)
{
  int i,j,k,dir;

#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )  //SUPERSASMARK atch
  //zero out the array with lower order fractions at every substep before c2e reconstructions are performed
  if(DOENOFLUX!=NOENOFLUX){
    COMPFULLLOOP {
      DIMENLOOP(dir) MACP1A0(weno_prim_lower_order_fraction,dir,i,j,k) = 0.0;
    }
  }
#endif


}






// whether to debug bound_flux()
#define DEBUGBOUNDFLUX 0

int vectorpot_useflux(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int flux_point2avg(int stage, int whichmaem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecother[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int Nvec[NDIM];
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  int i,j,k;


  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  fluxvec[1] = F1;
  fluxvec[2] = F2;
  fluxvec[3] = F3;
  
  //  FULLLOOP{
  //    dualfprintf(fail_file,"BEFORE i=%d j=%d k=%d : F1[B2]=%21.15g F2[B1]=%21.15g\n",i,j,k,MACP0A1(F1,i,j,k,B2),MACP0A1(F2,i,j,k,B1));
  //  }
  
  // when called, will treat EMF (i.e. F[B1,B2,B3]) terms as normal flux as required since the called functions know about FLUXB
  flux_point2avg(STAGEM1, ISEMONLY, pr, Nvec, fluxvec,NULL);

  //  FULLLOOP{
  //    dualfprintf(fail_file,"AFTER i=%d j=%d k=%d : F1[B2]=%21.15g F2[B1]=%21.15g\n",i,j,k,MACP0A1(F1,i,j,k,B2),MACP0A1(F2,i,j,k,B1));
  //  }
  

  return(0);
}


// integrate point fluxes AT face to get average fluxes OVER entire face surface
int flux_point2avg(int stage, int whichmaem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecother[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int flux_deaverage_fluxrecon(int stage, int whichmaem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecother[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int bound_flux_fluxrecon(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int debug_boundfluxinitial(FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int debug_boundfluxfinal(FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);

  int flux_integrate_fluxsplit(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);

  int flux_integrate_finitevolume(int stage, int whichmaem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);


  if(NSPECIAL!=0){
    dualfprintf(fail_file,"flux_point2avg() not setup for NSPECIAL!=0\n");
    myexit(29855252);
  }


  if(DOENOFLUX==ENOFLUXRECON){

    // need to bound flux before de-averaging since need extra cells to get active value of de-averaged flux
    // normally bound all NPR=NFLUXBOUND
    bound_flux_fluxrecon(stage,pr,Nvec,fluxvec);

    if(DEBUGBOUNDFLUX){
      debug_boundfluxinitial(fluxvec);
    }
    
    flux_deaverage_fluxrecon(stage, whichmaem, pr, Nvec, fluxvec,fluxvecother);

    if(DEBUGBOUNDFLUX){
      debug_boundfluxfinal(fluxvec);
    }
       
  }
  else if(DOENOFLUX==ENOFLUXSPLIT){

    // not setup
    flux_integrate_fluxsplit(stage, pr, Nvec, fluxvec);

  }
  else if(DOENOFLUX==ENOFINITEVOLUME){

    // setup so no bounding of anything except primitives needed (extra ghost+active layer used)
    flux_integrate_finitevolume(stage, whichmaem, pr, Nvec, fluxvec);
    
  }


#if(PRODUCTION==0)
  trifprintf( "e");
#endif
    

  return(0);

}




// ideas: GODMARK:
// 1) Could use fluxsplitting for MA terms without problem of field at face
//    Then can use a2e_i(F_i) [average to edge].  Problem is that interpolating flux implies interpolating different parts of Lorentz factor and other stiff parameters.
//

// quasi-deaverage fluxes along their direction
// assumes only gets fluxvec on active region using ghost+active fluxes
// need to change ENOFLUXRECONTYPE to ENOFLUXRECONTYPEGHOSTACTIVE if have fluxes everywhere and want in ghost+active + active regions
int flux_deaverage_fluxrecon(int stage, int whichmaem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*fluxvecother[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int odir1,odir2;
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  extern void flux_interp(int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);
  int simple_a2c_limit_onepl(int dotransverse, int dir, int pl, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  FTYPE limit_fluxc2a_prim_change( 
                                  int dir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                                  FTYPE (*fluxvec_point)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                  FTYPE (*fluxvec_avg)[NSTORE2][NSTORE3][NPR+NSPECIAL]);   //atch
  int pl,pliter,i,j,k;
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int dotrans,dotransem,dotransma;
  int ftemp;
  int emforvectorpot_fluxrecon(int stage, int isemf, int *fluxdirlist, int *pldirlist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  FTYPE signforfluxlist[NDIM];
  int treatemfspecial;
  int dotranslocal[NPR];
  int whichquantity,interporflux;
  int extralimiting;
  int weightsplittype;
  FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR];
  FTYPE (*otherin)[NSTORE2][NSTORE3][NPR];
  int odirindex, theotherdir, odirarray[2];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;






#if(0)
  LOOPF2 dualfprintf(fail_file,"BEFORE: emf(0)=%21.15g %21.15g emf(N1)=%21.15g %21.15g diff0-N1 = %21.15g\n",MACP1A1(fluxvec,1,0,j,0,B2),-MACP1A1(fluxvec,2,0,j,0,B1),MACP1A1(fluxvec,1,N1,j,0,B2),-MACP1A1(fluxvec,2,N1,j,0,B1),MACP1A1(fluxvec,1,0,j,0,B2)-MACP1A1(fluxvec,1,N1,j,0,B2));
#endif

  // convert quasi-averaged fluxes to quasi-point fluxes
  // applies to fields the same as other quantities even for FLUXCTSTAG (although flux at slightly different position --GODMARK)


  //  bound_prim(STAGEM1,pr);


  PLOOP(pliter,pl) dotranslocal[pl]=do_transverse_flux_integration[pl];
  // now modify local given whichmaem setting
  if(whichmaem==ISMAANDEM){
    // then no modification
  }
  else if(whichmaem==ISMAONLY){
    // then don't need to do dB/dt parts -- never will MA-based flux be in dB/dt's EMF flux
    PLOOPBONLY(pl) dotranslocal[pl]=0;
  }
  else if(whichmaem==ISEMONLY){
    // then just do everything in case field there
  }


  ////////////////////
  //
  // choose interpolation and quantity types and choose to do extra limiting or not
  //
  ////////////////////

  interporflux=ENOFLUXRECONTYPE;
    
  if(whichmaem==ISMAANDEM){
    whichquantity=ENOFLUX;
  }
  else if(whichmaem==ISMAONLY){
    whichquantity=ENOMAFLUX;
  }
  else if(whichmaem==ISEMONLY){
    whichquantity=EMFLUXFIELDTERMTYPE;
  }

  higherorder_set(whichquantity, CVT_A2C, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_AC_FRAC_CHANGE) extralimiting=1;
    else extralimiting=0;
  }


  // any trans
  dotrans=0;
  PLOOP(pliter,pl) dotrans+=dotranslocal[pl];
  dotrans*=(DOENOFLUX == ENOFLUXRECON);



  dotransma=(DOENOFLUX == ENOFLUXRECON)&&
    (dotranslocal[RHO]
     ||dotranslocal[UU]
     ||dotranslocal[U1]
     ||dotranslocal[U2]
     ||dotranslocal[U3]);
  
  dotransem=(DOENOFLUX == ENOFLUXRECON)&&
    (dotranslocal[B1]
     ||dotranslocal[B2]
     ||dotranslocal[B3]);
  

  ////////////////////////////
  //
  // only use special FV emf method for staggered field method.
  // Assume all other methods use centered field and regular flux for field such that EMF is not treated specially
  // Do not treat emf as special when updating point fields
  if(FLUXB==FLUXCTSTAG && extrazones4emf){
    treatemfspecial=1;
  }
  else treatemfspecial=0;




  // CHANGINGMARK: want to avoid going into flux_interp() if not necessary (just check if any variables left?)
  //  if(whichmaem==ISMAANDEM || whichmaem==ISMAONLY || (whichmaem==ISEMONLY&&treatemfspecial==0) ){
  if(dotrans){

    // then simpler to have separate "emf" section, so remove fields for all calculations
    // default is to average all fluxes
    PALLLOOP(pl) whichpltoavg[pl]=dotranslocal[pl];// default
    PALLLOOP(pl) ifnotavgthencopy[pl]=0; // no need to copy since already exists

    if(treatemfspecial){
      // never do field since always do as EMF separately below
      PLOOPBONLY(pl) whichpltoavg[pl]=0;
      PLOOPBONLY(pl) ifnotavgthencopy[pl]=0;
    }
    else{
      PLOOPBONLY(pl) whichpltoavg[pl]=dotranslocal[pl];
      PLOOPBONLY(pl) ifnotavgthencopy[pl]=0;

    }


    ///////////////////////////
    //
    // loop over flux directions and obtain quasi-deaveraged flux
    //
    // F1 input and F1 as output is ok since each full 1D line is operated on, and each 1D line never used information from other lines, and only 1 direction per F1/F2/F3, so completely done with that dataset
    //
    ///////////////////////////
    DIMENLOOP(dir){
    
      if(Nvec[dir]!=1){


        // other dimensions
        odir1=dir%3+1;
        odir2=(dir+1)%3+1;

        odirarray[0]=odir1;
        odirarray[1]=odir2;



        if(whichmaem==ISMAONLY && splitmaem==1){
          // override above and deal with flux[FLUXSPLITPMA(dir)] term
          // this term really contains diagonal gas pressure
          pl=FLUXSPLITPMA(dir);
          whichpltoavg[pl]=dotranslocal[UU+dir];
          ifnotavgthencopy[pl]=0;
        }
    
        //////////////////////////////////
        // go ahead and control nprlist
        //
        // needed for Sasha's limit_fluxc2a (flux_interp internally would remove unwanted pl's based upon whichpltoavg/ifnotavgthencopy)
        //
        //////////////////////////////////
        addremovefromnpr(REMOVEFROMNPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);

        ////////////////////////////////////////////
        //
        // if no variables to do, then continue to next direction or if no more directions just done
        //
        ////////////////////////////////////////////

        if(nprlocalend>-1){

          ////////////////////////////////////////////
          //
          // Setup idel,jdel,kdel
          //
          ////////////////////////////////////////////
          idel = fluxloop[dir][FIDEL];
          jdel = fluxloop[dir][FJDEL];
          kdel = fluxloop[dir][FKDEL];


          /////////////
          //
          // store original lower-order flux
          // used for any extra limiting after higher-order calculation is done
          //
          /////////////
          if(extralimiting){
            COMPFULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(fluxvectemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl);
          }


 
          /////////////
          //
          // Determine variable to be used to determine stencil for higher-order calculation
          //
          // ok to use stencilvartemp since the flux_intepr() function used below never calls flux_interp_multiple() or avg2cen() that use stencilvartemp
          //
          /////////////
          // if fluxvecother!=NULL then assume "stencil" chosen by sum of original and other fluxes
          if(fluxvecother==NULL){
            COMPFULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl);
          }
          else{ // then using fluxvectother[] for some part of stencilvartemp[]
            // using 1|| above since this is less accurate than separating fluxes for stencil when using SPLITMAEM
            if(whichmaem==ISMAANDEM || (whichmaem==ISMAONLY && splitmaem==0) || whichmaem==ISEMONLY){
              COMPFULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl) + MACP1A1(fluxvecother,dir,i,j,k,pl);
            }
            else if(whichmaem==ISMAONLY && splitmaem==1){
              // have to be more careful now with MA+EM split and pdiag term moved to FLUXSPLITPMA(dir)
              // for no obvious reason, for test=106 the below is different than the above (this code leads to no failure for no obvious reasons -- probably minimization code or something is wrong)


              COMPFULLLOOP{
                PLOOPNOB1(pl)  GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl) + MACP1A1(fluxvecother,dir,i,j,k,pl);
                PLOOPNOB2(pl)  GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl) + MACP1A1(fluxvecother,dir,i,j,k,pl);

                // GODMARK: In stiff flows need to couple pressure and velocity, but otherwise pressure can be independent
                // By "couple" I mean when minimizing weights need to include pressure in stiff regimes
                // Must keep pressure out when velocity small so can have static equilibrium, for example, at high order.
                // so create stiffness that is d\gamma/\gamma and dv/u almost like already have, but not static stiff things like b^2/rho or u/rho...just want to compare velocity and pressure/internal energy
 
#if(SPLITPRESSURETERMINFLUXMA)
                // assume pdiag can't be mixed with EM term, but for pressure term add back in ram pressure term
                // don't want gas pressure term inside ram pressure because pressure can be constant with variations in velocity and so won't resolve variations in velocity if pressure nearly constant and noisy
                pl=UU+dir; GLOBALMACP0A1(stencilvartemp,i,j,k,pl) += MACP1A1(fluxvec,dir,i,j,k,FLUXSPLITPMA(dir));
                // to the gas pressure term, add full magnetic term instead of just some ad-hoc pressure term
                // This is done to ensure higher-order force-balance when magnetic and gas fluids spatially meet
                // go ahead and add full MA and full EM term as stencil for pressure term since if v->0 still good and want higher order if this flux is continuous if other terms (when minimized) still give mass-energy-momentum flux that doesn't reduce -- so now UU+dir and FLUXSPLITPMA terms both have same stencil weight choice at this level (then minimization later processes further)
                pl=FLUXSPLITPMA(dir); GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl) + MACP1A1(fluxvec,dir,i,j,k,UU+dir) + MACP1A1(fluxvecother,dir,i,j,k,UU+dir);


                /////////////////////
                // BEGIN DEBUG:
#if(0)
                // assume pdiag can't be mixed with EM term, but for pressure term add back in ram pressure term
                // CHANGINGMARK:
                pl=FLUXSPLITPMA(dir); GLOBALMACP0A1(stencilvartemp,i,j,k,pl) += MACP1A1(fluxvec,dir,i,j,k,UU+dir);
                pl=UU+dir; GLOBALMACP0A1(stencilvartemp,i,j,k,pl) += MACP1A1(fluxvec,dir,i,j,k,FLUXSPLITPMA(dir));
#endif
#if(0)
                // CHANGINGMARK:
                pl=UU+dir;  MACP1A1(fluxvec,dir,i,j,k,pl) += MACP1A1(fluxvec,dir,i,j,k,FLUXSPLITPMA(dir));
                pl=FLUXSPLITPMA(dir); MACP1A1(fluxvec,dir,i,j,k,pl) =0.0;
#endif
                // END  DEBUG:
                /////////////////////
#endif


              }// end COMPFULLLOOP
            }// end else if splitmaem==1 and ismaonly==1
            else{
              dualfprintf(fail_file,"No such stencilvartemp[] condition: %d %d\n",whichmaem,splitmaem);
              myexit(1751536);
            }


            if(weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS){
              ///////////////
              //
              // modify strencilvar so all quantities have comparable values (should be able to use this method under any case of splitmaem or SPLITPRESSUREFLUXMA/EM)
              //
              ///////////////
              COMPFULLLOOP{
                get_geometry(i,j,k,FACE1+dir-1,ptrgeom); // location of flux is always on face

                // modify stress-energy terms that are orthogonal to dir-flux along dir
                // energy stencilvar not modified since this method combines RHO,UU,UU+dir weights without direct comparison of weight values and then UU+dir version used as basis for comparison
                for(odirindex=0;odirindex<=1;odirindex++){
                  theotherdir=odirarray[odirindex];
                  pl=UU+theotherdir;
                  GLOBALMACP0A1(stencilvartemp,i,j,k,pl) *= sqrt(fabs(ptrgeom->gcon[GIND(theotherdir,theotherdir)])/(fabs(ptrgeom->gcon[GIND(dir,dir)])+SMALL));
                }

                PLOOP(pliter,pl){
                  if(pl==RHO || (pl>=B3+1 && pl<NPR)){
                    // KORALTODO: Needs different scaling as with vectors above.
                    // GODMARK: assume all other are advective and get into form so like pressure but with advective component weighting
                    // Then have \rho u^i u_i Y_e for term that was \rho Y_e u^i
                    // no fabs() around first stencilvartemp term so picks up sign of u^i (e.g. u^i u_i is constant in reflecting wall type shock)
                    // as currently written has same signature as original term but scales like pressure now
                    GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = fabs(MACP1A1(fluxvec,dir,i,j,k,UU+dir))*GLOBALMACP0A1(stencilvartemp,i,j,k,pl)/(fabs(MACP1A1(fluxvec,dir,i,j,k,RHO))+SMALL);
                  }
                }
              }// end COMPFULLLOOP
            }// end if special weightsplittype such that compare weights across quantities


          }// end else if using separate stencilvar


          ////////////////////
          //
          // Set stencilvar for higher-order stencil calculation
          //
          ////////////////////

          if(fluxvecother!=NULL || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS){
            stencilvar=GLOBALPOINT(stencilvartemp); // use modified specialized stencil variable
          }
          else stencilvar=NULL; // just use normal variables

          ////////////////////
          //
          // only single quasi-deaverage along flux direction if to be used to update point conserved quantities
          //
          ////////////////////

          //            COMPFULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(stencilvartemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl);

          flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dir, idel, jdel, kdel, pr, stencilvar, fluxvec[dir], NULL, fluxvec[dir], NULL); 


          if(extralimiting){
            // 0 below means limit along direction
            // limit this deaveragingn process since similar to differentiating that creates more structure
            // CHANGINGMARK: Must couple fluxes as coupled to determine stencil for a2c itself
            // CHANGINGMARK: so need to split "weight" and calculation
            PLOOP(pliter,pl) simple_a2c_limit_onepl(0, dir, pl, pr, GLOBALPOINT(fluxvectemp), fluxvec[dir]);
          }




        } // end if something to do

        // restore list
        addremovefromnpr(RESTORENPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


      }// end if dimension exists
    }// end over dimensions


  }// end if dotrans



  if(dotransem){
    if(treatemfspecial){
      // quasi-deaverage EMF
      // 1 means is an emf, second/third NULL isn't used for emf
      // don't really use plforfluxlist or signforfluxlist since just operate generically on one valid fluxvec for a given EMF and then assign other fluxvec an opposite sign when associated with the same EMF
      DIMENLOOP(dir) get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir],&signforfluxlist[dir]);
      emforvectorpot_fluxrecon(stage, 1, fluxdirlist, pldirlist, pr, Nvec, fluxvec);
    }
  }



#if(0)
  LOOPF2 dualfprintf(fail_file,"AFTEr: emf(0)=%21.15g %21.15g emf(N1)=%21.15g %21.15g diff0-N1 = %21.15g\n",MACP1A1(fluxvec,1,0,j,0,B2),-MACP1A1(fluxvec,2,0,j,0,B1),MACP1A1(fluxvec,1,N1,j,0,B2),-MACP1A1(fluxvec,2,N1,j,0,B1),MACP1A1(fluxvec,1,0,j,0,B2)-MACP1A1(fluxvec,1,N1,j,0,B2));
#endif



  return(0);

}











// applies to conserved quantity form with \detg of quasifield and pointfield, not without \detg
// reintegrate quasi-field to point field
// arrays quasifield and pointfield should be different pointers as required anyways for mathematical model to make sense
// need to reaverage the quasi-conserved magnetic flux
// c2a re-averages quasi-conserved magnetic flux to obtain point value of field
// notice that for FLUXRECON, here the field is quasi-deaveraged already as some fake conserved differenced field, which we now convert to a point value.  See flux_point2avg.c.  Hence ordering of utoinvert and myupoint looks funny given ENOCENT2AVGTYPE used
//
// convert quasi-averaged fluxes to quasi-point fluxes
int field_integrate_fluxrecon(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR])
{
  int odir1,odir2;
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  extern void flux_interp(int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);
  FTYPE limit_fluxc2a_prim_change( 
                                  int dir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                                  FTYPE (*fluxvec_point)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                  FTYPE (*fluxvec_avg)[NSTORE2][NSTORE3][NPR+NSPECIAL]);   //atch
  int pl,pliter,i,j,k;
  int pldir;
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int dotransem;
  int ftemp;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  int Nvec[NDIM];
  int dodir[NDIM];
  int whichmaem;
  int interporflux, whichquantity;




  dotransem=(DOENOFLUX == ENOFLUXRECON)&&
    (do_transverse_flux_integration[B1]
     ||do_transverse_flux_integration[B2]
     ||do_transverse_flux_integration[B3]);


  // setup behavior to be EM only
  whichmaem=ISEMONLY;
  interporflux=ENOQUASIFIELDFLUXRECONTYPE;

  if(whichmaem==ISMAANDEM){
    whichquantity=ENOFLUX;
  }
  else if(whichmaem==ISMAONLY){
    whichquantity=ENOMAFLUX;
  }
  else if(whichmaem==ISEMONLY){
    whichquantity=FIELDTERMTYPE;
  }


  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;



  // determine if doing a certain field direction
  DIMENLOOP(dir){
    // other dimensions
    odir1=dir%3+1;
    odir2=(dir+1)%3+1;

    // only do if dimension exists and not true that both other dimensions don't exist since then field along dir must be constant (optimization for 1D problems)
    dodir[dir] = (Nvec[dir]!=1 &&  (!(Nvec[odir1]==1 && Nvec[odir2]==1)));
    if(dotransem==0) dodir[dir]=0;
  }




  if(dotransem){

    ///////////////////////////
    //
    // loop over field directions (B1,B2,B3) for dir=1,2,3
    //
    ///////////////////////////
    DIMENLOOP(dir){

      if(dodir[dir]){

        pldir=B1+dir-1;
        // other dimensions
        odir1=dir%3+1;
        odir2=(dir+1)%3+1;


        // then simpler to have separate "emf" section, so remove fields for all calculations
        // default is to average all fluxes
        PALLLOOP(pl) whichpltoavg[pl]=0;
        PALLLOOP(pl) ifnotavgthencopy[pl]=0; // no need to copy since already exists
        PLOOPBONLY(pl) ifnotavgthencopy[pl]=1; // need to copy if not doing
        PLOOPBONLY(pl) whichpltoavg[pl]=0; // default
        // now choose which field to operate on
        whichpltoavg[pldir]=do_conserved_integration[pldir];
 
        // go ahead and control nprlist
        // needed for Sasha's limit_fluxc2a (flux_interp internally would remove unwanted pl's based upon whichpltoavg/ifnotavgthencopy)
        addremovefromnpr(REMOVEFROMNPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);

        ////////////////////////////////////////////
        //
        // if no variables to do, then continue to next direction or if no more directions just done
        //
        ////////////////////////////////////////////

        if(nprlocalend>-1){



          idel = fluxloop[dir][FIDEL];
          jdel = fluxloop[dir][FJDEL];
          kdel = fluxloop[dir][FKDEL];


          ////////////////////
          //
          // reintegrate quasifield to point field
          // ENOQUASIFIELDFLUXRECONTYPE defines an interporflux type that means c2a but along direction of flux and on face
          //
          ////////////////////
          flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dir, idel, jdel, kdel, pr, NULL, quasifield, NULL, pointfield, NULL); 
      

        } // end if something to do

        // restore list
        addremovefromnpr(RESTORENPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


      }// end if dimension exists
    }// end over dimensions

  }// end if dotrans



  ////////////////////////////
  //
  // Need to copy the field not created
  //
  ///////////////////////////
  DIMENLOOP(dir){

    // other dimensions
    odir1=dir%3+1;
    odir2=(dir+1)%3+1;
    pldir=B1+dir-1;
    
    if(dodir[dir]==0 ){
      // then need to copy over
      COMPFULLLOOP MACP0A1(pointfield,i,j,k,pldir)=MACP0A1(quasifield,i,j,k,pldir);
    }
  }



  return(0);

}




// initial Upoint -> Uavg for FLUXRECON
int initial_averageu_fluxrecon(int *fieldfrompotential, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*Upoint)[NSTORE2][NSTORE3][NPR], FTYPE (*Uavg)[NSTORE2][NSTORE3][NPR])
{
  int field_deaverage_fluxrecon(int *fieldfrompotential, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR]);
  int i,j,k,pl,pliter;

  
  // non-field U stays as point, but field Upoint needs to be quasi-deaveraged if field didn't come from potential in order to obtain correct "Uavg"
  // don't overwrite Uavg for field since may be from vector potential and so correct divb=0 form
  COMPFULLLOOP{
    PLOOPNOB1(pl) MACP0A1(Uavg,i,j,k,pl)=MACP0A1(Upoint,i,j,k,pl);
    PLOOPNOB2(pl) MACP0A1(Uavg,i,j,k,pl)=MACP0A1(Upoint,i,j,k,pl);
  }

  // now deal with conserved (divb=0) fields
  // using original Upoint and new Uavg to determine if to use Uavg
  if(extrazones4emf) field_deaverage_fluxrecon(fieldfrompotential, prim, Upoint, Uavg);
  else{
    // only points ever for all conserved quantities
    PLOOPBONLY(pl) MACP0A1(Uavg,i,j,k,pl)=MACP0A1(Upoint,i,j,k,pl);
  }

  //  FULLLOOP PLOOP(pliter,pl){
  //    dualfprintf(fail_file,"Uavg[%d][%d][%d][%d]=%21.15g\n",i,j,k,pl,MACP0A1(Uavg,i,j,k,pl));
  //  }



  return(0);

}





// quasi-deaverage the point conserved field (applied to \detg based field)
// assume generally that quasifield needs to be filled with pointfield if not doing higher order on some field
// pointfield must exist *everywhere* so can reconstruct quasifield
int field_deaverage_fluxrecon(int *fieldfrompotential, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR])
{
  int odir1,odir2;
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  extern void flux_interp(int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);
  int simple_a2c_limit_onepl(int dotransverse, int dir, int pl, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  int pl,pliter,i,j,k;
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int dotransem;
  int ftemp;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  int Nvec[NDIM];
  int dodir[NDIM];
  int pldir;
  FTYPE (*intemp)[NSTORE2][NSTORE3][NPR];
  int whichmaem;
  int whichquantity,interporflux;
  int weightsplittype;
  int extralimiting;


  // convert quasi-averaged fluxes to quasi-point fluxes
  // applies to fields the same as other quantities even for FLUXCTSTAG (although flux at slightly different position --GODMARK)


  dotransem=(DOENOFLUX == ENOFLUXRECON)&&
    (do_transverse_flux_integration[B1]
     ||do_transverse_flux_integration[B2]
     ||do_transverse_flux_integration[B3]);


  // setup behavior to be EM only
  whichmaem=ISEMONLY;
  interporflux=ENOFLUXRECONTYPEGHOSTACTIVE;

  if(whichmaem==ISMAANDEM){
    whichquantity=ENOFLUX;
  }
  else if(whichmaem==ISMAONLY){
    whichquantity=ENOMAFLUX;
  }
  else if(whichmaem==ISEMONLY){
    whichquantity=FIELDTERMTYPE;
  }

  higherorder_set(whichquantity, CVT_A2C, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_AC_FRAC_CHANGE) extralimiting=1;
    else extralimiting=0;
  }

  

  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;



  // determine if doing a certain field direction
  DIMENLOOP(dir){
    // other dimensions
    odir1=dir%3+1;
    odir2=(dir+1)%3+1;

    dodir[dir] = (Nvec[dir]!=1 &&  (!(Nvec[odir1]==1 && Nvec[odir2]==1)));
    if(dotransem==0) dodir[dir]=0;
    if(fieldfrompotential[dir]) dodir[dir]=0;

  }



  if(dotransem){

    ///////////////////////////
    //
    // loop over field directions (B1,B2,B3) for dir=1,2,3
    //
    ///////////////////////////
    DIMENLOOP(dir){

      if(dodir[dir]){ // only do if dimension exists and not true that both other dimensions don't exist since then field along dir must be constant (optimization for 1D problems)


        pldir=B1+dir-1;
        // other dimensions
        odir1=dir%3+1;
        odir2=(dir+1)%3+1;


        // DEBUG:
        //      dualfprintf(fail_file,"DOING dir=%d pldir=%d fieldfrompotential[%d]=%d\n",dir,pldir,dir,fieldfrompotential[dir]);




        // then simpler to have separate "emf" section, so remove fields for all calculations
        // default is to average all fluxes
        PALLLOOP(pl) whichpltoavg[pl]=0;
        PALLLOOP(pl) ifnotavgthencopy[pl]=0; // no need to copy since already exists
        PLOOPBONLY(pl) ifnotavgthencopy[pl]=1-fieldfrompotential[dir]; // don't copy if not doing since assume output is correct from vector potential
        PLOOPBONLY(pl) whichpltoavg[pl]=0; // default
        // now choose which field to operate on
        whichpltoavg[pldir]=do_conserved_integration[pldir];
 
        // go ahead and control nprlist
        // needed for Sasha's limit_fluxc2a (flux_interp internally would remove unwanted pl's based upon whichpltoavg/ifnotavgthencopy)
        addremovefromnpr(REMOVEFROMNPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


        ////////////////////////////////////////////
        //
        // if no variables to do, then continue to next direction or if no more directions just done
        //
        ////////////////////////////////////////////

        if(nprlocalend>-1){



          idel = fluxloop[dir][FIDEL];
          jdel = fluxloop[dir][FJDEL];
          kdel = fluxloop[dir][FKDEL];

          if(extralimiting){
            if(pointfield==quasifield){
              intemp=GLOBALPOINT(fluxvectemp);
              // below not needed if input and output arrays are different
              COMPFULLLOOP MACP0A1(intemp,i,j,k,pldir) = MACP0A1(pointfield,i,j,k,pldir);
            }
            else intemp=pointfield; // then just access pointfield directly since unchanged by flux_interp operation
          }// otherwise don't have to define intemp

          ////////////////////
          //
          // reintegrate quasifield to point field
          // ENOQUASIFIELDFLUXRECONTYPE defines an interporflux type that means c2a but along direction of flux and on face
          //
          ////////////////////
          flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dir, idel, jdel, kdel, pr, NULL, pointfield, NULL, quasifield, NULL); 
      

          if(extralimiting){
            // 0 below means limit along direction
            simple_a2c_limit_onepl(0, dir, pldir, pr, intemp, quasifield);
          }

        } // end if something to do


        // restore list
        addremovefromnpr(RESTORENPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


      }// end if dimension exists
    }// end over dimensions


  }// end if dotrans


  ////////////////////////////
  //
  // Need to copy the field not created and if not already set by vector potential
  //
  ///////////////////////////
  DIMENLOOP(dir){

    // other dimensions
    odir1=dir%3+1;
    odir2=(dir+1)%3+1;
    pldir=B1+dir-1;
    
    if(fieldfrompotential[dir]==0){
      if(dodir[dir]==0 ){
        // then need to copy over
        // dualfprintf(fail_file,"LASTCOPY dir=%d pldir=%d\n",dir,pldir);
        COMPFULLLOOP MACP0A1(quasifield,i,j,k,pldir)=MACP0A1(pointfield,i,j,k,pldir);
      }
    }
  }




  return(0);

}




// get Bhat that is actually only needed when computing divb in setfdivb()
int field_Bhat_fluxrecon(FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*pointfield)[NSTORE2][NSTORE3][NPR], FTYPE (*quasifield)[NSTORE2][NSTORE3][NPR])
{
  int fakefieldfrompotential[NDIM];
  int jj;

  // then get Bhat that is actually only needed when computing divb in setfdivb() in initbase.c
  fakefieldfrompotential[TT]=0;
  SLOOPA(jj) fakefieldfrompotential[jj]=0; // 0 means will modify and always want to obtain new quasifield
  field_deaverage_fluxrecon(fakefieldfrompotential,pr,pointfield,quasifield);
 

  return(0);

}








// Performs fluxrecon (quasi-deaverage) on emf or vector potential that has been placed in fluxvec
// does double-transverse quasi-deaveraging
// fluxdir[dir=1,2,3] gives which flux is used to store emf/vpot in dir-direction
// pldir[dir=1,2,3] gives the "field" B1+pldir-1 direction for the emf/vpot in dir-direction
// single input is treated as output, so no need to "copy" over results as some other functions do
int emforvectorpot_fluxrecon(int stage, int isemf, int *fluxdirlist, int *pldirlist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int plforflux,fluxdir,pldir;
  int myintdir;
  int orthodirs;
  int odir1,odir2;
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  extern void flux_interp_multiple(int *whichpltoavg, int *ifnotavgthencopy, int numdirs, int *whichquantitylist, int *interporfluxlist, int *dirmethodlist, int *Nvec, int *intdirlist, int *fluxdirlist, int *idellist, int *jdellist, int *kdellist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  FTYPE limit_fluxc2a_prim_change( 
                                  int dir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                                  FTYPE (*fluxvec_point)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                  FTYPE (*fluxvec_avg)[NSTORE2][NSTORE3][NPR+NSPECIAL]);   //atch
  int pl,pliter,i,j,k;
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int dotrans,dotransem,dotransma;
  int ftemp;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int numdirs, intdirlist[NDIM], idellist[NDIM],jdellist[NDIM],kdellist[NDIM];
  int whichquantitylist[NDIM],interporfluxlist[NDIM],dirmethodlist[NDIM];
  int fluxdirlistnew[NDIM];
  int otherdir;
  int simple_a2c_limit_onepl(int dotransverse, int dir, int pl, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  int whichmaem;
  int whichquantity;
  int interporflux,weightsplittype;
  int extralimiting;
  

  // convert quasi-averaged fluxes to quasi-point fluxes
  // applies to fields the same as other quantities even for FLUXCTSTAG (although flux at slightly different position --GODMARK)


  // setup behavior to be EM only
  whichmaem=ISEMONLY;

  // need ghost+active if [not emf] or [EMF but evolving Bhat instead of Bpoint]
  if(!isemf || extrazones4emf) interporflux=ENOFLUXRECONTYPEGHOSTACTIVE;
  else interporflux=ENOFLUXRECONTYPE;

  if(whichmaem==ISMAANDEM){
    whichquantity=ENOFLUX;
  }
  else if(whichmaem==ISMAONLY){
    whichquantity=ENOMAFLUX;
  }
  else if(whichmaem==ISEMONLY){
    whichquantity=EMFLUXFIELDTERMTYPE;
  }

  higherorder_set(whichquantity, CVT_A2C, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_AC_FRAC_CHANGE) extralimiting=1;
    else extralimiting=0;
  }
  


  dotransem=(DOENOFLUX == ENOFLUXRECON)&&
    (do_transverse_flux_integration[B1]
     ||do_transverse_flux_integration[B2]
     ||do_transverse_flux_integration[B3]);



  

  if(dotransem){

    ///////////////////////////
    //
    // loop over emf/vpot directions
    //
    ///////////////////////////
    DIMENLOOP(dir){


      // other dimensions
      odir1=dir%3+1;
      odir2=(dir+1)%3+1;


      // do emf[dir]
      if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){ // then emf[dir] is unimportant (cancels itself always), so assume original emf[dir] was set to 0 already

        // then doing emf[dir]
        fluxdir=fluxdirlist[dir];
        pldir=pldirlist[dir];
        plforflux = B1+pldir-1;
 

        if(Nvec[fluxdir]>1){


          // first set which quantities to deal with (fields only)
          PALLLOOP(pl) whichpltoavg[pl]=0;
          PALLLOOP(pl) ifnotavgthencopy[pl]=0;
 
          // now choose which pl to use (only 1 and SAME for a given dir)
          whichpltoavg[plforflux]=do_transverse_flux_integration[plforflux];
          ifnotavgthencopy[plforflux]=0; // already exists as some value, so no need to copy over into self


          // go ahead and control nprlist (changes for each dir)
          // needed for Sasha's limit_fluxc2a
          addremovefromnpr(REMOVEFROMNPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);

          ////////////////////////////////////////////
          //
          // if no variables to do, then continue to next direction or if no more directions just done
          //
          ////////////////////////////////////////////
   
          if(nprlocalend>-1){


            // get loop details
            idel1 = fluxloop[fluxdir][FIDEL];
            jdel1 = fluxloop[fluxdir][FJDEL];
            kdel1 = fluxloop[fluxdir][FKDEL];

            idel2 = fluxloop[pldir][FIDEL];
            jdel2 = fluxloop[pldir][FJDEL];
            kdel2 = fluxloop[pldir][FKDEL];


            // store original flux
            // limit this de-averaging process
            if(extralimiting){
              COMPFULLLOOP GLOBALMACP0A1(fluxvectemp,i,j,k,plforflux) = MACP1A1(fluxvec,fluxdir,i,j,k,plforflux);
            }

            ///////////////////////////
            //
            // Obtain quasi-deaveraged EMFs or vector potential
            //
            // The E_{dir} or A_{dir} are quasi-deaveraged in both orthogonal directions
            //
            // In this case divb=0 is maintained when using EMF update for quasi-deaveraged fields and vector potential defines quasi-deaveraged fields
            //
            // Then the field is some quasi "conserved" magnetic flux
            //
            // To maintain higher-order fields, they need to be re-averaged (c2a) along the field direction to obtain point value at that FACE to be used for flux calculation
            // GODMARK: Re-deaveraging of the flux (a2c) should lead to the same quasi-conserved field in some limits, but for WENO c2a(a2c(f))!=f, so this is a problem as in FV method.
            //
            ///////////////////////////

            //    dualfprintf(fail_file,"dir=%d fluxdir=%d pldir=%d plforflux=%d :: idel1=%d jdel1=%d kdel1=%d ::  idel2=%d jdel2=%d kdel2=%d :: \n",dir,fluxdir,pldir,plforflux,idel1,jdel1,kdel1,idel2,jdel2,kdel2);

            numdirs=2;
            // fluxdirlist is already set correctly upon input
            // dirmethodlist and intdirlist should be same for FLUXRECON

   
            otherdir=1; dirmethodlist[otherdir]=fluxdir; whichquantitylist[otherdir]=whichquantity; interporfluxlist[otherdir]=interporflux; fluxdirlistnew[otherdir]=fluxdir; intdirlist[otherdir]=fluxdir; idellist[otherdir]=idel1; jdellist[otherdir]=jdel1; kdellist[otherdir]=kdel1;
            otherdir=2; dirmethodlist[otherdir]=pldir;   whichquantitylist[otherdir]=whichquantity; interporfluxlist[otherdir]=interporflux; fluxdirlistnew[otherdir]=fluxdir; intdirlist[otherdir]=pldir;   idellist[otherdir]=idel2; jdellist[otherdir]=jdel2; kdellist[otherdir]=kdel2;
   
            flux_interp_multiple(whichpltoavg, ifnotavgthencopy, numdirs, whichquantitylist, interporfluxlist, dirmethodlist, Nvec, intdirlist, fluxdirlistnew, idellist, jdellist, kdellist, pr, fluxvec);


            // limit changes between original flux and new (quasi-deaveraged) flux
            // 1 below means doing both transverse directions (if they exist)
            if(extralimiting){
              simple_a2c_limit_onepl(1,fluxdir, plforflux, pr, GLOBALPOINT(fluxvectemp), fluxvec[fluxdir]);
            }
 
          }// end if something to do
 
          // restore list
          addremovefromnpr(RESTORENPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


          // now copy over the other flux for same emf as required
          // assume emf's that should be zero are zero already
          if(isemf) if(Nvec[fluxdir]!=1 && Nvec[pldir]!=1) COMPFULLLOOP MACP1A1(fluxvec,pldir,i,j,k,B1+fluxdir-1) = - MACP1A1(fluxvec,fluxdir,i,j,k,B1+pldir-1);


        }// end if fluxdir exists, as required for fluxvec[fluxdir] to be used

      }// end if orthogonal odir1 or odir2 dimensions exist
    }// end over emf/vpot directions


  }// end if dotransem


  return(0);

}






int debug_boundfluxinitial(FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int i,j,k,pl,pliter;

  // GODMARK : test of bound_flux()
  if(nstep==30){
    // test of bound_flux().  Appears to work correctly.
    //
    // read in SM:
    // first add header from dump file and make size of grid as FULLLOOP zones
    //
    // jrdpheader3d 0_fail.out
    // da dumps/0_fail.out lines 2 1000000
    // read '%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g' {i j k f10 f20 f11 f21 f12 f22 f13 f23 f14 f24 f15 f25 f16 f26 f17 f27}
    // set r=4+i set h=j set ph=k gsetup gammienew
    //
    FULLLOOP{
      dualfprintf(fail_file,"%d %d %d ",i,j,k);
      PLOOP(pliter,pl) dualfprintf(fail_file,"%g %g ",MACP1A1(fluxvec,1,i,j,k,pl),MACP1A1(fluxvec,2,i,j,k,pl));
      dualfprintf(fail_file,"\n");
    }
    //  myexit(0);
  }

  return(0);

}

int debug_boundfluxfinal(FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int i,j,k,pl,pliter;


  // GODMARK : test of bound_flux()
  if(nstep==30){
    // test of bound_flux().  Appears to work correctly.
    //
    // read in SM:
    // first add header from dump file and make size of grid as FULLLOOP zones
    //
    // jrdpheader3d 0_fail.out
    // da dumps/0_fail.out lines 2 1000000
    // read '%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g' {i j k f10 f20 f11 f21 f12 f22 f13 f23 f14 f24 f15 f25 f16 f26 f17 f27}
    // set r=4+i set h=j set ph=k gsetup gammienew
    //
    FULLLOOP{
      logfprintf("%d %d %d ",i,j,k);
      PLOOP(pliter,pl) logfprintf("%g %g ",MACP1A1(fluxvec,1,i,j,k,pl),MACP1A1(fluxvec,2,i,j,k,pl));
      logfprintf("\n");
    }
    myexit(0);
  }


  return(0);

}










int bound_flux_fluxrecon(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pl,pliter;



  if(DOENOFLUX==ENOFLUXRECON && BOUNDFLUXRECON==1){

  
    ////////////////////////////////////////////
    //
    // save choice for PLOOPMPI
    nprlocalstart=nprboundstart;
    nprlocalend=nprboundend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=nprboundlist[pl];
    
    // bound the fluxes (only needed by this method)
    // set range of PLOOPMPI
    //  nprboundstart=0;
    //  nprboundend=NFLUXBOUND-1;
    //  for(pl=nprboundstart;pl<=nprboundend;pl++) nprboundlist[pl]=pl;

    /////////////////////
    //
    // now just assume same range as PLOOP and assume set right range of npr before calling bound so don't duplicate bound
    //
    ////////////////////
    nprboundstart=nprstart;
    nprboundend=nprend;
    for(pl=nprboundstart;pl<=nprboundend;pl++) nprboundlist[pl]=nprlist[pl];


    // ghost cells only need *additional* flux for FLUXRECON integration along direction itself
    // otherwise fluxes exist (for example) to update volume average quantities in ghost+active region
    // so this is not a good function to use to diagnose boundary issues if within ghost+active layer
    int finalstep=0;
    bound_flux(STAGEM1,finalstep,t,BOUNDFLUXTYPE,fluxvec[1],fluxvec[2],fluxvec[3], USEMPI);
    
#if(0)
    // trying to diagnose ghost+active issue with test=1102 and EVOLVEVPOT
    if(N1>1) bound_prim(-1,fluxvec[1]);
    if(N2>1) bound_prim(-1,fluxvec[2]);
    if(N3>1) bound_prim(-1,fluxvec[3]);
#endif

    
    ////////////////////////////////////////////
    //
    // restore PLOOPMPI
    nprboundstart=nprlocalstart;
    nprboundend=nprlocalend;
    PMAXNPRLOOP(pl) nprboundlist[pl]=nprlocallist[pl];

  }
  else{
    // then assuming using ghost+active region and then only have to set primitives
  }

  return(0);

}




// de-average point emf/vector potential at EDGE (CORN1,2,3)
// wrapper for emforvectorpot_fluxrecon()
// single input is treated as output, so no need to "copy" over results as some other functions do
// This function (for FV or FLUXRECON) is such that A has unique colocation with Flux and not multi-valued at end so only 1 flux is needed for temporary space
// So strictly the signforflux isn't needed since whatever put in cancels in end
int vectorpot_fluxreconorfvavg(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*F1)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F2)[NSTORE2][NSTORE3][NPR+NSPECIAL], FTYPE (*F3)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int emforvectorpot_fluxrecon(int stage, int isemf, int *fluxdirlist, int *pldirlist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int emforvectorpot_fvavg(int stage, int isemf, int *fluxdir, int *pldir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);

  int odir1,odir2;
  FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  int Nvec[NDIM];
  int fluxdir,pldir,plforflux;
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  FTYPE signforfluxlist[NDIM];
  FTYPE signforflux;
  int i,j,k;
  int dotransem;
  int dir;
  int pl,pliter;



  dotransem=
    do_transverse_flux_integration[B1]
    ||do_transverse_flux_integration[B2]
    ||do_transverse_flux_integration[B3];
  

  if(dotransem){


    Nvec[0]=0;
    Nvec[1]=N1;
    Nvec[2]=N2;
    Nvec[3]=N3;

    fluxvec[1] = F1;
    fluxvec[2] = F2;
    fluxvec[3] = F3;


    DIMENLOOP(dir) get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir], &signforfluxlist[dir]);

      
    // now copy over vector potential into flux
    // looping over vector potential directions (dir)
    DIMENLOOP(dir){
      fluxdir=fluxdirlist[dir];
      pldir=pldirlist[dir];
      plforflux=plforfluxlist[dir];
      signforflux=signforfluxlist[dir];

      // e.g. A3 will be treated as if F1[B2] or F2[B1]
      if(fluxdir!=0 && Nvec[fluxdir]>1){
        COMPFULLLOOP MACP1A1(fluxvec,fluxdir,i,j,k,plforflux)=signforflux*MACP1A0(A,dir,i,j,k);
      }

    }




#if(0)
    ///////////////////////////////
    // DEBUG:
    DIMENLOOP(dir){
      fluxdir=fluxdirlist[dir];
      pldir=pldirlist[dir];
      plforflux=plforfluxlist[dir];
      signforflux=signforfluxlist[dir];
    
      // e.g. A3 will be treated as if F1[B2] or F2[B1]
      if(fluxdir!=0 && Nvec[fluxdir]>1){
        bound_prim(STAGEM1,fluxvec[fluxdir]);
      
      }
    }
#endif



#if(FLUXDUMP==1)
    // this accounts for final flux
    FULLLOOP{
      DIMENLOOP(dir){
        if(Nvec[dir]>1){
          PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=MACP1A1(fluxvec,dir,i,j,k,pl);
        }
        else{
          PLOOP(pliter,pl) GLOBALMACP0A1(fluxdump,i,j,k,4*NPR + (dir-1)*NPR*5 + NPR*0 + pl)=0.0L;
        }
      }
    }
#endif



    // these 2 emforvectorpot functions don't care about plforfluxlist or signforfluxlist since just treat other flux with opposite sign (i.e. no direct relationship between sign of A_i and flux needed to be set inside function)
    if(DOENOFLUX==ENOFLUXRECON){
      // quasi-deaverage A_i
      emforvectorpot_fluxrecon(stage, 0, fluxdirlist, pldirlist, pr, Nvec, fluxvec);
    }
    else if(DOENOFLUX==ENOFINITEVOLUME){
      // 0 means not emf but vector potential and with fluxdir direction chosen for fluxvec
      emforvectorpot_fvavg(stage, 0, fluxdirlist, pldirlist, pr, Nvec, fluxvec);
    }

    // need to bound these fluxes so exist everywhere as requied for vector potential so it can get field in ghost zones.
    // --GODMARK not really needed since ghost zone values are determined by boundary conditions
    //bound_flux_fluxrecon(stage, pr, Nvec, fluxvec);



    // copy over vector potential that has been quasi-deaveraged
    DIMENLOOP(dir){
      fluxdir=fluxdirlist[dir];
      pldir=pldirlist[dir];
      plforflux=plforfluxlist[dir];
      signforflux=signforfluxlist[dir];

      // e.g. A3 was treated as if F1[B2] or F2[B1]
      if(fluxdir!=0 && Nvec[fluxdir]>1) COMPFULLLOOP MACP1A0(A,dir,i,j,k)=signforflux*MACP1A1(fluxvec,fluxdir,i,j,k,plforflux);

    }

  }



  return(0);

}












int flux_integrate_fluxsplit(int stage, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

  // GODMARK: NOTICE: this will probably go inside fluxcalc_fluxspliteno() because the procedure is
    
  // -1) no need to bound fluxes since SHOULD exist everywhere (should be taken care of by new fluxcalc_fluxspliteno() function)
  // 0) assume primitives everywhere defined (all grid including ghost zones)
  // 1) Compute split F's (F^+ F^-) for all grid including ghost zones, store in F1m/p F2m/p F3m/p
  // 2) use flux_interp(ENOFLUXSPLITTYPE,...) for each dir and get F1l/r F2l/r F3l/r
  // 3) add together fluxes and store in F1,F2,F3 (can be same memory space as used for F1m/p F2m/p F3m/p since done with it, but can't be same memory space as F123l/r
    
  //    dir=1; flux_interp(ENOFLUXSPLITTYPE, dir, idel, jdel, kdel, F1m, F1p, F1l, F1r);
  //    dir=2; flux_interp(ENOFLUXSPLITTYPE, dir, idel, jdel, kdel, F2m, F2p, F2l, F2r);
  //    dir=3; flux_interp(ENOFLUXSPLITTYPE, dir, idel, jdel, kdel, F3m, F3p, F3l, F3r);
  // these should be read in correct form as Fr Fl in fluxcalc_fluxspliteno() such that for i-th index:
  // dir=1
  // Fr(i) = F1l(i)
  // Fl(i) = F1r(i-1)

  return(0);

}












// FV method to surface (trnasverse to flux direction) integrate fluxes
// single input is treated as output, so no need to "copy" over results as some other functions do
int flux_integrate_finitevolume(int stage, int whichmaem, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int avgdirtype[NDIM],odir1,odir2;
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  extern void flux_interp_multiple(int *whichpltoavg, int *ifnotavgthencopy, int numdirs, int *whichquantitylist, int *interporfluxlist, int *dirmethodlist, int *Nvec, int *intdirlist, int *fluxdirlist, int *idellist, int *jdellist, int *kdellist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  FTYPE limit_fluxc2a_prim_change( 
                                  int dir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                                  FTYPE (*fluxvec_point)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                  FTYPE (*fluxvec_avg)[NSTORE2][NSTORE3][NPR+NSPECIAL]);   //atch
  FTYPE (*fluxveca[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (*fluxvecb[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  int pl,pliter,i,j,k;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int emforvectorpot_fvavg(int stage, int isemf, int *fluxdirlist, int *pldirlist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL]);
  int dotrans,dotransma,dotransem;
  int fluxdirlist[NDIM],pldirlist[NDIM],plforfluxlist[NDIM];
  FTYPE signforfluxlist[NDIM];
  int numdirs, intdirlist[NDIM], idellist[NDIM],jdellist[NDIM],kdellist[NDIM];
  int whichquantitylist[NDIM],interporfluxlist[NDIM],dirmethodlist[NDIM];
  int otherdir;
  int fluxdir;
  int treatemfspecial;
  int whichquantity, weightsplittype,extralimiting;


  // any trans
  dotrans=0;
  PLOOP(pliter,pl) dotrans+=do_transverse_flux_integration[pl];
  dotrans*=(DOENOFLUX == ENOFINITEVOLUME);


  if(whichmaem==ISMAANDEM){
    whichquantity=ENOFLUX;
  }
  else if(whichmaem==ISMAONLY){
    whichquantity=ENOMAFLUX;
  }
  else if(whichmaem==ISEMONLY){
    whichquantity=EMFLUXFIELDTERMTYPE;
  }

  higherorder_set(whichquantity, CVT_C2A, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_FLUXC2A_PRIM_CHANGE) extralimiting=1;
    else extralimiting=0;
  }




  if(dotrans){


    dotransma=DOENOFLUX == ENOFINITEVOLUME&&
      (do_transverse_flux_integration[RHO]
       ||do_transverse_flux_integration[UU]
       ||do_transverse_flux_integration[U1]
       ||do_transverse_flux_integration[U2]
       ||do_transverse_flux_integration[U3]);

    dotransem=DOENOFLUX == ENOFINITEVOLUME&&
      (do_transverse_flux_integration[B1]
       ||do_transverse_flux_integration[B2]
       ||do_transverse_flux_integration[B3]);


    avgdirtype[1]=ENOFLUXAVG1TYPE;
    avgdirtype[2]=ENOFLUXAVG2TYPE;
    avgdirtype[3]=ENOFLUXAVG3TYPE;


    // only use special FV emf method for staggered field method.
    // Assume all other methods use centered field and regular flux for field such that EMF is not treated specially
    if(FLUXB==FLUXCTSTAG){
      treatemfspecial=1;
    }
    else treatemfspecial=0;
    


    ///////////////////////////////////////////////////////////
    //
    // integrate point fluxes to surface averaged fluxes
    //
    ///////////////////////////////////////////////////////////

    if(dotrans){


      // then simpler to have separate "emf" section, so remove fields for all calculations
      // default is to average all fluxes
      PALLLOOP(pl) whichpltoavg[pl]=do_transverse_flux_integration[pl];// default
      PALLLOOP(pl) ifnotavgthencopy[pl]=0; // already exists as some value so no need to copy over to self
      if(treatemfspecial){
        // never do field since always do field with specially located emf below
        PLOOPBONLY(pl) whichpltoavg[pl]=0;
        PLOOPBONLY(pl) ifnotavgthencopy[pl]=0;
      }
      else{
        PLOOPBONLY(pl) whichpltoavg[pl]=do_transverse_flux_integration[pl];
        PLOOPBONLY(pl) ifnotavgthencopy[pl]=0;
      }
    
      // go ahead and control nprlist
      // needed for Sasha's limit_fluxc2a (flux_interp internally would remove unwanted pl's based upon whichpltoavg/ifnotavgthencopy)
      addremovefromnpr(REMOVEFROMNPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);
    

      ////////////////////////////////////////////
      //
      // if no variables to do, then continue to next direction or if no more directions just done
      //
      ////////////////////////////////////////////
      
      if(nprlocalend>-1){

    
        ///////////////////////////////////////////////////////////
        //
        // LOOP OVER FACES
        //
        ///////////////////////////////////////////////////////////
        DIMENLOOP(dir){
      
          // other dimensions
          odir1=dir%3+1;
          odir2=(dir+1)%3+1;



      
          if(Nvec[dir]!=1){
            // skip to next dir if no such dimension

            if(extralimiting){
              //otherwise fluxvectemp points to space in memory set in set_arrays.c
              //intitialize it with the flux given by the riemann solver in the dir direction
              COMPFULLLOOP PLOOP(pliter,pl) GLOBALMACP0A1(fluxvectemp,i,j,k,pl) = MACP1A1(fluxvec,dir,i,j,k,pl);
            }

            // get loop details
            idel1 = fluxloop[odir1][FIDEL];
            jdel1 = fluxloop[odir1][FJDEL];
            kdel1 = fluxloop[odir1][FKDEL];

            idel2 = fluxloop[odir2][FIDEL];
            jdel2 = fluxloop[odir2][FJDEL];
            kdel2 = fluxloop[odir2][FKDEL];


            numdirs=2;
            fluxdir=dir;
            otherdir=1; intdirlist[otherdir]=odir1; dirmethodlist[otherdir]=dir; whichquantitylist[otherdir]=whichquantity; interporfluxlist[otherdir]=avgdirtype[odir1]; fluxdirlist[otherdir]=fluxdir; idellist[otherdir]=idel1; jdellist[otherdir]=jdel1; kdellist[otherdir]=kdel1;
            otherdir=2; intdirlist[otherdir]=odir2; dirmethodlist[otherdir]=dir; whichquantitylist[otherdir]=whichquantity; interporfluxlist[otherdir]=avgdirtype[odir2]; fluxdirlist[otherdir]=fluxdir; idellist[otherdir]=idel2; jdellist[otherdir]=jdel2; kdellist[otherdir]=kdel2;
   
            flux_interp_multiple(whichpltoavg, ifnotavgthencopy, numdirs, whichquantitylist, interporfluxlist, dirmethodlist, Nvec, intdirlist, fluxdirlist, idellist, jdellist, kdellist, pr, fluxvec); 


            if(extralimiting){
              // GODMARK: above parameter is off normally since Sasha says his method didn't work (caused more problems)
              //limits how different fluxvec[dir] (averaged fluxes) is from fluxvectemp (point fluxes) are based on how much rho and gamma would
              //change due to such a difference in flux during the timestep, dt
              limit_fluxc2a_prim_change( dir, pr, GLOBALPOINT(fluxvectemp), fluxvec[dir] );
            }


          }// end if F_{dir} does not exist since not doing dimension in dir-direction
        } // end over DIMENLOOP
      } // end if something to do

      // restore list
      addremovefromnpr(RESTORENPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


    }
      


    if(dotransem && treatemfspecial){
      // line-integrate EMF
      // 1 means is an emf
      DIMENLOOP(dir) get_fluxpldirs(Nvec, dir, &fluxdirlist[dir], &pldirlist[dir], &plforfluxlist[dir],&signforfluxlist[dir]);
      // emforvectorpot_fvavg() doesn't need to set signature between A_i and fluxvec
      emforvectorpot_fvavg(stage, 1, fluxdirlist, pldirlist, pr, Nvec, fluxvec);
    }

  }// end if dotrans


  return(0);

}







// integrate point emf/vector potential at EDGE (CORN1,2,3) to get line-integrated value
// input in fluxvec form
// single input is treated as output, so no need to "copy" over results as some other functions do
int emforvectorpot_fvavg(int stage, int isemf, int *fluxdirlist, int *pldirlist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int avgdirtype[NDIM],odir1,odir2;
  int dir;
  int idel, jdel, kdel, face, idel1, idel2, jdel1, jdel2, kdel1, kdel2;
  int is, ie, js, je, ks, ke;
  extern void flux_interp(int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);
  FTYPE limit_fluxc2a_prim_change( 
                                  int dir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                                  FTYPE (*fluxvec_point)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                  FTYPE (*fluxvec_avg)[NSTORE2][NSTORE3][NPR+NSPECIAL]);   //atch
  int pl,pliter,i,j,k;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int fluxdir,pldir, plforflux;
  int orthogonaldir;
  int dotransem;
  int whichquantity,whichmaem,weightsplittype,extralimiting;
  



  dotransem=
    do_transverse_flux_integration[B1]
    ||do_transverse_flux_integration[B2]
    ||do_transverse_flux_integration[B3];


  whichmaem=ISEMONLY;
  whichquantity=EMFLUXFIELDTERMTYPE;

  higherorder_set(whichquantity, CVT_C2A, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_FLUXC2A_PRIM_CHANGE) extralimiting=1;
    else extralimiting=0;
  }





  if(dotransem){


    // definitions:
    avgdirtype[1]=ENOFLUXAVG1TYPE;
    avgdirtype[2]=ENOFLUXAVG2TYPE;
    avgdirtype[3]=ENOFLUXAVG3TYPE;



    ////////////////////////
    //
    // now deal with FLUXCTSTAG (no need for split or non-split integration since only 1 direction to integrate per emf
    //
    ////////////////////////

    // now do emf's assuming dir means emf direction (EDGES)
    // dir here is direction of averaging
    // fluxdir is which flux direction
    // plforflux is which quantity to use for that fluxdir
    DIMENLOOP(dir){


      // other dimensions
      odir1=dir%3+1;
      odir2=(dir+1)%3+1;
 

      // do emf[dir]
      if(!(Nvec[odir1]==1 && Nvec[odir2]==1)){ // then emf[dir] is unimportant (cancels itself always), so assume original emf[dir] was set to 0 already
   
        // then doing emf[dir]
        fluxdir=fluxdirlist[dir];
        pldir=pldirlist[dir];
        plforflux=B1+pldir-1;
        orthogonaldir=fluxdir;


        // first set which quantities to deal with (fields only)
        PALLLOOP(pl) whichpltoavg[pl]=0;
        PALLLOOP(pl) ifnotavgthencopy[pl]=0;
        // now choose which pl to use
        whichpltoavg[plforflux]=do_transverse_flux_integration[plforflux];
        ifnotavgthencopy[plforflux]=0; // already exists as some value so no need to copy over to self

        // go ahead and control nprlist (changes for each dir)
        // needed for Sasha's limit_fluxc2a
        addremovefromnpr(REMOVEFROMNPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);

        ////////////////////////////////////////////
        //
        // if no variables to do, then continue to next direction or if no more directions just done
        //
        ////////////////////////////////////////////

        if(nprlocalend>-1){


          // save the point flux for comparison
          if(isemf){ // otherwise vectorpotential which is not flux for conserved quantity
            if(extralimiting){
              // only need 1 flux value per dir
              COMPFULLLOOP GLOBALMACP0A1(fluxvectemp,i,j,k,plforflux) = MACP1A1(fluxvec,fluxdir,i,j,k,plforflux);
            }
          }


          // as consistent with above normal fluxes, integration directly is controlled via these "del"'s, and integration direction is in dir-direction for emf
          idel = fluxloop[dir][FIDEL];
          jdel = fluxloop[dir][FJDEL];
          kdel = fluxloop[dir][FKDEL];

          //      dualfprintf(fail_file,"isemf: %d %d %d %d %d %d\n",isemf,vpotfluxdir,dir,avgdirtype[dir],orthogonaldir,fluxdir);

#if(0) // INCOMPLETE
          // GODMARK: Should be able to use Taylor series expansion or something to be able to generalize this to 3D
          // Here we subtract off the emf term related to divb=0 constraint assuming 2D problem
          if(isemf){
            if(dir==2 && Nvec[3]==1){ // then E2 with term $\int_y v^z B^x dy$ can have vz0 subtracted out
              if(fluxdir==3 && pldir==1){ // then E2=F3[B1] has B^1 v^3 term, so remove it
                //     COMPFULLLOOP MACP1A1(fluxvec,fluxdir,i,j,k,plforflux) -= (
              }
            }
          }

#endif



          // 3rd argument "dir" corresponds to direction of averaging (i.e. along dir)
          // 4th argument "fluxdir" corresponds to orthogonal direction
          if(Nvec[dir]>1 && Nvec[fluxdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, avgdirtype[dir], orthogonaldir, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxvec[fluxdir], NULL);



          if(isemf){
            if(extralimiting){
              // below controlled by PLOOP, so only use over range changed fluxes
              limit_fluxc2a_prim_change( fluxdir, pr, GLOBALPOINT(fluxvectemp), fluxvec[fluxdir] );
            }
          }

        }//end if something to do

        // go ahead and control nprlist
        // needed for Sasha's limit_fluxc2a
        addremovefromnpr(RESTORENPR,whichpltoavg,ifnotavgthencopy,&nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);


        // now copy over the other flux for same emf as required
        // assume emf's that should be zero are zero already
        if(isemf) if(Nvec[fluxdir]!=1 && Nvec[pldir]!=1) COMPFULLLOOP MACP1A1(fluxvec,pldir,i,j,k,B1+fluxdir-1) = - MACP1A1(fluxvec,fluxdir,i,j,k,B1+pldir-1);


      }// end if doing this emf[dir]
    } // end over DIMENLOOP


  }// end if transem      
      


  return(0);

}





// FV method to deaverage surface-averaged fields to get point field at FACE (used for IC)
// GODMARK: this current doesn't choose limiting since 
int deaverage_fields_fv(FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR])
{
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int doconsem;
  int locpl[NPR];
  extern int avg2cen_interp(int *locpl, int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int whichavg2cen, FTYPE (*prims_from_avg_cons)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  int simple_a2c_limit_field(FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  int pl,pliter;
  int i,j,k;
  int whichmaem,whichquantity,weightsplittype,extralimiting;
  int interporflux;


  doconsem=(do_conserved_integration[B1]||do_conserved_integration[B2]||do_conserved_integration[B3]);

  whichmaem=ISEMONLY;
  whichquantity=FIELDFVTYPE;
  interporflux=ENOAVG2CENTTYPE;

  higherorder_set(whichquantity, CVT_A2C, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_AC_FRAC_CHANGE) extralimiting=1;
    else extralimiting=0;
  }



    
  // de-average initial fields like in advance.c finitevolume method
  if(doconsem){
    PALLLOOP(pl) whichpltoavg[pl]=0;// default
    PALLLOOP(pl) ifnotavgthencopy[pl]=0;// default
    PLOOPBONLY(pl) whichpltoavg[pl]=do_conserved_integration[pl];
    PLOOPBONLY(pl) ifnotavgthencopy[pl]=1-do_conserved_integration[pl]; //  need to copy over if not operating on in->out
    if(FLUXB==FLUXCTSTAG){
      locpl[B1]=FACE1;
      locpl[B2]=FACE2;
      locpl[B3]=FACE3;
    }
    else{
      locpl[B1]=CENT;
      locpl[B2]=CENT;
      locpl[B3]=CENT;
    }


    avg2cen_interp(locpl, whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, primreal, in, out);


    //////////////////////
    //
    // assume if user wants to limit change in field, then out and in are separate memory
    //
    /////////////////////

    if(in!=out && extralimiting){
      simple_a2c_limit_field(primreal,in,out);
    }


  }   
  else{
    COMPFULLLOOP PLOOP(pliter,pl) MACP0A1(out,i,j,k,pl)=MACP0A1(in,i,j,k,pl);
    // ulast now has point value of field at staggered position
  }

  
  return(0);
 
}




// limit a2c for field only
// assume de-averaged in transverse direction to direction of field
int simple_a2c_limit_field(FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;
  int dir;
  int simple_a2c_limit_onepl(int dotransverse, int dir, int pl, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  int dotransverse;


  dotransverse=1;
  
  // then limit change since deaveraging can create noise
  PLOOPBONLY(pl){
    
    // direction of field
    dir=pl-B1+1;
    
    simple_a2c_limit_onepl(dotransverse, dir, pl, primreal, in, out);

  }


  return(0);


}





#define OUTEND (MAX_AC_FRAC_CHANGE)
#define MIXEND (2.0*MAX_AC_FRAC_CHANGE)

// small section loop
#define LOOPSURF(i,j,k) for(ii=i+myim1;ii<=i+myip1;ii++) for(jj=j+myjm1;jj<=j+myjp1;jj++) for(kk=k+mykm1;kk<=k+mykp1;kk++)

// limiting for transverse de-averaging to dir-direction
// limit a2c for one quantity
int simple_a2c_limit_onepl(int dotransverse, int dir, int pl, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR])
{
  FTYPE xfrac,yfrac;
  int i,j,k;
  int myim1,myip1,myjm1,myjp1,mykm1,mykp1;
  FTYPE ref;
  int ii,jj,kk;
  int numref;
  FTYPE diff;
  FTYPE myin;


  if(dotransverse){
    if(dir==1){
      myip1=0;
      myim1=0;
      myjp1=SHIFT2;
      myjm1=-SHIFT2;
      mykp1=SHIFT3;
      mykm1=-SHIFT3;
    }
    else if(dir==2){
      myip1=SHIFT1;
      myim1=-SHIFT1;
      myjp1=0;
      myjm1=0;
      mykp1=SHIFT3;
      mykm1=-SHIFT3;
    }
    else if(dir==3){
      myip1=SHIFT1;
      myim1=-SHIFT1;
      myjp1=SHIFT2;
      myjm1=-SHIFT2;
      mykp1=0;
      mykm1=0;
    }
    else{
      dualfprintf(fail_file,"No such dir=%d in simple_a2c_limit_onepl\n",dir);
      myexit(18619626);
    }
  }
  else{ // then assume just along dir
    if(dir==1){
      myip1=SHIFT1;
      myim1=-SHIFT1;
      myjp1=0;
      myjm1=0;
      mykp1=0;
      mykm1=0;
    }
    else if(dir==2){
      myip1=0;
      myim1=0;
      myjp1=SHIFT2;
      myjm1=-SHIFT2;
      mykp1=0;
      mykm1=0;
    }
    else if(dir==3){
      myip1=0;
      myim1=0;
      myjp1=0;
      myjm1=0;
      mykp1=SHIFT3;
      mykm1=-SHIFT3;
    }
    else{
      dualfprintf(fail_file,"No such dir=%d in simple_a2c_limit_onepl\n",dir);
      myexit(18619627);
    }
  }



    
  // now loop over positions
  COMPFULLLOOP{

    // get average of averaged value around that was used to do a2c
    ref = 0.0;
    numref = 0;
    LOOPSURF(i,j,k){
      ref+=fabs(MACP0A1(in,ii,jj,kk,pl));
      numref++;
    }
    ref/=( (FTYPE)numref); // now ref is absolute average value over surface used to do a2c

    myin=MACP0A1(in,i,j,k,pl);

    // difference between average and point value
    diff=MACP0A1(out,i,j,k,pl)-myin;

    // see if diff is comparable to reference value
    xfrac = fabs(diff)/(fabs(ref)+SMALL);
    if(xfrac>1.0) xfrac=1.0;
    if(xfrac<0.0) xfrac=0.0;

    // now use interpolation between average and point depending upon how different from reference value
    if(xfrac<OUTEND){
      yfrac=1.0;
    }
    else if(xfrac<MIXEND){
      yfrac = 1.0 - (xfrac-OUTEND)/(MIXEND-OUTEND);
    }
    else{
      yfrac=0.0;
    }

    // DEBUG:
    //dualfprintf(fail_file,"i=%d j=%d k=%d :: xfrac=%21.15g yfrac=%21.15g\n",i,j,k,xfrac,yfrac);
      
    // obtain final result per point
    MACP0A1(out,i,j,k,pl) = yfrac*MACP0A1(out,i,j,k,pl) + (1.0-yfrac)*myin;


  }// loop over positions i,j,k
    



  return(0);


}










// FV method for initial averaging of point conserved quantities to get volume-averaged conserved quantities
int initial_averageu_fv(int *fieldfrompotential, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*Upoint)[NSTORE2][NSTORE3][NPR], FTYPE (*Uavg)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  int pl,pliter;
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int docons;
  int locpl[NPR];
  extern int avg2cen_interp(int *locpl, int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux,FTYPE (*prims_from_avg_cons)[NSTORE2][NSTORE3][NPR], FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR]);
  int whichmaem,whichquantity,weightsplittype,extralimiting;
  int interporflux;


  // any cons
  docons=0;
  PLOOP(pliter,pl) docons+=do_conserved_integration[pl];
  docons*=(DOENOFLUX == ENOFINITEVOLUME);

  whichmaem=ISMAANDEM;
  whichquantity=ENOCONSERVED;
  interporflux=ENOCENT2AVGTYPE;

  higherorder_set(whichquantity, CVT_C2A, &weightsplittype);
  if(weightsplittype==CONSTANT_ALL_WEIGHTS || weightsplittype==MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS) extralimiting=0;
  else{
    if(LIMIT_AC_FRAC_CHANGE) extralimiting=1;
    else extralimiting=0;
  }



  // volume integrate initial Upoint
  // GODMARK: in reality should have analytical solution for this
  // GODMARK: causes smoothing of initial condition in incorrect way in general
  // 00000000000 && because CENT2AVG is inaccurate in general -- although required for convergence study -- SUPERGODMARK

  //if(0000000000 && doconsma ) {
  if(docons) {
    // c2a_1 c2a_2 c2a_3
    PALLLOOP(pl) whichpltoavg[pl]=do_conserved_integration[pl];// default
    PALLLOOP(pl) ifnotavgthencopy[pl]=1-do_conserved_integration[pl];// default
    PALLLOOP(pl) locpl[pl]=CENT;
    PLOOPBONLY(pl) if(fieldfrompotential[pl-B1+1]) whichpltoavg[pl]=0;
    PLOOPBONLY(pl) if(fieldfrompotential[pl-B1+1]) ifnotavgthencopy[pl]=0;

    PALLLOOP(pl){
      trifprintf("initial averaging: whichpltoavg[%d]=%d ifnotavgthencopy[%d]=%d\n",pl,whichpltoavg[pl],pl,ifnotavgthencopy[pl]);
    }

    avg2cen_interp(locpl, whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, prim, Upoint, Uavg);


  }
  else {
    if(1 || FLUXB==FLUXCTSTAG || FLUXB==FLUXCTTOTH){ //(1 || because currently all FLUXB methods set Uavg before pi2Uavg either through vpot or in init.c's  init_conservatives() before the pi2Uavg call)
      // can't overwrite staggered conserved field that MUST be set
      COMPFULLLOOP{
        PLOOPNOB1(pl) MACP0A1(Uavg,i,j,k,pl) = MACP0A1(Upoint,i,j,k,pl);
        PLOOPNOB2(pl) MACP0A1(Uavg,i,j,k,pl) = MACP0A1(Upoint,i,j,k,pl);
      }
    }
    else{
      COMPFULLLOOP PLOOP(pliter,pl) {
        MACP0A1(Uavg,i,j,k,pl) = MACP0A1(Upoint,i,j,k,pl);
      }
    }
  }


  return(0);


}














// Get primitive for a2c/c2a operations
int get_primitive_centerlocation(int *locpl, int *whichpltoavg, int interporflux, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR], FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;
  //Utoprimgen is already declared in global.h
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1; 


#if(0)
  // presently only thing that uses this is FV method for a2c for getting point field from average field
  // assume ok to use CENT for all -- GODMARK
  PALLLOOP(pl){
    if(whichpltoavg[pl]){
      if(locpl[pl]!=CENT){
        dualfprintf(fail_file,"avg2cen_interp() not setup for non-CENT locations\n");
        myexit(1659562);
      }
    }
  }
#endif
 

  // 1|| since should always be doing since always need primitive
  if(1|| avgscheme[1] == WENO5FLAT ||avgscheme[2] == WENO5FLAT ||avgscheme[3] == WENO5FLAT || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ) {

    if(USEPRIMITIVEFROMAVGCONSERVED){
      int eomtype=EOMDEFAULT;

      COMPFULLLOOP{
        get_geometry(i, j, k, CENT, ptrgeom);
        // set guess
        PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl)=MACP0A1(primreal,i,j,k,pl);
        FTYPE dissmeasure=-1.0; // assume energy try ok
        int whichcap=CAPTYPEBASIC;
        int whichmethod=MODEDEFAULT;
        int modprim=0;
        int checkoninversiongas=CHECKONINVERSION;
        int checkoninversionrad=CHECKONINVERSIONRAD;
        MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM, UEVOLVE, MAC(U,i,j,k), NULL, ptrgeom, dissmeasure, MAC(prim_goodlocation,i,j,k), MAC(prim_goodlocation,i,j,k),&newtonstats),"interpline.c:avg2cen_interp()", "Utoprimgen", 1);
        // if problem with inversion, then reduce to using primreal
        if(GLOBALMACP0A1(pflag,i,j,k,FLAGUTOPRIMFAIL)) PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl)=MACP0A1(primreal,i,j,k,pl);
        // pflag will be reset by real inversion routine in advance.c before used elsewhere
      }

    }
    else{
      // then just use primreal, which is offset in time when used on updated average conserved quantity, but spatially in the right place
      // probably very good indicator of the future behavior

      if( primreal == NULL ) {
        //no primitive quantities values supplied -> obtain them by treating input as conserved quantities and inverting them
        dualfprintf( fail_file, "interpline.c: avg2cen_interp: WENO5FLAT requires supplying the primitive quantities\n" );
        myexit( 1 );
      }

      // GODMARK: for FLUXCTSTAG offset in transverse direction for non-dir fields that are de-averaged
      COMPFULLLOOP  PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl) = MACP0A1(primreal,i,j,k,pl);
      //      prim_goodlocation = primreal; // could just assign pointer
    }
  }
  else {
    //do not compute the primitive quantities that correspond to conserved ones because they will not be used; the case with NULL should be set up to be handled properly
    prim_goodlocation = NULL;
  }



  return(0);

}




// get primitive for c2a flux-type operations
int get_primitive_fluxlocation(int dir, int interporflux, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k,pl,pliter;

  // GODMARK -- assume for all locations fluxlocation is FACE_dir


  // 1|| because should always have this even if not doing WENOFLAT
  if(1|| avgscheme[dir] == WENO5FLAT || CONTACTINDICATOR || COMPUTEDRHODP|| SHOCKINDICATOR ) {

    // 1|| because should always have this even if not doing WENOFLAT
    if(1||USEAVGPRIMITIVEFORWENOFLAT && interporflux!=ENOFLUXSPLITTYPE){
      if(dir==1){
        // need to average primitive in dir=1 direction so in same location
        // GODMARK: loop is too big, but data isn't accessed later in bad region
        COMPLOOPINFP1dir23full PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl) = 0.5*(MACP0A1(primreal,i,j,k,pl)+MACP0A1(primreal,im1mac(i),j,k,pl));
      }
      else if(dir==2){
        COMPLOOPINFP1dir13full PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl) = 0.5*(MACP0A1(primreal,i,j,k,pl)+MACP0A1(primreal,i,jm1mac(j),k,pl));
      }
      else if(dir==3){
        COMPLOOPINFP1dir12full PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl) = 0.5*(MACP0A1(primreal,i,j,k,pl)+MACP0A1(primreal,i,j,km1mac(k),pl));
      }
      else{
        dualfprintf(fail_file,"No such dir=%d in flux_interp()\n",dir);
        myexit(100);
      }
    }
    else{
      // flux split is zone centered
      //      prim_goodlocation=pr; // could also loop as above if want to protect pr
      COMPFULLLOOP PALLLOOP(pl) MACP0A1(prim_goodlocation,i,j,k,pl) = MACP0A1(primreal,i,j,k,pl);
    }
  }
  else {
    //do not compute the primitive quantities that correspond to conserved ones because they will not be used; the case with NULL should be set up to be handled properly
    prim_goodlocation = NULL;
  }




  return(0);

}













// flux_interp() is provided p2interp and returns pleft/pright
//
// |=interface
// i=zone center of ith zone
//
// |              |p2interpm(i) if interporflux==ENOFLUXSPLITTYPE
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//

// dir here is each direction for each flux (F1,F2,F3)
// does a 1D integration or differentiation
void flux_interp(int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR])
{

  void slope_lim_linetype(int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prims_guess)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);

  int Nvec[NDIM];
  int i, j, k, pl, pliter;

  //use ptemparray as a temp array
  FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR] = GLOBALPOINT(ptemparray);
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int get_primitive_fluxlocation(int dir, int interporflux, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR]);




  // get reference primitive but in correct location
  get_primitive_fluxlocation(dir, interporflux, primreal, prim_goodlocation);


  //////////////////////////////
  //
  // Perform interpolation
  //
  //////////////////////////////

  // remove fields should not average or deaverage
  // assume never copying if not averaging
  addremovefromnpr(REMOVEFROMNPR, whichpltoavg, ifnotavgthencopy, &nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);
  
  ////////////////////////////////////////////
  //
  // if no variables to do, then continue to next direction or if no more directions just done
  //
  ////////////////////////////////////////////
  
  if(nprlocalend>-1){
    slope_lim_linetype(whichquantity, interporflux, dir, idel, jdel, kdel, prim_goodlocation, stencilvar, p2interpm, p2interpp, pleft, pright);
  }

  // restore dir-field
  addremovefromnpr(RESTORENPR, whichpltoavg, ifnotavgthencopy, &nprlocalstart,&nprlocalend,nprlocallist,NULL,NULL);



}







// dir here is each direction for each flux (F1,F2,F3)
// does a 1D integration or differentiation for multiple directions in a way that preserves commutivity of directions by using a single original quantity to compute weights
void flux_interp_multiple(int *whichpltoavg, int *ifnotavgthencopy, int numdirs, int *whichquantitylist, int *interporfluxlist, int *dirmethodlist, int *Nvec, int *intdirlist, int *fluxdirlist, int *idellist, int *jdellist, int *kdellist, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int i, j, k, pl, pliter;
  int intdir,fluxdir,idel,jdel,kdel,interporflux,whichquantity,dirmethod;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  //use ptemparray as a temp array
  int otherdir;
  FTYPE (*fluxveca[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (*fluxvecb[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL];
  FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR];
  void flux_interp(int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);
  int dodouble;



  // only used for non-PERFECTUNSPLIT methods that assume only doing 2 directions
  dodouble=(Nvec[intdirlist[1]]>1) && (Nvec[intdirlist[2]]>1);


  if(FLUXDIMENSPLIT==PERFECTUNSPLIT){
    // assumes stuff in interpline is setup to use 1 quantity for weights for both directions

    // storage for quantity to use for determining stencil
    stencilvar=GLOBALPOINT(stencilvartemp);

    for(otherdir=1;otherdir<=numdirs;otherdir++){

      whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];

      if(numdirs>1){
        // then need to store original to use for weights
        COMPFULLLOOP PLOOP(pliter,pl) MACP0A1(stencilvar,i,j,k,pl) = MACP1A1(fluxvec,fluxdir,i,j,k,pl); // assumes not flux splitting method since only one input and one output
      }
      else stencilvar=NULL; // says use normal input for stencil

      // FLUXRECON intdir=dirmethod=fluxdir/pldir
      if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, stencilvar, fluxvec[fluxdir], NULL, fluxvec[fluxdir], NULL); 
    }


  }
  else if( FLUXDIMENSPLIT==QUASISTRANG || (FLUXDIMENSPLIT==UNSPLIT && dodouble==0) ){
    // rest assume always 2 directions

    // Quasi-Strang Splitting (alternate order of dimension) of ordering of c2a integration of flux over 2 orthogonal directions
    // shouldn't start with pldir if pldir direction is dimensionless and so not done
    if((nstep*(long long)numstepparts+(long long)steppart)%2 || dodouble==0){


      //       dualfprintf(fail_file,"nstep=%ld steppart=%d :: dir=%d fluxdir=%d pldir=%d plforflux=%d\n",nstep,steppart,dir,fluxdir,pldir,plforflux);
      //       dualfprintf(fail_file,"nprstart=%d nprend=%d\n",nprstart,nprend);
      //       PALLLOOP(pl) dualfprintf(fail_file,"nprlist[%d]=%d\n",pl,nprlist[pl]);

      // over at most 2 orthogonal directions
      // GODMARK: ordering issue.  If a2c operations commuted, then would not matter.  How to ensure this?
       
      // first do integration along flux dir, then along pldir
      // only change integration direction, not which quantity (e.g. emf[3]=+-F2[B1] or F1[B2] -- only do 1 both times)
      otherdir=1;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
      if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxvec[fluxdir], NULL); 
       
      otherdir=2;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
      if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxvec[fluxdir], NULL); 
    }
    else{
      otherdir=2;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
      if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxvec[fluxdir], NULL); 
       
      otherdir=1;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
      if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxvec[fluxdir], NULL); 
    }
  }
  else if(FLUXDIMENSPLIT==UNSPLIT){ // not expensive compared to matter terms
    // only here if doing both directions

    // then both cross directions exist and need to store in-between
    // storage for UNSPLIT flux_interp method
    fluxveca[1]=GLOBALPOINT(Fa);
    fluxveca[2]=GLOBALPOINT(Fa);
    fluxveca[3]=GLOBALPOINT(Fa);
    fluxvecb[1]=GLOBALPOINT(Fb);
    fluxvecb[2]=GLOBALPOINT(Fb);
    fluxvecb[3]=GLOBALPOINT(Fb);

    // if here, always will do fluxdir since fluxdir was defined to necessarily be an existing dimension
    otherdir=1;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
    if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxveca[fluxdir], NULL); 
     
    otherdir=2;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
    if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxveca[fluxdir], NULL, fluxveca[fluxdir], NULL); 
   

    otherdir=2;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
    if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvec[fluxdir], NULL, fluxvecb[fluxdir], NULL); 
     
    otherdir=1;  whichquantity=whichquantitylist[otherdir]; interporflux=interporfluxlist[otherdir]; dirmethod=dirmethodlist[otherdir]; intdir=intdirlist[otherdir]; fluxdir=fluxdirlist[otherdir]; idel=idellist[otherdir]; jdel=jdellist[otherdir];  kdel=kdellist[otherdir];
    if(Nvec[intdir]>1) flux_interp(whichpltoavg, ifnotavgthencopy, whichquantity, interporflux, dirmethod, idel, jdel, kdel, pr, NULL, fluxvecb[fluxdir], NULL, fluxvecb[fluxdir], NULL); 
     
    // symmetrize flux
    // assumes pl is chosen correctly (i.e. plforflux for FLUXRECON)
    // GODMARK: assumes fluxdir is same for both dimensions, which is normal
    COMPFULLLOOP PLOOP(pliter,pl) MACP1A1(fluxvec,fluxdir,i,j,k,pl)=0.5*(MACP1A1(fluxveca,fluxdir,i,j,k,pl)+MACP1A1(fluxvecb,fluxdir,i,j,k,pl));

  }





}











// convert avg->cent or cent->avg
// whichavg2cen = ENOCENT2AVGTYPE or ENOAVG2CENTTYPE
// does a full 3D deaverage or average for all directions
int avg2cen_interp(int *locpl, int *whichpltoavg, int *ifnotavgthencopy, int whichquantity, int whichavg2cen, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*origin)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR])
{
  int get_primitive_centerlocation(int *locpl, int *whichpltoavg, int interporflux, FTYPE (*primreal)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR], FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR]);
  void slope_lim_linetype(int whichquantity, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR], FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR], FTYPE (*p2interpm)[NSTORE2][NSTORE3][NPR],FTYPE (*p2interpp)[NSTORE2][NSTORE3][NPR], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR], FTYPE (*pright)[NSTORE2][NSTORE3][NPR]);
  int dir;
  int idel,jdel,kdel;
  int Nvec[NDIM];
  int i, j, k, pl, pliter;
  FTYPE (*current_in)[NSTORE2][NSTORE3][NPR]; //atch
  FTYPE (*current_out)[NSTORE2][NSTORE3][NPR]; //atch
  int numdirs,dirlist[NDIM],dirlisti;
  int numperms;
  long long itemp,permi;
  int get_dirs_list(int *numdirs, int *dirlist);
  int get_compdimen(int *numdirs, long long *itemp);
  int dirlist3d(long long itemp, int *dirlist);
  int dirlist2d(long long itemp, int *dirlist);
  int dirlist1d(long long itemp, int *dirlist);
  //use ptemparray as a temp array
  FTYPE (*prim_goodlocation)[NSTORE2][NSTORE3][NPR] = GLOBALPOINT(ptemparray);
  void (*multidir_pre_slope_lim_linetype)(void);
  void (*multidir_post_slope_lim_linetype)(void);
  extern void multidir_pre_slope_lim_linetype_weno(void);
  extern void multidir_post_slope_lim_linetype_weno(void);
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  FTYPE (*stencilvar)[NSTORE2][NSTORE3][NPR];
  FTYPE (*in)[NSTORE2][NSTORE3][NPR];








  // see if need to be here at all
  if(interporder[avgscheme[1]]<3 && interporder[avgscheme[2]]<3 && interporder[avgscheme[3]]<3){
    dualfprintf(fail_file,"Shouldn't be here in avg2cen_interp\n");
    myexit(333);
  }


  if(WENOINTERPTYPE(avgscheme[1]) || WENOINTERPTYPE(avgscheme[2]) || WENOINTERPTYPE(avgscheme[3])){
    // then assume avg2cen using WENO, so use weno version of multidir_pre and multidir_post functions
    multidir_pre_slope_lim_linetype=&multidir_pre_slope_lim_linetype_weno;
    multidir_post_slope_lim_linetype=&multidir_post_slope_lim_linetype_weno;
  }
  else{
    dualfprintf(fail_file,"No such defined multidir_pre/post_slope_lim_linetype() for avgscheme[1,2,3]=%d %d %d\n",avgscheme[1],avgscheme[2],avgscheme[3]);
    myexit(91758726);
  }
 


  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  //////////////////
  //
  // get primitive at CENT (or locpl) location where origin[] is located
  //
  //////////////////

  get_primitive_centerlocation(locpl, whichpltoavg, whichavg2cen, primreal, origin, prim_goodlocation);




  //////////////////////////////
  //
  // Perform interpolation (don't touch all outputs!)
  //
  //////////////////////////////


  //////////////////////////
  //
  // memory used for input so input isn't changed
  //
  //////////////////////////
  in=GLOBALPOINT(a2cin);
  

  ////////////
  //
  // get number of dimensions (ignore itemp)
  //
  ////////////
  get_compdimen(&numdirs,&itemp);


  // can't change input origin[], so copy over to in
  // from here on use in/out and not a2cin
  COMPFULLLOOP PLOOP(pliter,pl){
    // now copy over input to output and process only output (don't want to change input!)
    MACP0A1(in,i,j,k,pl) = MACP0A1(origin,i,j,k,pl);
    // don't touch all outputs!
    if(whichpltoavg[pl] || ifnotavgthencopy[pl]){
      MACP0A1(out,i,j,k,pl)=MACP0A1(in,i,j,k,pl);
    }
  }
  

  // choose quasistrang or do this process if only 1 direction since no splitting
  if(A2CDIMENSPLIT==PERFECTUNSPLIT){
    // assumes out[] is filled with original in[]

    //////////////
    //
    // get # of dimension and the directions to loop over
    //
    //////////////
    get_dirs_list(&numdirs, dirlist);


    ///////////////
    //
    // Get stencil var
    // don't necessarily touch all outputs
    // setup input array so don't have to assume user sends in=out
    //
    ///////////////

    if(numdirs>1){

      // storage for quantity to use for determining stencil
      stencilvar=GLOBALPOINT(stencilvartemp);

      // then need to store original to use for weights
      COMPFULLLOOP PLOOP(pliter,pl) MACP0A1(stencilvar,i,j,k,pl) = MACP0A1(in,i,j,k,pl);
    }
    else stencilvar=NULL; // says use normal input for stencil



    ////////////////////////////////
    //
    // loop over dimensions/directions
    //
    ////////////////////////////////
    for(dirlisti=1;dirlisti<=numdirs;dirlisti++){

      // get direction to integrate
      dir=dirlist[dirlisti];
      //    dualfprintf(fail_file,"steppart=%d nstep=%d : dirlisti=%d dir=%d\n",steppart,nstep,dirlisti,dir);


      // get loop details
      // GODMARK: determine idel,jdel,kdel
      idel = fluxloop[dir][FIDEL];
      jdel = fluxloop[dir][FJDEL];
      kdel = fluxloop[dir][FKDEL];


      // remove fields should not average or deaverage
      addremovefieldfromnpr(REMOVEFROMNPR, whichpltoavg, ifnotavgthencopy, whichavg2cen, dir, &nprlocalstart,&nprlocalend,nprlocallist,out,out);

      if(nprlocalend>-1){
        // a2c or c2a
        slope_lim_linetype(whichquantity, whichavg2cen, dir, idel, jdel, kdel, prim_goodlocation, stencilvar, out, NULL, out, NULL);
      }

      // restore dir-field
      addremovefieldfromnpr(RESTORENPR, whichpltoavg, ifnotavgthencopy, whichavg2cen, dir, &nprlocalstart,&nprlocalend,nprlocallist,out,out);

      // only processing out[] as final output
    }


  }
  else{
    if( (A2CDIMENSPLIT==QUASISTRANG)||(numdirs==1) ){


      //      dualfprintf(fail_file,"origin=%d a2cin=%d in=%d out=%d\n",origin,a2cin,in,out);

      // choose quasistrang or do this process if only 1 direction since no splitting
      // Quasi-strang splitting
      // Only applicable for 2D and 3D, where compared to UNSPLIT method this is 2X faster in 2D and 6X faster in 3D

      // get # of dimension and the directions to loop over
      get_dirs_list(&numdirs, dirlist);

      // setup input array so don't have to assume user sends in=out
      current_in = in; 
      current_out = out;

      // loop over dimensions/directions
      for(dirlisti=1;dirlisti<=numdirs;dirlisti++){

        // get direction to integrate
        dir=dirlist[dirlisti];
        // dualfprintf(fail_file,"steppart=%d nstep=%d : dirlisti=%d dir=%d\n",steppart,nstep,dirlisti,dir);

        // get loop details
        // GODMARK: determine idel,jdel,kdel
        idel = fluxloop[dir][FIDEL];
        jdel = fluxloop[dir][FJDEL];
        kdel = fluxloop[dir][FKDEL];


        // remove fields should not average or deaverage
        addremovefieldfromnpr(REMOVEFROMNPR, whichpltoavg, ifnotavgthencopy, whichavg2cen, dir, &nprlocalstart,&nprlocalend,nprlocallist,current_in,current_out);

        if(nprlocalend>-1){
          // a2c or c2a
          slope_lim_linetype(whichquantity, whichavg2cen, dir, idel, jdel, kdel, prim_goodlocation, NULL, current_in, NULL, current_out, NULL);
        }

        // restore dir-field
        addremovefieldfromnpr(RESTORENPR, whichpltoavg, ifnotavgthencopy, whichavg2cen, dir, &nprlocalstart,&nprlocalend,nprlocallist,current_in,current_out);


        // apply deaveraging one after another on the same array (out)
        current_in = current_out;

      }


    }
    else if(A2CDIMENSPLIT==UNSPLIT){


      // UNSPLIT a2c method (2X more expensive in 2D and 6X more expensive in 3D)

      if(numdirs==3){
        // then have 6 permutations to average out
        numperms=6;
      }
      else if(numdirs==2){
        // then have 2 permutations to average out
        numperms=2;
      }
      else numperms=1;




      // don't necessarily touch all outputs
      PLOOP(pliter,pl){
        if(whichpltoavg[pl] || ifnotavgthencopy[pl]){
          COMPFULLLOOP{
            // inialize output (if in=out, then ok since just above stored in)
            MACP0A1(out,i,j,k,pl) = 0.0;
          }
        }
      }




      // loop over permutations
      for(permi=0;permi<numperms;permi++){

        // get permutation
        if(numdirs==3) dirlist3d(permi,dirlist);
        else if(numdirs==2) dirlist2d(permi,dirlist);
        else if(numdirs==1) dirlist1d(permi,dirlist);

        // setup input array so don't have to assume user sends in=out
        current_in = in;
        // output array is the temporary storage for this permutation
        current_out = GLOBALPOINT(a2cout);


        // things to do before all directions for slope_lim_linetype()
        multidir_pre_slope_lim_linetype();


        // loop over dimensions/directions
        for(dirlisti=1;dirlisti<=numdirs;dirlisti++){

          // get direction to integrate
          dir=dirlist[dirlisti];
          //    dualfprintf(fail_file,"steppart=%d nstep=%d : dirlisti=%d dir=%d\n",steppart,nstep,dirlisti,dir);

          // get loop details
          // GODMARK: determine idel,jdel,kdel
          idel = fluxloop[dir][FIDEL];
          jdel = fluxloop[dir][FJDEL];
          kdel = fluxloop[dir][FKDEL];


          // remove fields should not average or deaverage
          addremovefieldfromnpr(REMOVEFROMNPR, whichpltoavg, ifnotavgthencopy, whichavg2cen, dir, &nprlocalstart,&nprlocalend,nprlocallist,current_in,current_out);


          if(nprlocalend>-1){
            // a2c or c2a
            //atch here every next interpolation will use the global value of whether we have reduced while interpolating in 
            //the other direction, and if we did reduce then, then it will reduce now in the same way
            slope_lim_linetype(whichquantity, whichavg2cen, dir, idel, jdel, kdel, prim_goodlocation, NULL, current_in, NULL, current_out, NULL);
          }

          // restore dir-field
          addremovefieldfromnpr(RESTORENPR, whichpltoavg, ifnotavgthencopy, whichavg2cen, dir, &nprlocalstart,&nprlocalend,nprlocallist,current_in,current_out);



          // apply de/averaging one after another on the a2cout array
          // and don't overwrite the a2cin array since needed for every permutation
          current_in = current_out;

        }// end loop over dirlisti


        multidir_post_slope_lim_linetype();


        // add this permutation to out (current_out=a2cout)
        PLOOP(pliter,pl) {
          if(whichpltoavg[pl] || ifnotavgthencopy[pl]){
            COMPFULLLOOP MACP0A1(out,i,j,k,pl) += MACP0A1(current_out,i,j,k,pl);
          }
        }
    

      }// end loop over permutations

      // normalize the permutations
      PLOOP(pliter,pl){
        if(whichpltoavg[pl] || ifnotavgthencopy[pl]){
          COMPFULLLOOP MACP0A1(out,i,j,k,pl) /= ( (FTYPE)numperms ) ;
        }
      }



    }// end UNSPLIT method
  }// end else if not PERFECTUNSPLIT method






  return( 0 );
}






// acts on globals, assumes static internals that get recalled upon reentering
int addremovefieldfromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int interptype, int dir, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[NSTORE2][NSTORE3][NPR], FTYPE (*out)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;
  int pl2,pl3;
  int i,j,k;


  if(doadd==REMOVEFROMNPR){
    ////////////////////////////////////////////
    //
    // save choice for interpolations
    *nprlocalstart=nprstart;
    *nprlocalend=nprend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=nprlist[pl];



    // now remove any other pl's not wanting to average for whatever reason
    for(pl3=0;pl3<NPR;pl3++){ // as above, but loop over all undesired quantities
      for(pl=nprstart;pl<=nprend;pl++){
        if(whichpltoavg[pl3]==0 && nprlist[pl]==pl3){
          // need to copy over unchanged quantity
          if(ifnotavgthencopy[pl3] && out!=NULL && in!=NULL) COMPFULLLOOP MACP0A1(out,i,j,k,pl3)=MACP0A1(in,i,j,k,pl3);
          for(pl2=pl+1;pl2<=nprend;pl2++) nprlist[pl2-1]=nprlist[pl2]; // moving upper to lower index
          nprend-=1; // removed dir-field
          break;
        }
      }
    }



    // avoid averaging or de-averaging field in dir-direction since already per-point in that direction for FLUXCTSTAG
    // field is never de-averaged along itself in staggered method
    if(interptype==ENOAVG2CENTTYPE && FLUXB==FLUXCTSTAG){


      // just remove dir-field
      pl3=B1+dir-1;
      for(pl=nprstart;pl<=nprend;pl++){
        if(nprlist[pl]==pl3){
          // need to copy over unchanged quantity
          // dir-field for avg2cent always gets copied!
          if(out!=NULL && in!=NULL) COMPFULLLOOP MACP0A1(out,i,j,k,pl3)=MACP0A1(in,i,j,k,pl3);
          for(pl2=pl+1;pl2<=nprend;pl2++) nprlist[pl2-1]=nprlist[pl2]; // moving upper to lower index
          nprend-=1; // removed dir-field
          break;
        }
      }




    }

  }
  else if(doadd==RESTORENPR){


    ////////////////////////////////////////////
    //
    // restore choice for interpolations
    nprstart= *nprlocalstart;
    nprend= *nprlocalend;
    PMAXNPRLOOP(pl) nprlist[pl]=nprlocallist[pl];
  }


  //  PALLLOOP(pl){
  //    dualfprintf(fail_file,"dir=%d interptype=%d nprstart=%d nprend=%d nprlist[%d]=%d\n",dir,interptype,nprstart,nprend,pl,nprlist[pl]);
  //  }
  

  return(0);

}




// Quasi-Strang Splitting (alternate order of dimension) of ordering of c2a integration of flux over 2 orthogonal directions
int get_dirs_list(int *numdirs, int *dirlist)
{
  long long itemp;
  int get_compdimen(int *numdirs, long long *itemp);
  int dirlist3d(long long itemp, int *dirlist);
  int dirlist2d(long long itemp, int *dirlist);
  int dirlist1d(long long itemp, int *dirlist);

  get_compdimen(numdirs,&itemp);

  if(*numdirs==3) dirlist3d(itemp,dirlist);
  else if(*numdirs==2) dirlist2d(itemp,dirlist);
  else if(*numdirs==1) dirlist1d(itemp,dirlist);

  return(0);
}



// get number of computational directions and the iterator unique to that dimension for choosing number of permutations
int get_compdimen(int *numdirs, long long *itemp)
{

  if(N1NOT1&&N2NOT1&&N3NOT1){
    *itemp=(nstep*(long long)numstepparts+(long long)steppart)%6;
    *numdirs=3;
  }
  else if( (N1NOT1&&N2NOT1)||(N1NOT1&&N3NOT1)||(N2NOT1&&N3NOT1) ){
    *itemp=(nstep*(long long)numstepparts+(long long)steppart)%2;
    *numdirs=2;
  }
  else{
    *itemp=(nstep*(long long)numstepparts+(long long)steppart)%1;
    *numdirs=1;
  }

  return(0);

}




// sequence of directions to interpolate in order of 1,2,3
int dirlist3d(long long itemp, int *dirlist)
{
  // 3D, number of permutations per call is: 2,2,1,2,2,1
  // normal ZEUS code has 1,2,1,2,1,1, so we are better mixed
  if(itemp==0){
    dirlist[1]=1;
    dirlist[2]=2;
    dirlist[3]=3;
  }
  else if(itemp==1){
    dirlist[1]=3;
    dirlist[2]=1;
    dirlist[3]=2;
  }
  else if(itemp==2){
    dirlist[1]=2;
    dirlist[2]=3;
    dirlist[3]=1;
  }
  else if(itemp==3){
    dirlist[1]=3;
    dirlist[2]=2;
    dirlist[3]=1;
  }
  else if(itemp==4){
    dirlist[1]=1;
    dirlist[2]=3;
    dirlist[3]=2;
  }
  else if(itemp==5){
    dirlist[1]=2;
    dirlist[2]=1;
    dirlist[3]=3;
  }

  return(0);
}


// sequence of directions to interpolate in order of 1,2
int dirlist2d(long long itemp, int *dirlist)
{
  if(N1NOT1&&N2NOT1){
    // 1 2
    if(itemp==0){
      dirlist[1]=1;
      dirlist[2]=2;
    }
    else if(itemp==1){
      dirlist[1]=2;
      dirlist[2]=1;
    }
  }
  else if(N1NOT1&&N3NOT1){
    // 1 3
    if(itemp==0){
      dirlist[1]=1;
      dirlist[2]=3;
    }
    else if(itemp==1){
      dirlist[1]=3;
      dirlist[2]=1;
    }
  }
  else if(N2NOT1&&N3NOT1){
    // 2 3
    if(itemp==0){
      dirlist[1]=2;
      dirlist[2]=3;
    }
    else if(itemp==1){
      dirlist[1]=3;
      dirlist[2]=2;
    }
  }

  return(0);
}



// sequence of directions to interpolate in order of 1
int dirlist1d(long long itemp, int *dirlist)
{

  if(N1NOT1){
    // 1
    dirlist[1]=1;
  }
  else if(N2NOT1){
    // 2
    dirlist[1]=2;
  }
  else if(N3NOT1){
    // 3
    dirlist[1]=3;
  }

  return(0);

}

















// GODMARK: Sasha says the limit_fluxc2a_prim_change() didn't work
// this code seems to have incorrect position for flux (uses CENT for pr, but doesn't at least interpolate/average to flux position)
// dir here is always direction of flux as it would be differenced
FTYPE limit_fluxc2a_prim_change( 
                                int dir, FTYPE (*pr)[NSTORE2][NSTORE3][NPR],
                                FTYPE (*fluxvec_point)[NSTORE2][NSTORE3][NPR+NSPECIAL],
                                FTYPE (*fluxvec_avg)[NSTORE2][NSTORE3][NPR+NSPECIAL])
{
  int is, ie, js, je, ks, ke;
  int i, j, k;
  int i1, j1, k1;
  int pl,pliter;
  int index;
  int idel, jdel, kdel;
  FTYPE delta_u_left[NPR];
  FTYPE delta_u_right[NPR];
  FTYPE delta_u[2][NPR];
  FTYPE Upoint[NPR], Upoint_updated[NPR];
  FTYPE pr_updated[NPR];
  PFTYPE pflag_current;
  PFTYPE pflag_backup;
  FTYPE frac_point_flux_used_array[2];
  FTYPE frac_point_flux_used;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_state q; //atch
  extern FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout );
  struct of_newtonstats newtonstats; setnewtonstatsdefault(&newtonstats);
  int showmessages=1;
  int allowlocalfailurefixandnoreport=1;

  //unit vector in the flux direction
  idel = (dir == 1);
  jdel = (dir == 2);
  kdel = (dir == 3);

  //start and end indices for the loop over the fluxes in the dir direction that will be used for evolving the conserved quantities
  //since the conserved quantites are confined to the range of indices given by Uconsloop and there is one more flux
  //that conserved quantities, we expand the loop by one grid cell in the direction of the flux
  // Uconsloop is correct instead of Uconsevolveloop that includes larger region where averaged fluxes are updated but there is no point value there
  is = Uconsloop[FIS]; 
  ie = Uconsloop[FIE] + idel;
  js = Uconsloop[FJS];
  je = Uconsloop[FJE] + jdel;
  ks = Uconsloop[FKS];
  ke = Uconsloop[FKE] + kdel;

  COMPZSLOOP( is, ie, js, je, ks, ke ) {
    PLOOP(pliter,pl) { //loop over the interfaces where the fluxes that determine the evolution of conserved quantities are
      //
      // !! Insert c2a limiting code here !! SASMARK atch
      //

      //
      // compute the changes that would be made to conserved quantities due to the c2a corrections for the fluxes at all interfaces during dt
      // fluxvec_avg contains the averaged fluxes
      //


      //pseudo-update to conserved quantity located to the left of the current interface, i.e. located at (i - idel, j - jdel, k - kdel),
      //due to c2a correction to the flux at (i,j,k) interface in the dir direction 
      delta_u_left[pl] = - dt * (MACP0A1(fluxvec_avg,i,j,k,pl) - MACP0A1(fluxvec_point,i,j,k,pl)) / dx[dir];

      //to the right:
      delta_u_right[pl] = - delta_u_left[pl];

      //unite them in one array
      delta_u[0][pl] = delta_u_left[pl];
      delta_u[1][pl] = delta_u_right[pl];
    }

    //loop through the two neighbouring cells that are affected by the MAC(fluxvec_avg,i,j,k);
    //index = 0 is the left and index = 1 is the right cells
    for( index = 0; index <= 1; index++ ) {

      i1 = i + (index-1) * idel;
      j1 = j + (index-1) * jdel;
      k1 = k + (index-1) * kdel;

      //1. get u_point_adjacent[0..1] from pr's -- inconsistent but this way avoid bounding

      // set geometry for zone to be updated
      get_geometry(i1, j1, k1, CENT, ptrgeom);

      // find U(pr)
      MYFUN(get_state(MAC(pr,i1,j1,k1), ptrgeom, &q),"flux.c:fluxcalc()", "get_state()", 1);
      MYFUN(primtoU(UEVOLVE,MAC(pr,i1,j1,k1), &q, ptrgeom, Upoint, NULL),"step_ch.c:advance()", "primtoU()", 1);


      //2. add delta_u[0..1] to it
      //that's what Upoint would become were it evolve only due to c2a flux correction on one interface
      PLOOP(pliter,pl) Upoint_updated[pl] = Upoint[pl] + delta_u[index][pl]; 
      PLOOP(pliter,pl) pr_updated[pl] = MACP0A1(pr,i1,j1,k1,pl);  //fill in the initial guess for inversion

      //3. invert & check if the difference between pr_adjacent_updated[0..1] and pr_adjacent_original[0..1]= pr is too large
      // invert point Upoint_updated-> point pr_updated
      pflag_backup = GLOBALMACP0A1(pflag,i1,j1,k1,FLAGUTOPRIMFAIL); //back up the old inversion flag, just in case, probably not needed anyway
   
      int eomtype=EOMDEFAULT;
      FTYPE dissmeasure=-1.0; // assume energy try ok
      int whichcap=CAPTYPEBASIC;
      int whichmethod=MODEDEFAULT;
      int modprim=0;
      int checkoninversiongas=CHECKONINVERSION;
      int checkoninversionrad=CHECKONINVERSIONRAD;
      MYFUN(Utoprimgen(showmessages,checkoninversiongas,checkoninversionrad,allowlocalfailurefixandnoreport, 0, &eomtype,whichcap,whichmethod,modprim,EVOLVEUTOPRIM, UEVOLVE, Upoint_updated, &q, ptrgeom, dissmeasure, pr_updated, pr_updated,&newtonstats),"flux.c:fluxcalc()", "Utoprimgen", 1);

      pflag_current = GLOBALMACP0A1(pflag,i1,j1,k1,FLAGUTOPRIMFAIL);  //backup the new inversion flag
      GLOBALMACP0A1(pflag,i1,j1,k1,FLAGUTOPRIMFAIL) = pflag_backup;   //restore the inversion flag

      if( 
         //(pflag_backup==UTOPRIMFAILUNEG || IFUTOPRIMNOFAILORFIXED(pflag_backup)) && //not sure if can use this since it will always fix up the failure in the end
         (pflag_current==UTOPRIMFAILUNEG || IFUTOPRIMNOFAILORFIXED(pflag_current))   //if only u < 0 or no failure it is OK (since only check rho & gamma)
          ) {

        //4. limit the flux c2a correction (difference btw. pr & pr_updated based on this difference
        frac_point_flux_used_array[index] = limit_prim_correction(MAX_AC_PRIM_FRAC_CHANGE, ptrgeom, MAC(pr,i1,j1,k1), pr_updated);

      }
      else {
        frac_point_flux_used_array[index] = 1.; //if inversion from new value failed, let there be no c2a correction done to flux
      }
    }

    frac_point_flux_used = MAX( frac_point_flux_used_array[0],  frac_point_flux_used_array[1] );

    if( 0 < frac_point_flux_used ) { 
      dualfprintf( fail_file, "Limited flux integration, dir = %d, i = %d, j = %d\n", dir, i, j );
    }

    PLOOP(pliter,pl) MACP0A1(fluxvec_avg,i1,j1,k1,pl) = frac_point_flux_used * MACP0A1(fluxvec_point,i1,j1,k1,pl) + (1.0 - frac_point_flux_used) * MACP0A1(fluxvec_avg,i1,j1,k1,pl);
  
  } //end COMPZSLOOP

  return( 0 );

}



// choose how to handle weights for higher order interpolations
// generally can be user-defined, but chosen as best right now
int higherorder_set(int whichquantity, int recontype, int*weightsplittype)
{

  if(whichquantity==ENOPRIMITIVE){
    if(recontype==CVT_C2E) *weightsplittype=DO_SPLITC2E;
    // no alternative at the moment
    else *weightsplittype=DO_SPLITC2E;
  }
  else if(whichquantity==ENOCONSERVED){
    if(recontype==CVT_A2C) *weightsplittype=DO_SPLITA2C;
    else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITA2C;
    else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITA2C;
  }
  else if(whichquantity==ENOSMOOTHCONSERVED){
    if(recontype==CVT_A2C) *weightsplittype=DO_SPLITA2CSMOOTH;
    else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITA2CSMOOTH;
    else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITA2CSMOOTH;
  }
  else if(whichquantity==ENOFLUX){
    if(recontype==CVT_A2C) *weightsplittype=DO_SPLITA2C4FLUX;
    else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITA2C4FLUX;
    else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITA2C;
  }
  else if(whichquantity==ENOMAFLUX){ // only MA flux
    if(splitmaem){ // only choose special weight method if really doing splitmaem
      if(recontype==CVT_A2C) *weightsplittype=DO_SPLITA2C4MAFLUX;
      else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITA2C4MAFLUX;
      else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITA2C;
    }
    else{
      // as for ENOFLUX
      if(recontype==CVT_A2C) *weightsplittype=DO_SPLITA2C4FLUX;
      else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITA2C4FLUX;
      else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITA2C;
    }
  }
  else if(whichquantity==ENOSMOOTHFLUX){ // probably only EM flux
    if(recontype==CVT_A2C) *weightsplittype=DO_SPLITA2C4SMOOTHFLUX;
    else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITA2C4SMOOTHFLUX;
    else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITA2C;
  }
  else if(whichquantity==ENOSOURCETERM){
    if(recontype==CVT_A2C) *weightsplittype=DO_SPLITSOURCE;
    else if(recontype==CVT_C2A) *weightsplittype=DO_SPLITSOURCE;
    else if(recontype==CVT_C2E) *weightsplittype=DO_SPLITSOURCE;
  }

  return(0);


}


// which quantity to set as reference pl for weight calculation
// generally can be user-defined, but chosen as best right now
int plstart_set(int whichquantity, int dir, int recontype, int *plstart)
{
  int Batemp,Bbtemp;
  // GODMARK: Very dependent upon what quantities dealing with -- so why depends upon what's in nprlist
      
  if(dir==1){
    Batemp=B2;
    Bbtemp=B3;
  }
  else if(dir==2){
    Batemp=B1;
    Bbtemp=B3;
  }
  else if(dir==3){
    Batemp=B1;
    Bbtemp=B2;
  }
    
  // If a2c or c2a and flux, use F_{dir}, if conserved quantity use E
  // GODMARK: Very dependent upon what quantities dealing with -- so why depends upon what's in nprlist

  if(
     // full conserved
     (nprstart==0 && nprend==NPR-1 && nprlist[nprstart]==RHO && nprlist[nprend]==NPR-1)
     // non-field
     ||(nprstart==0 && nprend==U3 && nprlist[nprstart]==RHO && nprlist[nprend]==U3)
     // conserved with field along dir removed as for FLUXCTSTAG
     ||(nprstart==0 && nprend==NPR-2 && nprlist[0]==RHO && nprlist[1]==UU && nprlist[2]==U1 && nprlist[3]==U2 && nprlist[4]==U3 && nprlist[5]==Batemp && nprlist[6]==Bbtemp)
     ){
    // then doing all or MA quantities so can use Sasha's method for controlling weights
    

    if( whichquantity == ENOFLUX || whichquantity == ENOSMOOTHFLUX || whichquantity == ENOMAFLUX) {
      //for fluxes, use common weights computed for the flux of dir-momentum in the direction of dir (because that flux does not vanish even if the velocities are zero)
      // That is, pressure term remains
      *plstart = dir + U1 - 1; 
    }
    else if( whichquantity == ENOCONSERVED || whichquantity == ENOSMOOTHCONSERVED){
      //for conserved quantities, use common weights computed for the energy (because in general total energy is unlikely to vanish)
      *plstart = UU;
    }
    else {
      dualfprintf( fail_file, "plstart_set not setup for %d %d\n",whichquantity,recontype);
      myexit( 138629762 );
    }

  }
  else if(
          (nprstart==0 && nprend==2 && nprlist[nprstart]==B1 && nprlist[nprend]==B3)||
          (nprstart==0 && nprend==0 && (nprlist[nprstart]==B1 && nprlist[nprend]==B1 || nprlist[nprstart]==B2 && nprlist[nprend]==B2 || nprlist[nprstart]==B3 && nprlist[nprend]==B3) )
          ){
 
    // then doing fraction of normal quantities (e.g. field or emf)
    *plstart=nprlist[nprstart];
  }
  else{
    // force user to define this setup if didn't exist
    dualfprintf(fail_file,"This list of variables not setup in plstart_set()\n");
    int pl,pliter;
    PMAXNPRLOOP(pl) dualfprintf(fail_file,"nprstart=%d nprend=%d nprlist[%d]=%d\n",nprstart,nprend,pl,nprlist[pl]);
    myexit(834672346);
  }


  return(0);


}

