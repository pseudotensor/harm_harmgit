/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define SLOWFAC 1.0  /* reduce u_phi by this amount */

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;
double tabrho[512][128];
double tabom[512][128];

FTYPE rin;

int fieldfrompotential[NDIM];


int prepre_init_specific_init(void)
{

  // choice// GODMARK: not convenient location, but needed for init_mpi()
  periodicx1=0;
  periodicx2=0;
#if(USEMPI&&N3!=1)
  periodicx3=1;// GODMARK: periodic in \phi for 3D spherical polar
#else
  periodicx3=0;
#endif


  return(0);

}


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.1; // REBECCAMARK
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  return(0);
}


int init_conservatives(FTYPE (*p)[NSTORE2][NSTORE3][NPR], FTYPE (*Utemp)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;
  
  PLOOPBONLY(pl) trifprintf("fieldfrompotential[%d]=%d\n",pl-B1+1,fieldfrompotential[pl-B1+1]);


  trifprintf("begin init_conservatives\n");
  pi2Uavg(fieldfrompotential, p, Utemp, U);
  trifprintf("end init_conservatives\n");

  return(0);

}


int post_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;

  //  UTOPRIMVERSION =   UTOPRIM5D1;
  //    UTOPRIMVERSION =   UTOPRIM2DFINAL;

  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  return(0);

}


#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1

#define WHICHPROBLEM 0

int init_grid(void)
{
  SFTYPE rh;
  SFTYPE rminus;


  // metric stuff first
  a = 0.7 ;
  rminus=1.-pow(1.-a*a, 0.5);


#if(WHICHPROBLEM==NORMALTORUS)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.3;  // REBECCAMARK
  Rout = 35.;
  // define coordinate type
  //defcoord=0;
  defcoord = REBECCAGRID;
  
#elif(WHICHPROBLEM==GRBJET)
  R0 = -3.0;
  Rout = 1E5;
  // define coordinate type
  defcoord = JET4COORDS;

  // defcoord=0.;
#endif

 
  Rhor=rhor_calc(0);
  // Rin=setRin(setihor());  // REBECCAMARK
  Rin = 0.95* Rhor;
  
  hslope = 0.15;  // REBECCAMARK




  return(0);

}

int init_global(void)
{
  int pl,pliter;

  DODIAGEVERYSUBSTEP = 0;


  SAFE=1.3;
  //  cour = 0.9;
  cour=0.8;
  //  cour = 0.5;

  ///////////////////////
  //
  // ENO-RELATED STUFF
  //
  ///////////////////////
  
  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;

  LIMIT_AC_FRAC_CHANGE=0; // CHANGINGMARK: avoiding complication for now
  // have to make consistent with weno-minimization for fluxes
  MAX_AC_FRAC_CHANGE=0.2;

  // need MAXBND==17 if not evolving point field and doing full WENO5BND
  // test1103 N1=8 is smallest tried and new simple_ limiting works at 10% very well
  //dofluxreconevolvepointfield=0;
  // only need MAXBND==11 like normal higher-order method (like FV method)
  dofluxreconevolvepointfield=1;



#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)


  //avgscheme[1]=avgscheme[2]=avgscheme[3]=WENO5BND;
  lim[1] = lim[2] = lim[3] = PARALINE;

  avgscheme[1]=avgscheme[2]=avgscheme[3]=DONOR; // CHANGINGMARK
  //lim[1] = lim[2] = lim[3] = MC; // CHANGINGMARK


  DOENOFLUX = NOENOFLUX;
  //DOENOFLUX = ENOFLUXRECON; // CHANGINGMARK
  //DOENOFLUX = ENOFINITEVOLUME;

  if(DOENOFLUX == ENOFLUXRECON){
    // below applies to all fluxes
    PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
    PLOOPBONLY(pl) do_transverse_flux_integration[pl] = 1;
    // below used for re-averaging of field in advance.c for dBhat/dt method
    PALLLOOP(pl) do_conserved_integration[pl] = 1;
    PLOOPBONLY(pl) do_conserved_integration[pl] = 1;
  }

  if(DOENOFLUX == ENOFINITEVOLUME){
    PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
    PLOOPBONLY(pl) do_transverse_flux_integration[pl] = 1;
    PALLLOOP(pl) do_source_integration[pl] = 0;
    PLOOPBONLY(pl) do_source_integration[pl] = 0;
    PALLLOOP(pl) do_conserved_integration[pl] = 1;
    PLOOPBONLY(pl) do_conserved_integration[pl] = 1;
    //    do_conserved_integration[B1-1+DIRZ] = 1;
  }



  //FLUXB = FLUXCTSTAG;  // CHANGINGMARK
  //  FLUXB = FLUXCTHLL;
  FLUXB = FLUXCTTOTH;
  TIMEORDER=2;
  //TIMEORDER=4;
  //  TIMEORDER=3;
  fluxmethod= HLLFLUX;
  //fluxmethod= LAXFFLUX; // generally more robust than HLLFLUX for various reasons
  

  //  UTOPRIMVERSION=UTOPRIM5D1;
  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;
  //  UTOPRIMVERSION = UTOPRIM1DFINAL;


#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim[1] =lim[2]=lim[3]= MC;
  TIMEORDER=2;


  // below applies to all fluxes
  PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
  PLOOPBONLY(pl) do_transverse_flux_integration[pl] = 1;
  // below used for re-averaging of field in advance.c for dBhat/dt method
  PALLLOOP(pl) do_conserved_integration[pl] = 1;
  PLOOPBONLY(pl) do_conserved_integration[pl] = 1;



  fluxmethod=LAXFFLUX; // generally more robust than HLLFLUX
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
  // whether/which ENO used to interpolate fluxes
  DOENOFLUX = ENOFINITEVOLUME;
  //  DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
#endif



  ranc(1,7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=COOLGAMMIETHINDISK;


#if(WHICHPROBLEM==NORMALTORUS)
  BCtype[X1UP]=OUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  //BSQORHOLIMIT=1E3; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  //BSQOULIMIT=1E4; // was 1E3 but latest BC test had 1E4
  RHOMIN=1.e-6;
  UUMIN =1.e-6;
  RHOMINLIMIT=1.e-16;
  UUMINLIMIT =1.e-16;
  BSQORHOLIMIT=1E1;
  BSQOULIMIT=1E1;
  UORHOLIMIT=1E1;


#elif(WHICHPROBLEM==GRBJET)
  BCtype[X1UP]=FIXEDOUTFLOW;
  rescaletype=4;
  BSQORHOLIMIT=1E3;
  BSQOULIMIT=1E4;
  RHOMIN = 23.0;
  UUMIN = 1.7;
#endif


  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;



  // below floor model is only used if rescaletype!=4
  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prfloorcoef[RHO]=RHOMIN/100.0;
    prfloorcoef[UU]=UUMIN/100.0;
  }
  else{
    prfloorcoef[RHO]=RHOMIN;
    prfloorcoef[UU]=UUMIN;
  }


#if(WHICHPROBLEM==NORMALTORUS)
  /* output choices */
  tf = 8000.0;

  DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 50.;   /* dumping frequency, in units of M */
  DTdumpgen[DISSDUMPTYPE]=DTdumpgen[MAINDUMPTYPE]/5.0;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 50.0;
  DTdumpgen[FIELDLINEDUMPTYPE]=DTdumpgen[ENERDUMPTYPE] = 2.0;   /* logfile frequency, in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 2.0;   /* image file frequ., in units of M */
  DTdumpgen[DEBUGDUMPTYPE] = 1000.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 1000;   /* restart file period in steps */
#elif(WHICHPROBLEM==GRBJET)
  /* output choices */
  tf = 5E5;
  
  DTd = 250.;                 /* dumping frequency, in units of M */
  DTavg = 250.0;
  DTener = 2.0;                       /* logfile frequency, in units of M */
  DTi = 10.0;                 /* image file frequ., in units of M */
  DTdebug = 250.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;                  /* restart file period in steps */
#endif



  return(0);

}


int read_data(void)
{


  int i,j, k;
  SFTYPE tabr;
  SFTYPE tabo;
  FILE * fpr;
  FILE * fpo;
  int ti, tj , tk;
  int procno;
  int nscanned1, nscanned2;
  long lineno, numused;
  int myk;




  myk=2;

#if( USEMPI )
  for( procno = 0; procno < numprocs; procno++ ) {
    if( myid == procno ){
#endif
      //////
      //
      //  Open file for reading

      trifprintf("myid is %d\n",myid);

      //      fp = fopen( fname, "rb" );

      fpr=fopen("./mytest", "r");
      


      if( NULL != fpr ) {
        //trifprintf( "done.\npi2Uavg_specific(): reading data from %s... ", fpr );

        lineno = 0;
        numused = 0;

        while( (nscanned1 = fscanf(fpr, "%d %d %d %lg %lg \r\n", &ti, &tj, &tk , &tabr, &tabo )) == 5 ) {



          ///////
          //
          //  Check if within range and if within, save it to the array
          //  i, j - local grid indices for this processor, ti, tj - global

          // trifprintf("startpos1=%d, startpos2=%d \n", startpos[1], startpos[2]);

          i = ti -1 - startpos[1];
          j = tj-1 - startpos[2];
          k = tk-1-startpos[3];

          if( i >= INFULL1 && i <= OUTFULL1 &&
              j >= INFULL2 && j <= OUTFULL2 && k>=INFULL3 && k<=OUTFULL3 ) {


            MACP0A1(p,i,j,k,RHO) = tabr;
            MACP0A1(p,i,j,k,U3) = tabo;

            trifprintf(" k = %d , i=%d, j=%d, rho=%lg\n", k,i, j, tabr);


            //   trifprintf("i=%d, j=%d, k=%d, ti=%d, tj=%d, rho=%lg, om=%lg\n", i, j, k, ti, tj, tabr, tabo);

            //  MACP0A1(Uavg,i,j,k,U1) = p1a / coordparams.timescalefactor;
            // MACP0A1(Uavg,i,j,k,U2) = p2a / coordparams.timescalefactor;
            // MACP0A1(Uavg,i,j,k,U3) = 0.0;
            // MACP0A1(Uavg,i,j,k,B1) = 0.0;
            // MACP0A1(Uavg,i,j,k,B2) = 0.0;
            // MACP0A1(Uavg,i,j,k,B3) = 0.0;



            //              get_geometry(i, j, k, CENT, &geom);

            // PLOOP(pliter,pl) {
            //        MACP0A1(Uavg,i,j,k,pl) *= geom.e[pl];  //multiply conserved quantities by the gdet to have the same representation as inside the code
            /// }



            // FTYPE   coord(i, j, k, CENT, X);
            //FTYPE   bl_coord(X, V);
            //FTYPE   get_geometry(i,j,k,CENT,&geom) ;


            //convert the conserved momenta (covariant quantities) from orthonormal basis to coordinate one.
            // MACP0A1(Uavg,i,j,k,U1) *= dxdxp[1][1];
            //      MACP0A1(Uavg,i,j,k,U2) *= dxdxp[2][2];

            //*whichvel=VEL3;
            //*whichcoord=BLCOORD;

            numused++;
          }
          lineno++;
        }
        // Close file
        fclose( fpr );
        //
        ///////
      }
      else{
        // if null report and fail!
        dualfprintf(fail_file,"No Rebecca file 1\n");
        myexit(246346);
      }
#if(USEMPI)
    }
    MPI_Barrier(MPI_COMM_GRMHD);
  }
#endif


  if(numprocs==1) {
    
    
    fpr=fopen("./mytest", "r");
    

    if(fpr!=NULL){
    
      for(j=0; j<N2; j++){
 
        for(i=0; i<N1; i++){
   
          for(k=0; k<N3; k++){
     
     
            //      dualfprintf(fail_file,"i=%d j=%d\n",i,j);
     
            fscanf(fpr, "%d %d %d %lg %lg \r\n", &ti, &tj, &tk , &tabr, &tabo );
     
            MACP0A1(p,i,j,k,RHO) = tabr;
            MACP0A1(p,i,j,k,U3) = tabo;
     
     
            //      printf("tabom %lf  ", tabom[i][j]);

     
            //      dualfprintf(fail_file,"i=%d j=%d tabr=%21.15g\n",i,j,tabr);
          }
        }
      }
    }
    else{
      // if null report and fail!
      dualfprintf(fail_file,"No Rebecca file 2\n");
      myexit(246346);
    }
    //}
    printf("files read");
    fclose(fpr);
    // fclose(fpo);

  }

  return(0);
}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl,pliter;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  set_atmosphere(-1,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  if(pr[RHO]<pratm[RHO]){
    PLOOP(pliter,pl) pr[pl]=pratm[pl];
  }
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}

int init_grid_post_set_grid(void)
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);




  //SASMARK restart: need to populate panalytic with IC's
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives
    MYFUN(init_primitives(panalytic),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }




  
  // diagnostic
  // determine nature of inner radial edge (assumes myid==0 is always there)
  if(myid==0){
    i=INFULL1;
    j=k=0;
    coord(i,j,k, FACE1, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    trifprintf("rmin: %21.15g\n", r);
    trifprintf("rmin(i=%d,X=%21.15g): %21.15g\n", i,X[1],r);
    //    trifprintf("rmin/rsing: %21.15g\n", r / (a+SMALL));
    if(r/Rhor<=1.0){
      trifprintf("inner grid is inside horizon\n");
    }
    else{
      trifprintf("inner grid is outside horizon\n");
    }
  }

  // check that singularities are properly represented by code
  check_spc_singularities_user();

  
  return(0);

}



int init_primitives(FTYPE (*p)[NSTORE2][NSTORE3][NPR])
{
  int whichvel, whichcoord;
  int initreturn;
  int i = 0, j = 0, k = 0, l;
  struct of_geom geom;
  FTYPE r,th,X[NDIM],V[NDIM];
  int normalize_densities(FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
  int normalize_field(FTYPE (*p)[NSTORE2][NSTORE3][NPR]);
  int init_dsandvels(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *p, FTYPE *pstag);
  int init_atmosphere(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *pr);
  int pl,pliter;


  ///////////////////////////////////
  //
  // Assign primitive variables
  //
  ///////////////////////////////////
  trifprintf("Assign primitives\n");

  // assume we start in bl coords and convert to KSprim
  COMPFULLLOOP{
    PLOOP(pliter,pl) MACP0A1(p,i,j,k,pl)=0.0; // so field defined when get to floor model (fixup)
  }

  ////////////////////////////////////
  //
  // read data from Rebecca's file
  //
  ////////////////////////////////////
  read_data();


  // assume we start in bl coords and convert to KSprim
  COMPFULLLOOP{
    initreturn=init_dsandvels(&whichvel, &whichcoord,i,j,k,MAC(p,i,j,k),MAC(pstagscratch,i,j,k)); // request densities for all computational centers
    if(initreturn>0) return(1);
    else MYFUN(transform_primitive_vB(whichvel, whichcoord, i,j,k, p, pstagscratch),"init.c:init_primitives","transform_primitive_vB()",0);
  }


  /////////////////////////////
  //
  // normalize density if wanted
  //
  /////////////////////////////// 
  // at this point densities are still standard, so just send "p"
  trifprintf("Normalize densities\n");
  normalize_densities(p);


  /////////////////////////////
  //
  // Define an atmosphere if wanted
  //
  /////////////////////////////// 

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  // normalized atmosphere
  trifprintf("Add atmosphere\n");
  COMPZLOOP{
    initreturn=init_atmosphere(&whichvel, &whichcoord,i,j,k,MAC(p,i,j,k));
    if(initreturn>0) return(1);
    else{
      // transform from whichcoord to MCOORD
      if (bl2met2metp2v(whichvel, whichcoord,MAC(p,i,j,k), i,j,k) >= 1)
        FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
    }
    //  trifprintf("done with atmo\n");
  }
#endif


  // copy over initial solution as analytic solution
  // NEEDED FOR BOUND in case uses panalytic
  COMPFULLLOOP{
    PLOOP(pliter,pl) MACP0A1(panalytic,i,j,k,pl)=MACP0A1(p,i,j,k,pl);
  }


  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #1\n");

#if(EOMTYPE!=EOMFFDE)
  // trifprintf("doing something\n");
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  if (bound_prim(STAGEM1,0.0,p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_prim()", 1);

  if(pre_fixup(STAGEM1,p)>=1)
    FAILSTATEMENT("init.c:init()", "pre_fixup()", 1);
  //trifprintf("donefix\n");
#endif
  
  
  // trifprintf("done fixup\n");
  /////////////////////////////
  //
  // Initialize field from vector potential
  //
  ///////////////////////////////
  
 
#if(1)
  init_vpot(p);
  normalize_field(p); // normalizes p and pstagscratch and unew and vpotarray if tracked
#else
  // no field
  init_zero_field(p);
#endif
  



  // copy over initial solution as analytic solution
  // NEEDED FOR BOUND in case uses panalytic                                                                                                                                                            
  COMPFULLLOOP{
    PLOOP(pliter,pl) MACP0A1(panalytic,i,j,k,pl)=MACP0A1(p,i,j,k,pl);
  }


  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  //
  // BOUND AGAIN IN CASE USING PANALYTIC TO BOUND!
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #2\n");

#if(EOMTYPE!=EOMFFDE)
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  if (bound_allprim(STAGEM1,0.0,p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_allprim()", 1);

  if(pre_fixup(STAGEM1,p)>=1)
    FAILSTATEMENT("init.c:init()", "pre_fixup()", 1);
#endif




  return(0);


}



// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, mrho;
  FTYPE X[NDIM],V[NDIM],r,th,R;
  struct of_geom realgeom,geom;
  
  
  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  int pl,pliter;


  coord(i, j,k,  CENT, X);
  bl_coord(X, V);
 
  r=V[1];
  th=V[2];

  //  trifprintf("this is i, j, k ,r, th %d, %d, %d , %g, %g\n", i, j, k,r, th);
  //  fprintf(stderr,"proc: %d : %d %d\n",myid,i,j); fflush(stderr);
  
  //  beta=1.e2;
  beta=100;
  randfact=0.1;
  //rin = (1. + h_over_r)*Risco;
  rin = 8.5 ;
  
  
  
  /* region outside disk */
  R = r*sin(th) ;
  
  // pr[RHO]=MAC(tabrho,i,j,0);
  //pr[U3]=MAC(tabom,i,j,0);   


  if(r<rin || pr[0] < 10e-24)
    {
      trifprintf("density still zero r is %g \n", r);
      pr[RHO]=1E-30;
      // small density will indicate to atmosphere to change it

      *whichvel=WHICHVEL;
      *whichcoord=PRIMECOORDS;


    }

  else{
    mrho=pr[0];
    u = 0.01*(pow(mrho,gam))/(gam - 1.) ;
    //    u=0;
    ur = 0. ;
    uh = 0. ;

    //  up = 1./(pow(r,1.5) + a) ;
    // solution for 3-vel




    //  MAC(p,i,j,UU) = u* (1. + randfact * (ranc(0,0) - 0.5));
    //  MAC(p,i,j,UU)= u* (1. + randfact * (ranc(0,0) - 0.5));
    pr[UU]=u* (1. + randfact * (ranc(0,0) - 0.5));

    pr[U1] = ur ;
    pr[U2] = uh ;    

    //pr[U3]=1./pow(r, 1.5);
    trifprintf("r=%lf, rho=%lf, i=%d, j=%d\n", r, mrho, i, j);

    PLOOPBONLY(pl) pr[pl] = pstag[pl] = 0.0;


    *whichvel=VEL3;
    *whichcoord=BLCOORDS;


  }
  return(0);
 
}


#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2

#define FIELDTYPE DISKFIELD

// assumes normal field in pr
int init_vpot_user(int *whichcoord, int l, int i, int j, int k, int loc, FTYPE (*pr)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, u_av, q, fieldhor;
  FTYPE r,th;

  *A=0.0;

  if(l==3){// A_\phi

    r=V[1];
    th=V[2];

    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      *A += 0.5*r*sin(th) ;
    }
    /* field-in-disk version */

    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)){
      // average of density that lives on CORN3
#define STARTFIELD (1.1*rin)
      fieldhor=0.06;
      // since init_vpot() is called for all i,j,k, can't use
      // non-existence values, so limit averaging:
      if((i==-N1BND)&&(j==-N2BND)){
        rho_av = MACP0A1(pr,i,j,k,RHO);
      }
      else if(i==-N1BND){
        rho_av = AVGN_2(pr,i,j,k,RHO);
      }
      else if(j==-N2BND){
        rho_av = AVGN_1(pr,i,j,k,RHO);
      }
      else{ // normal cells
        rho_av = AVGN_for3(pr,i,j,k,RHO);
        u_av=AVGN_for3(pr,i,j,k,UU);
      }

      //q = rho_av / rhomax - 0.2;

      //if (q > 0.)      *A += q;

      if(r > STARTFIELD) q = ((u_av/umax) - 0.2)*pow(r,0.75) ;
      else q = 0. ;
      trifprintf("rhoav=%g q=%g\n", rho_av, q);

      //     if (q > 0.)      *A += q;
      //      trifprintf("%d %d rho_av=%g A=%g\n",i,j,rho_av, *A);

      if(q > 0.){

        //    *A+ = q*sin(log(r/STARTFIELD)/fieldhor)* (1. + 0.1 * (ranc(0,0) - 0.5));
        *A+=q*q*sin(log(r/STARTFIELD)/fieldhor)* (1. + 0.02 * (ranc(0,0) - 0.5))  ;

        //  *A+=q;
      }

    }
  }

  //////////////////////////////////
  //
  // finally assign what's returned
  //
  //////////////////////////////////
  *whichcoord = MCOORD;



  return(0);

}



int init_vpot2field_user(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  extern int vpot2field(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*pr)[NSTORE2][NSTORE3][NPR]);
  int i,j,k,pl,pliter;
  int toreturn;
 


  // obtain primitive magnetic field from vector potential
  toreturn=vpot2field(A,panalytic); // uses panalytic as temporary variable

  // default (assume all fields are from potential)
  PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=1;


  // Can override vector potential choice for some field components, like B3 in axisymmetry
  // see init.sasha.c

  ////////////////////
  //
  // don't override
  //
  ////////////////////
  COMPFULLLOOP{
    PLOOPBONLY(pl) MACP0A1(pr,i,j,k,pl)=MACP0A1(panalytic,i,j,k,pl);
  }

  return(toreturn);

}



// assumes we are fed the true densities
int normalize_densities(FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;


  rhomax=0;
  umax=0;
  COMPZLOOP{
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];

    if (MACP0A1(pr,i,j,k,RHO) > rhomax)   rhomax = MACP0A1(pr,i,j,k,RHO);
    if (MACP0A1(pr,i,j,k,UU) > umax && r > rin)    umax = MACP0A1(pr,i,j,k,UU);
  }

  mpimax(&rhomax);
  mpimax(&umax);
  trifprintf("rhomax: %21.15g umax: %21.15g\n", rhomax, umax);


  COMPZLOOP{
    MACP0A1(pr,i,j,k,RHO) /= rhomax;
    MACP0A1(pr,i,j,k,UU) /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

  trifprintf("now rhomax=%lg\n", rhomax);

  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE (*pr)[NSTORE2][NSTORE3][NPR])
{
  int i,j,k;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;

  bsq_max = 0.;
  COMPZLOOP {
    get_geometry(i, j, k, CENT, &geom);    

    if(FIELDTYPE==VERTFIELD){
      coord(i, j, k, CENT, X);
      bl_coord(X, V);
      r=V[1];
      th=V[2];
      
      if((r>rin)&&(fabs(th-M_PI*0.5)<4.0*M_PI*dx[2]*hslope)){
        if (bsq_calc(MAC(pr,i,j,k), &geom, &bsq_ij) >= 1)
          FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
 
        if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
      }
    }
    else{
      if (bsq_calc(MAC(pr,i,j,k), &geom, &bsq_ij) >= 1)
        FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      
      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %21.15g\n", bsq_max);

  /* finally, normalize to set field strength */
  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
  trifprintf("initial beta: %21.15g (should be %21.15g)\n", beta_act,beta);
  norm = sqrt(beta_act / beta);
  
  
  // not quite right since only correct static field energy, not moving field energy
  normalize_field_withnorm(norm);


  // check bsq_max again
  bsq_max = 0.;
  COMPZLOOP {
    get_geometry(i, j, k, CENT, &geom);
    if (bsq_calc(MAC(pr,i,j,k), &geom, &bsq_ij) >= 1)
      FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
    if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    
  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %21.15g\n", bsq_max);

  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

  trifprintf("new bsq_max: %21.15g\n", bsq_max);
  trifprintf("final beta: %21.15g (should be %21.15g)\n", beta_act,beta);

  return(0);
}



#undef SLOWFAC

SFTYPE lfish_calc(SFTYPE r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
           ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
            sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
            ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) +
                                          pow(a,
                                              2) * (2. + r))) / sqrt(1 +
                                                                     (2.
                                                                      *
                                                                      a)
                                                                     /
                                                                     pow
                                                                     (r,
                                                                      1.5)
                                                                     -
                                                                     3.
                                                                     /
                                                                     r)))
          / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
             (pow(a, 2) + (-2. + r) * r))
          );
}

// UUMIN/RHOMIN used for atmosphere

// for each WHICHVEL possibility, set atmosphere state for any coordinate system
// which=0 : initial condition
// which=1 : evolution condition (might also include a specific angular momentum or whatever)
// which==1 assumes pr set to something locally reasonable, and we adjust to that slowly

#define TAUADJUSTATM (10.0) // timescale for boundary to adjust to using preset inflow
int set_atmosphere(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE rho,u,ur,uh,up;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  int pl,pliter;

  // Bondi like initial atmosphere
  //    rho = RHOMIN * 1.E-14;
  //    u = UUMIN * 1.E-14;
  coord(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X);
  bl_coord(X,V);
  r=V[1];
  th=V[2];

  // default
  PLOOP(pliter,pl) prlocal[pl]=pr[pl];

#if((EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD))
  // Bondi-like atmosphere
  if(rescaletype==4){
#if(WHICHPROBLEM==NORMALTORUS)
    // couple rescaletype to atmosphere type
    //prlocal[RHO] = RHOMIN*pow(r,-2.0);
    prlocal[RHO] = RHOMIN*pow(r, -1.5);
    prlocal[UU]=0.00064*pow(prlocal[RHO], gam)/(gam-1.);
#elif(WHICHPROBLEM==GRBJET)
    // couple rescaletype to atmosphere type
    if(r>40.0) prlocal[RHO] = RHOMIN*pow(r,-2.0);
    else prlocal[RHO] = RHOMIN*pow(40.0,-2.0);
#endif
  }
  else{
    prlocal[RHO] = RHOMIN*pow(r,-1.5);
  }
#else
  prlocal[RHO] = 0;
#endif

#if(EOMTYPE==EOMGRMHD)
  // Bondi-like atmosphere
  //prlocal[UU]  = UUMIN*pow(r,-2.5);
  prlocal[UU]=0.00064*pow(prlocal[RHO], gam)/(gam-1.);

#else
  prlocal[UU]  = 0;
#endif

    
  // bl-normal observer (4-vel components)
  
  // normal observer velocity in atmosphere
  if(whichvel==VEL4){
    prlocal[U1] = -ptrgeom->gcon[GIND(0,1)]/sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
    prlocal[U2] = -ptrgeom->gcon[GIND(0,2)]/sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
    prlocal[U3] = -ptrgeom->gcon[GIND(0,3)]/sqrt(-ptrgeom->gcon[GIND(0,0)]) ;
  }
  else if(whichvel==VEL3){
    prlocal[U1] = ptrgeom->gcon[GIND(0,1)]/ptrgeom->gcon[GIND(0,0)] ;
    prlocal[U2] = ptrgeom->gcon[GIND(0,2)]/ptrgeom->gcon[GIND(0,0)] ;
    prlocal[U3] = ptrgeom->gcon[GIND(0,3)]/ptrgeom->gcon[GIND(0,0)] ;
    // GAMMIE
    //ur = -1./(r*r);
    //uh=up=0.0;
  }
  else if(whichvel==VELREL4){
    prlocal[U1] = 0.0;
    prlocal[U2] = 0.0;
    prlocal[U3] = 0.0;
  }
  
  if(whichcond==1){
    if(100.0*dt>TAUADJUSTATM){
      dualfprintf(fail_file,"dt=%21.15g and TAUADJUSTATM=%21.15g\n",dt,TAUADJUSTATM);
      myexit(1);
    }
    // TAUADJUSTATM must be >> dt always in order for this to make sense (i.e. critical damping to fixed solution)
    PLOOP(pliter,pl) pr[pl] = pr[pl]+(prlocal[pl]-pr[pl])*dt/TAUADJUSTATM;
  }
  else if(whichcond==0){ 
    PLOOP(pliter,pl) pr[pl] = prlocal[pl];
    // very specific
    // always outflow field
    //    pr[B1] = pr[B2] = pr[B3] = 0;
  }
  else if(whichcond==-1){ 
    // t=0, just sets as "floor"
    PLOOP(pliter,pl) pr[pl] = prlocal[pl];
    pr[B1] = pr[B2] = pr[B3] = 0;
  }


  return(0);

}


int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  int retval;

  // set defaults
  retval=set_density_floors_default(ptrgeom, pr, prfloor);


  // Now do Rebecca's modifications
  prfloor[UU]=MAX(prfloor[UU],0.00016*pow(prfloor[RHO],gam)/(gam-1.));

  return(retval);
}



