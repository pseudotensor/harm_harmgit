
/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */
#define MAXPASSPARMS 10

#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2
#define THINDISKFROMMATHEMATICA 3
#define THINTORUS 4

// For blandford problem also need to set:
// 0) WHICHPROBLEM 2
// 1) a=0.92
// 2) Rout=1E3
// 3) hslope=0.3
// 4) BSQORHOLIMIT=1E2;
// 5) BSQOULIMIT=1E3;
// 6) UORHOLIMIT=1E3;
// 7) tf=1E4; // and maybe DTd's
// 8) randfact = .2;
// 8.5) Choose fieldtype: FIELDTYPE BLANDFORDQUAD or DISKFIELD (SS dipole)
// 9) lim=PARALINE; FLUXB=FLUXCTSTAG; TIMEORDER=4;
// 10) N1,N2,N3 in init.h
// 11) MAXBND 4
// 12) PRODUCTION 1
// 13) USE<systemtype>=1
// 14) setup batch

//#define WHICHPROBLEM THINDISKFROMMATHEMATICA // choice
#define WHICHPROBLEM THINTORUS // choice

static FTYPE lfunc( FTYPE lin, FTYPE *parms );
static void compute_gu( FTYPE r, FTYPE th, FTYPE a, FTYPE *gutt, FTYPE *gutp, FTYPE *gupp );
static FTYPE compute_udt( FTYPE r, FTYPE th, FTYPE a, FTYPE l );
static FTYPE compute_omega( FTYPE r, FTYPE th, FTYPE a, FTYPE l );
static FTYPE compute_l_from_omega( FTYPE r, FTYPE th, FTYPE a, FTYPE omega );
static FTYPE thintorus_findl( FTYPE r, FTYPE th, FTYPE a, FTYPE k, FTYPE al );
static SFTYPE rhomax=0,umax=0,bsq_max=0; // OPENMPMARK: These are ok file globals since set using critical construct
static SFTYPE beta,randfact,rin; // OPENMPMARK: Ok file global since set as constant before used
static FTYPE rhodisk;


static FTYPE nz_func(FTYPE R) ;

static int read_data(FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR]);

int prepre_init_specific_init(void)
{
  int funreturn;
  
  dofull2pi = 0;

  funreturn=user1_prepre_init_specific_init();
  if(funreturn!=0) return(funreturn);

  return(0);

}


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
#if( WHICHPROBLEM == THINDISKFROMMATHEMATICA )
  h_over_r=0.1;
#elif( WHICHPROBLEM == THINTORUS )
  h_over_r=0.1;
#else
  h_over_r=0.3;
#endif
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  rhodisk=1.0;

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;

  return(0);
}

int set_fieldfrompotential(int *fieldfrompotential)
{
  int pl,pliter;

  // default (assume all fields are from potential)
  PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=1;


  return(0);
}


int init_conservatives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*Utemp)[NSTORE2][NSTORE3][NPR], FTYPE (*U)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  int fieldfrompotential[NDIM];

  set_fieldfrompotential(fieldfrompotential);

  funreturn=user1_init_conservatives(fieldfrompotential, prim,pstag, Utemp, U);
  if(funreturn!=0) return(funreturn);


  return(0);
 
}


int post_init_specific_init(void)
{
  int funreturn;

  funreturn=user1_post_init_specific_init();

  TIMEORDER = 2;
  DTr = 2000;

  if(funreturn!=0) return(funreturn);

  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  return(0);

}


/*

Models to run:

Constant parameters:

1) Rout=1E3 and run for tf=1E4 (so will take 5X longer than compared to Orange run at 128x128x32)

2) BSQORHOLIMIT=1E3, etc.

3) PARALINE, FLUXCTSTAG, TO=4

4) Form of A_\phi fixed

Field parameter studies in 2D axisymmetry at 256^2:

1) H/R=0.3, a=0.9: LS quadrapole,  LS dipole, SS quadrapole, SS dipole

 

Spin parameter study in 2D axisymmetry at 256^2:

 

1) H/R=0.3, LS quadrapole: a=-.999,-.99,-.9,-.5,-0.2,0,.2,.5,.9,.99,.999

H/R parameter study in 2D axisymmetry at 256^2:

1) a=0.9 LS quadrapole with H/R=0.1,0.3,0.9,1.5

2D Fiducial Models:

1) Using a=0.9, H/R=0.3, LS quad and LS dipole, do two 2D fudicial models at: 1024^2

3D Fiducial Models:

1) Using a=0.9, H/R=0.3, LS quadrapole and LS dipole, do two 3D fiducial models at 2 different resolutions: 128x128x32 and 256x256x64

Questions for Roger:

1) Choice for disk thickness?
2) Choice for field shape -- specifically?
3) Choice for flux threading disk vs. BH initially?
4) Ask about BZ77 and residual A_\phi at pole
5) 

*/



int init_defcoord(void)
{
  
  // define coordinate type
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  defcoord = SJETCOORDS;
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  defcoord = REBECCAGRID ;
#elif(WHICHPROBLEM==THINTORUS)
  defcoord = REBECCAGRID ;
#elif(WHICHPROBLEM==GRBJET)
  // define coordinate type
  defcoord = JET4COORDS;
#endif

  return(0);
}


int init_grid(void)
{
  
  // metric stuff first


#if(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  a = 0.;
#elif(WHICHPROBLEM==THINTORUS)
  a = 0.;
#else
  a = 0.95;   //so that Risco ~ 2
#endif

 
  Rhor=rhor_calc(0);

  hslope = 0.13;  //sas: use a constant slope as Jon suggests in the comments
  //hslope = 1.04*pow(h_over_r,2.0/3.0);


#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.8 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3 * Rin;
  Rout = 1.e4;
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.92 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3;
  Rout = 50.;
#elif(WHICHPROBLEM==THINTORUS)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  Rin = 0.92 * Rhor;  //to be chosen manually so that there are 5.5 cells inside horizon to guarantee stability
  R0 = 0.3;
  Rout = 50.;
#elif(WHICHPROBLEM==GRBJET)
	setRin_withchecks(&Rin);
	R0 = -3.0;
  Rout = 1E5;
#endif


  return(0);
}



int init_global(void)
{
  int pl,pliter;
  int funreturn;


  funreturn=user1_init_global();
  if(funreturn!=0) return(funreturn);


  //////////////////
  // overrides for more detailed problem dependence


  TIMEORDER=2; // no need for 4 unless higher-order or cold collapse problem.
  lim[1] = lim[2] = lim[3] = PARALINE; //sas: it's already set in init.tools.c but reset it here just to make sure
  //also need to ensure that in para_and_paraenohybrid.h JONPARASMOOTH is set to 0 (resolves disk best) or 1 (resolves jet best)
#if(  WHICHPROBLEM==THINDISKFROMMATHEMATICA )
  cooling = COOLREBECCATHINDISK; //do Rebecca-type cooling; make sure enk0 is set to the same value as p/rho^\Gamma in the initial conditions (as found in dump0000).
#elif( WHICHPROBLEM==THINTORUS )
  cooling = COOLREBECCATHINDISK; //do Rebecca-type cooling; make sure enk0 is set to the same value as p/rho^\Gamma in the initial conditions (as found in dump0000).
#else
  cooling = NOCOOLING; //no cooling
#endif

  //FLUXB = FLUXCTTOTH;
  FLUXB = FLUXCTSTAG;

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==THINDISKFROMMATHEMATICA \
  || WHICHPROBLEM == THINTORUS)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  //SASMARK: decrease magnetization by 2x to make it easier (still is around ~45>>1)
  BSQORHOLIMIT=0.5*1E2; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=0.5*1E3; // was 1E3 but latest BC test had 1E4
  UORHOLIMIT=0.5*1E3;
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
#elif(WHICHPROBLEM==GRBJET)
  BCtype[X1UP]=FIXEDOUTFLOW;
  BCtype[X1DN]=FREEOUTFLOW;
  rescaletype=4;
  BSQORHOLIMIT=1E3;
  BSQOULIMIT=1E4;
  RHOMIN = 23.0;
  UUMIN = 1.7;
#endif






#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  /* output choices */
  tf = 1e4;

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 50.;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 50.0;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = 2.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 2.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = 50.0;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 2000;

#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS)
  /* output choices */
  tf = 2000.;

  /* dumping frequency, in units of M */
  DTdumpgen[FAILFLOORDUDUMPTYPE]=DTdumpgen[RESTARTDUMPTYPE]=DTdumpgen[RESTARTMETRICDUMPTYPE]=DTdumpgen[GRIDDUMPTYPE]=DTdumpgen[DEBUGDUMPTYPE]=DTdumpgen[ENODEBUGDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[DISSDUMPTYPE]=DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE] = 100.;
  DTdumpgen[AVG1DUMPTYPE]=DTdumpgen[AVG2DUMPTYPE]= 100.0;
  // ener period
  DTdumpgen[ENERDUMPTYPE] = 10.0;
  /* image file frequ., in units of M */
  DTdumpgen[IMAGEDUMPTYPE] = 10.0;
  // fieldline locked to images so can overlay
  DTdumpgen[FIELDLINEDUMPTYPE] = DTdumpgen[IMAGEDUMPTYPE];

  /* debug file */  
  DTdumpgen[DEBUGDUMPTYPE] = 100.0;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 2000;

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

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int funreturn;

  funreturn=user1_init_atmosphere(whichvel, whichcoord,i, j, k, pr);
  if(funreturn!=0) return(funreturn);

  return(0);

}

int init_grid_post_set_grid(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);



  beta = 1.e2 ;
  randfact = 4.e-2; //sas: as Jon used for 3D runs but use it for 2D as well

#if(WHICHPROBLEM==NORMALTORUS)
  //rin = Risco;
  rin = 6. ;
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  rin = 20. ;
#elif(WHICHPROBLEM==THINTORUS)
  rin = 20. ;
#elif(WHICHPROBLEM==KEPDISK)
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;
#elif(WHICHPROBLEM==GRBJET)
  rin = Risco;
#endif

  


#if( ANALYTICMEMORY == 1 && WHICHPROBLEM != THINDISKFROMMATHEMATICA )
  //SASMARK restart: need to populate panalytic with IC's; DO NOT do this 
  //when reading the ICs in from a file since then need to carry the file around
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives (but use ulast so don't overwrite unew read-in from file)
    MYFUN(init_primitives(panalytic,pstaganalytic,GLOBALPOINT(utemparray),vpotanalytic,Bhatanalytic,panalytic,pstaganalytic,vpotanalytic,Bhatanalytic,F1,F2,F3,Atemp),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }
#endif
  // check rmin
  check_rmin();


  // check that singularities are properly represented by code
  check_spc_singularities_user();

  
  return(0);

}



int init_primitives(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR], FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*pstaganalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*vpotanalytic)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhatanalytic)[NSTORE2][NSTORE3][NPR], FTYPE (*F1)[NSTORE2][NSTORE3][NPR], FTYPE (*F2)[NSTORE2][NSTORE3][NPR], FTYPE (*F3)[NSTORE2][NSTORE3][NPR], FTYPE (*Atemp)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int funreturn;

#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA ) 
  //read initial conditions from input file for the Mathematica-generated thin-disk ICs
#if(ANALYTICMEMORY==0)
#error Cannot do THINDISKFROMMATHEMATICA problem with ANALYTICMEMORY==0.  Please set ANALYTICMEMORY = 1.
#endif
  read_data(panalytic);
#endif

  funreturn=user1_init_primitives(prim, pstag, ucons, vpot, Bhat, panalytic, pstaganalytic, vpotanalytic, Bhatanalytic, F1, F2, F3,Atemp);
  if(funreturn!=0) return(funreturn);

  return(0);


}


//reads in ICs from a file and stores them in panalytic array for future use in init_dsandvels
int read_data(FTYPE (*panalytic)[NSTORE2][NSTORE3][NPR])
{
  int i,j, k;
  SFTYPE tabr,tabo;
  SFTYPE kappa;
  FILE * fpr;
  int ti, tj , tk;
  int procno;
  int nscanned1;
  long lineno, numused;

  //BOBMARK: should be the value of KK in mathematica file
  kappa = 0.01 ;


#if( USEMPI )
  for( procno = 0; procno < numprocs; procno++ ) {
    if( myid == procno ){
#endif
      //////
      //
      //  Open file for reading
      trifprintf("myid is %d\n",myid);
      fpr=fopen("./mytest", "r");
      

      if( NULL != fpr ) {
        lineno = 0;
        numused = 0;

	while( (nscanned1 = fscanf(fpr, "%d %d %d %lg %lg \r\n", &ti, &tj, &tk , &tabr, &tabo )) == 5 ) {

          ///////
          //
          //  Check if within range and if within, save it to the array
          //  i, j - local grid indices for this processor, ti, tj - global

          // trifprintf("startpos1=%d, startpos2=%d \n", startpos[1], startpos[2]);

          i = ti - 1 - startpos[1];
          j = tj - 1 - startpos[2];
          k = tk - 1 - startpos[3];

          if( i >= INFULL1 && i <= OUTFULL1 &&
              j >= INFULL2 && j <= OUTFULL2 && k>=INFULL3 && k<=OUTFULL3 ) {


	    //BOBMARK: converted to new, macrofy-d format
	    //	    p[i][j][k][RHO] = tabr;
	    //	    p[i][j][k][U3] = tabo;
	    MACP0A1(panalytic,i,j,k,RHO) = tabr;
	    MACP0A1(panalytic,i,j,k,UU) = kappa * pow(tabr, gam) / (gam - 1.);

	    MACP0A1(panalytic,i,j,k,U1) = 0.0;
	    MACP0A1(panalytic,i,j,k,U2) = 0.0;
	    MACP0A1(panalytic,i,j,k,U3) = tabo;
	    
	    MACP0A1(panalytic,i,j,k,B1) = 0.0;
	    MACP0A1(panalytic,i,j,k,B2) = 0.0;
	    MACP0A1(panalytic,i,j,k,B3) = 0.0;


	    trifprintf(" k = %d , i=%d, j=%d, rho=%lg\n", k, i, j, tabr);

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

	    
	    //BOBMARK: converted to new, macrofy-d format
	    //	    p[i][j][k][RHO] = tabr;
	    //	    p[i][j][k][U3] = tabo;
	    MACP0A1(panalytic,i,j,k,RHO) = tabr;
	    MACP0A1(panalytic,i,j,k,UU) = kappa * pow(tabr, gam) / (gam - 1.);

	    MACP0A1(panalytic,i,j,k,U1) = 0.0;
	    MACP0A1(panalytic,i,j,k,U2) = 0.0;
	    MACP0A1(panalytic,i,j,k,U3) = tabo;
	    
	    MACP0A1(panalytic,i,j,k,B1) = 0.0;
	    MACP0A1(panalytic,i,j,k,B2) = 0.0;
	    MACP0A1(panalytic,i,j,k,B3) = 0.0;
	    
	    
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

  }
  return(0);
}








int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindiskfrommathematica(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

#if(WHICHPROBLEM==NORMALTORUS)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==THINDISKFROMMATHEMATICA)
  return(init_dsandvels_thindiskfrommathematica(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#endif

}


// unnormalized density
int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeomdontuse;
  struct of_geom *ptrrealgeom=&realgeomdontuse;
  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];
  int pl,pliter;






  kappa = 1.e-3 ;
  rmax = 12. ;
  l = lfish_calc(rmax) ;
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];



  sth = sin(th);
  cth = cos(th);

  /* calculate lnh */
  DD = r * r - 2. * r + a * a;
  AA = (r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth;
  SS = r * r + a * a * cth * cth;
  
  thin = M_PI / 2.;
  sthin = sin(thin);
  cthin = cos(thin);
  DDin = rin * rin - 2. * rin + a * a;
  AAin = (rin * rin + a * a) * (rin * rin + a * a)
    - DDin * a * a * sthin * sthin;
  SSin = rin * rin + a * a * cthin * cthin;
  
  if (r >= rin) {
    lnh = 0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
			       (AA * sth * AA * sth))) / (SS * DD /
							  AA))
      - 0.5 * sqrt(1. +
		   4. * (l * l * SS * SS) * DD / (AA * AA * sth *
						  sth))
      - 2. * a * r * l / AA -
      (0.5 *
       log((1. +
	    sqrt(1. +
		 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin *
						      sthin *
						      sthin))) /
	   (SSin * DDin / AAin))
       - 0.5 * sqrt(1. +
		    4. * (l * l * SSin * SSin) * DDin / (AAin *
							 AAin *
							 sthin *
							 sthin))
       - 2. * a * rin * l / AAin);
  } else
    lnh = 1.;
  

  
  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {


    get_geometry(i, j, k, CENT, ptrrealgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrrealgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  /* region inside magnetized torus; u^i is calculated in
     Boyer-Lindquist coordinates, as per Fishbone & Moncrief, so it
     needs to be transformed at the end */
  else {
    hm1 = exp(lnh) - 1.;
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));
    u = kappa * pow(rho, gam) / (gam - 1.);
    ur = 0.;
    uh = 0.;
    
    /* calculate u^phi */
    expm2chi = SS * SS * DD / (AA * AA * sth * sth);
    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5));
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }
}



int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th,ph;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;
  int pl,pliter;


  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  ph=V[3];




  /* region outside disk */
  R = r*sin(th) ;

  if(R < rin) {

    get_geometry(i, j, k, CENT, ptrgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  else {

    H = h_over_r*R ;
    nz = nz_func(R) ;
    z = r*cos(th) ;
    S = 1./(H*H*nz) ;
    cs = H*nz ;

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H)) * taper_func(R,rin) ;
    u = rho*cs*cs/(gam - 1.) ;
    ur = 0. ;
    uh = 0. ;
    up = 1./(pow(r,1.5) + a) ;
    // solution for 3-vel


    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5));

    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    if(FLUXB==FLUXCTSTAG){
      PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
    }

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
  }
}


// unnormalized density
int init_dsandvels_thindiskfrommathematica(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int pl;

    PALLLOOP(pl) pr[pl] = MACP0A1(GLOBALPOINT(panalytic),i,j,k,pl);

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
}


FTYPE lfunc( FTYPE lin, FTYPE *parms )
{    
   FTYPE gutt, gutp, gupp, al, k;
   FTYPE ans;
   
   gutt = parms[0];
   gutp = parms[1];
   gupp = parms[2];
   al = parms[3]; 
   k = parms[4];
   
   ans = (gutp - lin * gupp)/( gutt - lin * gutp) - k *pow(lin,al);
   
   return(ans);
}
    
void compute_gu( FTYPE r, FTYPE th, FTYPE a, FTYPE *gutt, FTYPE *gutp, FTYPE *gupp )
{
  //metric (expressions taken from eqtorus_c.nb):
  *gutt = -1 - 4*r*(pow(a,2) + pow(r,2))*
      pow((-2 + r)*r + pow(a,2),-1)*
      pow(pow(a,2) + cos(2*th)*pow(a,2) + 2*pow(r,2),-1);
      
  *gutp = -4*a*r*pow((-2 + r)*r + pow(a,2),-1)*
      pow(pow(a,2) + cos(2*th)*pow(a,2) + 2*pow(r,2),-1);
    
  *gupp = 2*((-2 + r)*r + pow(a,2)*pow(cos(th),2))*
      pow(sin(th),-2)*pow((-2 + r)*r + pow(a,2),-1)*
      pow(pow(a,2) + cos(2*th)*pow(a,2) + 2*pow(r,2),-1);
}  
  
FTYPE thintorus_findl( FTYPE r, FTYPE th, FTYPE a, FTYPE k, FTYPE al )
{
  FTYPE gutt, gutp, gupp;
  FTYPE parms[5];
  FTYPE l;
  
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  //store params in an array before function call
  parms[0] = gutt;
  parms[1] = gutp;
  parms[2] = gupp;
  parms[3] = al; 
  parms[4] = k;
  
  //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3) 
  //demand accuracy 5x machine prec.
  //in non-rel limit l_K = sqrt(r), use 10x that as the upper limit:
  l = rtbis( lfunc, parms, 1e-7, 10.*sqrt(Rout), 5.*DBL_EPSILON );
  
  return( l );
}

FTYPE compute_udt( FTYPE r, FTYPE th, FTYPE a, FTYPE l )
{
  FTYPE gutt, gutp, gupp;
  FTYPE udt;
   
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  udt = -sqrt(- 1/(gutt - 2 * l * gutp + l * l * gupp));
  
  return( udt );
}
    
    
FTYPE compute_omega( FTYPE r, FTYPE th, FTYPE a, FTYPE l )
{
  FTYPE gutt, gutp, gupp;
  FTYPE omega;
   
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  omega = (gutp - gupp*l)*pow(gutt - gutp*l,-1);
  
  return( omega );
}
    
FTYPE compute_l_from_omega( FTYPE r, FTYPE th, FTYPE a, FTYPE omega )
{
  FTYPE gutt, gutp, gupp;
  FTYPE l;
   
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  l = (gutp - omega * gutt)/(gupp - omega * gutp);
  
  return( l );
}
    
int init_dsandvels_thintorus(int *whichvel, int*whichcoord, int ti, int tj, int tk, FTYPE *pr, FTYPE *pstag)
{
  FTYPE kappa, n, rmax; //kappa <-> KK in mathematica file
  FTYPE r, th, omk, lk, c;
  FTYPE k, al;
  FTYPE X[NDIM], V[NDIM];
  FTYPE lin, l;
  FTYPE utin, flin;
  FTYPE udt, hh, f, eps;
  FTYPE rho, u, om, ur, uh, up;
  int pl;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;

  coord(ti, tj, tk, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  
  /* region outside disk? */
  R = r*sin(th) ;

  //avoid computations way outside the torus: fill the region with atmosphere
  if( R < rin ) {

    get_geometry(ti, tj, tk, CENT, ptrgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }

  ///
  /// Parameters
  ///
 
  kappa = 0.01;
  n = 2. - 1.65; 
  rmax = 15.;

  
  ///
  /// Computations at pressure max
  ///
  
  r = rmax;
  th = M_PI_2l;
  omk = 1./(pow(rmax,1.5) + a);
  lk = compute_l_from_omega( r, th, a, omk );
  //log(omk) == (2/n) log(c) + (1 - 2./n) log(lk) <-- solve for c:
  c = pow( lk, 1 - n/2. ) * pow( omk, n/2. );
  
  k = pow(c, 2./n);
  al = (n - 2.)/n;
  
  
  ///
  /// Computations at torus inner edge
  ///
  
  //l = lin at inner edge, r = rin
  r = rin;
  th = M_PI_2l;
  lin = thintorus_findl( r, th, a, k, al );
      
  //finding DHK03 lin, utin, f (lin)
  utin = compute_udt( r, th, a, lin );
  flin = pow(fabs(1 - k*pow(lin,1 + al)),pow(1 + al,-1));
  

  ///
  /// Computations at current point: r, th
  ///
  
  coord(ti, tj, tk, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  
  //l at current r, th
  l = thintorus_findl( r, th, a, k, al );
  
  
  udt = compute_udt( r, th, a, l );
  f = pow(fabs(1 - k*pow(l,1 + al)),pow(1 + al,-1));
  hh = flin*utin*pow(f,-1)*pow(udt,-1);
  eps = (-1 + hh)*pow(gam,-1);
  
  //avoid computations outside the torus: fill the region with atmosphere
  if( eps < 0 ) {
    get_geometry(ti, tj, tk, CENT, ptrgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
   
  rho = pow((-1 + gam)*eps*pow(kappa,-1),pow(-1 + gam,-1));
  u = kappa * pow(rho, gam) / (gam - 1.);
  om = compute_omega( r, th, a, l );
  
  ur = 0.;
  uh = 0.;
  up = om; // = u^phi/u^t
 
  
  pr[RHO] = rho ;
  pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5) );

  pr[U1] = ur ;
  pr[U2] = uh ;    
  pr[U3] = up;

  if(FLUXB==FLUXCTSTAG){
    PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;
  }

  *whichvel=VEL3;
  *whichcoord=BLCOORDS;
  return(0);

}




#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define BLANDFORDQUAD 3

//#define FIELDTYPE BLANDFORDQUAD
#define FIELDTYPE DISKFIELD

FTYPE setgpara(FTYPE myr, FTYPE th, FTYPE thpower)
{
  FTYPE fneg,fpos;
  FTYPE gpara;

  fneg=1.0-pow(cos(th),thpower);
  fpos=1.0+pow(cos(th),thpower);
  gpara=0.5*(myr*fneg + 2.0*fpos*(1.0-log(fpos)));
  // remove BZ77 Paraboloidal divb!=0 at pole
  gpara=gpara-2.0*(1.0-log(2.0));

  return(gpara);


}

FTYPE setblandfordfield(FTYPE r, FTYPE th)
{
  FTYPE setgpara(FTYPE myr, FTYPE th, FTYPE thpower);
  FTYPE rshift,myr,rpower,myz,myR,myvert;
  FTYPE thother,thpower,gparalow,gparahigh,mygpara;
  FTYPE aphi;


  rshift=4.0;
  rpower=0.75;
  thpower=4.0;
  

  myr=pow(r+rshift,rpower);
  myz=myr*cos(th);
  myR=myr*sin(th);
  myvert = (th>M_PI*0.5) ? (myr*sin(th)) : (myr*sin(-th));

  thother=M_PI-th;
  gparalow=setgpara(myr,th,thpower);
  gparahigh=setgpara(myr,thother,thpower);
  mygpara=(th<0.5*M_PI) ? gparalow : gparahigh;

  // GOOD:
  // aphi=mygpara;
  // aphi=mygpara*cos(th); // B1 diverges at pole
  //  aphi=mygpara*cos(th)*sin(th); // doesn't diverge as much
  //  aphi=mygpara*cos(th)*sin(th)*sin(th); // old choice before subtracted original BZ77 problem
  aphi=mygpara*cos(th); // latest choice
  //aphi=myvert*cos(th); // vert with quad
  //aphi=myR*cos(th);

  // BAD:
  // aphi=myvert;
  
  

  return(aphi);


}

// assumes normal field in pr
int init_vpot_user(int *whichcoord, int l, int i, int j, int k, int loc, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, u_av, q;
  FTYPE r,th;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);
#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
  FTYPE fieldhor;
#endif



  vpot=0.0;


  if(l==3){// A_\phi

    r=V[1];
    th=V[2];



    // Blandford quadrapole field version
    if(FIELDTYPE==BLANDFORDQUAD){
      vpot += setblandfordfield(r,th);
    }

    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      FTYPE rpow;
      rpow=3.0/4.0; // Using rpow=1 leads to quite strong field at large radius, and for standard atmosphere will lead to \sigma large at all radii, which is very difficult to deal with -- especially with grid sectioning where outer moving wall keeps opening up highly magnetized region
      vpot += 0.5*pow(r,rpow)*sin(th)*sin(th) ;
    }


    /* field-in-disk version */
    
    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)){

#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
#define STARTFIELD (1.1*rin)
      fieldhor=0.28;
#endif
      // average of density that lives on CORN3
      // since init_vpot() is called for all i,j,k, can't use
      // non-existence values, so limit averaging:
      if((i==-N1BND)&&(j==-N2BND)){
	rho_av = MACP0A1(prim,i,j,k,RHO);
        u_av = MACP0A1(prim,i,j,k,UU);
      }
      else if(i==-N1BND){
	rho_av = AVGN_2(prim,i,j,k,RHO);
        u_av = AVGN_2(prim,i,j,k,UU);
      }
      else if(j==-N2BND){
	rho_av = AVGN_1(prim,i,j,k,RHO);
        u_av = AVGN_1(prim,i,j,k,UU);
      }
      else{ // normal cells
	rho_av = AVGN_for3(prim,i,j,k,RHO);
	u_av=AVGN_for3(prim,i,j,k,UU);
      }

#if( WHICHPROBLEM==THINDISKFROMMATHEMATICA || WHICHPROBLEM == THINTORUS ) 
      //SASMARK: since u was randomly perturbed, may need to sync the u across tiles to avoid monopoles
      if(r > STARTFIELD) q = ((u_av/umax) - 0.2)*pow(r,0.75) ;
      else q = 0. ;
      trifprintf("rhoav=%g q=%g\n", rho_av, q);

      if(q > 0.){
       vpot += q*q*sin(log(r/STARTFIELD)/fieldhor)* (1. + 0.02 * (ranc(0,0) - 0.5))  ;
      }
#else
      q = rho_av / rhomax - 0.2;
      if (q > 0.)      vpot += q;
#endif
    
    }
  }

  //////////////////////////////////
  //
  // finally assign what's returned
  //
  //////////////////////////////////
  *A = vpot;
  *whichcoord = MCOORD;



  return(0);

}



int init_vpot2field_user(FTYPE (*A)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3],FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  int fieldfrompotential[NDIM];

  funreturn=user1_init_vpot2field_user(fieldfrompotential, A, prim, pstag, ucons, Bhat);
  if(funreturn!=0) return(funreturn);
 

  return(0);


}



// assumes we are fed the true densities
int normalize_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;
  FTYPE parms[MAXPASSPARMS];
  int eqline;

  eqline=1;
  parms[0]=rin;
  parms[1]=rhodisk;

  funreturn=user1_normalize_densities(eqline, parms, prim, &rhomax, &umax);
  if(funreturn!=0) return(funreturn);
 

  return(0);
}



// assumes we are fed the true densities
int getmax_densities(FTYPE (*prim)[NSTORE2][NSTORE3][NPR],SFTYPE *rhomax, SFTYPE *umax)
{
  int funreturn;

  funreturn=user1_getmax_densities(prim,rhomax, umax);
  if(funreturn!=0) return(funreturn);
 
  return(0);
}



// get maximum b^2 and p_g
int get_maxes(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE *bsq_max, FTYPE *pg_max, FTYPE *beta_min)
{
  int funreturn;
  int eqslice;
  FTYPE parms[MAXPASSPARMS];
  
  if(FIELDTYPE==VERTFIELD || FIELDTYPE==BLANDFORDQUAD){
    eqslice=1;
  }
  else{
    eqslice=0;
  }

  parms[0]=rin;

  funreturn=user1_get_maxes(eqslice, parms,prim, bsq_max, pg_max, beta_min);
  if(funreturn!=0) return(funreturn);
 
  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], FTYPE (*pstag)[NSTORE2][NSTORE3][NPR], FTYPE (*ucons)[NSTORE2][NSTORE3][NPR], FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3], FTYPE (*Bhat)[NSTORE2][NSTORE3][NPR])
{
  int funreturn;

 
  funreturn=user1_normalize_field(beta, prim, pstag, ucons, vpot, Bhat);
  if(funreturn!=0) return(funreturn);
 
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
								     pow(r,
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
  int funreturn;
  int atmospheretype;

  if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK || WHICHPROBLEM==THINDISKFROMMATHEMATICA){
    atmospheretype=1;
  }
  else if(WHICHPROBLEM==GRBJET){
    atmospheretype=2;
  }
  else {
    atmospheretype=1; // default
  }
 
  funreturn=user1_set_atmosphere(atmospheretype, whichcond, whichvel, ptrgeom, pr);
  if(funreturn!=0) return(funreturn);
 
  return(0);

}



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  int funreturn;
  
  funreturn=set_density_floors_default(ptrgeom, pr, prfloor);
  if(funreturn!=0) return(funreturn);

  return(0);
}




FTYPE nz_func(FTYPE R)
{
  return(
	 sqrt(
	      (3.*a*a - 4.*a*sqrt(R) + R*R)/
	      pow(R*(a + pow(R,1.5)),2)
	      )
	 ) ;


}



// Setup problem-dependent grid sectioning
int theproblem_set_enerregiondef(int forceupdate, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int (*enerregiondef)[NDIM] )
{

  // Torus problem
  //  torus_set_enerregiondef(forceupdate, timeorder, numtimeorders, thenstep, thetime, enerregiondef);
  //jet_set_enerregiondef(forceupdate, timeorder, numtimeorders, thenstep, thetime, enerregiondef);

#if(1)
  enerregiondef[POINTDOWN][1]=0;
  enerregiondef[POINTUP][1]=totalsize[1]-1;
  enerregiondef[POINTDOWN][2]=0;
  enerregiondef[POINTUP][2]=totalsize[2]-1;
  enerregiondef[POINTDOWN][3]=0;
  enerregiondef[POINTUP][3]=totalsize[3]-1;
#endif

  return(0);
}


int theproblem_set_enerregionupdate(int forceupdate, int timeorder, int numtimeorders, long int thenstep, FTYPE thetime, int *updateeverynumsteps, int *everynumsteps)
{

  ////////////////////
  //
  // Setup update period
  //
  ///////////////////

  ////////
  //
  // number of steps after which position/size of active section is updated
  //
  ////////
#if(N3==1)
  //  *updateeverynumsteps=100;
  *updateeverynumsteps=1; // update every step since otherwise flow runs into wall at outer boundary
#else
  //  *updateeverynumsteps=10;
  *updateeverynumsteps=1; // update every step since otherwise flow runs into wall at outer boundary
#endif

  ////////
  //
  //number of steps after which position/size of active section is reported to file
  //
  ////////
  *everynumsteps = *updateeverynumsteps*100;

  return(0);
}


// specify MPI task rank ordering
// example user-dependent code
int theproblem_set_myid(void)
{
  int retval;
 
  // default is to do nothing
  //  retval=jet_set_myid();
  retval=0;

  // do other things?

  return(retval);

}


