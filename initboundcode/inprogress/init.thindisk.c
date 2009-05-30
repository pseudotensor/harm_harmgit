
/* 
 *
 thin disk solution
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"

#define SLOWFAC 1.0		/* reduce u_phi by this amount */

FTYPE nz_func(FTYPE R) ;


FTYPE rin ;
SFTYPE rhomaxold,umaxold,rhomax=0,umax=0,bsq_max=0,beta;

int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.05;
  h_over_r_jet=h_over_r*2.0;

  return(0);
}

int post_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations

  return(0);
}

int init_grid(void)
{
  // set metric parameters first
  a = 0.9375 ;

  // choose coordinate type
  //  defcoord = 5;
  defcoord = 0;



  // make changes to primary coord parms
  R0 = 0.0; // R0=-1.0 is bad if it leads to unresolved horizon
  Rhor=rhor_calc(0); // needed here
  Rout = 40.;

  Rin=setRin(setihor());
  //  Rin = 0.98 * Rhor;

  hslope = 1.04*pow(h_over_r,2.0/3.0);

  return(0);
}

int init_global(void)
{

  ranc(0); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=1;

  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

  /* output choices */
  tf = 2000.0;

  DTd = DTavg = 50.;			/* dumping frequency, in units of M */
  DTener = 2.0;			/* logfile frequency, in units of M */
  DTi = 2.0;			/* image file frequ., in units of M */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;			/* restart file period in steps */
  DTdebug = 50.0;
  return(0);

}


int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, FTYPE *pr)
{
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;

  

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

  //  beta=1.e2;
  beta=1.e2/3.0;
  randfact=1E-3;
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;



  /* region outside disk */
  R = r*sin(th) ;

  if(R < rin) {

    get_geometry(i, j, CENT, &realgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

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

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H))
      * taper_func(R,rin) ;
    u = rho*cs*cs/(gam - 1.) ;
    ur = 0. ;
    uh = 0. ;
    up = 1./(pow(r,1.5) + a) ;
    // solution for 3-vel


    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0) - 0.5));

    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;
    

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
  }
}

#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2

#define FIELDTYPE DISKFIELD
// 0: disk only
// 1: vert only
// 2: vert+disk field

int init_vpot(int i, int j,FTYPE *A)
{
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE u_av, u_ref,ftemp, q;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;
  FTYPE fieldhor;


  fieldhor=3.0*h_over_r;


  coord(i, j, CORN, X);
  bl_coord(X, &r, &th);

  *A=0;

  // vertical field
  if((FIELDTYPE==1)||(FIELDTYPE==2)){
   *A += 0.5*r*sin(th) ;
  }

  if((FIELDTYPE==0)||(FIELDTYPE==2)){
    R = r*sin(th) ;

    u_av = 0.25*(
		 p[i][j][UU] +
		 p[i-1][j][UU] +
		 p[i][j-1][UU] +
		 p[i-1][j-1][UU]) ;

    // in MPI, don't have "equatorial value" for all cpus so easily, so just compute it
    // just replaced R with r and remaining th with Pi/2
    H = h_over_r*r ;
    nz = nz_func(r) ;
    z = 0.0 ;
    S = 1./(H*H*nz) ;
    cs = H*nz ;

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H))
      * taper_func(r,rin) ;
    u_ref = rho*cs*cs/(gam - 1.) ;
    u_ref/=rhomaxold;
    /*
      u_ref = 0.25*(
      p[i][N2/2][UU] +
      p[i-1][N2/2][UU] +
      p[i][N2/2-1][UU] +
      p[i-1][N2/2-1][UU]) ;
    */
#define STARTFIELD (1.1*rin)

    if(r > STARTFIELD) q = ((u_av/u_ref) - 0.2)*pow(r,0.25) ;
    else q = 0. ;

    if(q > 0.){
      *A += q*q*sin(log(r/STARTFIELD)/fieldhor)*taper_func(r,STARTFIELD);
      //trifprintf("%d %d u_ref=%21.15g A=%21.15g\n",i,j,u_ref,A[i][j]);
    }
  }

  return(0);

}




int init_vpot2field(FTYPE A[][N2M],FTYPE pr[][N2M][NPR])
{
  extern int vpot2field(FTYPE A[][N2M],FTYPE p[][N2M][NPR]);


  return(vpot2field(A,pr));

}


// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  set_atmosphere(0,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  if(pr[RHO]<pratm[RHO]){
    PLOOP(pl) pr[pl]=pratm[pl];
  }
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);



}




int normalize_densities(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE X[NDIM],r,th;


  rhomax=0;
  umax=0;
  ZLOOP{
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);

    if (p[i][j][RHO] > rhomax)   rhomax = p[i][j][RHO];
    if (p[i][j][UU] > umax && r > rin)    umax = p[i][j][UU];
  }

  mpimax(&rhomax);
  mpimax(&umax);
  trifprintf("rhomax: %21.15g umax: %21.15g\n", rhomax, umax);

  rhomaxold=rhomax;
  umaxold=umax;

  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    p[i][j][RHO] /= rhomax;
    p[i][j][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

  return(0);
}


int normalize_field(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE r,th,X[NDIM];

  bsq_max = 0.;
  ZLOOP {
    get_geometry(i, j, CENT, &geom);    

    if(FIELDTYPE==VERTFIELD){
      coord(i, j, CENT, X);
      bl_coord(X, &r, &th);
      
      if((r>rin)&&(fabs(th-M_PI*0.5)<4.0*M_PI*dx[2]*hslope)){
	if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
	  FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
	
	if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
      }
    }
    else{
      if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
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
  
  bsq_max = 0.;
  ZLOOP {
    p[i][j][B1] *= norm;
    p[i][j][B2] *= norm;
    p[i][j][B3] *= norm;

    get_geometry(i, j, CENT, &geom);
    if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
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



FTYPE nz_func(FTYPE R)
{
  return(
	 sqrt(
	      (3.*a*a - 4.*a*sqrt(R) + R*R)/
	      pow(R*(a + pow(R,1.5)),2)
	      )
	 ) ;


}
