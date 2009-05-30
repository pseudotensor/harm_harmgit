
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

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin;

int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.2;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  return(0);
}

int post_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations

  return(0);
}

int init_grid(void)
{
  SFTYPE rh;
  
  // metric stuff first
  a = 0.9375 ;
  
  // define coordinate type
  defcoord = 9;

  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = -3.0;
  Rhor=rhor_calc(0);
  Rout = 1E3;
 
  Rin=setRin(setihor());
  //  Rin = 0.98 * Rhor;
  
  hslope = 0.3;

  return(0);
}

int init_global(void)
{

  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;

  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

  /* output choices */
  tf = 1E4;

  DTd = 25;			/* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = 2;			/* logfile frequency, in units of M */
  DTi = 10;			/* image file frequ., in units of M */
  DTdebug = DTd; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;			/* restart file period in steps */

  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr)
{
  int k;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];

  
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}


// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, FTYPE *pr)
{
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],r,th;
  struct of_geom realgeom,geom;
  

  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax;
  SFTYPE rh;
  //  FTYPE pratm[NPR];
  int k;


  rin = 6. ;
  rmax = 12. ;
  kappa = 1.e-3 ;
  beta = 1.e-15 ;
  randfact = 4.e-2;
  

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

  sth = sin(th);
  cth = cos(th);

  get_geometry(i, j, CENT, &realgeom); // true coordinate system

  pr[U1]=pr[U2]=pr[U3]=0.0;
  //  if(th<M_PI*0.5)  pr[B1]=1.0*gdet[horizoni][j][CENT]/(gdet[i][j][CENT]);
  //  else pr[B1]=-1.0*gdet[horizoni][j][CENT]/(gdet[i][j][CENT]);
  pr[B2]=pr[B3]=0;
  pr[B1]=1.0*gdet[horizoni][N2-1-j][CENT]/(gdet[i][N2-1-j][CENT]);

  // Ruben's talk says they set $\dF^{tr} = C\sin{\theta}/\detg$.

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);
}


#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2

#define FIELDTYPE VERTFIELD

// assumes normal field in pr
int init_vpot(int i, int j,FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;


  return(0);


  coord(i, j, CORN, X);
  bl_coord(X, &r, &th);


  *A=0;

  /* vertical field version*/
  if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
    get_geometry(i,j,CORN,&geom);
    *A = -1.0*cos(th) ;
  }
  /* field-in-disk version */

  if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)){
    rho_av = 0.25 * (p[i][j][RHO] +
                     p[i - 1][j][RHO] +
                     p[i][j - 1][RHO] + p[i - 1][j - 1][RHO]);

    q = rho_av / rhomax - 0.2;

    if (q > 0.)      *A += q;
  }

  return(0);

}


int init_vpot2field(SFTYPE A[][N2+1],FTYPE pr[][N2M][NPR])
{
  extern int vpot2field(SFTYPE A[][N2+1],FTYPE p[][N2M][NPR]);

  return(vpot2field(A,pr));
}

// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE X[NDIM],r,th;

  return(0);

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


  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    p[i][j][RHO] /= rhomax;
    p[i][j][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM];
  FTYPE r,th;

  return(0);

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

  /* this is a little test calculation with a radial field, designed to 
     make the scheme fail */

  /*     Br0 = 1.0 ; ZLOOP { GSET(i,j,CENT) p[i][j][B1] = Br0/(rcurr*rcurr) 
     ; p[i][j][B2] = 0. ; p[i][j][B3] = 0. ; } */


#undef SLOWFAC
