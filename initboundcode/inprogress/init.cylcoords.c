
/* 
 *
 NS disk or not and various NS field setups
 *
 */

#include "decs.h"


FTYPE nz_func(FTYPE R) ;
int init_donut(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr);
int init_kepdisk(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr);
int init_nodisk(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr);


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

SFTYPE rin,rhomaxold,umaxold,rhomax=0,umax=0,bsq_max=0,beta;

int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=1.0;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=h_over_r;

  // see avery_NS_GRB.nb

  // polar field for NS in units of $\sqrt{\rho_{disk} c^2}$

  // weak field
  //Bpole=2.1E-6;
  // strong field
  Bpole=1.0; 

  Omegastar=0.0216; // in units of per GM/c^3

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

  // define metric parms
  a = 0.0 ;

  // define coord type
  defcoord = 101;


  // make changes to coord. parms
  Rin=1E-3; // in GM/c^2, NS surface

  //R0 = 0.0;
  //  R0=0.0;
  // see ns-boundarylayer.nb
  R0 = 0.903*Rin; // accounts for boundary layer with 20 zones above surface for N1=256 Rout=100GM/c^2

  Rhor=rhor_calc(0);
  Rout = 400.;
 

  //  Rin=setRin(setihor());
  //  Rin = 0.98 * Rhor;
  
  //  hslope = 0.5;
  hslope = 1.04*pow(h_over_r,2.0/3.0);


  return(0);
}

int init_global(void)
{

  // for NS:


  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;

  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=OUTFLOW;
  BCtype[X2DN]=OUTFLOW;

  /* output choices */
  tf = 2000.0;

  DTd = 1;			/* dumping frequency, in units of M */
  DTavg = DTd;
  DTener = DTd/1E5;			/* logfile frequency, in units of M */
  DTi = DTd;			/* image file frequ., in units of M */
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

#define KEPDISK 0
#define DONUTDISK 1
#define NODISK 2

#define DISKTYPE NODISK

// unnormalized density
int init_dsandvels(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr)
{

  if(DISKTYPE==KEPDISK){
    return(init_kepdisk(whichvel, whichcoord,i,j,pr));
  }
  else if(DISKTYPE==DONUTDISK){
    return(init_donut(whichvel, whichcoord,i,j,pr));
  }
  else if(DISKTYPE==NODISK){
    return(init_nodisk(whichvel, whichcoord,i,j,pr));
  }

  return(1); // shouldn't reach here
}



int init_donut(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr)
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
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];
  int k;


  rin = 6. ;
  rmax = 20. ;
  l = lfish_calc(rmax) ;
  kappa = 1.e-3 ;
  beta = 1.e2 ;
  randfact = 4.e-2;
  

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

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
  

  get_geometry(i, j, CENT, &realgeom); // true coordinate system
  //  set_atmosphere(0,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system
  
  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {

    pr[RHO]=1E-30;
    // small density will indicate to atmosphere to change it

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
    pr[UU] = u* (1. + randfact * (ranc(0) - 0.5));
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;
 
    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }


}



int init_kepdisk(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr)
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

  beta=1.e2;
  randfact=1E-3;
  //rin = (1. + h_over_r)*Risco;
  //rin=Risco;
  rin = 6;



  /* region outside disk */
  R = r*sin(th) ;

  if(R < rin) {

    rho = 1.e-7*RHOMIN ;
    u = 1.e-7*UUMIN ;

    get_geometry(i, j, CENT, &geom); // true coordinate system

    // normal observer velocity
    if(WHICHVEL==VEL4){
      ur = -geom.gcon[0][1]/sqrt(-geom.gcon[0][0]) ;
      uh = -geom.gcon[0][2]/sqrt(-geom.gcon[0][0]) ;
      up = -geom.gcon[0][3]/sqrt(-geom.gcon[0][0]) ;
    }
    else if(WHICHVEL==VEL3){
      ur = geom.gcon[0][1]/geom.gcon[0][0] ;
      uh = geom.gcon[0][2]/geom.gcon[0][0] ;
      up = geom.gcon[0][3]/geom.gcon[0][0] ;
    }
    else if(WHICHVEL==VELREL4){
      ur = 0.0;
      uh = 0.0;
      up = 0.0;
    }

    pr[RHO] = rho;
    pr[UU] = u;
    pr[U1] = ur;
    pr[U2] = uh;
    pr[U3] = up;

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



int init_nodisk(int *whichvel, int*whichcoord,int i, int j, FTYPE *pr)
{
  SFTYPE randfact;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],r,th;
  struct of_geom realgeom,geom;
  


  get_geometry(i, j, CENT, &realgeom); // true coordinate system
  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

  // choose something reasonable
  //  set_atmosphere(0,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system
  
  // small density will indicate to atmosphere to change it to atmosphere model in fixup.c
  pr[RHO]=1E-20;
  pr[UU]=1E-22;
  //  pr[U1]=-.01*r;
  pr[U1]=0.0;
  pr[U2]=0.0;
  pr[U3]=0.00;
  
  
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}





// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][NPR])
{
  int i,j;
  FTYPE X[NDIM],r,th;

#if(DISKTYPE!=NODISK)

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

#endif

  return(0);
}




#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define DIPOLEFIELD 3 // i.e. no disk field with dipole field (like pulsar model)
#define DISKDIPOLE 4


#define FIELDTYPE DIPOLEFIELD

// assumes normal field in pr
int init_vpot(int i, int j,FTYPE *A)
{
  SFTYPE rho_av, u_av, q;
  FTYPE X[NDIM],r,th;
  struct of_geom geom;
  
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE u_ref,ftemp,rho;
  SFTYPE sth, cth;
  FTYPE fieldhor;


  coord(i, j, CORN, X);
  bl_coord(X, &r, &th);


  *A=0;
  //  *A += 1E-5*r ;

  return(0);

}

int init_vpot2field(SFTYPE A[][N2+1],FTYPE pr[][N2M][NPR])
{
  extern int vpot2field(SFTYPE A[][N2+1],FTYPE p[][N2M][NPR]);

  return(vpot2field(A,pr));
}


// assume already normalized how like with vector potential
int normalize_field(FTYPE p[][N2M][NPR])
{

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





FTYPE nz_func(FTYPE R)
{
  return(
	 sqrt(
	      (3.*a*a - 4.*a*sqrt(R) + R*R)/
	      pow(R*(a + pow(R,1.5)),2)
	      )
	 ) ;


}
