
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

int init()
{
  int i = 0, j = 0, k = 0;
  SFTYPE randfact;
  FTYPE r, th;
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  SFTYPE Br0;
  FTYPE X[NDIM];
  struct of_geom geom,ksgeom;

  /* for disk interior */
  SFTYPE l, rin, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;

  /* for magnetic field */
  SFTYPE A[N1 + 1][N2 + 1];
  FTYPE bsq_ij;
  SFTYPE rho_av, rhomax, umax, beta,  bsq_max, norm, q, beta_act;
  SFTYPE rmax, lfish_calc(double rmax);
  SFTYPE rh;

  //  ranc(7);
  ranc(0); // no MPI method yet, so just pure randomization
  randfact = 4.e-2;
  defcoord = 0;
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;

  /* disk parameters (use fishbone.m to select new solutions) */
  /* 
     a = 0.999 ; l= lfish_calc(rmax); */

  /*
  a = 0.0;
  l = 5.;
  rin = 5.7;

  kappa = 1.e-3;
  beta = 1.e2;
  */  
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

  a = 1.0 ;
  rin = 6. ;
  rmax = 12. ;
  l = lfish_calc(rmax) ;
  kappa = 1.e-3 ;
  beta = 1.e2 ;
  

  /* 
     a = 0.5 ; l = 4.4 ; rin = 5.0 ; kappa = 1.e-3 ; beta = 4.e2 ; */

  /* 
     a = 0. ; l = 5. ; rin = 5.7 ; kappa = 1.e-3 ; beta = 4.e2 ; */

  /* some numerical parameters */
  failuremode = 0;		// clean start
  failed = 0;
  cour = 0.9;
  lim = MC;
  // initial dt
  dt = 1.0e-5;
  R0 = 0.5;
  rh=(1. + sqrt(1. - a * a));
  Rin = 0.98 * rh;
  Rout = 40.;
  hslope = 0.3;

  set_grid();

  // determine nature of inner radial edge (assumes myid==0 is always there)
  if(myid==0){
    coord(-2, 0, FACE1, X);
    bl_coord(X, &r, &th);
    trifprintf("rmin: %g\n", r);
    trifprintf("rmin/rh: %g\n", r / rh );
    trifprintf("rmin/rsing: %g\n", r / a);
    if(r/rh<=1.0){
      trifprintf("inner grid is inside horizon\n");
    }
    else{
      trifprintf("inner grid is outside horizon\n");
    }
    if(r/a<=1.0){
      dualfprintf(fail_file,"inner grid is inside singularity\n");
      return(1);
    }
  }
  


  t = 0.;
  /* output choices */
  tf = 4000.0;


  DTd = 50.;			/* dumping frequency, in units of M */
  DTener = 2.0;			/* logfile frequency, in units of M */
  DTi = 2.0;			/* image file frequ., in units of M */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;			/* restart file period in steps */

  rhomax = 0.;
  umax = 0.;
  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    blgset(i,j,&geom) ;
    coord(geom.i, geom.j, CENT, X);
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


    /* regions outside torus */
    // this region is already in Kerr Schild prime
    if (lnh < 0. || r < rin) {
      rho = RHOMIN * 1.E-14;
      u = UUMIN * 1.E-14;
      

      // bl-normal observer (4-vel components)
      get_geometry(i, j, CENT, &geom);
      //ur = -geom.gcon[0][1]/sqrt(-geom.gcon[0][0]) ;
      //uh = -geom.gcon[0][2]/sqrt(-geom.gcon[0][0]) ;
      //up = -geom.gcon[0][3]/sqrt(-geom.gcon[0][0]) ;
      ur = geom.gcon[0][1]/geom.gcon[0][0] ;
      uh = geom.gcon[0][2]/geom.gcon[0][0] ;
      up = geom.gcon[0][3]/geom.gcon[0][0] ;

      p[i][j][RHO] = rho;
      p[i][j][UU] = u;
      p[i][j][U1] = ur;
      p[i][j][U2] = uh;
      p[i][j][U3] = up;
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


      p[i][j][RHO] = rho ;
      if (rho > rhomax)
	rhomax = rho;
      p[i][j][UU] = u* (1. + randfact * (ranc(0) - 0.5));
      if (u > umax && r > rin)
	umax = u;
      p[i][j][U1] = ur ;
      p[i][j][U2] = uh ;

      p[i][j][U3] = SLOWFAC * up;

      /* transform four-velocity to kerr-schild */
      if (bl2met2metp2v(0,p[i][j], i,j) >= 1)
	FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);

    }




    p[i][j][B1] = 0.;
    p[i][j][B2] = 0.;
    p[i][j][B3] = 0.;
  }

  mpimax(&rhomax);
  mpimax(&umax);
  trifprintf("rhomax: %g umax: %g\n", rhomax, umax);


  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    p[i][j][RHO] /= rhomax;
    p[i][j][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;
  if(fixup(p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
  if (bound_prim(p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_prim()", 1);
  if(postbc_fixup(p,0)>=1)
    FAILSTATEMENT("init.c:init()", "postbc_fixup()", 1);

  /* first find corner-centered vector potential */
  ZSLOOP(0, N1, 0, N2) A[i][j] = 0.;
  ZSLOOP(0, N1, 0, N2) {
    /* vertical field version */
    /* 
       coord(i,j,CORN,X) ; bl_coord(X,&r,&th) ;

       A[i][j] = 0.5*r*sin(th) ; */

    /* field-in-disk version */
    rho_av = 0.25 * (p[i][j][RHO] +
		     p[i - 1][j][RHO] +
		     p[i][j - 1][RHO] + p[i - 1][j - 1][RHO]);

    q = rho_av / rhomax - 0.2;
    if (q > 0.)      A[i][j] = q;
  }

  /* now differentiate to find cell-centered B, and begin normalization 
   */
  bsq_max = 0.;
  ZLOOP {
    get_geometry(i, j, CENT, &geom);

    /* flux-ct */
    p[i][j][B1] = (A[i][j] - A[i][j + 1]
		   + A[i + 1][j] - A[i + 1][j +
					    1]) / (2. * dx[2] * geom.g);
    p[i][j][B2] = -(A[i][j] + A[i][j + 1]
		    - A[i + 1][j] - A[i + 1][j +
					     1]) / (2. * dx[1] *
						    geom.g);

    p[i][j][B3] = 0.;

    if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
      FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

    if (bsq_ij > bsq_max)
      bsq_max = bsq_ij;
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %g\n", bsq_max);

  /* finally, normalize to set field strength */
  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
  trifprintf("initial beta: %g (should be %g)\n", beta_act,beta);
  norm = sqrt(beta_act / beta);

  bsq_max = 0.;
  ZLOOP {
    p[i][j][B1] *= norm;
    p[i][j][B2] *= norm;

    get_geometry(i, j, CENT, &geom);
    if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
      FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
    if (bsq_ij > bsq_max)
      bsq_max = bsq_ij;

  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %g\n", bsq_max);

  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

  trifprintf("new bsq_max: %g\n", bsq_max);
  trifprintf("final beta: %g (should be %g)\n", beta_act,beta);

  /* this is a little test calculation with a radial field, designed to 
     make the scheme fail */
  /* 
     Br0 = 1.0 ; ZLOOP { GSET(i,j,CENT) p[i][j][B1] = Br0/(rcurr*rcurr) 
     ; p[i][j][B2] = 0. ; p[i][j][B3] = 0. ; } */

  /* enforce boundary conditions */
  if(fixup(p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 2);
  if (bound_prim(p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_prim()", 2);
  if(postbc_fixup(p,0)>=1)
    FAILSTATEMENT("init.c:init()", "postbc_fixup()", 2);


  trifprintf("end init.c\n");
  return (0);

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
