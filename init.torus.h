//
//
// Description: 
//
// Common initialization routines that can be shared between multiple init files
//
// Author: Alexander Tchekhovskoy <atchekho@cfa.harvard.edu>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

static FTYPE lfunc( FTYPE lin, FTYPE *parms );
static void compute_gu( FTYPE r, FTYPE th, FTYPE a, FTYPE *gutt, FTYPE *gutp, FTYPE *gupp );
static FTYPE compute_udt( FTYPE r, FTYPE th, FTYPE a, FTYPE l );
static FTYPE compute_omega( FTYPE r, FTYPE th, FTYPE a, FTYPE l );
static FTYPE compute_l_from_omega( FTYPE r, FTYPE th, FTYPE a, FTYPE omega );
static FTYPE thintorus_findl( FTYPE r, FTYPE th, FTYPE a, FTYPE c, FTYPE al );
static FTYPE nz_func(FTYPE R) ;


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

  // AKMARK: entropy constant
  //BOBMARK: should be the value of KK in mathematica file
  kappa = toruskappa; //now set at top of init.xxxtorus.c


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
            //     p[i][j][k][RHO] = tabr;
            //     p[i][j][k][U3] = tabo;
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
            //     p[i][j][k][RHO] = tabr;
            //     p[i][j][k][U3] = tabo;
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






  kappa = toruskappa ;
  rmax = torusrmax ;   // AKMARK: torus pressure max
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
  FTYPE r, th, R, rho, u;
  FTYPE X[NDIM], V[NDIM];
  int pl;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;

  
  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  
  /* region outside disk? */
  R = r*sin(th) ;

  //get atmospheric values
  get_geometry(i, j, k, CENT, ptrgeom); // true coordinate system
  set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system
  
  //pull out analytic density
  rho = MACP0A1(GLOBALPOINT(panalytic),i,j,k,RHO);
  u = MACP0A1(GLOBALPOINT(panalytic),i,j,k,UU);

  //avoid computations way outside the torus: fill the region with atmosphere
  if( R < rin || rho < pr[RHO] ) {
    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
    
    
  PALLLOOP(pl) {
    if( pl == UU && u < pr[UU] ) {
      //leave u at floor value if u wants to be lower than floor value
      continue;
    }
    pr[pl] = MACP0A1(GLOBALPOINT(panalytic),i,j,k,pl);
  }


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
  FTYPE gutt, gutp, gupp, al, c;
  FTYPE ans;
   
  gutt = parms[0];
  gutp = parms[1];
  gupp = parms[2];
  al = parms[3]; 
  c = parms[4];
   
  ans = (gutp - lin * gupp)/( gutt - lin * gutp) - c *pow(lin/c,al); // (lin/c) form avoids catastrophic cancellation due to al = 2/n - 1 >> 1 for 2-n << 1
   
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
  
FTYPE thintorus_findl( FTYPE r, FTYPE th, FTYPE a, FTYPE c, FTYPE al )
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
  parms[4] = c;
  
  //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3) 
  //demand accuracy 5x machine prec.
  //in non-rel limit l_K = sqrt(r), use 10x that as the upper limit:
  l = rtbis( &lfunc, parms, 1, 10*sqrt(r), 5.*DBL_EPSILON );
  
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
#if(TORUSHASBREAKS == 1)
  FTYPE rbreak1, rbreak2, lbreak1, lbreak2;
#endif
  
    
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
 
  kappa = toruskappa;   // AKMARK: entropy constant KK from mathematica file
  n = torusn;   // AKMARK: n from mathematica file (power of lambda in DHK03)
  rmax = torusrmax;   // AKMARK: torus pressure max
#if(TORUSHASBREAKS == 1)   //AKMARK: midplane radii at which break in angular momentum profile occurs
  rbreak1 = 25.;
  rbreak2 = 75.;
#endif
  
  
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
  lin = thintorus_findl( r, th, a, c, al );
#if(TORUSHASBREAKS == 1)
  lbreak1 = thintorus_findl( rbreak1, th, a, c, al);
  lbreak2 = thintorus_findl( rbreak2, th, a, c, al);
#endif
  
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
  l = thintorus_findl( r, th, a, c, al );
#if(TORUSHASBREAKS == 1)
  if (l < lbreak1) l = lbreak1;
  if (l > lbreak2) l = lbreak2;
#endif
  
  
  udt = compute_udt( r, th, a, l );
  f = pow(fabs(1 - k*pow(l,1 + al)),pow(1 + al,-1));
  hh = flin*utin*pow(f,-1)*pow(udt,-1);
  eps = (-1 + hh)*pow(gam,-1);
  rho = pow((-1 + gam)*eps*pow(kappa,-1),pow(-1 + gam,-1));
  
  //compute atmospheric values
  get_geometry(ti, tj, tk, CENT, ptrgeom); // true coordinate system
  set_atmosphere(-1,WHICHVEL,ptrgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system
  
  //avoid computations outside the torus (either imaginary or negative density): fill the region with atmosphere
  if( eps < 0 || rho < pr[RHO] ) {
    //then use atmospheric values 
    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
   
  u = kappa * pow(rho, gam) / (gam - 1.);
  om = compute_omega( r, th, a, l );
  
  ur = 0.;
  uh = 0.;
  up = om; // = u^phi/u^t
 
  
  pr[RHO] = rho ;
  if( u > pr[UU] ) {
    pr[UU] = u* (1. + randfact * (ranc(0,0) - 0.5) );
  }

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

FTYPE nz_func(FTYPE R)
{
  return(
         sqrt(
              (3.*a*a - 4.*a*sqrt(R) + R*R)/
              pow(R*(a + pow(R,1.5)),2)
              )
         ) ;


}

// specify MPI task rank ordering
// example user-dependent code
// groups MPI tasks in n1xn2xn3 "tiles" (n# -- user adjustable)
// [taken from jet_set_myid()]
int group_by_node_set_myid(int n1tile, int n2tile, int n3tile)
{
  int ranki,rankj,rankk,origid,newid;
  int ranki0node, rankj0node, rankk0node, id0node;
  int n1 = n1tile, n2 = n2tile, n3 = n3tile;
  
  if(USEMPI==0){
    return(0); // nothing to do ever
  }

  // Can choose to rearrange MPI tasks

  // Reorder MPI tasks such that each node contains spatially-close cores in a 2x2x2 configuration
  
  //check if ncpux# are each even -- otherwise cannot decompose domain into 2x2x2 configurations
  if( 0 != (ncpux1 % n1) || 0 != (ncpux2 % n2) || 0 != (ncpux3 % n3) ) {
    dualfprintf( fail_file, "queenbee_set_myid(): no changes to MPI mapping: MPI domain size must be multiples of: ncpux1 = %d (%d), ncpux2 = %d (%d), ncpux3 = %d (%d)\n", ncpux1, n1, ncpux2, n2, ncpux3, n3 );
    return(0);
  }

  for(rankk=0;rankk<ncpux3;rankk++){
    for(rankj=0;rankj<ncpux2;rankj++){
      for(ranki=0;ranki<ncpux1;ranki++){
        origid=ranki + rankj*ncpux1 + rankk*ncpux1*ncpux2;
        ranki0node = ranki/n1;
        rankj0node = rankj/n2;
        rankk0node = rankk/n3;
        id0node=ranki0node + rankj0node*(ncpux1/n1) + rankk0node*(ncpux1/n1)*(ncpux2/n2);
        newid=(ranki%n1) + (rankj%n2)*n1 + (rankk%n3)*n1*n2;  //id on the node tile (size n1xn2xn3 cores)
        MPIid[origid]=id0node*n1*n2*n3+newid;
      }
    }
  }// end over rankk


  // Note that unlike TACC Ranger, TACC Lonestar has 2 sockets with 2 cores per socket, but sockets *share* main memory.  So no special socket association is required for optimal memory use.
  // http://services.tacc.utexas.edu/index.php/lonestar-user-guide

 

  return(0);
}

