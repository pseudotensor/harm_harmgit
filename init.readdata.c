

/*! \file init.readdata.c
    \brief GRB/KAZFULLEOS related reading of data from Fortran code's generated data
    
*/


/// initialize stellar structure
int init_star(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
  int set_phi_velocity_grb2(FTYPE *V, FTYPE *prstellar, FTYPE *pr);
  FTYPE X[NDIM],V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE r,th;
  void set_stellar_solution(int ii, int jj, int kk,FTYPE *pr, FTYPE *hcmsingle, FTYPE *ynu0single, FTYPE *ynusingle);
  FTYPE prstellar[NPR];
  FTYPE przamobl[NPR];
  FTYPE przamo[NPR];
  FTYPE pratm[NPR];
  struct of_geom geom;
  struct of_geom *ptrgeom;
  struct of_geom geombl;
  int pl,pliter;
  FTYPE hcmsingle,ynu0single,ynusingle;
  FTYPE parlist[MAXPARLIST];
  int numparms;








  //////////////////////////////////
  //
  // set free parameters
  //
  beta=1E6;
  Rbeta=1000.0*2.0*G*Mcgs/(C*C)/Lunit; // 1000R_S where Mcgs is the mass

  // DEBUG
  beta=1E30; // DEBUG
  // DEBUG


  //////////////////////////////////
  //
  // Interpolate read-in data to computational grid
  set_stellar_solution(i,j,k,prstellar,&hcmsingle,&ynu0single,&ynusingle);
  // prstellar has 3-velocity in prstellar[U1], need to convert


  /////////////////////////////
  //
  // Go ahead and set quantities that need no conversion of any kind
  //
  /////////////////////////////
  if(!isfinite(hcmsingle)){
    dualfprintf(fail_file,"read-in or interpolated bad hcmsingle: %d %d %d\n",i,j,k);
  }

  if(!isfinite(ynu0single)){
    dualfprintf(fail_file,"read-in or interpolated bad ynu0single: %d %d %d\n",i,j,k);
  }

  if(!isfinite(ynusingle)){
    dualfprintf(fail_file,"read-in or interpolated bad ynusingle: %d %d %d\n",i,j,k);
  }

#if(DOYL!=DONOYL)
  pr[YL] = prstellar[YL];
  if(!isfinite(pr[YL])){
    dualfprintf(fail_file,"read-in or interpolated bad YL: %d %d %d\n",i,j,k);
  }
#endif
#if(DOYNU!=DONOYNU)
  // already accounted for WHICHEVOLVEYNU
  pr[YNU] = prstellar[YNU];
  if(!isfinite(pr[YNU])){
    dualfprintf(fail_file,"read-in or interpolated bad YNU: %d %d %d\n",i,j,k);
  }
#endif

  //  dualfprintf(fail_file,"YLYNU: i=%d yl=%21.15g ynu=%21.15g\n",i,pr[YL],pr[YNU]);


  ////////////////////////////////////
  //
  // Setup EOSextra for kazfulleos.c
  // parlist starts its index at 0 with EOSextra["TDYNORYEGLOBAL"]
  //
  // why not just call to compute EOSglobal things? (only because of H)
  // only matters for Kaz EOS and should be in correct order and type of quantity
  parlist[TDYNORYEGLOBAL-FIRSTEOSGLOBAL]=pr[YE]; // now using YE as primitive and YL as conserved   ///pr[YL]-ynusingle; // Y_e for any WHICHEVOLVEYNU
  parlist[YNU0OLDGLOBAL-FIRSTEOSGLOBAL]=parlist[YNU0GLOBAL-FIRSTEOSGLOBAL]=ynu0single; // Y^0_\nu
  parlist[YNUOLDGLOBAL-FIRSTEOSGLOBAL]=ynusingle; // for WHICHEVOLVEYNU, no older yet, so indicate that by using same value


  int hi;
  for(hi=0;hi<NUMHDIRECTIONS;hi++){
    parlist[HGLOBAL-FIRSTEOSGLOBAL+hi]=hcmsingle; // H
  }
  // first guess is neutrinos have U=P=S=0
  parlist[UNUGLOBAL-FIRSTEOSGLOBAL]=0.0;
  parlist[PNUGLOBAL-FIRSTEOSGLOBAL]=0.0;
  parlist[SNUGLOBAL-FIRSTEOSGLOBAL]=0.0;

  parlist[IGLOBAL-FIRSTEOSGLOBAL]=i;
  parlist[JGLOBAL-FIRSTEOSGLOBAL]=j;
  parlist[KGLOBAL-FIRSTEOSGLOBAL]=k;

  store_EOS_parms(WHICHEOS,NUMEOSGLOBALS,GLOBALMAC(EOSextraglobal,i,j,k),parlist);



  //////////////////////////////////
  //
  // Set other aspects of stellar model, such as rotational velocity
  //
  ///////////////////////////////////
  ptrgeom=&geom;
  get_geometry(i,j,k,CENT,ptrgeom);
  coord_ijk(i, j, k, CENT, X);
  bl_coord_ijk(i, j, k, CENT,V);
  dxdxprim_ijk(i, j, k, CENT,dxdxp);
  r=V[1];
  th=V[2];



  //  if(i==1 && j==0 && k==0) dualfprintf(fail_file,"BANG1\n");
  // PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d j=%d k=%d prstellar[%d]=%21.15g\n",i,j,k,pl,prstellar[pl]);




  ////////////////////////////
  // assume in BL-coords
  // also adds up stellar 3-velocity to my additional atmosphere 3-velocity
  PLOOP(pliter,pl) pratm[pl]=prstellar[pl]; // default value (this copies densities and fields)
  set_phi_velocity_grb2(V,prstellar,pratm);


  //  dualfprintf(fail_file,"prstellar[UU]=%21.15g pratm[UU]=%21.15g\n",prstellar[UU],pratm[UU]);


  //  dualfprintf(fail_file,"%d %d :: star_rho=%g star_u=%g star_v1=%g star_v3=%g\n",i,j,pratm[RHO],pratm[UU],pratm[U1],pratm[U3]);




  ////////////////////////////////////////////////////////////
  //
  // Set primitives for different radial regions
  //
  ////////////////////////////////////////////////////////////


  // MBH is in length units
  if(fabs(r)<4.0*MBH){
    // fabs(r) is because r can be less than 0 and if far from BH in negative sense, then don't want to still use this
    // chose 4MBH since then smoothly matches onto freely-falling frame from stellar model
    // if this close to BH, then set velocity in PRIMECOORDS since can't use BL coords inside horizon and dubiously set inside ergosphere
    // even if using VEL4 we can't use BL-coords inside horizon
    PLOOP(pliter,pl) przamo[pl]=0.0;
    set_zamo_velocity(WHICHVEL,ptrgeom, przamo); // in PRIMECOORDS/WHICHVEL

    pr[RHO]=pratm[RHO];
    pr[UU]=pratm[UU];
    for(pl=U1;pl<=U3;pl++) pr[pl]=przamo[pl];
    for(pl=B1;pl<=B3;pl++) pr[pl]=pratm[pl]; // although really need vector potential+B3 or just B3 for 2D

  }
  else{

    // GODMARK: at the moment the initial-value problem is setup only with 1 step iteration (essentially ignore star initially)
    // GODMARK: velocity can be quite bad if arbitrarily chosen and put in black hole that requires certain velocity near it.
    // so add in ZAMO observer in BL-coords here instead of afterwards

    // this gets BL-geometry in native BL spc coordinates (not PRIMECOORDS)
    gset(0,BLCOORDS,i,j,k,&geombl);
  
    przamobl[U1] = (geombl.gcon[GIND(0,1)])/(geombl.gcon[GIND(0,0)]) ;
    przamobl[U2] = (geombl.gcon[GIND(0,2)])/(geombl.gcon[GIND(0,0)]) ;
    przamobl[U3] = (geombl.gcon[GIND(0,3)])/(geombl.gcon[GIND(0,0)]) ;

    // can avoid doing below if using VEL4
    //    if(pratm[U1]<-0.1) pratm[U1]=-0.1; // near horizon time slows down so that 3-velocity actually goes to 0
    for(pl=U1;pl<=U3;pl++)  pratm[pl] += przamobl[pl];

    //    if(i==1 && j==0 && k==0) dualfprintf(fail_file,"BANG\n");
    //PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d j=%d k=%d prstellar[%d]=%21.15g pratm[%d]=%21.15g przamobl[%d]=%21.15g\n",i,j,k,pl,prstellar[pl],pl,pratm[pl],pl,przamobl[pl]);


    // convert BL-coordinate velocity to PRIMECOORDS
    // now conversion should be safe since have at least ZAMO + some small modification due to the star
    // GODMARK: converting to KSCOORDS since don't at the moment have TOV solution since haven't yet setup fluid!
    //    *whichvel=VEL3;
    *whichvel=VEL4; // use 4-velocity so near BH velocity can be like in KS-coords and be -0.5 as in stellar model
    *whichcoord=BLCOORDS;
    if (bl2met2metp2v_gen(*whichvel,*whichcoord, WHICHVEL, KSCOORDS, pratm, i,j,k) >= 1)
      FAILSTATEMENT("init.readdata.c:get_data()", "bl2ks2ksp2v()", 1);

    
    ///////////////////////
    //
    // add up velocities in PRIMECOORDS
    // zamo is used in case near black hole so solution still good, where other term is not expected to account for black hole
    // zamo will be small if black hole starts off with small mass compared to self-gravity mass
    pr[RHO]=pratm[RHO];
    pr[UU]=pratm[UU];
    //for(pl=U1;pl<=U3;pl++) pr[pl]=przamo[pl]+pratm[pl]; // already accounted for ZAMO term in BL-coords above
    for(pl=U1;pl<=U3;pl++) pr[pl]=pratm[pl];

    // These field components are overwritten by vector potential solution
    // If want B3, should provide A_r or A_\theta in init.c
    for(pl=B1;pl<=B3;pl++) pr[pl]=pratm[pl]; // although really need vector potential+B3 or just B3 for 2D
    
    // DEBUG
    pr[U3]=0.0;// DEBUG
    // DEBUG
  }






  //  if(i==256 && j==0 && k==0) dualfprintf(fail_file,"BANG\n");
  //  PLOOP(pliter,pl) dualfprintf(fail_file,"i=%d j=%d k=%d prstellar[%d]=%21.15g przamo=%21.15g pratmpost=%21.15g\n",i,j,k,pl,prstellar[pl],przamo[pl],pratm[pl]);


  //  dualfprintf(fail_file,"prstellar[UU]=%21.15g pratm[UU]=%21.15g pr[UU]=%21.15g\n",prstellar[UU],pratm[UU],pr[UU]);



  // assume same for now GODMARK (only non-field set so far anyways)
  PLOOP(pliter,pl){
    pstag[pl]=pr[pl];
  }


  //////////////////////////////////
  //
  // Choose conversion of velocity
  //

  // assume already converted everything to PRIMECOORDS/WHICHVEL

  // use if setting in PRIMECOORDS
  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}



// set in BL-coords
int set_phi_velocity_grb2(FTYPE *V, FTYPE *prstellar, FTYPE *pr)
{
  FTYPE roR0;
  FTYPE R0;
  FTYPE r,th;
  FTYPE omega,omegak,llimit,lspecific,lspeclimit;
  FTYPE Rhecgs,Rhe;
  FTYPE rotcoef;

  r=V[1];
  th=V[2];



  Rhecgs=2.11E8; // in cm (may be wrong, see grbmodel.m)
  Rhe=Rhecgs/Lunit;


#if( (ROTATIONPROFILE==ROT_FUJFAST)||(ROTATIONPROFILE==ROT_FUJSLOW))


#if(ROTATIONPROFILE==ROT_FUJFAST)
  // 10 radians per second
  omega0 = 10.0*Tunit;
  // 1000km
  R0 = 1000.0*1E3*1E2/Lunit;
#elif(ROTATIONPROFILE==ROT_FUJSLOW)
  // 0.5 radians per second
  omega0 = 0.5*Tunit;
  // 5000km
  R0 = 5000.0*1E3*1E2/Lunit;  
#endif

  roR0=r/R0;
  pr[U1] = prstellar[U1];
  pr[U2] = prstellar[U2];
  pr[U3] = omega0/(1+roR0);

#elif(ROTATIONPROFILE==ROT_MW99) // Similar to Nakataki (2006)
  // Similar to ROT_PROGA, except here rotation is specified for all angles and radii with no limit to iron core

  rotcoef=0.05;
  lspeclimit=1E17; // in cm^2/s

  if(r<2.0*MBH){ // recall MBH is in length (Lunit) units
    // else not rotation anymore than from ZAMO frame that is added later
    pr[U1]=prstellar[U1];
    pr[U2]=prstellar[U2];
    pr[U3]=prstellar[U3];
  }
  else if(r<Rhe){ // Inner radius of Helium envelope
    pr[U1]=prstellar[U1];
    pr[U2]=prstellar[U2];
    // can't be Keplerian with Mcgs close to GMcgs/C^2, so make iron core rotate uniformly

    omegak=sqrt(G*Mcgs/pow(Rhecgs,3.0))*Tunit; // dimensionless now
    omega=sqrt(rotcoef)*omegak; // rotcoef = centrifugal/(R-cylindrical gravity) 
    lspecific=pow(Rhe*sin(th),2.0)*omega; // u_\phi 
    llimit=lspeclimit/(Lunit*Lunit)*Tunit; // upper limit on specific angular momentum
    if(lspecific>llimit) lspecific=llimit;
    pr[U3] = lspecific/pow(Rhe*sin(th),2.0); // non-relativistic - diagonal approximation for spherical polar coordinates
  }
  else{
    
    // here they use radial velocity from stellar model
    pr[U1]=prstellar[U1];
    pr[U2]=prstellar[U2];
    // Fcentrifugal/Fgrav(perp rot) = rotcoef
    
    omegak=sqrt(G*Mcgs/pow(r*Lunit,3.0))*Tunit; // dimensionless now
    omega=sqrt(rotcoef)*omegak; // rotcoef = centrifugal/(R-cylindrical gravity) 
    lspecific=pow(r*sin(th),2.0)*omega; // u_\phi 
    llimit=lspeclimit/(Lunit*Lunit)*Tunit; // upper limit on specific angular momentum
    if(lspecific>llimit) lspecific=llimit;
    pr[U3] = lspecific/pow(r*sin(th),2.0); // non-relativistic - diagonal approximation for spherical polar coordinates
    
  }

#elif(ROTATIONPROFILE==ROT_PROGA)
    
  rotcoef=0.02;
  lspeclimit=1E17; // in cm^2/s

  if(r<2.0*MBH){ // recall MBH is in length (Lunit) units
    // else not rotation anymore than from ZAMO frame that is added later
    pr[U1]=prstellar[U1];
    pr[U2]=prstellar[U2];
    pr[U3]=prstellar[U3];
  }
  else if(r<Rhe/Lunit){ // Inner radius of Helium envelope
    // this actually only includes the Iron core
    // non-relativistic free-fall onto a gravitational object of mass Mcgs
    //pr[U1]=sqrt(2.0*G*Mcgs/(r*Lunit))/Vunit; // overrides stellar model
    pr[U1]=prstellar[U1]; // stellar model now has free-fall velocity, and set using enclosed mass
    pr[U2]=prstellar[U2];
    pr[U3]=prstellar[U3];
  }
  else{ // outside R_{He}, the inner radius of the helium envelope (beyond iron core)

    // here they use radial velocity from stellar model
    pr[U1]=prstellar[U1];
    pr[U2]=prstellar[U2];
    // Fcentrifugal/Fgrav(perp rot) = rotcoef

    omegak=sqrt(G*Mcgs/pow(r*Lunit,3.0))*Tunit; // dimensionless now

    omega=sqrt(rotcoef)*omegak; // rotcoef = centrifugal/(R-cylindrical gravity) 
    lspecific=pow(r*sin(th),2.0)*omega; // u_\phi 
    llimit=lspeclimit/(Lunit*Lunit)*Tunit; // upper limit on specific angular momentum
    if(lspecific>llimit) lspecific=llimit;
    pr[U3] = lspecific/pow(r*sin(th),2.0); // non-relativistic - diagonal approximation for spherical polar coordinates
  }
#endif


  return(0);


}





// stellar models are read in from file created by SM macro "grbmodel.m"

// GODMARK:
// if doing evolving metric, then don't need "truncaterho" in grbmodel.m when doing rdstar
// if doing self-gravity, then don't need "setfreefall vr" in grbmodel.m when doing rdstar
//
void get_stellar_data(void)
{
  int i;
  FILE * indata;
  FTYPE ftemp1;
  FTYPE ye;



  if(myid==0){

    // report read-in
    trifprintf("Begin reading in stellar model\n");
    
    
    
    ///////////////////
    //
    // open file
    //
#if(STELLARTYPE==0)
   
    // for all models the "Helium core" has been kept and the Hydrogen envelope replaced with a stellar wind
#if(DOSELFGRAVVSR==0 && DOEVOLVEMETRIC==1) 
    // star is given free-fall velocity so don't have to follow self-gravity (although still apparently needed for outer region)
#define FILENAME "stellar0.txt"
#elif(DOSELFGRAVVSR==1 && DOEVOLVEMETRIC==1) 
    // vr and rho are left as is so can follow collapse
#define FILENAME "stellar1.txt"
#elif(DOSELFGRAVVSR==0 && DOEVOLVEMETRIC==0) 
    // star is given free-fall velocity and mass is lost from iron core
#define FILENAME "stellar2.txt"
#endif


    if( (indata=fopen(FILENAME,"rt"))==NULL){
      dualfprintf(fail_file,"Cannot open %s\n",FILENAME);
      myexit(100);
    }
    
    
#elif(STELLARTYPE==1)
    dualfprintf(fail_file,"No such STELLARTYPE==%d\n",STELLARTYPE);
    myexit(266);
#endif 
    
    
    
    ///////////////////////////////////////////
    //
    // READ IN DATA FILE
    //
    ////////////////////////////////////////////
    
    // skip first blank line
    //while(fgetc(indata)!='\n'); // skip first line, a comment
    
    // read in data file
    
    // read in first row (parameters)
    // NOTE: a0 not right since no rotation in stellar model.  Need to determine that later
    fscanf(indata,"%d %lf %lf %lf %lf %lf",&NRADIAL,&MBH0,&a0,&QBH0,&EP30,&THETAROT0);
    
    if(NRADIAL>=NRADIALMAX){
      dualfprintf(fail_file,"NRADIAL=%d and NRADIALMAX=%d\n",NRADIAL,NRADIALMAX);
      myexit(7);
    }
    
    
    // convert MBH0, a0, and QBH0
    
    MBH0 = (G*MBH0/(C*C))/Lunit; // mass in length units
    // no conversion yet for a0 and QBH0 since not set right in stellar model


    // from SM macro "outputmodel" in grbmodel.m
    //  print $filename '%d %21.15g %21.15g %21.15g\n' {numlines MBH0 aBH0 QBH0}
    //  #
    //  #
    //  print + $filename '%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n' \
    //           {r rho utot vr temp ye Ynu Ynu0 hcm nucleons helium carbon oxygen neon magnesium si iron j omega inertia}

    
    // report parameters
    trifprintf("Stellar parameters: %d :: %g %g %g %g %g\n",NRADIAL,MBH0,a0,QBH0,EP30,THETAROT0);
    
    // 18 values
    // read in 1-D data
    for(i=0;i<NRADIAL;i++){
      fscanf(indata,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(radius[i]),&(rho[i]),&(ie[i]),&(vr[i]),&ftemp1,&ye,&(ynu[i]),&(ynu0[i]),&(hcm[i]),&ftemp1,&ftemp1,&ftemp1,&ftemp1,&ftemp1,&ftemp1,&ftemp1,&ftemp1,&ftemp1,&(omega3[i]),&ftemp1);

      
      // some assignments
      // generally true for any WHICHEVOLVEYNU
      yl[i] = ye + ynu[i];


      //////////////////////////////////
      //
      // Convert units for rho,u,orhonormal v^r

      //      dualfprintf(fail_file,"i=%d ie=%21.15g Pressureunit=%21.15g\n",i,ie[i],Pressureunit);


      radius[i]/=Lunit;
      if(rho0unittype==0) rho[i]/=rhounit;
      else rho[i]/=rhomassunit; // read-in g/cc, while rhounit is erg/cc (code will intreprety rho as in erg/cc just like pressure)
      ie[i]/=Pressureunit;
      vr[i]/=Vunit;
      omega3[i]/=(1.0/Tunit);
      hcm[i]/=(Lunit);
      // ye and ynu (and so yl) are dimensionless

      
    }


    // report read-in
    trifprintf("Done reading in stellar model\n");
    
    
  } // end if myid==0



  // send data to all CPUs
  // so all CPUs have full 1-D data
#if(USEMPI)
  MPI_Bcast(&NRADIAL,1,MPI_INT,MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&MBH0,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&a0,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&QBH0,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&EP30,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&THETAROT0,1,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&(radius[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(rho[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(ie[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(vr[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(omega3[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(yl[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(ynu[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(ynu0[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&(hcm[0]),NRADIAL,MPI_FTYPE,MPIid[0], MPI_COMM_GRMHD);
#endif


}






////////////////////////////////////////////////////
//
// Interpolate radial stellar data onto grid
//
///////////////////////////////////////////////////
// OPENMPMARK: Assume set_stellar_solution() not called by multiple threads, so ok to  have static and firsttime.
void set_stellar_solution(int ii, int jj, int kk,FTYPE *pr, FTYPE *hcmsingle, FTYPE *ynu0single, FTYPE *ynusingle)
{
  FTYPE X[NDIM],V[NDIM],r;
  struct of_geom geom;
  FTYPE x2;
  int i;
  static int firsttime=1;
  FTYPE myrho,myie,myvr,myomega3,myyl,myynu;


  ///////////////////////////////////////////
  //
  // Already have data, now interpolate data to HARM grid into "pr"
  //
  ////////////////////////////////////////////


  coord_ijk(ii, jj, kk, CENT, X);
  bl_coord_ijk(ii, jj, kk, CENT, V);
  //get_geometry(ii,jj,kk,CENT,&geom);
  r=V[1];


  // radius is not uniform, so loop to find first radius
  for(i=0;i<NRADIAL;i++){
    if(radius[i]>r) break;
  }

  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, rho, &myrho);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, ie, &myie);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, vr, &myvr);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, omega3, &myomega3);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, yl, &myyl);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, ynu, &myynu);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, hcm, hcmsingle);
  interpfun(QUADRATICTYPE,NRADIAL,i, r, radius, ynu0, ynu0single);
    
  // units already converted
  pr[RHO]=myrho;
  pr[UU]=myie;
  pr[U1]=myvr; // in BL-coords treated as 4-velocity...just a place holder for 3-velocity, will convert-add to later
  pr[U2]=0.0;
  pr[U3]=myomega3; // treated as with pr[U1]
  pr[B1]=0.0;
  pr[B2]=0.0;
  pr[B3]=0.0;
#if(DOYL!=DONOYL)
  //  pr[YL]=myyl;
  // now using YE as primitive and YL as conserved
  pr[YE]=myyl-myynu;
#endif
#if(DOYNU!=DONOYNU)
  *ynusingle=myynu;

  if(WHICHEVOLVEYNU==EVOLVEYNUNOTRAD){
    // evolving Ynu0, which should be primitive
    pr[YNU]=*ynu0single;
  }
  else{
    // evolving Ynu, which should be primitive
    pr[YNU]=myynu;
  }
#endif
}








