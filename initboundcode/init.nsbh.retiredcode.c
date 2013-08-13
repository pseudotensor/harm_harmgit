 

// NOT DONE -- see below
// based upon do_reduce_near_cusps() in reconstructeno.c
// detect cusp
#define MAXCUSP 10
#define SQRT_CUSP_EPSILON (500.0L*NUMEPSILON)
#define CUSP_EPSILON (SQRT_CUSP_EPSILON*SQRT_CUSP_EPSILON)
int checkcusp(int firstpos, int lastpos, FTYPE *xpos, FTYPE *yfun)
{
  int do_reduce = 0;
  FTYPE a_first_der_arr[4];
  FTYPE a_second_der_arr[3];
  FTYPE *first_der_arr;
  FTYPE *second_der_arr;
  FTYPE normu;
  FTYPE epsilon;
  int i;
  FTYPE a_u[MAXCUSP],a_pos[MAXCUSP];
  FTYPE *u,*pos;
  int min_index,max_index;


  // setup stored x and y
  min_index=firstpos-2;
  max_index=lastpos-2;
  pos=&a_pos[MAXCUSP/2];
  u=&a_u[MAXCUSP/2];
  first_der_arr=&a_first_der_arr[2];
  second_der_arr=&a_second_der_arr[1];

  for( i = -2; i <= 2; i++ ){
    pos[i] = xpos[firstpos+2+i];
    u[i] = yfun[firstpos+2+i];
  }


  normu=0.0;
  for( i = -2; i <= 2; i++ ){
    normu += fabs(u[i]);
  }

  for( i = -2; i <= 1; i++ ){
    first_der_arr[i]=(u[i+1]-u[i])/(pos[i+1]-pos[i]);
  }
  for( i = -1; i <= 1; i++ ){
    second_der_arr[i]=(first_der_arr[i]-u[i-1])/(0.5*(pos[i+1]+pos[i]) - 0.5*(pos[i]+pos[i-1]));
  }


  epsilon = normu * CUSP_EPSILON;


  // NOT DONE -- nothing below is correct/checked/etc.


  if( min_index <= -2 && max_index >= 2 ) {  //only proceed if there is enough data around
    second_der_arr[0] = (u[0] + u[2]) - 2 * u[1];
    second_der_arr[1] = (u[-1] + u[1]) - 2 * u[0];
    second_der_arr[2] = (u[0] + u[-2]) - 2 * u[-1];

    if( second_der_arr[2] * ( u[0] - u[-1] ) > -epsilon &&
        second_der_arr[1] * ( u[0] - u[-1] ) < epsilon)
      {
        if( ( u[2] - u[1] ) * ( u[0] - u[-1] ) < epsilon  ) {
          do_reduce |= 2;  //to 2nd order  //atch correct reduction
        }
        else if(  ( u[1] - u[0] ) * ( u[0] - u[-1] ) < epsilon ) {
          do_reduce |= 1;  //to 1st order  //atch correct reduction
        }
      }

    if( second_der_arr[0] * ( u[0] - u[1] ) > -epsilon &&
        second_der_arr[1] * ( u[0] - u[1] ) < epsilon
        )
      {
        if( ( u[-2] - u[-1] ) * ( u[0] - u[1] ) < epsilon ) {
          do_reduce |= 2;  //to 2nd order  //atch correct reduction
        }
        else if( ( u[-1] - u[0] ) * ( u[0] - u[1] ) < epsilon ) {
          do_reduce |= 1;  //to 1st order  //atch correct reduction
        }
      }
  }

  return( do_reduce );
}




// NOT DONE
// detect kink
#define BADKINK 1.0
int checkkink(int firstpos, int lastpos, FTYPE *xpos, FTYPE *yfun)
{
  int kink;
  FTYPE der[2],ders;

  // xpos[0],yfun[0] hold extrapolated x,y
  // xpos[firstpos],ypos[firstpos] holds next x,y
  // xpos[firstpos+1],ypos[firstpos+1] holds next x,y

  if(firstpos+1<lastpos){
    dualfprintf(fail_file,"Not enough points in checkink\n");
    myexit(189482345);
  }

  // get first derivatives
  der[0]=(yfun[firstpos]-yfun[0])/(xpos[firstpos]-xpos[0]);
  der[1]=(yfun[firstpos+1]-yfun[firstpos])/(xpos[firstpos+1]-xpos[firstpos]);
  ders=(yfun[firstpos+1]-yfun[0])/(xpos[firstpos+1]-xpos[0]);

  if(sign(der[0])!=sign(der[1]) && (fabs(der[0])/(fabs(ders)+SMALL)>BADKINK || fabs(der[1])/(fabs(ders)+SMALL)>BADKINK ) ){
    kink=1;
  }
  else kink=0;


  return(kink);
}




// special NS boundary code for centerered primitives
// averages all points within NS, but low order
//int bound_prim_user_dir_nsbh(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
int bound_prim_user_dir_nsbh_old1(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;
  static int firstwhichdir;
  static int firsttime=1;
  FTYPE prcum[NPR];
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  struct of_geom geomdontuseii;
  struct of_geom *ptrgeomii=&geomdontuseii;
  int rescale_Bcon_to_Bcovphi(FTYPE *pr,struct of_geom *ptrrgeom, FTYPE *Bd3);
  int rescale_Bconp_and_Bcovphi_to_Bconphi(FTYPE *pr, struct of_geom *ptrgeom, FTYPE Bd3);





  //////////////////
  //
  // Only need to do this once, not per-dimension, so setup which dimension can always trigger on (assumes dimensionality of problem doesn't change)
  //
  //////////////////
  if(firsttime==1){
    firsttime=0;
    firstwhichdir=whichdir;
  }





  if(WHICHPROBLEM==NSBH){
    //////////////////
    //
    // check if NS has sufficient grid depth for boundary conditions to be used for interpolations on smoothish functions
    //
    // if NS is moving, then must do this every substep and before boundary call since boundary call corresponds to at new time.
    //
    //////////////////
    //    check_nsdepth(boundtime);    // GODMARK: TODO: In 3D will need this to be time-dependent, but now that code for getting distances is too slow.
  }







  //////////////////////
  //
  // Set (time-dependent) primitives except pstag and field centered
  // copied from init.tools.c from init_primitives user1 function
  //
  //////////////////////


  if(whichdir==firstwhichdir){ // only need to set NS boundary conditions once -- not per dimension

    // DEBUG:
    dualfprintf(fail_file,"BOUNDTYPE: inittypeglobal=%d\n",inittypeglobal);


    int whichvel,whichcoord,i,j,k;
    int initreturn;
    int inittype=0; // evolve type
    //  don't set inittypeglobal, set by init_primitives
    int ii,jj,kk;

#pragma omp parallel private(i,j,k,initreturn,whichvel,whichcoord) OPENMPGLOBALPRIVATEFULL
    {
      OPENMP3DLOOPVARSDEFINE;
      ////////  COMPFULLLOOP{
      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);


        if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1){
          // if inside NS, then use shell values to form average (distance weighted) answer for boundary values



          // temp things
          FTYPE weight,iweight;
          FTYPE totalweight;
          FTYPE maxdist;
          FTYPE Bd3;


          // get geometry for i,j,k
          get_geometry(i, j, k, CENT, ptrgeom);

   
          // initialize
          PLOOP(pliter,pl) prcum[pl]=0.0;
          totalweight=0.0;
          // avoid using completely opposite side of NS, where field values flip, so only go 2-3 cells beyond minimum distance to shell.
          // assume (int) on dist=sqrt() would potentially lead to 0.9->0, so minimum value is 1
          maxdist=GLOBALMACP0A1(nsmask,i,j,k,NSMASKDISTTOSHELL)+2.0;
   
   
          // LOOP over possible ii,jj,kk for prim's to add to i,j,k prim
          COMPLOOPNiijjkk(ii,jj,kk){

            int iim=ii-N1NOT1;
            int jjm=jj-N2NOT1;
            int kkm=kk-N3NOT1;

            int iip=ii+N1NOT1;
            int jjp=jj+N2NOT1;
            int kkp=kk+N3NOT1;

            if(
               (ii>=-N1BND && ii<=N1M-1 &&  jj>=-N2BND && jj<=N2M-1 &&  kk>=-N3BND && kk<=N3M-1) && // has to be on grid
               (GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKINSIDE)==0) && // can't itself be inside NS
               (
                (GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==1)|| // if a shell grid
                // grab a few more aligned with grid (ala 1D interpolation)
                (iim>=-N1BND && iim<=N1M-1 && GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==0 && GLOBALMACP0A1(nsmask,iim,jj,kk,NSMASKSHELL)==1)|| // not on shell, but shell next to this cell
                (jjm>=-N2BND && jjm<=N2M-1 && GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==0 && GLOBALMACP0A1(nsmask,ii,jjm,kk,NSMASKSHELL)==1)|| // not on shell, but shell next to this cell
                (kkm>=-N3BND && kkm<=N3M-1 && GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==0 && GLOBALMACP0A1(nsmask,ii,jj,kkm,NSMASKSHELL)==1)|| // not on shell, but shell next to this cell
                (iip>=-N1BND && iip<=N1M-1 && GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==0 && GLOBALMACP0A1(nsmask,iip,jj,kk,NSMASKSHELL)==1)|| // not on shell, but shell next to this cell
                (jjp>=-N2BND && jjp<=N2M-1 && GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==0 && GLOBALMACP0A1(nsmask,ii,jjp,kk,NSMASKSHELL)==1)|| // not on shell, but shell next to this cell
                (kkp>=-N3BND && kkp<=N3M-1 && GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==0 && GLOBALMACP0A1(nsmask,ii,jj,kkp,NSMASKSHELL)==1) // not on shell, but shell next to this cell
                )
               ){
              // then found shell value, add to the i,j,k value

              //////////////////
              //
              // first get \gamma to ensure not going to copy-in bad (i.e. limited due to bsq low and gamma\sim GAMMAMAX) values.
              //
              //////////////////
              get_geometry(ii, jj, kk, CENT, ptrgeomii);
              FTYPE Bcon[NDIM],ucon[NDIM],vcon[NDIM],others[NUMOTHERSTATERESULTS];
              ucon_calc(MAC(prim,ii,jj,kk),ptrgeomii,ucon,others);

              // only use point if not out of wack when doing FFDE (or lower GAMMAMAX if doing MHD)
              if(ucon[TT]<0.5*GAMMAMAX){

                // weight is grid distance (never will be zero since never inside and shell at same position)
                FTYPE dist;
                dist=sqrt((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k));


                iweight=((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k));
                weight=1.0/iweight;


                if(dist<=maxdist){

                  // add this ii,jj,kk value to i,j,k with weight
                  PLOOP(pliter,pl){
                    if(pl==B3){
                      // get B_\phi and use it as thing to interpolate
                      rescale_Bcon_to_Bcovphi(MAC(prim,ii,jj,kk),ptrgeomii, &Bd3);
                      prcum[pl] += Bd3*weight;
                    }
                    else{
                      prcum[pl] += MACP0A1(prim,ii,jj,kk,pl)*weight;
                    }
                  }
                  totalweight+=weight;
                }
              }

              // DEBUG:
              //       dualfprintf(fail_file,"WEIGHTCENT: %d %d %d : %d %d %d : %21.15g %21.15g : maxdist=%21.15g\n",i,j,k,ii,jj,kk,iweight,totalweight,maxdist);

            }// end if found grid shell
          }// end over ii,jj,kk

          // normalize
          if(totalweight>0.0){
            PLOOP(pliter,pl){
              if(pl==B3){
                // final Bd3
                Bd3=prcum[pl]/totalweight;
                // puts Bconphi into prim[B3]
                rescale_Bconp_and_Bcovphi_to_Bconphi(MAC(prim,i,j,k), ptrgeom, Bd3);
              }
              else{
                MACP0A1(prim,i,j,k,pl)=prcum[pl]/totalweight;
              }
            }


          }
          else{
            dualfprintf(fail_file,"Never found shell to use for: %d %d %d.  Means inside NS, but no nearby shell (can occur for MPI).  Assume if so, then value doesn't matter (not used), so revert to fixed initial value for NS.\n",i,j,k);
            dualfprintf(fail_file,"maxdist=%d\n",maxdist);

            initreturn=init_dsandvels(inittype, CENT, &whichvel, &whichcoord,i,j,k,boundtime,MAC(prim,i,j,k),NULL);

            if(initreturn>0){
              FAILSTATEMENT("init.c:init_primitives()", "init_dsandvels()", 1);
            }
            else if(initreturn==0) MYFUN(transform_primitive_vB(whichvel, whichcoord, i,j,k, prim, NULL),"init.c:init_primitives","transform_primitive_vB()",0);
          }

        }// end if inside NS

      }// end loop block
    }// end parallel region


  }

  return(0);

}













// special NS boundary code for centerered primitives
// only averages if near surface (so bad to higher order)
int bound_prim_user_dir_nsbh_old2(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int pl,pliter;
  static int firstwhichdir;
  static int firsttime=1;
  FTYPE prext[NPR];
  FTYPE *prextptr;

  //////////////////
  //
  // Only need to do this once, not per-dimension, so setup which dimension can always trigger on (assumes dimensionality of problem doesn't change)
  //
  //////////////////
  if(firsttime==1){
    firsttime=0;
    firstwhichdir=whichdir;
  }




  if(WHICHPROBLEM==NSBH){
    //////////////////
    //
    // check if NS has sufficient grid depth for boundary conditions to be used for interpolations on smoothish functions
    //
    // if NS is moving, then must do this every substep and before boundary call since boundary call corresponds to at new time.
    //
    //////////////////
    //check_nsdepth(boundtime);    
  }







  //////////////////////
  //
  // Set (time-dependent) primitives except pstag and field centered
  // copied from init.tools.c from init_primitives user1 function
  //
  //////////////////////


  if(whichdir==firstwhichdir){ // only need to set NS boundary conditions once -- not per dimension

    // DEBUG:
    dualfprintf(fail_file,"BOUNDTYPE: inittypeglobal=%d\n",inittypeglobal);


    int whichvel,whichcoord,i,j,k;
    int initreturn;
    int inittype=0; // evolve type
    //  don't set inittypeglobal, set by init_primitives

#pragma omp parallel private(i,j,k,initreturn,whichvel,whichcoord) OPENMPGLOBALPRIVATEFULL
    {
      OPENMP3DLOOPVARSDEFINE;
      ////////  COMPFULLLOOP{
      OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMP3DLOOPBLOCK{
        OPENMP3DLOOPBLOCK2IJK(i,j,k);

        // initialize
        PLOOP(pliter,pl) prext[pl]=0.0;
        FTYPE count=0.0;

        int dir;
        for(dir=1;dir<=3;dir++){

          // fill pstag with primitive that is nearest external from NS surface
          // see if FACEdir is NS boundary
          int im,jm,km;
          im=i-(dir==1)*N1NOT1;
          jm=j-(dir==2)*N2NOT1;
          km=k-(dir==3)*N3NOT1;

          int im2,jm2,km2;
          im2=i-2*(dir==1)*N1NOT1;
          jm2=j-2*(dir==2)*N2NOT1;
          km2=k-2*(dir==3)*N3NOT1;
   
          int ip,jp,kp;
          ip=i+(dir==1)*N1NOT1;
          jp=j+(dir==2)*N2NOT1;
          kp=k+(dir==3)*N3NOT1;

          int ip2,jp2,kp2;
          ip2=i+2*(dir==1)*N1NOT1;
          jp2=j+2*(dir==2)*N2NOT1;
          kp2=k+2*(dir==3)*N3NOT1;
   

          int faceleft,faceright,faceleft2,faceright2;

          // |    shell    |face_i    i inNS    |
          faceleft=(im>=-N1BND && jm>=-N2BND && km>=-N3BND && GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1    && GLOBALMACP0A1(nsmask,im,jm,km,NSMASKSHELL)==1);
          faceleft2=(im2>=-N1BND && jm2>=-N2BND && km2>=-N3BND && GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1    && GLOBALMACP0A1(nsmask,im2,jm2,km2,NSMASKSHELL)==1);

          // |   i inNS     |face_i   shell   |
          faceright =(ip<=N1M-1 && jp<=N2M-1 && kp<=N3M-1 && GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,ip,jp,kp,NSMASKSHELL)==1);
          faceright2 =(ip2<=N1M-1 && jp2<=N2M-1 && kp2<=N3M-1 && GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,ip2,jp2,kp2,NSMASKSHELL)==1);
   
   
          if(faceleft){
            PLOOP(pliter,pl) prext[pl] += MACP0A1(prim,im,jm,km,pl);
            count++;
          }
          if(faceleft2){
            PLOOP(pliter,pl) prext[pl] += MACP0A1(prim,im2,jm2,km2,pl);
            count++;
          }
          if(faceright){
            PLOOP(pliter,pl) prext[pl] += MACP0A1(prim,ip,jp,kp,pl);
            count++;
          }
          if(faceright2){
            PLOOP(pliter,pl) prext[pl] += MACP0A1(prim,ip2,jp2,kp2,pl);
            count++;
          }
        }// over dir

 
        // normalize
        if(count>=0.5){
          PLOOP(pliter,pl) prext[pl]/=count;
          prextptr=prext;
        }
        else{
          // default is same value as local
          PLOOP(pliter,pl) prext[pl] = MACP0A1(prim,i,j,k,pl);
          prextptr=NULL; // avoid special case, since otherwise reuses same values
        }

        // DEBUG:
        // dualfprintf(fail_file,"BCCOUNT: %d %d %d : %21.15g\n",i,j,k,count);


#if(0)
        // inputted (to function) prim[] is inputted here since want to change that CENT value of prim
        // prext is fed as pstag but really contains nearest primitive external to NS
        // request densities for all computational centers
        initreturn=init_dsandvels(inittype, CENT, &whichvel, &whichcoord,boundtime,i,j,k,MAC(prim,i,j,k),prextptr);
        // initreturn=init_dsandvels(inittype, CENT, &whichvel, &whichcoord,boundtime,i,j,k,MAC(prim,i,j,k),NULL);

        if(initreturn>0){
          FAILSTATEMENT("init.c:init_primitives()", "init_dsandvels()", 1);
        }
        else if(initreturn==0) MYFUN(transform_primitive_vB(whichvel, whichcoord, i,j,k, prim, NULL),"init.c:init_primitives","transform_primitive_vB()",0);

        // else initreturn<0 then no change (i.e. nothing to set for boundary condition)
#else
        // BCs are just proxy for on-grid behavior so that interpolation doesn't lead to jumps and so dissipative fluxes.  Since field is staggered, can even mess with centered field.
        // So, for CENT BCs's that are used for FACE fluxes (not on NS surface that's set directly) and for EMF velocities (except on NS surface that's set directly), this will be used for them.
 

#endif

      }
    }// end parallel region


  }

  return(0);

}


// For those A_i components that are constrained to be fixed in time, fix A_i to the init_vpot_user() version fo A_i within the NS and on the NS surface
// GODMARK TODO: Not right for 3D.
void adjust_fluxctstag_vpot_dosetfix_new(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int insideNS;
  int i,j,k;
  int l;
  int loc;
  FTYPE vpotlocal[NDIM];
  int ii,jj,kk;


  //////////////
  //
  // First set interior and surface of NS to analytical solution
  //
  // Different loops for each A_i
  //
  // don't include actual or MPI boundaries in count since they don't form real domain of NS that must be resolved
  //////////////


  // GODMARK TODO SUPERMARK: TO FIX FOR CURVATURE VERSION OF ON SURFACE (NO LONGER USING THIS FUNCTION, so can ignore)

  COMPFULLLOOP{ // First go over all non-surface + surface cells for simplicity


    if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1){

      // Set A_i in volume (lower CORNi of cell)
      // get_vpot_fluxctstag_primecoords gets A_i in PRIMECOORDS coords as required for vpot[]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
      if(Nvec[1]==1)      MACP1A0(vpot,1,i,j,k)       =vpotlocal[1];
      if(Nvec[2]==1)      MACP1A0(vpot,2,i,j,k)       =vpotlocal[2];
      if(Nvec[3]==1)      MACP1A0(vpot,3,i,j,k)       =vpotlocal[3];

      // Set A_i on surfaces completely (other CORNi's of cell)

      // A_x A_y [up k]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
      if(Nvec[1]==1)      MACP1A0(vpot,1,i,j,k+N3NOT1)       =vpotlocal[1];
      if(Nvec[2]==1)      MACP1A0(vpot,2,i,j,k+N3NOT1)       =vpotlocal[2];

      // A_x A_z [up j]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
      if(Nvec[1]==1)      MACP1A0(vpot,1,i,j+N2NOT1,k)       =vpotlocal[1];
      if(Nvec[3]==1)      MACP1A0(vpot,3,i,j+N2NOT1,k)       =vpotlocal[3];

      // A_y A_z [up i]
      get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
      if(Nvec[2]==1)      MACP1A0(vpot,2,i+N1NOT1,j,k)       =vpotlocal[2];
      if(Nvec[3]==1)      MACP1A0(vpot,3,i+N1NOT1,j,k)       =vpotlocal[3];


      // A_x [up j-k]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k+N3NOT1,prim,vpotlocal);
      if(Nvec[1]==1)      MACP1A0(vpot,1,i,j+N2NOT1,k+N3NOT1)       =vpotlocal[1];

      // A_y [up i-k]
      get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k+N3NOT1,prim,vpotlocal);
      if(Nvec[2]==1)      MACP1A0(vpot,2,i+N1NOT1,j,k+N3NOT1)       =vpotlocal[2];

      // A_z [up i-j]
      get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j+N2NOT1,k,prim,vpotlocal);
      if(Nvec[3]==1)      MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k)       =vpotlocal[3];

    }
  }
}




// Fix A_i to the init_vpot_user() version fo A_i within the NS and on the NS surface
// NOTE: this can't know future evolved A_i (e.g. non-trivial) except *within* NS can set as a beginning to boundary condition.  But must avoid surface A_i!
// But, with new deeppara, this dosetfix is not needed or wanted or correct.
void adjust_fluxctstag_vpot_dosetfix(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int insideNS;
  int i,j,k;
  int l;
  int loc;
  FTYPE vpotlocal[NDIM];
  int ii,jj,kk;


  //////////////
  //
  // First set interior and surface of NS to analytical solution
  //
  // Different loops for each A_i
  //
  // don't include actual or MPI boundaries in count since they don't form real domain of NS that must be resolved
  //////////////


  COMPFULLLOOP{ // First go over all non-surface + surface cells for simplicity


    if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1){

      // Set A_i in volume (lower CORNi of cell)
      // get_vpot_fluxctstag_primecoords gets A_i in PRIMECOORDS coords as required for vpot[]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
      MACP1A0(vpot,1,i,j,k)       =vpotlocal[1];
      MACP1A0(vpot,2,i,j,k)       =vpotlocal[2];
      MACP1A0(vpot,3,i,j,k)       =vpotlocal[3];

      // Set A_i on surfaces completely (other CORNi's of cell)

      // A_x A_y [up k]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
      MACP1A0(vpot,1,i,j,k+N3NOT1)       =vpotlocal[1];
      MACP1A0(vpot,2,i,j,k+N3NOT1)       =vpotlocal[2];

      // A_x A_z [up j]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
      MACP1A0(vpot,1,i,j+N2NOT1,k)       =vpotlocal[1];
      MACP1A0(vpot,3,i,j+N2NOT1,k)       =vpotlocal[3];

      // A_y A_z [up i]
      get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
      MACP1A0(vpot,2,i+N1NOT1,j,k)       =vpotlocal[2];
      MACP1A0(vpot,3,i+N1NOT1,j,k)       =vpotlocal[3];


      // A_x [up j-k]
      get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k+N3NOT1,prim,vpotlocal);
      MACP1A0(vpot,1,i,j+N2NOT1,k+N3NOT1)       =vpotlocal[1];

      // A_y [up i-k]
      get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k+N3NOT1,prim,vpotlocal);
      MACP1A0(vpot,2,i+N1NOT1,j,k+N3NOT1)       =vpotlocal[2];

      // A_z [up i-j]
      get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j+N2NOT1,k,prim,vpotlocal);
      MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k)       =vpotlocal[3];

    }
  }
}








// copy/extrapolate A_i to inside NS by 1 layer
void adjust_fluxctstag_vpot_dosetextrapdirect(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int insideNS;
  int i,j,k;
  int l;
  int loc;
  FTYPE vpotlocal[NDIM];
  int ii,jj,kk;




  //////////////
  //
  // Next go over all sets that can be extrapolated/copied so B^i is smooth for interpolation routines
  //
  // Do this first so some corner-ish positions will be first extrapolated/copied and then in the next COMPLOOPN they will be set correctly if need to be fixed
  // That is, some 2-cell zones (if viewed per boundary-shell-active 2-cell connection) will be extrapolted to even if later determined to be on shell if using another 2-cell pair.
  // Easier than checking if a given point position for A_i is on a corner-ish cell that would require checking multiple ii,jj,kk  for each A_i position.
  // Easier to do per 2-cell and just use fact that setting overrides copy/extrapolation.
  //
  // For B^n = B^x : Copy A_x into BC from active grid
  //                 Linearly interpolate A_z into BC
  //                 Linearly interpolate A_y into BC
  // etc.
  //
  //
  // For copy/extrpolate, just use grid-directed copy/extrapolation.  No need to use corners since interpolation (which will use this information) is anyways 1D per dimension.
  //
  // NOTE: Only 1 layer deep
  //
  // Need to account for possibly multiple interps/extraps to same point -- so add-up values and divide by number of them -- but then have to keep track of #
  //
  //////////////

  int offc,off1,off2,nooff;
  FTYPE y1,y2,mm,xx,bb;
  COMPFULLLOOP{


    int dir;
    for(dir=1;dir<=3;dir++){

      int pos;
      // set pos for computing geometry and primitive
      if(dir==1) pos=CORN1;
      else if(dir==2) pos=CORN2;
      else if(dir==3) pos=CORN3;

      int odir1,odir2;
      // get odir's
      get_odirs(dir,&odir1,&odir2);

#if(0)
      // DEBUG:
      dualfprintf(fail_file,"Gotoutside %d %d %d :: dir=%d odir1=%d odir2=%d\n",i,j,k,dir,odir1,odir2);
      dualfprintf(fail_file,"outside: %d %d %d %d\n",    (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKINSIDE))                                     ,
                  (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKINSIDE))                   ,
                  (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKINSIDE))                   ,
                  (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKINSIDE))
                  );
#endif
 

      FTYPE count=0.0;
      FTYPE vpotlocalsingle=0.0;

      // Only bound if A_i  is such that a boundary cell, which occurs if *all* cells touching A_i are boundary cells.  Otherwise, should be set as surface value.
      if(
         (i>=-N1BND+N1NOT1 && j>=-N2BND+N2NOT1 && k>=-N3BND+N3NOT1)         && // avoid out of bounds
         (GLOBALMACP0A1(nsmask,i                                    ,j                                    ,k,NSMASKINSIDE)==1)                                     && 
         (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                  ,j-(odir1==2)*N2NOT1                  ,k-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)                   && 
         (GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                  ,j-(odir2==2)*N2NOT1                  ,k-(odir2==3)*N3NOT1,NSMASKINSIDE)==1)                   && 
         (GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1-(odir2==1)*N1NOT1,j-(odir1==2)*N2NOT1-(odir2==2)*N2NOT1,k-(odir2==3)*N3NOT1-(odir1==3)*N3NOT1,NSMASKINSIDE)==1)
         ){



#if(0)
        // DEBUG:
        dualfprintf(fail_file,"Gotinside %d %d %d\n",i,j,k);
        dualfprintf(fail_file,"inside: %d %d %d %d %d %d\n",
                    GLOBALMACP0A1(nsmask,i+(dir==1)*N1NOT1                                    ,j+(dir==2)*N2NOT1                                    ,k+(dir==3)*N3NOT1,NSMASKSHELL) ,
                    GLOBALMACP0A1(nsmask,i-(dir==1)*N1NOT1                                    ,j-(dir==2)*N2NOT1                                    ,k-(dir==3)*N3NOT1,NSMASKSHELL) ,
                    GLOBALMACP0A1(nsmask,i+(odir2==1)*N1NOT1                                    ,j+(odir2==2)*N2NOT1                                    ,k+(odir2==3)*N3NOT1,NSMASKSHELL),
                    GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                                    ,j-(odir2==2)*N2NOT1                                    ,k-(odir2==3)*N3NOT1,NSMASKSHELL),
                    GLOBALMACP0A1(nsmask,i+(odir1==1)*N1NOT1                                    ,j+(odir1==2)*N2NOT1                                    ,k+(odir1==3)*N3NOT1,NSMASKSHELL),
                    GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                                    ,j-(odir1==2)*N2NOT1                                    ,k-(odir1==3)*N3NOT1,NSMASKSHELL)
                    );
#endif


        ///////////////
        // cases
        ///////////////

        if(GLOBALMACP0A1(nsmask,i+(dir==1)*N1NOT1                                    ,j+(dir==2)*N2NOT1                                    ,k+(dir==3)*N3NOT1,NSMASKSHELL)==1){
          // copy from +dir
          count++;
          // add
          vpotlocalsingle += MACP1A0(vpot,dir,i+(dir==1)*N1NOT1                                    ,j+(dir==2)*N2NOT1                                    ,k+(dir==3)*N3NOT1);
        }
        if(GLOBALMACP0A1(nsmask,i-(dir==1)*N1NOT1                                    ,j-(dir==2)*N2NOT1                                    ,k-(dir==3)*N3NOT1,NSMASKSHELL)==1){
          // copy from -dir
          count++;
          // add
          vpotlocalsingle += MACP1A0(vpot,dir,i-(dir==1)*N1NOT1                                    ,j-(dir==2)*N2NOT1                                    ,k-(dir==3)*N3NOT1);
        }
        if(GLOBALMACP0A1(nsmask,i+(odir2==1)*N1NOT1                                    ,j+(odir2==2)*N2NOT1                                    ,k+(odir2==3)*N3NOT1,NSMASKSHELL)==1){
          // interp from up in odir2
          count++;

          y1=MACP1A0(vpot,dir,i+1*(odir2==1)*N1NOT1                                    ,j+1*(odir2==2)*N2NOT1                                    ,k+1*(odir2==3)*N3NOT1);
          y2=MACP1A0(vpot,dir,i+2*(odir2==1)*N1NOT1                                    ,j+2*(odir2==2)*N2NOT1                                    ,k+2*(odir2==3)*N3NOT1);
          mm=(y2-y1);
          bb=y1;
          xx=-1.0;
          // add
          vpotlocalsingle+=mm*xx + bb;
        }
        if(GLOBALMACP0A1(nsmask,i-(odir2==1)*N1NOT1                                    ,j-(odir2==2)*N2NOT1                                    ,k-(odir2==3)*N3NOT1,NSMASKSHELL)==1){
          // interp from down in odir2
          count++;

          y1=MACP1A0(vpot,dir,i-1*(odir2==1)*N1NOT1                                    ,j-1*(odir2==2)*N2NOT1                                    ,k-1*(odir2==3)*N3NOT1);
          y2=MACP1A0(vpot,dir,i-2*(odir2==1)*N1NOT1                                    ,j-2*(odir2==2)*N2NOT1                                    ,k-2*(odir2==3)*N3NOT1);
          mm=(y2-y1);
          bb=y1;
          xx=-1.0;
          // add
          vpotlocalsingle+=mm*xx + bb;
        }
        if(GLOBALMACP0A1(nsmask,i+(odir1==1)*N1NOT1                                    ,j+(odir1==2)*N2NOT1                                    ,k+(odir1==3)*N3NOT1,NSMASKSHELL)==1){
          // interp from up in odir1
          count++;

          y1=MACP1A0(vpot,dir,i+1*(odir1==1)*N1NOT1                                    ,j+1*(odir1==2)*N2NOT1                                    ,k+1*(odir1==3)*N3NOT1);
          y2=MACP1A0(vpot,dir,i+2*(odir1==1)*N1NOT1                                    ,j+2*(odir1==2)*N2NOT1                                    ,k+2*(odir1==3)*N3NOT1);
          mm=(y2-y1);
          bb=y1;
          xx=-1.0;
          // add
          vpotlocalsingle+=mm*xx + bb;
        }
        if(GLOBALMACP0A1(nsmask,i-(odir1==1)*N1NOT1                                    ,j-(odir1==2)*N2NOT1                                    ,k-(odir1==3)*N3NOT1,NSMASKSHELL)==1){
          // interp from down in odir1
          count++;

          y1=MACP1A0(vpot,dir,i-1*(odir1==1)*N1NOT1                                    ,j-1*(odir1==2)*N2NOT1                                    ,k-1*(odir1==3)*N3NOT1);
          y2=MACP1A0(vpot,dir,i-2*(odir1==1)*N1NOT1                                    ,j-2*(odir1==2)*N2NOT1                                    ,k-2*(odir1==3)*N3NOT1);
          mm=(y2-y1);
          bb=y1;
          xx=-1.0;
          // add
          vpotlocalsingle+=mm*xx + bb;
        }

      }// done if A_i @ i,j,k at loc=CORN_i is a boundary cell of the kind where one just copies other cell's active value (i.e. interpolation along i-direction for A_i)


      // average all possibilities
      if(count>0.5){
        MACP1A0(vpot,dir,i,j,k)=vpotlocalsingle/count;

        // DEBUG:
        //   dualfprintf(fail_file,"DIRECT: %d %d %d : %21.15g %21.15g\n",i,j,k,count,vpotlocalsingle);

      }
      else{
        // then either part of active domain, deep inside NS, or surface of NS where A_i is fixed
      }


    }// end over dirs (i.e. A_{dir})

  }// end over COMPFULLLOOP
}// end DOSETEXTRAPDIRECT function
// copy/extrapolate A_i to inside NS
// super loop copied/emulated from boundary condition code for CENT quantitites
// GODMARK: NOTE: TODO: To do averaging, this uses NS surface values that are fixed (in time), while normal bound code only uses CENT and not fixed FACE values.
void adjust_fluxctstag_vpot_dosetextrapdirect_deep(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{




  //////////////
  //
  // Next go over all sets that can be extrapolated/copied so B^i is smooth for interpolation routines
  //
  // Do this first so some corner-ish positions will be first extrapolated/copied and then in the next COMPLOOPN they will be set correctly if need to be fixed
  // That is, some 2-cell zones (if viewed per boundary-shell-active 2-cell connection) will be extrapolted to even if later determined to be on shell if using another 2-cell pair.
  // Easier than checking if a given point position for A_i is on a corner-ish cell that would require checking multiple ii,jj,kk  for each A_i position.
  // Easier to do per 2-cell and just use fact that setting overrides copy/extrapolation.
  //
  // For B^n = B^x : Copy A_x into BC from active grid
  //                 Linearly interpolate A_z into BC
  //                 Linearly interpolate A_y into BC
  // etc.
  //
  //////////////


#pragma omp parallel OPENMPGLOBALPRIVATEFULL
  {
    int i,j,k;
    int ii,jj,kk;
    int l;
    FTYPE vpotlocal[NDIM];
    int dir;
    int pos;
    int odir1,odir2;
    FTYPE vpotlocalsingle;
    struct of_geom geomdontuse;
    struct of_geom *ptrgeom=&geomdontuse;
    struct of_geom geomdontuseii;
    struct of_geom *ptrgeomii=&geomdontuseii;
    int isinsideNSijk;    
    int isinsideNSiijjkk;    
    int hasmaskijk,hasinsideijk,reallyonsurfaceijk,cancopyfromijk;
    FTYPE weight4[4],function4[4];
    int ijk4[3][4];


    
    OPENMP3DLOOPVARSDEFINE;
    ////////  COMPFULLLOOP{
    OPENMP3DLOOPSETUPFULL;
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);



      for(dir=1;dir<=3;dir++){

        // set pos for computing geometry and primitive (not done so far)
        if(dir==1) pos=CORN1;
        else if(dir==2) pos=CORN2;
        else if(dir==3) pos=CORN3;

        // get odir's
        get_odirs(dir,&odir1,&odir2);

        // Only bound if A_i  is such that a boundary cell, which occurs if *all* cells directly *touching* A_i are boundary cells.  Otherwise, should be set as surface value.
        isinsideNSijk=is_dir_insideNS(dir,i,j,k, &hasmaskijk, &hasinsideijk, &reallyonsurfaceijk, &cancopyfromijk);


        if(isinsideNSijk){
          // if inside NS, then use shell values to form average (distance weighted) answer for boundary values
     
 
          // even for a single A_i at a single ii,jj,kk, have parallel and perpendicular contributions:
          // add copy contribution for differences in i-direction for A_i
          // add linear interpolation/extrapolation contribution for other directions

          FTYPE rescalefunc=rescale_A(fluxtime,dir,i,j,k);


          // temp things
          FTYPE weight,iweight;
          FTYPE totalweight;
          FTYPE maxdist;
   
          // initialize
          vpotlocalsingle=0.0;
          totalweight=0.0;
          // avoid using completely opposite side of NS, where field values flip, so only go 2-3 cells beyond minimum distance to shell.
          // assume (int) on dist=sqrt() would potentially lead to 0.9->0, so minimum value is 1
          maxdist=GLOBALMACP0A1(nsmask,i,j,k,NSMASKDISTTOSHELL+dir)+2.0;

   
   
          // LOOP over possible ii,jj,kk for prim's to add to i,j,k prim
          COMPLOOPNiijjkk(ii,jj,kk){

            int iim=ii-N1NOT1;
            int jjm=jj-N2NOT1;
            int kkm=kk-N3NOT1;

            int iip=ii+N1NOT1;
            int jjp=jj+N2NOT1;
            int kkp=kk+N3NOT1;
            int hasmaskiijjkk,hasinsideiijjkk,reallyonsurfaceiijjkk,cancopyfromiijjkk;

            // see if ii,jj,kk point is inside for A_i
            int isongrid=(ii>=-N1BND && ii<=N1M-1 &&  jj>=-N2BND && jj<=N2M-1 &&  kk>=-N3BND && kk<=N3M-1);
            isinsideNSiijjkk=is_dir_insideNS(dir,ii,jj,kk, &hasmaskiijjkk, &hasinsideiijjkk, &reallyonsurfaceiijjkk, &cancopyfromiijjkk);

     
            // 0 : average out
            // 1 : other average out
            // 2 : bilinear (then only get 4 nearest values)
#define WHICHINTERPDEEP 0

            // has to be on grid
            // if a shell grid (then either fixed or evolved point, which can be used for interpolation to i,j,k)
            // can't itself be inside NS, but needs to have a shell
            //     if(isongrid && isinsideNSiijjkk==0 && hasmaskiijjkk==1 && hasinsideiijjkk==1){
            if(isongrid && cancopyfromiijjkk==1){
              // then found shell value, add to the i,j,k value
              FTYPE dist;
              dist=sqrt((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k));

              // weight is grid distance (never will be zero since never inside and shell at same position)
#if(WHICHINTERPDEEP==0)
              // seems more generally nice, but more speculative than when dealing with B^i directly, so may be more unstable than copying A_dir
              //       iweight=sqrt((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k));
              iweight=((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k));
              weight=1.0/iweight;
#elif(WHICHINTERPDEEP==1)
              // assume want lowest order to get B^i, so keep like non-deep version and copy independent of distance along dir for A_dir
              // not sure this makes sense to be so non-local in general
              if(dir==1) iweight=sqrt((jj-j)*(jj-j) + (kk-k)*(kk-k));
              else if(dir==2) iweight=sqrt((ii-i)*(ii-i) + (kk-k)*(kk-k));
              else if(dir==3) iweight=sqrt((ii-i)*(ii-i) + (jj-j)*(jj-j));
              weight=1.0/iweight;
#elif(WHICHINTERPDEEP==2)
              // bilinear
              weight=(1.0-(i-ii))*(1.0-(j-jj))*(1.0-(k-kk));
#endif


              if(dist<=maxdist){
  
                // add this ii,jj,kk value to i,j,k with weight

                // get gdet for Bp?
                //  get_geometry(ii, jj, kk, CENT, ptrgeomii);

                FTYPE rescalefunciijjkk;
                rescalefunciijjkk=rescale_A(fluxtime,dir,ii,jj,kk);
   

                vpotlocalsingle += (MACP1A0(vpot,dir,ii,jj,kk)*rescalefunciijjkk)*weight;
                totalweight+=weight;
              }
       
              // DEBUG:
              //       dualfprintf(fail_file,"WEIGHT_deep: %d %d %d : %d %d %d : %21.15g %21.15g\n",i,j,k,ii,jj,kk,iweight,totalweight);

            }// end if found grid shell
          }// end over ii,jj,kk

          // normalize
          if(totalweight>0.0){
            // unrescale using local scale factor
            MACP1A0(vpot,dir,i,j,k)=(vpotlocalsingle/rescalefunc)/totalweight;
          }
          else{
            dualfprintf(fail_file,"adjust_fluxctstag_vpot_dosetextrapdirect_deep(): Never found shell to use for: %d %d %d.  Means inside NS, but no nearby shell (can occur for MPI).  Assume if so, then value doesn't matter (not used), so revert to fixed initial value for NS.\n",i,j,k);

            get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
            MACP1A0(vpot,dir,i,j,k)=vpotlocal[dir];

          }

        }// end if inside NS

      }// end over dirs (i.e. A_{dir})
      
    }// end loop block
  }// end parallel region

      
}// end DOSETEXTRAPDIRECTDEEP function



// generic abusive extrapolation of A_i into NS (no longer overwrites fixed A_i, but not complete copy/extrapolation when multiple cells involved)
void adjust_fluxctstag_vpot_dosetextrap(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int insideNS;
  int i,j,k;
  int l;
  int loc;
  FTYPE vpotlocal[NDIM];
  int ii,jj,kk;


  //////////////
  //
  // Next go over all sets that can be extrapolated/copied so B^i is smooth for interpolation routines
  //
  // Do this first so some corner-ish positions will be first extrapolated/copied and then in the next COMPLOOPN they will be set correctly if need to be fixed
  // That is, some 2-cell zones (if viewed per boundary-shell-active 2-cell connection) will be extrapolted to even if later determined to be on shell if using another 2-cell pair.
  // Easier than checking if a given point position for A_i is on a corner-ish cell that would require checking multiple ii,jj,kk  for each A_i position.
  // Easier to do per 2-cell and just use fact that setting overrides copy/extrapolation.
  //
  // For B^n = B^x : Copy A_x into BC from active grid
  //                 Linearly interpolate A_z into BC
  //                 Linearly interpolate A_y into BC
  // etc.
  //
  //
  // For copy/extrpolate, just use grid-directed copy/extrapolation.  No need to use corners since interpolation (which will use this information) is anyways 1D per dimension.
  //
  // NOTE: Only 1 layer deep
  //
  // need to account for possibly multiple interps/extraps to same point -- so add-up values and divide by number of them -- but then have to keep track of #
  //
  //////////////

  int offc,off1,off2,nooff;
  FTYPE y1,y2,mm,xx,bb;
  COMPFULLLOOP{

    if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1){


      // first 
      for(kk=k-N3NOT1;kk<=k+N3NOT1;kk++){ // TODO: avoid boundary cells?
        for(jj=j-N2NOT1;jj<=j+N2NOT1;jj++){
          for(ii=i-N1NOT1;ii<=i+N1NOT1;ii++){

            if(GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==1){
              // this should avoid trigger on ii==i && jj==j && kk==k since no cell can be both a mask and shell

              // +-x
              if(k==kk){
                if(j==jj){
                  if(i+N1NOT1==ii || i-N1NOT1==ii){// |   i inside   |   shell   ||    shell   |    i inside   |

                    // offc: centered offset
                    // off1,off2,nooff: face offsets
                    // xx is -1.0 in both cases because of how y1,y2 chosen as y2 always with larger (- or +) offset
                    if     (i+N1NOT1==ii){ offc=+N1NOT1; off1=+N1NOT1; off2=+2*N1NOT1; nooff=0;       xx=-1.0; } // shell is to right
                    else if(i-N1NOT1==ii){ offc=-N1NOT1; off1=0;       off2=-N1NOT1;   nooff=+N1NOT1; xx=-1.0; }    // shell is to left

                    // A_x for B^{n=x}: Copy A_x @ i+offc,j&j+1,k&k+1 -> i,j&j+1,k&k+1
                    if(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k-N3NOT1,NSMASKINSIDE)==1){
                      MACP1A0(vpot,1,i,j,k)=MACP1A0(vpot,1,i+offc,j,k);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1){
                      MACP1A0(vpot,1,i,j+N2NOT1,k)=MACP1A0(vpot,1,i+offc,j+N2NOT1,k);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1){
                      MACP1A0(vpot,1,i,j,k+N3NOT1)=MACP1A0(vpot,1,i+offc,j,k+N3NOT1);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1){
                      MACP1A0(vpot,1,i,j+N2NOT1,k+N3NOT1)=MACP1A0(vpot,1,i+offc,j+N2NOT1,k+N3NOT1);
                    }

                    if(
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N1NOT1 && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1)
                       ){
                      // A_y for B^{n=x}: Extrapolate @ i+off1/2,j,k -> i+nooff,j,k
                      y1=MACP1A0(vpot,2,i+off1,j,k); y2=MACP1A0(vpot,2,i+off2,j,k);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,2,i+nooff,j,k)=mm*xx + bb;
                    }

                    if(
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N1NOT1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1)
                       ){
                      // A_y for B^{n=x}: Extrapolate @ i+off1/2,j,k+1 -> i+nooff,j,k+1
                      y1=MACP1A0(vpot,2,i+off1,j,k+N3NOT1); y2=MACP1A0(vpot,2,i+off2,j,k+N3NOT1);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,2,i+nooff,j,k+N3NOT1)=mm*xx + bb;
                    }

                    if(
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N1NOT1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1)
                       ){
                      // A_z for B^{n=x}: Extrapolate @ i+off1/2,j,k -> i+nooff,j,k
                      y1=MACP1A0(vpot,3,i+off1,j,k); y2=MACP1A0(vpot,3,i+off2,j,k);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,3,i+nooff,j,k)=mm*xx + bb;
                    }

                    if(
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N1NOT1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1)
                       ){
                      // A_z for B^{n=x}: Extrapolate @ i+off1/2,j+1,k -> i+nooff,j+1,k
                      y1=MACP1A0(vpot,3,i+off1,j+N2NOT1,k); y2=MACP1A0(vpot,3,i+off2,j+N2NOT1,k);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,3,i+nooff,j+N2NOT1,k)=mm*xx + bb;
                    }

                  }
                }
              }


              // +-y
              if(k==kk){
                if(i==ii){
                  if(j+N2NOT1==jj || j-N2NOT1==jj){

                    // offc: centered offset
                    // off1,off2,nooff: face offsets
                    // xx is -1.0 in both cases because of how y1,y2 chosen as y2 always with larger (- or +) offset
                    if     (j+N2NOT1==jj){ offc=+N2NOT1; off1=+N2NOT1; off2=+2*N2NOT1; nooff=0;       xx=-1.0; }
                    else if(j-N2NOT1==jj){ offc=-N2NOT1; off1=0;       off2=-N2NOT1;   nooff=+N2NOT1; xx=-1.0; }

                    // A_y for B^{n=y}: Copy A_y @ i&i+1,j+offc,k&k+1 -> i&i+1,j,k&k+1
                    if(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1){ // dup. of 1st above nooff=0 from +-x
                      MACP1A0(vpot,2,i,j,k)=MACP1A0(vpot,2,i,j+offc,k);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1){ // dup. of 1st above nooff=+N1NOT1 from +-x
                      MACP1A0(vpot,2,i+N1NOT1,j,k)=MACP1A0(vpot,2,i+N1NOT1,j+offc,k);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1){ // dup. of 2nd above nooff=0 from +-x
                      MACP1A0(vpot,2,i,j,k+N3NOT1)=MACP1A0(vpot,2,i,j+offc,k+N3NOT1);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1){ // dup. of 2nd above nooff=+N1NOT1 from +-x
                      MACP1A0(vpot,2,i+N1NOT1,j,k+N3NOT1)=MACP1A0(vpot,2,i+N1NOT1,j+offc,k+N3NOT1);
                    }

                    if(// dup of a couple A_x terms from +-x and N1NOT1->N2NOT1
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k-N3NOT1,NSMASKINSIDE)==1)||
                       (nooff==+N2NOT1 && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1)
                       ){
                      // A_x for B^{n=y}: Extrapolate @ i,j+off1/2,k -> i,j+nooff,k
                      y1=MACP1A0(vpot,1,i,j+off1,k); y2=MACP1A0(vpot,1,i,j+off2,k);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,1,i,j+nooff,k)=mm*xx + bb;
                    }

                    if(// dup of a couple A_x terms from +-x and N1NOT1->N2NOT1
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1)||
                       (nooff==+N2NOT1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1)
                       ){
                      // A_x for B^{n=y}: Extrapolate @ i,j+off1/2,k+1 -> i,j+nooff,k+1
                      y1=MACP1A0(vpot,1,i,j+off1,k+N3NOT1); y2=MACP1A0(vpot,1,i,j+off2,k+N3NOT1);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,1,i,j+nooff,k+N3NOT1)=mm*xx + bb;
                    }

                    if(// dup of A_z from +-x but nooff==+N2NOT1 case modified and N1NOT1->N2NOT1
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N2NOT1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)
                       ){
                      // A_z for B^{n=y}: Extrapolate @ i,j+off1/2,k -> i,j+nooff,k
                      y1=MACP1A0(vpot,3,i,j+off1,k); y2=MACP1A0(vpot,3,i,j+off2,k);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,3,i,j+nooff,k)=mm*xx + bb;
                    }

                    if( // dup of A_z from +-x but nooff==0 & nooff==+N2NOT1 cases modified and N1NOT1->N2NOT1
                       (nooff==0       && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1)||
                       (nooff==+N2NOT1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1)
                        ){
                      // A_z for B^{n=y}: Extrapolate @ i+1,j+off1/2,k -> i+1,j+nooff,k
                      y1=MACP1A0(vpot,3,i+N1NOT1,j+off1,k); y2=MACP1A0(vpot,3,i+N1NOT1,j+off2,k);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,3,i+N1NOT1,j+nooff,k)=mm*xx + bb;
                    }

                  }
                }
              }

              // +-z
              if(i==ii){
                if(j==jj){
                  if(k+N3NOT1==kk || k-N3NOT1==kk){

                    // offc: centered offset
                    // off1,off2,nooff: face offsets
                    // xx is -1.0 in both cases because of how y1,y2 chosen as y2 always with larger (- or +) offset
                    if     (k+N3NOT1==kk){ offc=+N3NOT1; off1=+N3NOT1; off2=+2*N3NOT1; nooff=0;       xx=-1.0; }
                    else if(k-N3NOT1==kk){ offc=-N3NOT1; off1=0;       off2=-N3NOT1;   nooff=+N3NOT1; xx=-1.0; }

                    // A_z for B^{n=z}: Copy A_z @ i&i+1,j&j+1,k+offc -> i&i+1,j&j+1,k
                    if(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1){ // perfect dup
                      MACP1A0(vpot,3,i,j,k)=MACP1A0(vpot,3,i,j,k+offc);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1){ // perfect dup
                      MACP1A0(vpot,3,i+N1NOT1,j,k)=MACP1A0(vpot,3,i+N1NOT1,j,k+offc);
                    }
                    if(GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1){ // perfect dup
                      MACP1A0(vpot,3,i,j+N2NOT1,k)=MACP1A0(vpot,3,i,j+N2NOT1,k+offc);
                    }
                    if(GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j+N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1){ // perfect dup
                      MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k)=MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k+offc);
                    }

                    // A_y for B^{n=z}: Extrapolate @ i,j,k+off1/2 -> i,j,k+nooff
                    if( // perfect dup from 2 places
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N3NOT1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i-N1NOT1,j,k,NSMASKINSIDE)==1)
                        ){
                      y1=MACP1A0(vpot,2,i,j,k+off1); y2=MACP1A0(vpot,2,i,j,k+off2);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,2,i,j,k+nooff)=mm*xx + bb;
                    }

                    // A_y for B^{n=z}: Extrapolate @ i+1,j,k+off1/2 -> i+1,j,k+nooff
                    if( // perfect dup from 2 places
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1)||
                       (nooff==+N3NOT1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i+N1NOT1,j,k,NSMASKINSIDE)==1)
                        ){
                      y1=MACP1A0(vpot,2,i+N1NOT1,j,k+off1); y2=MACP1A0(vpot,2,i+N1NOT1,j,k+off2);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,2,i+N1NOT1,j,k+nooff)=mm*xx + bb;
                    }
        
                    // A_x for B^{n=z}: Extrapolate @ i,j,k+off1/2 -> i,j,k+nooff
                    if( // perfect dup from 2 places
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k-N3NOT1,NSMASKINSIDE)==1)||
                       (nooff==+N3NOT1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j-N2NOT1,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1)
                        ){
                      y1=MACP1A0(vpot,1,i,j,k+off1); y2=MACP1A0(vpot,1,i,j,k+off2);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,1,i,j,k+nooff)=mm*xx + bb;
                    }

                    // A_x for B^{n=z}: Extrapolate @ i,j+1,k+off1/2 -> i,j+1,k+nooff
                    if( // perfect dup from 2 places
                       (nooff==0       && GLOBALMACP0A1(nsmask,i,j,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k-N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1)||
                       (nooff==+N3NOT1 && GLOBALMACP0A1(nsmask,i,j,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k+N3NOT1,NSMASKINSIDE)==1 && GLOBALMACP0A1(nsmask,i,j+N2NOT1,k,NSMASKINSIDE)==1)
                        ){
                      y1=MACP1A0(vpot,1,i,j+N2NOT1,k+off1); y2=MACP1A0(vpot,1,i,j+N2NOT1,k+off2);
                      mm=(y2-y1); bb=y1;  MACP1A0(vpot,1,i,j+N2NOT1,k+nooff)=mm*xx + bb;
                    }


                  }
                }
              }






            }// end if ii,jj,kk is active domain

          }
        }
      }
      
    }// if i,j,k is inside boundary

  }// end over COMPLOOPN
}// end DOSETEXTRAP function





// used with DOSETEXTRAP function so correctly ensures A_i's that should be fixed are fixed
void adjust_fluxctstag_vpot_dosetfinalforce(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*vpot)[NSTORE1+SHIFTSTORE1][NSTORE2+SHIFTSTORE2][NSTORE3+SHIFTSTORE3])
{
  int insideNS;
  int i,j,k;
  int l;
  int loc;
  FTYPE vpotlocal[NDIM];
  int ii,jj,kk;


  //////////////
  //
  // Now go over places to set (force) so B^n remains fixed to desired value as if perfect conductor
  //
  // This is not just redundant with DOSETFIX because DOSETEXTRAP is kept simple by being lazy about how extrapolate.
  // DOSETFINALFORCE ensures really fix perpendicular components of A_i to surface normal for B^n accounting exactly for situation of cell and boundary
  //
  //////////////


  COMPFULLLOOP{


    if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1){

      // 

      // masks based upon loc=CENT, and only processing A_i on boundary of that cell, so only go i and i+1
      // so 27-1=26 positions of which 8 are not connected by a surface or line with each other, leaving only 18 positions that have sharing line or surface between inside and outside
      for(kk=k-N3NOT1;kk<=k+N3NOT1;kk++){ // TODO: avoid boundary cells?
        for(jj=j-N2NOT1;jj<=j+N2NOT1;jj++){
          for(ii=i-N1NOT1;ii<=i+N1NOT1;ii++){

            if(GLOBALMACP0A1(nsmask,ii,jj,kk,NSMASKSHELL)==1){
              // this should avoid trigger on ii==i && jj==j && kk==k since no cell can be both a mask and shell
              // then normal points from inside to outside through shell, and so need to fix Aperp1 and Aperp2 on surface between CENT inside and CENT outside
       
              if(i==ii){ // 8 unique cells
                if(j==jj){
                  if(k+N3NOT1==kk){
                    // so set: A_x @ k+1 && (j || j+1) and  A_y @ k+1 (i || i+1)

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j,k+N3NOT1)       =vpotlocal[1];
                    MACP1A0(vpot,2,i,j,k+N3NOT1)       =vpotlocal[2];

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j+N2NOT1,k+N3NOT1)=vpotlocal[1];

                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,2,i+N1NOT1,j,k+N3NOT1)=vpotlocal[2];

                    // extrapolate/copy A_z @ 

                  }
                  else if(k-N3NOT1==kk){
                    // so set: A_x @ k && (j || j+1) and  A_y @ k (i || i+1)

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j,k)       =vpotlocal[1];
                    MACP1A0(vpot,2,i,j,k)       =vpotlocal[2];

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j+N2NOT1,k)=vpotlocal[1];

                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
                    MACP1A0(vpot,2,i+N1NOT1,j,k)=vpotlocal[2];
                  }
                }// end if j==jj
                else if(j+N2NOT1==jj){
                  if(k+N3NOT1==kk){
                    // so set: A_x @ j+1 k+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j+N2NOT1,k+N3NOT1)=vpotlocal[1];
                  }
                  else if(k-N3NOT1==kk){
                    // so set: A_x @ j+1 k
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j+N2NOT1,k)=vpotlocal[1];
                  }
                  else if(k==kk){
                    // so set: A_x @ j+1,k AND j+1,k+1
                    //    set: A_z @ j+1,i AND j+1,i+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j+N2NOT1,k)=vpotlocal[1];
                    MACP1A0(vpot,3,i,j+N2NOT1,k)=vpotlocal[3];

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j+N2NOT1,k+N3NOT1)=vpotlocal[1];

                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k)=vpotlocal[3];
                  }
                }
                else if(j-N2NOT1==jj){
                  if(k+N3NOT1==kk){
                    // so set: A_x @ j k+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j,k+N3NOT1)=vpotlocal[1];
                  }
                  else if(k-N3NOT1==kk){
                    // so set: A_x @ j k
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j,k)=vpotlocal[1];
                  }
                  else if(k==kk){
                    // so set: A_x @ j,k AND j,k+1
                    //    set: A_z @ j,i AND j,i+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j,k)=vpotlocal[1];
                    MACP1A0(vpot,3,i,j,k)=vpotlocal[3];

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,1,i,j,k+N3NOT1)=vpotlocal[1];

                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i+N1NOT1,j,k)=vpotlocal[3];
                  }
                }
              }
              else if(i-N1NOT1==ii){ // only new cells are over 5 cells since that's all that shares line boundary.  Others share CORNT only.
                if(j==jj){ // 3 positions
                  if(k==kk){
                    // so set A_y @ k,j AND A_y @ k+1,j AND A_z @ k,j AND A_z @ k,j+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
                    MACP1A0(vpot,2,i,j,k)       =vpotlocal[2];
                    MACP1A0(vpot,3,i,j,k)       =vpotlocal[3];

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,2,i,j,k+N3NOT1)       =vpotlocal[2];

                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i,j+N2NOT1,k)       =vpotlocal[3];
                  }
                  else if(k+N3NOT1==kk){
                    // so set: A_y @ k+1 && j
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,2,i,j,k+N3NOT1)=vpotlocal[2];
                  }
                  else if(k-N3NOT1==kk){
                    // so set: A_y @ k && j
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
                    MACP1A0(vpot,2,i,j,k)=vpotlocal[2];
                  }
                }// end if j==jj
                else if(j+N2NOT1==jj){
                  if(k==kk){
                    // so set: A_z @ j+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i,j+N2NOT1,k)=vpotlocal[3];
                  }
                }
                else if(j-N2NOT1==jj){
                  if(k==kk){
                    // so set: A_z @ j
                    get_vpot_fluxctstag_primecoords(fluxtime,i,j,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i,j,k)=vpotlocal[3];
                  }
                }
              }// end if i-N1NOT1==ii
              else if(i+N1NOT1==ii){ // only new cells are over 5 cells since that's all that shares line boundary.  Others share CORNT only.
                // same as i-N1NOT1==ii but i-> i+N1NOT1 for A_\mu
                if(j==jj){ // 3 positions
                  if(k==kk){
                    // so set A_y @ k,j AND A_y @ k+1,j AND A_z @ k,j AND A_z @ k,j+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
                    MACP1A0(vpot,2,i+N1NOT1,j,k)       =vpotlocal[2];
                    MACP1A0(vpot,3,i+N1NOT1,j,k)       =vpotlocal[3];

                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,2,i+N1NOT1,j,k+N3NOT1)       =vpotlocal[2];

                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k)       =vpotlocal[3];
                  }
                  else if(k+N3NOT1==kk){
                    // so set: A_y @ k+1 && j
                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k+N3NOT1,prim,vpotlocal);
                    MACP1A0(vpot,2,i+N1NOT1,j,k+N3NOT1)=vpotlocal[2];
                  }
                  else if(k-N3NOT1==kk){
                    // so set: A_y @ k && j
                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
                    MACP1A0(vpot,2,i+N1NOT1,j,k)=vpotlocal[2];
                  }
                }// end if j==jj
                else if(j+N2NOT1==jj){
                  if(k==kk){
                    // so set: A_z @ j+1
                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j+N2NOT1,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i+N1NOT1,j+N2NOT1,k)=vpotlocal[3];
                  }
                }
                else if(j-N2NOT1==jj){
                  if(k==kk){
                    // so set: A_z @ j
                    get_vpot_fluxctstag_primecoords(fluxtime,i+N1NOT1,j,k,prim,vpotlocal);
                    MACP1A0(vpot,3,i+N1NOT1,j,k)=vpotlocal[3];
                  }
                }
              }// end if i+N1NOT1==ii


            }// end if other cell is shell

          }// end ii
        }// end jj
      }//end kk


    }// end if insideNS


  } // end COMPLOOP
}// end DOSETFINALFORCE function


// special NS boundary code for EMFs for staggered fields
void adjust_fluxctstag_emfs_old1(SFTYPE fluxtime, FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
// GODMARK TODO DEBUG
//void adjust_fluxctstag_emfs(FTYPE (*prim)[NSTORE2][NSTORE3][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[NSTORE2][NSTORE3][NPR+NSPECIAL])
{

  // v1) force EMFperp1=EMFperp2=0 after computed EMF so that B^n stays fixed.
  // only good if object is stationary since EMF at midpoint in time, but still ok to do since A_i BC will fix-up the final field on the NS surface and interior

  // v2) No, don't want to set EMF to zero since EMF is at midpoint in time, so would force no update to that grid cell where NS was, while NS will be in new position invalidating exact EMF=0.  In essence, the presence of the NS creates a modified EMF (due to the conductive surface) that forces A_i on surface to be what NS wants.  So no need to force EMF itself, just force A_i and it will correspond to some effectively modified EMF (compared to flux update) for a given cell.


  // v3) No, erroneous EMF on NS boundary due to dissipative terms can lead to large fields just above surface.
  // Next timestep will take care of proper field evolution near NS given (e.g.) motion of NS through grid.



  int insideNS;
  int i,j,k,pl,pliter;
  struct of_gdetgeom gdetgeomdontuse;
  struct of_gdetgeom *ptrgdetgeom=&gdetgeomdontuse;
  int inittype;
  int ii,jj,kk;
  int dir,pos,odir1,odir2;
  int initreturn;
  int whichvel,whichcoord;
  FTYPE emf,gdet;
  FTYPE pr[NPR];
  int lll;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  FTYPE Bcon[NDIM],ucon[NDIM],vcon[NDIM],others[NUMOTHERSTATERESULTS];




  
  if(DOSETFIXEMF==0) return;


#if(0)
  // DEBUG:
  COMPLOOPP1{
    dualfprintf(fail_file,"IJKBAD1: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,1,i,j,k,B2));
    dualfprintf(fail_file,"IJKBAD2: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,2,i,j,k,B1));
    dualfprintf(fail_file,"IJKBAD3: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,1,i,j,k,B3));
    dualfprintf(fail_file,"IJKBAD4: %d %d %d : %21.15g\n",i,j,k,MACP1A1(fluxvec,2,i,j,k,B3));
  }   
#endif 




  // -1 means evolve type and set conditions as if inside NS, but don't have good interpolation for corner currently, so don't assume input is at loc==pos
  //  inittype=-1;

  // -3 means evolve type and set conditions as if inside NS, but only assume input primitive is good for B orthogonal to EMF_{dir} located at loc==pos.
  inittype=-3;



  // LOOP over i,j,k
  COMPFULLLOOP{



    // only operate if cell is inside NS
    if(GLOBALMACP0A1(nsmask,i,j,k,NSMASKINSIDE)==1){

      // Loop over possible EMF positions relevant for this cell at i,j,k
      // redundant for interior of NS, but easier to manage memory and whether inside NS or not for that surface
      for(kk=k;kk<=k+N3NOT1;kk++){
        for(jj=j;jj<=j+N2NOT1;jj++){
          for(ii=i;ii<=i+N1NOT1;ii++){

            // over dir (i.e. EMF's corners)
            for(dir=1;dir<=3;dir++){

              // condition for when EMF @ corner is on cell boundary for cell CORN_{dir} @ ii,jj,kk
              if(dir==1 && i==ii || dir==2 && j==jj || dir==3 & k==kk){

                // set pos for computing geometry and primitive
                if(dir==1) pos=CORN1;
                else if(dir==2) pos=CORN2;
                else if(dir==3) pos=CORN3;

                // get odir's
                get_odirs(dir,&odir1,&odir2);

                // skip if no such dimension relevant for non-zero EMF (assumes other dimension has all homogeneous terms (including geometry!)
                if(Nvec[odir1]==1 && Nvec[odir2]==1) continue;

                ////////////////////////
                //
                // Setup guess for primitive at CORN_{dir}
                //
                ////////////////////////
                PLOOP(pliter,pl) pr[pl]=MACP0A1(prim,ii,jj,kk,pl);


                ////////////////////////////////
                //
                // use actual field that was interpolated to CORN_i
                // Would have just computed these if doing evolution, which is only time EMF modified as in this function
                // don't worry about VEL since is forced within NS and on surface -- i.e. it's fixed while only originally have field at FACEs.
                // See fluxctstag.c for how pvbcorn accessed:
                // pvbcorninterp[dir = EMF_dir][ii,jj,kk][which v or B component: odir1 or odir2][0,1 for v +- in odir1, NUMCS for B][0,1 for v +- in odir2 or B jump that is perp to both dir and component chosen);
                // 
                FTYPE Bconl[NDIM],Bconr[NDIM];
                FTYPE Bcontouse[NDIM];
                int localisinside,localhasmask,localhasinside,localreallyonsurface,localcancopyfrom;

                Bconl[odir1]=GLOBALMACP1A3(pvbcorninterp,dir,ii,jj,kk,odir1,NUMCS,0); // jump in odir2
                Bconr[odir1]=GLOBALMACP1A3(pvbcorninterp,dir,ii,jj,kk,odir1,NUMCS,1); // jump in odir2
                //  GLOBALMACP1A3(pvbcorninterp,dir,ii,jj,kk,odir2,m,l); // vel
                Bconl[odir2]=GLOBALMACP1A3(pvbcorninterp,dir,ii,jj,kk,odir2,NUMCS,0); // jump in odir1
                Bconr[odir2]=GLOBALMACP1A3(pvbcorninterp,dir,ii,jj,kk,odir2,NUMCS,1); // jump in odir1
                //  GLOBALMACP1A3(pvbcorninterp,dir,ii,jj,kk,odir1,m,l); // vel

                // TODO: decide which corner is more inside NS (so more representative of desired surface field).  Might be better (more stable) than just averaging.
                // Also, could compute normal field analytically at CORN_i and use that -- most correct.
                Bcontouse[odir1]=0.5*(Bconl[odir1]+Bconr[odir1]);
                Bcontouse[odir2]=0.5*(Bconl[odir2]+Bconr[odir2]);


                // get "gdet" factor (really, whatever factor also used in fluxctstag.c:fluxcalc_fluxctstag_emf_1d() )
                // used to get true B^i since needed, and used later to get flux from emf
                get_geometry_geomeonly(ii, jj, kk, pos, ptrgdetgeom); // get geometry at CORN[dir] where emf is located
                gdet=ptrgdetgeom->EOMFUNCMAC(B1-1+odir1); // which field ptrgeom->e doesn't matter as mentioned below


#if(CORNGDETVERSION==1)
                // then pvbcorninterp does *not* contain any extra gdet factors, just B^i
                pr[B1+odir1-1]=Bcontouse[odir1];
                pr[B1+odir2-1]=Bcontouse[odir2];
#else
                // then pvbcorninterp contains gdet factor, so \detg B^i  needs to be converted to B^i
                FTYPE igdetnosing=sign(gdet)/(fabs(gdet)+SMALL);
                pr[B1+odir1-1]=Bcontouse[odir1]*igdetnosing;
                pr[B1+odir2-1]=Bcontouse[odir2]*igdetnosing;
#endif
                // still need to get pr[B1+dir-1], which will be from interpolation in init_dsandvels_nsbh()  

    
                ///////////////////////
                //
                // now constrain velocity (and possibly also field)
                //
                // get raw primitive in whichvel/whichcoord coordinates
                // Note: Don't need inputted (to function) prim[] since fixing values
                // prim is default field if interpolation-extrapolation to loc=pos fails to be contained entirely within NS
                // ok that don't fill "pstag" entry with pext[] type quantity, since not setting densities here and don't care about v||B
                ///////////////////////
                initreturn=init_dsandvels(inittype, pos, &whichvel, &whichcoord, fluxtime, ii, jj, kk, pr, NULL);

                // if successfully got raw primitive, then transform raw primitive to WHICHVEL/PRIMECOORD primitive and compute EMF
                if(initreturn>0){
                  FAILSTATEMENT("init.c:set_plpr_nsbh()", "init_dsandvels()", 1);
                }
                else if(initreturn==0){   
                  // general location transformation for v and B from whichvel/whichcoord to WHICHVEL/(MCOORD->PRIMECOORDS)
                  bl2met2metp2v_genloc(whichvel, whichcoord, pr, ii, jj, kk, pos);


                  // set B^\mu [internally uses pglobal and pstagglobal to set B, so uses previous time-step values for B -- problem for RK? GODMARK]
                  Bcon[TT]=0.0;
                  SLOOPA(lll) Bcon[lll] = pr[B1+lll-1];
      
                  // get u^\mu
                  get_geometry(ii, jj, kk, pos, ptrgeom);
                  ucon_calc(pr, ptrgeom, ucon, others);

                  // set v^i
                  vcon[TT]=0.0;
                  SLOOPA(lll) vcon[lll] = ucon[lll]/ucon[TT];


                  // now set EMF_i
                  // see fluxct.c for signature of emf[v,B] as related to fluxes
                  //
                  // so:
                  // emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]
                  // emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]
                  // emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]
                  //
                  // get_odirs(): dir=3 gives odir1=1 and odir2=2
                  // so cyclic in: dir,odir1,odir2
                  emf = Bcon[odir2] * vcon[odir1] - Bcon[odir1] * vcon[odir2];


                  // DEBUG:
                  if(Nvec[odir1]>1) dualfprintf(fail_file,"1EMF(%d): ijk=%d %d %d : old: %21.15g  new: %21.15g\n",dir,ii,jj,kk,MACP1A1(fluxvec,odir1,ii,jj,kk,B1-1+odir2),emf*gdet);
                  if(Nvec[odir2]>1) dualfprintf(fail_file,"2EMF(%d): ijk=%d %d %d : old: %21.15g  new: %21.15g\n",dir,ii,jj,kk,MACP1A1(fluxvec,odir2,ii,jj,kk,B1-1+odir1),emf*gdet);


#if(1)
                  if(dir==1 || dir==2){
                    // DEBUG: CHECK that \Omega_F is as expected
                    // get NS parameters
                    FTYPE rns,omegak,omegaf,Rns,Rsoft,v0;
                    setNSparms(fluxtime, &rns, &omegak, &omegaf,  &Rns, &Rsoft, &v0);

                    FTYPE dxdxp[NDIM][NDIM];
                    dxdxprim_ijk(ii,jj,kk,pos,dxdxp);

                    FTYPE omegaftest;
                    if(dir==1){
                      omegaftest=(-emf/Bcon[2])*dxdxp[PH][PH];
                    }
                    if(dir==2){
                      omegaftest=(emf/Bcon[1])*dxdxp[PH][PH];
                    }

                    dualfprintf(fail_file,"OMEGA: %d %d %d : %21.15g : %21.15g =?= %21.15g\n",i,j,k,omegak,omegaf,omegaftest);
                    if(fabs(omegaftest-omegaf)>1E-5){
                      dualfprintf(fail_file,"OMEGADIFF: %d %d %d : %21.15g : %21.15g =?= %21.15g\n",i,j,k,omegak,omegaf,omegaftest);
                    }
                  }
#endif

    
                  // assign
                  if(Nvec[odir1]>1) MACP1A1(fluxvec,odir1,ii,jj,kk,B1-1+odir2) = + emf*gdet;
                  if(Nvec[odir2]>1) MACP1A1(fluxvec,odir2,ii,jj,kk,B1-1+odir1) = - emf*gdet;
                  if(Nvec[dir]>1) MACP1A1(fluxvec,dir,ii,jj,kk,B1-1+dir)     =   0.0;
                }//end if initreturn==0
                // else initreturn<0 then no change (i.e. nothing to set for boundary condition)

              }// condition if EMF on shell of cell
            }// end over dir

          }// over ii
        }// over jj
      }// over kk

    }// end if inside NS

  }// end i,j,k LOOP



  return;


}

