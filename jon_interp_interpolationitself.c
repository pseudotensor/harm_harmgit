#include "decs.h"

 

// local declarations of local functions
static void new_xycoord(int h, int i, int j, int k, FTYPE *tc, FTYPE *xc, FTYPE *yc, FTYPE *zc);
static void new_coord(int h, int i,int j, int k, FTYPE *t, FTYPE *r, FTYPE *th, FTYPE *ph) ;
static int old_coord(FTYPE t, FTYPE r,FTYPE th,FTYPE ph,FTYPE *X) ;
static FTYPE old_distance(FTYPE t1, FTYPE x1, FTYPE y1, FTYPE z1,  FTYPE t2, FTYPE x2, FTYPE y2, FTYPE z2);
static FTYPE new_distances(FTYPE t1, FTYPE r1, FTYPE th1, FTYPE ph1,   FTYPE t2, FTYPE r2, FTYPE th2, FTYPE ph2);
static void old_xyzcoord(FTYPE t, FTYPE r, FTYPE th, FTYPE ph, FTYPE *tc, FTYPE *xc, FTYPE*yc, FTYPE *zc);
static int plane_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref, int hold, int iold, int jold, int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
static int bilinear_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref, int hold, int iold, int jold, int kold, unsigned char*****oldimage, FTYPE*****olddata, FTYPE *newdatatemp);
static int bilinear_interp_ij(int hold, int iold, int jold,int kold, FTYPE dhold, FTYPE diold,FTYPE djold,FTYPE dkold,unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
static int bilinear_interp_simple(int hold, int iold, int jold, int kold, FTYPE *X, FTYPE *startx, FTYPE *dx, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
static void interpicoord(FTYPE *X,int loc, int *h, int *i, int *j, int *k);
static int is_atboundary(int *holdvar, int *ioldvar, int *joldvar, int *koldvar, int *atboundary, int *nearboundary);
static int nearest_interp(int hold, int iold,int jold, int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
static int bicubic_interp_wrap(int nt, int nx, int ny, int nz, int hold, int iold, int jold, int kold, FTYPE x0, FTYPE x1, FTYPE x2, FTYPE x3, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
static void X2spc(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_curr);

static int root_find1(FTYPE t, FTYPE r,FTYPE th, FTYPE ph, FTYPE*X) ;
static int root_find2(FTYPE t, FTYPE r, FTYPE th, FTYPE ph, FTYPE *X);  

static int nearest_interp_ij(int hold, int iold,int jold,int kold,unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
static void old_ijk2rthph(int hold, int iold, int jold, int kold, FTYPE *t, FTYPE*r, FTYPE*th, FTYPE *ph);
static void oldf_ijk2rthph(FTYPE hold, FTYPE iold, FTYPE jold, FTYPE kold, FTYPE *t, FTYPE*r, FTYPE*th, FTYPE *ph);
static void dervs_for_bicubic(int nx, int ny, FTYPE **ya, FTYPE **fy1a, FTYPE **fy2a, FTYPE **fy12a); // only 2D
static void oldf_ijk2x123(FTYPE hold, FTYPE iold, FTYPE jold, FTYPE kold, FTYPE*X);
static void old_ijk2x123(int hold, int iold, int jold, int kold, FTYPE*X);
static void old_x1232rthphcoord(FTYPE *X, FTYPE *t, FTYPE *r, FTYPE *th, FTYPE *ph);

static int root_find_parameter1(FTYPE *parms, FTYPE *returnvalue);
static void resid_parameter1(int n, FTYPE *parms, FTYPE *guess, FTYPE *diff);
static int usrfun_parameter1(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha);
static int usrfun_joninterp(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *beta, FTYPE **alpha, FTYPE*norm);
static int usrfun_joninterp2(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha);


static void low2high(int ntlow, int nxlow, int nylow, int nzlow,  int nthigh, int nxhigh, int nyhigh, int nzhigh, unsigned char*****oldimage,FTYPE*****olddata);
static void high2low(int nthigh, int nxhigh, int nyhigh, int nzhigh,  int ntlow, int nxlow, int nylow, int nzlow, unsigned char*****oldimage,FTYPE*****olddata);

static int lowericoord(FTYPE *X,int loc, int *h, int *i, int *j, int *k);







// just copy over data from old -> new
//unsigned char*****oldimage,FTYPE*****olddata,unsigned char*****newimage,FTYPE*****newdata
void copy_old2new(void)
{
  int coli,h,i,j,k;

  if(ALLOCATENEWIMAGEDATA){
    if(DATATYPE==0) for(coli=0;coli<numoutputcols;coli++) LOOPINTERP newimage[coli][h][i][j][k]=oldimage[coli][h][i][j][k];
    else for(coli=0;coli<numoutputcols;coli++) LOOPINTERP newdata[coli][h][i][j][k]=olddata[coli][h][i][j][k];
  }
  else{
    // otherwise write to file directly instead of storing in newimage or newdata
    OUTPUTLOOP{
      output2file_perpointcoli_postinterpolation(oldimage[coli][h][i][j][k], olddata[coli][h][i][j][k]);
    }
  }

}



// function that actually does interpolation
//unsigned char*****oldimage,FTYPE*****olddata,unsigned char*****newimage,FTYPE*****newdata
int compute_spatial_interpolation(void)
{
  //
  // used for image specific details
  int interptypetodo,getinterp;
  //
  int coli;
  //
  int hold,iold,jold,kold ;
  /////////////////////////////
  // coordinate stuff
  int h,i,j,k;
  int hh,ii,jj,kk;
  //  FTYPE Xmax[NDIM];
  FTYPE t,t2,r,r2,th,th2,ph,ph2;

  /////////////////////////////
  // special interpolation stuff
  int nosolution;
  FTYPE dhold,diold,djold,dkold;
  unsigned char uctemp;

  // setup temp column space
  FTYPE *newdatatemp=(FTYPE*)malloc((unsigned)(numoutputcols)*sizeof(FTYPE));
  if(newdatatemp==NULL){
    fprintf(stderr,"Couldn't allocate newdatatemp\n");
    exit(1);
  }

  // setup temp column space
  unsigned char *newimagetemp=(unsigned char*)malloc((unsigned)(numoutputcols)*sizeof(unsigned char));
  if(newimagetemp==NULL){
    fprintf(stderr,"Couldn't allocate newimagetemp\n");
    exit(1);
  }

  
  if(DATATYPE==0)   fprintf(stderr,"start to interpolate image (in 3D, each . is each k from %d to %d)\n",0,nN3);
  else fprintf(stderr,"start to interpolate data (in 3D, each . is each k from %d to %d)\n",0,nN3);
  fflush(stderr);


  //interpolate to new image


  if(DEBUGINTERP){
    fprintf(stderr,"%g %g ::%g %g :: %g %g :: %g %g\n",startx[0],Xmax[0],startx[1],Xmax[1],startx[2],Xmax[2],startx[3],Xmax[3]) ; fflush(stderr);
  }



  


  
  int kprogress=0;
  LOOPINTERP{ // no coli in LOOPINTERP!  COLIMARK: loop over columns (coli) is done for each position, so don't have to repeat physical space-time mapping for each data value

    // assume if 3D want indicator of progress every k
    if(k==kprogress){
      fprintf(stderr,".");
      kprogress++;
    }

    // DEBUG
    //    fprintf(stderr,"%d %d %d %d\n",h,i,j,k); fflush(stderr);


    if(DEBUGINTERP){
      if(j==0) fprintf(stderr,"%d %d %d %d\n",h,i,j,k); fflush(stderr);
    }

    ////////////////////////////
    // coordinate solver
    ///////////////////////////

    // find new real coordinates (e.g. real t,r,th,ph) for h,i,j,k
    new_coord(h,i,j,k,&t,&r,&th,&ph) ;

    //HACK
    //    if(ph>0 && ph<M_PI/2.0){
    //      r=100;
    //    }

    
    // DEBUG
    //    fprintf(stderr,"pre-old_coord(): %d %d %d %d :: %g %g %g %g\n",h,i,j,k,t,r,th,ph); fflush(stderr);


    // find old grid's X from that new grid's t,r,th,ph
    // input must be those coordinates (here t,r,th,ph) corresponding to metric's original coordinates after transforming back using dxdxp from PRIMECOORDS.  Otherwise, not seeking correct solution in root finder that obtains X
    // In cases where X must diverge to +-infinity to reach desired r,th,ph, assume this is out of bounds corresponding to nosolution==2
    nosolution=old_coord(t,r,th,ph,X) ;


    // DEBUG
    //    fprintf(stderr,"pre-interp: %d %d %d %d :: %g %g %g %g\n",h,i,j,k,t,r,th,ph); fflush(stderr);


    ////////////////////////
    // interpolate
    /////////////////////////

    
    // shift coordinate if periodic and went "out of bounds"
    if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
      if(X[3]-Xmax[3] >= 0){
        X[3] = X[3] - (Xmax[3]-startx[3]);
      }
      else if(X[3]-startx[3] < 0){
        X[3] = X[3] + (Xmax[3]-startx[3]);
      }
    }

    
    if(DATATYPE==0){
      // make black if no solution or if out of bounds
      if(nosolution){
        for(coli=0;coli<numoutputcols;coli++) newimagetemp[coli] = FLOAT2IMAGE(defaultvalue[coli]);
        getinterp=0;
      }
      else if(
              (oN0!=1)&&((X[0]-startx[0] < 0) || (X[0]-Xmax[0] >= 0)) ||
              (oN1!=1)&&((X[1]-startx[1] < 0) || (X[1]-Xmax[1] >= 0)) ||
              (oN2!=1)&&((X[2]-startx[2] < 0) || (X[2]-Xmax[2] >= 0)) ||
              (oN3!=1)&&((X[3]-startx[3] < 0) || (X[3]-Xmax[3] >= 0))
              ) {

        if(EXTRAPOLATE){
          getinterp=1;
          interptypetodo=0;
        }
        else{
          for(coli=0;coli<numoutputcols;coli++) newimagetemp[coli] = FLOAT2IMAGE(defaultvalue[coli]);
          getinterp=0;
          interptypetodo=0;
        }
      }
      else{
        getinterp=1;
        interptypetodo=INTERPTYPE;
      }
    }
    else{
      // make defaultvalue only if no solution
      if(nosolution) {
        for(coli=0;coli<numoutputcols;coli++) newdatatemp[coli]=defaultvalue[coli];
        getinterp=0;
      }
      else if(
              (oN0!=1)&&((X[0]-startx[0] < 0) || (X[0]-Xmax[0] >= 0)) ||
              (oN1!=1)&&((X[1]-startx[1] < 0) || (X[1]-Xmax[1] >= 0)) ||
              (oN2!=1)&&((X[2]-startx[2] < 0) || (X[2]-Xmax[2] >= 0)) ||
              (oN3!=1)&&((X[3]-startx[3] < 0) || (X[3]-Xmax[3] >= 0))
              ){
        if(EXTRAPOLATE){
          getinterp=1;
          interptypetodo=0;
        }
        else{
          for(coli=0;coli<numoutputcols;coli++) newdatatemp[coli]=defaultvalue[coli];
          getinterp=0;
          interptypetodo=0;
        }
      }
      else{
        // normal interpolation
        getinterp=1;
        interptypetodo=INTERPTYPE;
      }
    }

    // DEBUG
    if(DEBUGINTERP){
      fprintf(stderr,"get=%d nosol=%d :: %d %d %d %d :: %g %g %g %g :: %g %g %g %g\n",getinterp,nosolution, h,i,j,k, t,r,th,ph, X[0],X[1],X[2],X[3]) ; fflush(stderr);
    }




    if(getinterp){ // done for both DATATYPE==0 and 1.  Fills newdatatemp with result, so have to copy to newimagetemp at the end


      // get old grid's h,i,j,k from determined X
      interpicoord(X,CENT,&hold,&iold,&jold,&kold);

      //      fprintf(stderr,"DOG: %d %d %d %d\n",hold,iold,jold,kold);

      // GODMARK: for now only use nearest neighbor in k -> kold
      //      if(DEBUGINTERP){
      // if(kold==15 || kold==0){
      //   fprintf(stderr,"interp: %d %d %d\n",hold,iold,jold,kold);
      //   fflush(stderr);
      // }
      //      }

      // check if symmetric -- yes, it is
      //      if(j==nN2/3 && (i==nN1/4 || i==3*nN1/4-1)){
      //       fprintf(stderr,"i=%d j=%d r=%21.15g th=%21.15g X[1]=%21.15g X[2]=%21.15g iold=%d jold=%d\n",i,j,r,th,X[1],X[2],iold,jold);
      //      }

      // DEBUG
      if(DEBUGINTERP){
        fprintf(stderr,"doing hold=%d iold=%d jold=%d kold=%d h=%d i=%d j=%d k=%d\n",hold,iold,jold,kold,h,i,j,k); fflush(stderr);
      }

      //////////////////////////////
      // TRANSPLANT
      //////////////////////////////

      // DEBUG
      if(DEBUGINTERP){
        fprintf(stderr,"Done copy\n"); fflush(stderr);
      }

      if(interptypetodo==0){
        nearest_interp(hold,iold,jold,kold,oldimage,olddata, newdatatemp);
      }
      else if(interptypetodo==1){
        bilinear_interp_simple(hold,iold,jold,kold,X,startx,dx,oldimage,olddata,newdatatemp);
        // bilinear_interp(t,r,th,ph, hold, iold, jold, kold, oldimage,olddata,newdatatemp); // not working
      }
      else if(interptypetodo==2){
        plane_interp(t, r, th, ph, hold, iold, jold, kold, oldimage,olddata, newdatatemp);
      }
      else if(interptypetodo==3){
        bicubic_interp_wrap(oN0,oN1,oN2,oN3, hold,iold,jold,kold, X[0],X[1],X[2],X[3],oldimage,olddata, newdatatemp);
      }
      else{
        fprintf(stderr,"no such interptypetodo=%d\n",interptypetodo);
        exit(1);
      }

      // DEBUG
      if(DEBUGINTERP){
        fprintf(stderr,"Done interp h=%d i=%d j=%d k=%d\n",h,i,j,k); fflush(stderr);
      }

      // DEBUG
      if(DEBUGINTERP){
        fprintf(stderr,"Done set newdata(temp) or write to file the interpolation result\n"); fflush(stderr);
      }

      // test symmetry -- yes!
      //      if(j==nN2/3 && (i==nN1/4 || i==3*nN1/4-1)){
      //       fprintf(stderr,"i=%d j=%d data=%21.15g\n",i,j,newdatatemp[coli]);
      //      }


      // DEBUG
      //      if((r>95)&&(r<105)&&(fabs(th-M_PI*0.5)<0.1)){
      // fprintf(stderr,"%d %d : %g %g : %g %g : %d %d : %g\n",i,j,r,th,X[1],X[2],iold,jold,ftemp);
      //      }

      //      if(newdatatemp[coli]<=0.0){
      //       fprintf(stderr,"NEGATIVE FOUND interptypetodo=%d :: hold=%d iold=%d jold=%d kold=%d h=%d i=%d j=%d k=%d\n",interptypetodo,hold,iold,jold,kold,h,i,j,k); fflush(stderr);
      //      }

      // DEBUG
      if(DEBUGINTERP){
        fprintf(stderr,"%d %d %d %d: %g %g %g %g: %g %g %g %g: %d %d %d %d\n",h,i,j,k, t,r,th,ph, X[0],X[1],X[2],X[3], hold,iold,jold,kold);
        for(coli=0;coli<numoutputcols;coli++){ // over all independent columsn of data
          fprintf(stderr,"coli=%d newdatatemp=%g\n",coli,newdatatemp[coli]);
        }
        fflush(stderr);
      }


      // copy over newdatatemp to newimagetemp if doing DATATYPE==0
      if(DATATYPE==0){
        for(coli=0;coli<numoutputcols;coli++) newimagetemp[coli] = FLOAT2IMAGE(newdatatemp[coli]) ;
      }


    }// end if getinterp==1
    else{
      for(coli=0;coli<numoutputcols;coli++){ // over all independent columsn of data
        if(DATATYPE==0)      newimagetemp[coli] = FLOAT2IMAGE(defaultvalue[coli]) ;
        else newdatatemp[coli]=defaultvalue[coli];
      }
    }// end else


    ////////////////
    // OUTPUT or STORE result (either defaultvalue, interpoled value, or no value)
    int whichdatatype=DATATYPE; // 1 means using newdatatemp
    output2file_perpoint_postinterpolation(whichdatatype, h,i,j,k, newimagetemp, newdatatemp);


  }// end loop over new zones


  // free temp column space
  free(newdatatemp);
  free(newimagetemp);


  // report that done
  if(DATATYPE==0)   fprintf(stderr,"done interpolating image\n");
  else fprintf(stderr,"done interpolating data\n");
  fflush(stderr);


  return(0);

}





// set dtc,dxc,dyc,dyz only used by new_xycoord() and outputted to file as $dx1,2,3  ($dx0 not outputted)
void setup_newgrid(void)
{
  FTYPE ftemp;



  /////////////////////////////////
  //
  // new grid
  //
  /////////////////////////////////

  // assumes 1:1 aspect ratio and center of zoom is x=0,y=0
  // used by  new_xycoord() below

  if(newgridtype==GRIDTYPENOCHANGE){
    // here d?c, start?c, faked?c define internal (not real) coordinate system that is used (with bl_coord()) to get real coordinates

    //    dxc = (Rout-Rin)/(FTYPE)nN1 ;
    //    dyc = (M_PI-0)/(FTYPE)nN2 ;
    //    dzc = (2.0*M_PI-0)/(FTYPE)nN3 ;

    dtc=dX[0]*totalsize[0]/nN0;
    dxc=dX[1]*totalsize[1]/nN1;
    dyc=dX[2]*totalsize[2]/nN2;
    dzc=dX[3]*totalsize[3]/nN3;

    starttc=(startx[0]+starttdata); // starttdata added because set_points() forces startx[0]=0
    startxc=startx[1];
    startyc=startx[2];
    startzc=startx[3]; // was 0

    //      &t,&totalsize[1],&totalsize[2],&startx[1],&startx[2],&dX[1],&dX[2],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dt,&defcoord,&MBH,&QBH,&EP3,&THETAROT);

    fakedtc=dtc;
    fakedxc=dxc;
    fakedyc=dyc;
    fakedzc=dzc;

  }
  else if(newgridtype==GRIDTYPECART){
    dtc = (endtc-starttc)/(FTYPE)nN0 ;
    dxc = (endxc-startxc)/(FTYPE)nN1 ;
    dyc = (endyc-startyc)/(FTYPE)nN2 ;
    dzc = (endzc-startzc)/(FTYPE)nN3 ;

    fakedtc=dtc;
    fakedxc=dxc;
    fakedyc=dyc;
    fakedzc=dzc;
  }
  else if(newgridtype==GRIDTYPELOGLOGCYL){
    fakedtc= (endtc-starttc)/(FTYPE)nN0 ;
    fakedxc= (endxc-startxc)/(FTYPE)nN1 ;
    fakedyc= (endyc-startyc)/(FTYPE)nN2 ;
    fakedzc= (endzc-startzc)/(FTYPE)nN3 ;

    // don't do split log in time
    dtc=(endtc-starttc)/(FTYPE)nN0 ;

    // spatial coordinates
    // log or split log for each direction
    fakeRin=Rin;

    if(startxc>=0.0) dxc = log10(endxc/(fakeRin))/(FTYPE)nN1 ;
    else if((startxc<0.0)&&(endxc>0.0)){
      ftemp=MAX(-startxc,endxc);
      dxc = 2.0*log10(ftemp/(fakeRin))/(FTYPE)nN1 ;
    }
    else if(endxc<0.0){
      dxc = log10(-startxc/(fakeRin))/(FTYPE)nN1 ;
    }

    if(startyc>=0.0) dyc = log10(endyc/(fakeRin))/(FTYPE)nN2 ;
    else if((startyc<0.0)&&(endyc>0.0)){
      ftemp=MAX(-startyc,endyc);
      dyc = 2.0*log10(ftemp/(fakeRin))/(FTYPE)nN2 ;
    }
    else if(endyc<0.0){
      dyc = log10(-startyc/(fakeRin))/(FTYPE)nN2 ;
    }

    if(startzc>=0.0) dzc = log10(endzc/(fakeRin))/(FTYPE)nN3 ;
    else if((startzc<0.0)&&(endzc>0.0)){
      ftemp=MAX(-startzc,endzc);
      dzc = 2.0*log10(ftemp/(fakeRin))/(FTYPE)nN3 ;
    }
    else if(endzc<0.0){
      dzc = log10(-startzc/(fakeRin))/(FTYPE)nN3 ;
    }

    //    if(startyc>=0.0) dyc = log10(endyc/(fakeRin))/(FTYPE)nN2 ;
    //else dyc = 2.0*log10(endyc/(fakeRin))/(FTYPE)nN2 ;

    //    if(startzc>=0.0) dzc = log10(endzc/(fakeRin))/(FTYPE)nN3 ;
    //else dzc = 2.0*log10(endzc/(fakeRin))/(FTYPE)nN3 ;


  }
  else if(newgridtype==GRIDTYPELOGSPC){

    // GODMARK: Need to find coefficients for r(i)

    fakedtc= (endtc-starttc)/(FTYPE)nN0 ;
    fakedxc= (endxc-startxc)/(FTYPE)nN1 ;
    fakedyc= (endyc-startyc)/(FTYPE)nN2 ;
    fakedzc= (endzc-startzc)/(FTYPE)nN3 ;

    // normal in time
    dtc=fakedtc;

    // spatial coordinates
    // start and ending x',y',z'
    FTYPE startxp = log10(1.0+fabs(startxc))*sign(startxc);
    FTYPE endxp = log10(1.0+fabs(endxc))*sign(endxc);

    if(startxp<=endxp){
      dualfprintf(fail_file,"Cannot have startxp=%g = endxp=%g\n",startxp,endxp);
      myexit(348634);
    }

    FTYPE startyp = log10(1.0+fabs(startyc))*sign(startyc);
    FTYPE endyp = log10(1.0+fabs(endyc))*sign(endyc);

    if(startyp<=endyp){
      dualfprintf(fail_file,"Cannot have startyp=%g = endyp=%g\n",startyp,endyp);
      myexit(348634);
    }

    FTYPE startzp = log10(1.0+fabs(startzc))*sign(startzc);
    FTYPE endzp = log10(1.0+fabs(endzc))*sign(endzc);

    if(startzp<=endzp){
      dualfprintf(fail_file,"Cannot have startzp=%g = endzp=%g\n",startzp,endzp);
      myexit(348634);
    }

    // really dxc = dxp \equiv dx' and same for each x,y,z
    dxc = (endxp-startxp)/(FTYPE)nN1;
    dyc = (endyp-startyp)/(FTYPE)nN2;
    dzc = (endzp-startzp)/(FTYPE)nN3;

    // determine gridr0
    FTYPE parms[5];
    parms[1]=endxp;
    parms[2]=startxp;
    parms[3]=startxc;
    parms[4]=endxc;
    root_find_parameter1(parms, &gridr0global);
    gridAAglobal = startxp/log10(1.0+startxc/gridr0global);
    // only valid if startxc==startyc==startzc and endxc==endyc==endzc


  }
  else if(newgridtype==GRIDTYPECARTLIGHT){
    dtc = (endtc-starttc)/(FTYPE)nN0 ;  // assumed time range is linear over same range as original time range

    dxc = (endxc-startxc)/(FTYPE)nN1 ;
    dyc = (endyc-startyc)/(FTYPE)nN2 ;
    dzc = (endzc-startzc)/(FTYPE)nN3 ;

    fakedtc=dtc;

    fakedxc=dxc;
    fakedyc=dyc;
    fakedzc=dzc;
  }
  else{
    //  NOTE: newgridtype==GRIDTYPESPC not setup yet since haven't wanted to transform to SPC from some other coords
    dualfprintf(fail_file,"No such newgridtype=%d\n",newgridtype);
    myexit(24326264);
  }

  // DEBUG:
  //  N1=nN1;
  //N2=nN2;
  //N3=nN3;
}


// Get new real coordinates (often real Cartesian t,x,y,z) from new grid's h,i,j,k
static void new_xycoord(int h, int i, int j, int k, FTYPE *tc, FTYPE *xc, FTYPE *yc, FTYPE *zc)
{
  FTYPE tnew,xnew,ynew,znew;
  FTYPE newX0,newX1,newX2,newX3;
  FTYPE newX0u,newX0d;
  FTYPE newX1u,newX1d;
  FTYPE newX2u,newX2d;
  FTYPE newX3u,newX3d;
  FTYPE tt,xt,yt,zt;
  FTYPE X[NDIM], V[NDIM];

  /*
    xnew = (i+0.5)*dxc ;
    // *y = (j-nN2/2+0.5)*dyc ;
    ynew = (nN2/2-(j+0.5))*dyc ;

    *x= sqrt(xnew*xnew + ynew*ynew) ;
    *y = atan2(xnew,ynew) ; // deliberately interchanged args
    */

  if(newgridtype==GRIDTYPENOCHANGE){
    // GODMARK: Assumes input data at CENT

    X[0] = starttc+(i+0.5)*dtc;
    X[1] = startxc+(i+0.5)*dxc;
    X[2] = startyc+(j+0.5)*dyc ;
    X[3] = startzc+(k+0.5)*dzc ;
    
    // use bl_coord()
    interp_bl_coord(X, V);

    *tc = V[0];
    *xc = V[1];
    *yc = V[2];
    *zc = V[3];
  }
  else if(newgridtype==GRIDTYPELOGSPC){
    // for this grid, only outer radius corresponds to startxc/endxc, etc., while in between tc,xc,yc,zc are arbitrary
    *tc = starttc+(h+0.5)*dtc ;
    *xc = startxc+(i+0.5)*dxc ; // varies from startxc to endxc
    *yc = startyc+(j+0.5)*dyc ;
    *zc = startzc+(k+0.5)*dzc ;

    //    *tc = starttc+(endtc-starttc)*(h+0.5)/((FTYPE)nN0); // varies from starttc to endtc
    //    *xc = startxc+(endxc-startxc)*(i+0.5)/((FTYPE)nN1); // varies from startxc to endxc
    //    *yc = startyc+(endyc-startyc)*(j+0.5)/((FTYPE)nN2); // varies from startyc to endyc
    //    *zc = startzc+(endzc-startzc)*(k+0.5)/((FTYPE)nN3); // varies from startzc to endzc
  }
  else if(newgridtype==GRIDTYPECART){
    // GODMARK: Assumes input data at CENT
    *tc = starttc+(h+0.5)*dtc ;
    *xc = startxc+(i+0.5)*dxc ;
    *yc = startyc+(j+0.5)*dyc ;
    *zc = startzc+(k+0.5)*dzc ;
  }
  else if(newgridtype==GRIDTYPELOGLOGCYL){

    *tc = starttc+(h+0.5)*dtc ; // normal time coordinate

    if(startxc>=0.0){
      // assumes i=0 is x=0
      newX1=log10(fakeRin)+(i+0.5)*dxc;
      *xc= pow(10.0,newX1);
    }
    else if(endxc<=0.0){
      // assumes i=0 is x=0
      newX1=log10(fakeRin)+(i+0.5)*dxc;
      *xc= -pow(10.0,newX1);
    }
    else{ // assumes symmetric around x=0
      // assumes i=nN1/2 is x=0
      if(i>=nN1/2){
        newX1u=log10(fakeRin)+(i-nN1/2+0.5)*(dxc);
        *xc = pow(10.0,newX1u) ;
      }
      else if(i<nN1/2){
        newX1d=log10(fakeRin)+(nN1/2-i-0.5)*(dxc);
        *xc = -pow(10.0,newX1d) ;
      }
    }

    if(startyc>=0.0){
      // assumes i=0 is x=0
      newX2=log10(fakeRin)+(j+0.5)*dyc;
      *yc= pow(10.0,newX2);
    }
    else if(endyc<=0.0){
      newX2=log10(fakeRin)+(j+0.5)*dyc;
      *yc= -pow(10.0,newX2);
    }
    else{ // assumes symmetric around x=0
      // assumes j=nN2/2 is y=0
      if(j>=nN2/2){
        newX2u=log10(fakeRin)+(j-nN2/2+0.5)*(dyc);
        *yc = pow(10.0,newX2u) ;
      }
      else if(j<nN2/2){
        newX2d=log10(fakeRin)+(nN2/2-j-0.5)*(dyc);
        *yc = -pow(10.0,newX2d) ;
      }
    }
    if(startzc>=0.0){
      // assumes i=0 is x=0
      newX3=log10(fakeRin)+(k+0.5)*dzc;
      *zc= pow(10.0,newX3);
    }
    else if(endzc<=0.0){
      newX3=log10(fakeRin)+(k+0.5)*dzc;
      *zc= -pow(10.0,newX3);
    }
    else{ // assumes symmetric around x=0
      // assumes k=nN3/3 is z=0
      if(k>=nN3/2){
        newX3u=log10(fakeRin)+(k-nN3/2+0.5)*(dzc);
        *zc = pow(10.0,newX3u) ;
      }
      else if(k<nN3/2){
        newX3d=log10(fakeRin)+(nN3/2-k-0.5)*(dzc);
        *zc = -pow(10.0,newX3d) ;
      }
    }
  }
  else if(newgridtype==GRIDTYPECARTLIGHT){
    // use new grid's h,i,j,k to get new coordinates
    *tc = starttc+(h+0.5)*dtc ; // position of tobserver is uniformly spaced
    *xc = startxc+(i+0.5)*dxc ; // varies from startxc to endxc
    *yc = startyc+(j+0.5)*dyc ;
    *zc = startzc+(k+0.5)*dzc ;
  }
  else{
    //  NOTE: newgridtype==GRIDTYPESPC not setup yet since haven't wanted to transform to SPC from some other coords
    dualfprintf(fail_file,"No such newgridtype=%d\n",newgridtype);
    myexit(2467346);
  }

  // DEBUG:
  //  fprintf(stderr,"%d %d %g %g : %g %g %g : %g %g %g\n",i,j,*xc,*yc,newX1,log10(fakeRin),(i+0.5)*dxc, pow(newX1,10.0),pow(newX2u,10.0),-pow(newX2d,10.0)); fflush(stderr);
 
}


// Get original-type real coordinates (e.g. for spherical polar : new values of t,r,th,ph) from new grid's h,i,j,k (implicitly from tc,xc,yc,zc in Cartesian)
// This must really give those coordinates (here labelled t,r,th,ph) that correspond to original metric's coordinates that are obtained after transforming from PRIMECOORDS to normal coords using dxdxp matrix
// of course, input h,i,j,k or result from new_xycoord() can be tc,xc,yc,zc that are manipulated so that new h,i,j,k doesn't have to correspond to h,x,y,z but some other t',x',y',z'
static void new_coord(int h, int i,int j,int k, FTYPE *t, FTYPE *r,FTYPE *th,FTYPE *ph)
{
  FTYPE tc,xc,yc,zc;
  FTYPE Rc,rspc;



  // get new real coordiantes for new h,i,j,k
  new_xycoord(h,i,j,k,&tc,&xc,&yc,&zc);


  // for Cart:
  // tc is t
  // xc like R-cyl
  // yc like Z-cyl
  // zc like y-direction

  
  // spatial coordinates (again, t,r,th,ph should be whatever coordinates (post application of dxdxp[][]) the metric used via oldgridtype, where input tc,xc,yc,zc are new grid uniform labels


  if(newgridtype==GRIDTYPENOCHANGE){
    // then doesn't matter what input grid type is
    *t = tc;
    *r = xc;
    *th = yc;
    *ph = zc;
  }
  else if(newgridtype==GRIDTYPELOGSPC && oldgridtype==GRIDTYPESPC){
    // original grid must be SPC
 

    // here xc,yc,zc are related to r,th,ph via
    //
    // tc like t = t'
    // xc like R-cyl = x'
    // yc like Z-cyl = z'
    // zc like y-direction = y'
    //
    // i.e. tc = t' = t
    // i.e. xc = x' = log_{10}(1+r) sin\theta cos\phi
    // i.e. zc = y' = log_{10}(1+r) sin\theta sin\phi
    // i.e. yc = z' = log_{10}(1+r) cos\theta
    //
    // tan(\theta) =  x'/z' = xc/yc
    // r = -1 + 10^{\sqrt{x'^2 + y'^2 + z'^2}}
    // tan(\phi) = y'/x' = zc/xc
    //
    // assume startxc,etc. correspond to outermost radial or x,y,z coordinate in Cart coords

    // time coordinate
    *t = tc;

    // spatial coordinates:


    //    *r = -1.0 + pow(sqrt( xc*xc + yc*yc + zc*zc),10.0);
    FTYPE rfake = sqrt(xc*xc + yc*yc + zc*zc);
    *r = gridr0global*(pow(rfake/gridAAglobal,10.0)-1.0); // wrong unless range in startxc=startyc, etc.

    Rc = sqrt(xc*xc+zc*zc); // R' = \sqrt(x'^2 + y'^2)

    *th = atan2(Rc,yc) ; /* deliberately interchanged args */
    if(*th<0) *th+=M_PI;
    if(*th<0) *th+=M_PI;
    if(*th<0) *th+=M_PI;
    if(*th>M_PI) *th-=M_PI;
    if(*th>M_PI) *th-=M_PI;
    if(*th>M_PI) *th-=M_PI;

    *ph = atan2(zc,xc) ; // OLD: Note that since \hat{i}\times\hat{j} = -\hat{k}, then \hat{k}=-\hat{y} so we flip xc,zc = x',y' to get back normal r,\theta,\phi coordinates (or is this just consistent with how we setup vis5d+ files?)
    if(*ph<0) *ph+=2.0*M_PI;
    if(*ph<0) *ph+=2.0*M_PI;
    if(*ph<0) *ph+=2.0*M_PI;
    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
   

  }
  else if(newgridtype==GRIDTYPECART || newgridtype==GRIDTYPECARTLIGHT){// Cart output (i.e. tc,xc,yc,zc are Cartesian labels for new grid)


    //    if(newgridtype==GRIDTYPECARTLIGHT && ROTATECARTLIGHT==1){
    
    // now assume always have rotation
    // Rotate x,y,z before assigning spherical polar value
    // As here, xc,yc,zc are really new coordinates, so already have rotation.
    // Our job here is to get back the true Cartesian coordinates so the below can be computed as normal
    // recall: xc is Cart-x , yc is Cart-z, and zc is Cart-y
      
    // Keep things simple at first.
    // Since xc,yc,zc "already" rotated, then rotate back to get Cartesian
    // Back rotation means rotate -tnrdegrees degrees around (e.g.) x or y axis.  Let's stick to y-axis, so y-axis is unchanged.  This means zc is unchanged.

    // Needs to be consistent with generate_lambdacoord() such that formula between rotated and nonrotated coordinates is same
    // such that BACK-rotation with tnrdegrees=90deg takes : -xrot -> +znonrot  & +zrot -> +xnonrot  so BACK-rotation is from x-axis towards z-axis.  Or with y-axis pointed at you, the rotation is clockwise.

    FTYPE xcnew = +xc*cos(-tnrdegrees*M_PI/180.0) + +yc*sin(-tnrdegrees*M_PI/180.0); // x // then xcCartnonrot = -zcCartrot for tnrdegrees=90deg.
    FTYPE zcnew = -xc*sin(-tnrdegrees*M_PI/180.0) + +yc*cos(-tnrdegrees*M_PI/180.0); // z // then zcCartnonrot = xcCartrot for tnrdegrees=90deg.
    FTYPE ycnew  = zc; // y

    xc=xcnew;
    yc=ycnew;
    zc=zcnew;

    // now we have the back-rotation around the y-axis (i.e. zc)
      
    // now below use of xc,yc,zc will be correctly using true (unrotated) Cartesian coordinates.  This ends up so that new grid (rotated Cart) will show the rotation of the old grid data.
      
    //    }

    if(oldgridtype==GRIDTYPECART){ // Cart input
      *r=xc;
      *th=yc;
      *ph=zc;
      
      // get rspc
      rspc=sqrt(xc*xc+yc*yc+zc*zc);
    }
    else if(oldgridtype==GRIDTYPESPC){ // Spherical polar input (normal)

      *r = sqrt(xc*xc + yc*yc + zc*zc) ;
      rspc=*r;
      
      Rc=sqrt(xc*xc+zc*zc);
 
      // \theta = atan(R/z)
      *th = atan2(Rc,yc) ; /* deliberately interchanged args */
      if(*th<0) *th+=M_PI;
      if(*th<0) *th+=M_PI;
      if(*th<0) *th+=M_PI;
      if(*th>M_PI) *th-=M_PI;
      if(*th>M_PI) *th-=M_PI;
      if(*th>M_PI) *th-=M_PI;
      //     *th = atan(Rc/yc) ;
      //dualfprintf(fail_file," got here %d %d\n",newgridtype,oN3);


      // NOTE: \phi = atan(x/y) has x>0,y=0 plane as \phi=0, as in Griffiths.
      if(oN3==1){
        *ph = atan2(fabs(zc),fabs(xc));
        if(*ph<0) *ph+=M_PI;
        if(*ph>M_PI) *ph-=M_PI;
        //  dualfprintf(fail_file,"newcoord here ph=%g\n",*ph);
      }
      else{
        *ph = atan2(zc,xc) ;
        if(*ph<0) *ph+=2.0*M_PI;
        if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
      }

      
      if(defcoord==666){ // then override above
        *r=xc;
        *th=yc;
        *ph=zc;
        // get rspc
        rspc=sqrt(xc*xc+yc*yc+zc*zc); // GODMARK: Probably wrong, but no time to fix obscure 666 defcoord
      }
    }// end if oldgridtype==GRIDTYPESPC
    else if(oldgridtype==GRIDTYPECYL){ // Cylindrical coordinates input

      //     rspc = sqrt(xc*xc + yc*yc + zc*zc) ;
      Rc=sqrt(xc*xc+zc*zc);
     
      *r = Rc; // Cyl radius
      *th = yc; // Cyl height
      *ph = atan2(zc,xc) ; // Cyl angle
      if(*ph<0) *ph+=2.0*M_PI;
      if(*ph>2.0*M_PI) *ph-=2.0*M_PI;

      // get rspc
      rspc=sqrt(Rc*Rc+yc*yc);


    }
    else if(oldgridtype==GRIDTYPECARTLIGHT){ // CartLIGHT input

      if(ROTATECARTLIGHT==0){
        *r=xc;
        *th=yc;
        *ph=zc;
      }
      else{
        // recall: xc is Cart-x , yc is Cart-z, and zc is Cart-y
        fprintf(stderr,"NOT YET 295638667 oldgridtype=%d\n",oldgridtype);
        exit(1);
      }
      
      // get rspc
      rspc=sqrt(xc*xc+yc*yc+zc*zc);
    }
    else{
      dualfprintf(fail_file,"No such combination of newgridtype=%d and oldgridtype=%d\n",newgridtype,oldgridtype);
      myexit(249646);
    }



    // adjust time coordinate for the 2 different Cartesian newgrid types
    if(newgridtype==GRIDTYPECARTLIGHT){
      // then modify time coordinate (assume all other [old] coordinates use normal time)
      // t0 = t - nvec . rvec  = t - r*cos(theta between nhat and rhat)
      // tc : new grid label = t0
      // t : original time = t
      // assumes rspc computed above (note rspc doesn't change under the rotation if ROTATECARTLIGHT==1)
      *t = tc + rspc*cos(tnrdegrees*M_PI/180.0); // Note sign is + since this equation is for t = t0 + ... That is, *t is the original lab-frame time coordinate, which is one of the oldgridtype's above, and so normal lab-frame time
    }
    else if(oldgridtype==GRIDTYPECARTLIGHT){
      // assume all other grid use normal time, so specify tobs as *t such that tc is normal time
      *t = tc - rspc*cos(tnrdegrees*M_PI/180.0);
    }
    else{
      // standard time coordinate
      *t = tc;
    }

  }// end if newgrid is of Cartesian type
  else{
    dualfprintf(fail_file,"No such combination of newgridtype=%d and oldgridtype=%d\n",newgridtype,oldgridtype);
    myexit(978463);
  }




  return ;
}




// This function only used for distances for final interpolation, so not crucial how input and output related, but most correct gives best quality interpolation result
// notice that y and z are switched compared to in Griffiths book final page
// Cartesian t,x,y,z from spherical polar t,r,th,ph
// always the same for any coordinates unless setup t,x,y,z to not be Cartesian elsewhere also
// used for Cartesian distance, so really want Cartesian coordinates!
static void old_xyzcoord(FTYPE t, FTYPE r, FTYPE th, FTYPE ph, FTYPE *tc, FTYPE *xc, FTYPE*yc, FTYPE *zc)
{
  FTYPE rspc;

  if(oldgridtype==GRIDTYPECART){ // Cart input 2 Cart
    *tc = t;
    *xc = r;
    *yc = th;
    *zc = ph;
  }
  else if(oldgridtype==GRIDTYPECARTLIGHT){ // CartLIGHT input 2 Cart

    //    if(ROTATECARTLIGHT==0){
    //      *xc = r;
    //      *yc = th;
    //      *zc = ph;
    //    }
    //    else{
    ////NO: Assumes ROTATECARTLIGHT==1 when grid was created

    // won't change distance, but still do it
    // rotate newgrid back to Cart
    // xc: x
    // yc: z
    // zc: y
    // r : xold rotated
    // th : zold rotated
    // ph : yold rotated
    // should be same transformation as done in new_coord()
    *xc = +r*cos(-tnrdegrees*M_PI/180.0) + +th*sin(-tnrdegrees*M_PI/180.0);
    *zc = -r*sin(-tnrdegrees*M_PI/180.0) + +th*cos(-tnrdegrees*M_PI/180.0);
    *yc  = ph;
    //    }

    // get rspc
    rspc=sqrt((*xc)*(*xc) + (*yc)*(*yc) + (*zc)*(*zc));

    *tc = t + rspc*cos(tnrdegrees*M_PI/180.0);
  }
  else if(oldgridtype==GRIDTYPESPC){ // spc input 2 Cart
    //  if(defcoord!=666){


    
    if(newgridtype==GRIDTYPENOCHANGE){
      // test for fake distance
      *tc = t;
      *xc = r;
      *yc = th;
      *zc = ph;

    }
    else{// really used oldgridtype and going to new type:  for SPC 2 Cart


      if(newgridtype==GRIDTYPELOGSPC){
        // then fake Cartesian generated
        *tc = t;
        *xc = log(r+1.0)*sin(th)*cos(ph); // my x
        *zc = log(r+1.0)*cos(th); // my z
        *yc = log(r+1.0)*sin(th)*sin(ph); // my y
      }
      // below seems unnecessary and was even wrong before
      //      else if(newgridtype==GRIDTYPECART && oN3==1){
      // *tc = t;
      // *xc = 0; // not really equivalent to code before because before *xc was just unset, but unsure why that was or why doing this new version GODMARK
      // *zc = r*cos(th); // my z
      // *yc = r*sin(th)*sin(ph); // my y
      //      }
      else{// normal SPC 2 Cart
        *tc = t;
        *xc = fabs(r*sin(th)*cos(ph)); // my x
        *zc = r*cos(th); // my z
        *yc = r*sin(th)*sin(ph); // my y
      }
      

      if(newgridtype==GRIDTYPECARTLIGHT){ // override default tc
        // r is really spherical polar "r", so good to use directly
        *tc = t + r*cos(tnrdegrees*M_PI/180.0); // map back to true Cartesian (note sign is +, not -)

        // Note that if(ROTATECARTLIGHT==1), then we don't rotate or back rotate since we really want SPC -> True Cart here for distances!

      }

    }
  }
  else if(oldgridtype==GRIDTYPECYL){ // Cyl input 2 Cart
    *tc = t;
    *xc = r*sin(ph);
    *zc = r*cos(ph);
    *yc = th;
  }
  else if(oldgridtype==GRIDTYPECARTLIGHT){ // Cart input with light travel time t 2 Cart

    //    if(ROTATECARTLIGHT==0){
    //      *xc = r;
    //      *yc = th;
    //      *zc = ph;
    //    }
    //    else{

    // again, should be same transformation as in new_coord()
    *xc = +r*cos(-tnrdegrees*M_PI/180.0) + +th*sin(-tnrdegrees*M_PI/180.0);
    *zc = -r*sin(-tnrdegrees*M_PI/180.0) + +th*cos(-tnrdegrees*M_PI/180.0);
    *yc  = ph;

    //    }

    rspc = sqrt((*xc)*(*xc) + (*yc)*(*yc) + (*zc)*(*zc));
    *tc = t + rspc*cos(tnrdegrees*M_PI/180.0); // map back to true Cartesian (note sign is +, not -)
  }
  else{
    dualfprintf(fail_file,"No such oldgridtype=%d\n",oldgridtype);
    myexit(294763);
  }

}

// Cartesian distance between 2 spherical polar points
static FTYPE new_distance(FTYPE t1, FTYPE r1, FTYPE th1, FTYPE ph1, FTYPE t2, FTYPE r2, FTYPE th2, FTYPE ph2)
{
  FTYPE tc1,xc1,yc1,zc1,tc2,xc2,yc2,zc2;

  old_xyzcoord(t1,r1,th1,ph1, &tc1,&xc1,&yc1,&zc1);
  old_xyzcoord(t2,r2,th2,ph2, &tc2,&xc2,&yc2,&zc2);

  return(old_distance(tc1,xc1,yc1,zc1, tc2,xc2,yc2,zc2)); // can use same Cartesian function
  //return(sqrt((tc1-tc2)*(tc1-tc2)+(xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2)+(zc1-zc2)*(zc1-zc2)));
} 

// Cartesian distance between 2 Cartesian points
static FTYPE old_distance(FTYPE tc1, FTYPE xc1, FTYPE yc1, FTYPE zc1,  FTYPE tc2, FTYPE xc2, FTYPE yc2, FTYPE zc2)
{
  return(sqrt((tc1-tc2)*(tc1-tc2)+(xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2)+(zc1-zc2)*(zc1-zc2)));
} 

// Get t,r,th,ph from old grid's h,i,j,k
static void old_ijk2rthph(int hold, int iold, int jold, int kold, FTYPE *t, FTYPE*r, FTYPE*th, FTYPE *ph)
{
  FTYPE X[NDIM];

  old_ijk2x123(hold,iold,jold,kold,X);
  old_x1232rthphcoord(X,t,r,th,ph);

}

// Get t,r,th,ph from old grid's h,i,j,k
static void oldf_ijk2rthph(FTYPE hold, FTYPE iold, FTYPE jold, FTYPE kold,  FTYPE *t, FTYPE*r, FTYPE*th, FTYPE *ph)
{
  FTYPE X[NDIM];


  oldf_ijk2x123(hold,iold,jold,kold,X);
  old_x1232rthphcoord(X,t,r,th,ph);

}

// get t,r,th,ph for given x0,x1,x2,x3 for old grid
static void old_x1232rthphcoord(FTYPE *X,  FTYPE *t, FTYPE *r, FTYPE *th, FTYPE *ph)
{
  FTYPE V[NDIM];

  interp_bl_coord(X, V);

  *t=V[0];
  *r=V[1];
  *th=V[2];
  *ph=V[3];

}


// local version of bl_coord()
void interp_bl_coord(FTYPE *X, FTYPE *V)
{

  // this is only place normally need bl_coord() for all interpolation.  This gets old grid's V using new grid's determination of X so that iterative procedure can continue making errors smaller between input r,t,th,ph and V.

  // for time, V[0]=X[0] is set in bl_coord(), but set_points() forces startx[0]=0.  In order to avoid changing bl_coord() or set_points(), we change how X[0] is assinged in the interp code (see below old_ijk2x123(), etc.)
  // Then X[0] = told = datastartt + (hold + 0.5) datadt; where datadt = (dataendt-datastartt)/oN0;
  bl_coord(X,V);

  V[0] += starttdata;  // starttdata was added because set_points() assumes startx[0]=0 and bl_coord() assumes X[0]=V[0]
  // Note: since this is a constant offset, this doesn't affect derivatives part of root finding procedure.  Only affects what is meant by the solution in terms of an offset in time.

}

// get X0, X1, X2, X3 from h,i,j,k for old uniform X-grid
static void old_ijk2x123(int hold, int iold, int jold, int kold, FTYPE*X)
{
  coord(iold,jold,kold,CENT,X);

  X[0] = (startx[0]+ starttdata) + (hold + startpos[0] + 0.5) * dx[0] ; // starttdata was added because set_points() assumes startx[0]=0 and bl_coord() assumes X[0]=V[0] and can't adjust.  So X[0] must return back origgrid time.
  
}

// get X0, X1, X2, X3 from h,i,j,k for old uniform X-grid
static void oldf_ijk2x123(FTYPE hold, FTYPE iold, FTYPE jold, FTYPE kold, FTYPE*X)
{
  coordf(iold,jold,kold,CENT,X);

  X[0] = (startx[0]+ starttdata) + (hold + startpos[0] + 0.5) * dx[0] ; // starttdata was added because set_points() assumes startx[0]=0 and bl_coord() assumes X[0]=V[0] and can't adjust.  So X[0] must return back origgrid time.

}

// get h,i,j,k from X0, X1, X2, X3 for old uniform X-grid
static void old_x1232ijk(FTYPE *X, int *hold, int *iold, int *jold, int *kold)
{
  lowericoord(X,CENT,hold,iold,jold,kold);
}

// find old grid's X for given r,th,ph
static int old_coord(FTYPE t, FTYPE r, FTYPE th, FTYPE ph, FTYPE *X)
{

#if(ROOTMETHOD<=1)
  return(root_find1(t,r,th,ph,X)) ;
#else
  return(root_find2(t,r,th,ph,X)) ;
#endif

  fprintf(stderr,"shouldn't reach end of old_coord()\n");
  myexit(1);
  return(-1) ;
}










// finds X for a given t,r,theta,phi
static int root_find1(FTYPE t, FTYPE r, FTYPE th, FTYPE ph, FTYPE *X)
{
  int jj;
  int ntrial,mintrial;
  FTYPE tolx,tolf,tolxallowed,tolfallowed,tolxreport,tolfreport;
  FTYPE Xtrial[1+NDIM]; // to be used internally by newt() or broydn()
  int notfoundsolution;
  int ni; // newt type index

  // target t,r,th,ph
  ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
  ni++; spc_target[ni]=t;
  ni++; spc_target[ni]=r;
  ni++; spc_target[ni]=th;
  ni++; spc_target[ni]=ph;


  // mnewt parameters
  ntrial = 5;
  mintrial=5;
  tolx = 1.e-9;
  tolf = 1.e-9;
  tolxallowed=tolfallowed=1.e-6;
  tolxreport=1.e3*tolx;
  tolfreport=1.e3*tolf;

#if(GOODGUESS==0)
  coord(oN1/2-1,oN2/2-1,0,CENT,X);
#else
  // setup initial guess (can't just choose middle point, upper Pi/4's are lost
  if((th>M_PI*0.25)&&(th<0.75*M_PI)){
    coord(oN1/2-1,oN2/2-1,oN3/2-1,CENT,X);
  }
  else if(th<=M_PI*0.25){
    coord(oN1/2-1,oN2/4-1,oN3/2-1,CENT,X);
  }
  else if(th>=0.75*M_PI){
    coord(oN1/2-1,3*oN2/4-1,oN3/2-1,CENT,X);
  }
#endif

  // using different X because different ranks
  ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
  ni++; Xtrial[ni]=X[0];
  ni++; Xtrial[ni]=X[1];
  ni++; Xtrial[ni]=X[2];
  ni++; Xtrial[ni]=X[3];
    
  // iterate to find solution
  nrerrorflag=0;
  notfoundsolution=mnewt(ntrial, mintrial, Xtrial, 4, tolx, tolf, tolxallowed, tolfallowed, tolxreport, tolfreport, NULL, usrfun_joninterp);
  if(notfoundsolution) DLOOPA(jj) fprintf(stderr,"didn't find solution X[%d]=%21.15g\n",jj,X[jj]);
  if(nrerrorflag) DLOOPA(jj) fprintf(stderr,"nrerror: didn't find solution X[%d]=%21.15g\n",jj,X[jj]);
  if(notfoundsolution||nrerrorflag){
    return(1);
  }
  else{
    // assign answer
    ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
    ni++; X[0]=Xtrial[ni];
    ni++; X[1]=Xtrial[ni];
    ni++; X[2]=Xtrial[ni];
    ni++; X[3]=Xtrial[ni];

    return(0);
  }


}


#define NUMATTEMPTS 100

// finds X for a given t,r,theta,phi
static int root_find2(FTYPE t, FTYPE r, FTYPE th, FTYPE ph, FTYPE *X)
{
  int jj;
  int ntrial,mintrial;
  FTYPE tolx,tolf,tolxallowed,tolfallowed,tolxreport,tolfreport;
  FTYPE Xtrial[1+NDIM]; // to be used internally by newt() or broydn()
  int notfoundsolution;
  int attempt;
  int ni;



  // target t,r,th,ph
  ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
  ni++; spc_target[ni]=t;
  ni++; spc_target[ni]=r;
  ni++; spc_target[ni]=th;
  ni++; spc_target[ni]=ph;

#if(GOODGUESS==0)
  coord(oN1/2-1,oN2/2-1,0,CENT,X);
#else
  // setup initial guess (can't just choose middle point, upper Pi/4's are lost)
  if((th>M_PI*0.25)&&(th<0.75*M_PI)){
    coord(oN1/2-1,oN2/2-1,oN3/2-1,CENT,X);
  }
  else if(th<=M_PI*0.25){
    coord(oN1/2-1,0,oN3/2-1,CENT,X);
  }
  else if(th>=0.75*M_PI){
    coord(oN1/2-1,7*oN2/8-1,oN3/2-1,CENT,X);
  }
#endif




  // assume didn't find
  notfoundsolution=1;
  for(attempt=0;attempt<NUMATTEMPTS;attempt++){

    // using different X because different ranks
    ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
    ni++; Xtrial[ni]=X[0];
    ni++; Xtrial[ni]=X[1];
    ni++; Xtrial[ni]=X[2];
    ni++; Xtrial[ni]=X[3];
    
    // iterate to find solution (only need X2spc, newt computes Jacobian)
    //          newt(ANALYTICJAC,root guess, num dimension, &check,vecfun,usrfun)
    nrerrorflag=0;

    // DEBUG:
    //    for(jj=1;jj<=4;jj++) fprintf(stderr,"Xtrial[%d]=%21.15g\n",jj,Xtrial[jj]);




#if(ROOTMETHOD==2)
    newt  (ANALYTICJAC,NULL,Xtrial, 4, &notfoundsolution,X2spc,usrfun_joninterp2);
#elif(ROOTMETHOD==3)
    broydn(ANALYTICJAC,NULL,Xtrial, 4, &notfoundsolution,X2spc,usrfun_joninterp2);
#endif


    // DEBUG:
    //    fprintf(stderr,"notfoundsolution=%d\n",notfoundsolution); fflush(stderr);


    if(notfoundsolution==2){
      // then fjac bad -- probably out of bounds off grid for grid defined with singularity on interpolated grid.
      fprintf(stderr,"fjac singularity found -- treating as off grid and no solution\n");
      for(jj=1;jj<=4;jj++) fprintf(stderr,"fjab used: Xtrial[%d]=%21.15g\n",jj,Xtrial[jj]);
      return(2);
    }


    // quick message, and try again
    if(SIMPLEDEBUGINTERP){
      if(notfoundsolution) DLOOPA(jj) fprintf(stderr,"ns  X[%d]=%21.15g\n",jj,X[jj]);
      if(nrerrorflag) DLOOPA(jj) fprintf(stderr,"nr  X[%d]=%21.15g\n",jj,X[jj]);
    }

    if((notfoundsolution==0)&&(nrerrorflag==0)) break; // otherwise try until NUMATTEMPTS attempts have been tried
    else{
      // try a new starting position, which is ultimately why failed typically
      coord((int)((oN1/2-1)*ranc(0,0)),(int)((oN2/2-1)*ranc(0,0)),(int)((oN3/2-1)*ranc(0,0)),CENT,X);
    }
  }


  // DEBUG:
  //  fprintf(stderr,"attempt=%d : %d %d",attempt,notfoundsolution,nrerrorflag); fflush(stderr);



  if(notfoundsolution||nrerrorflag){
    if(1||SIMPLEDEBUGINTERP){
      fprintf(stderr,"Couldn't find solution after %d attempts\n",attempt); fflush(stderr);
    }
    return(1);
  }
  else{

    // assign answer
    ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
    ni++; X[0]=Xtrial[ni];
    ni++; X[1]=Xtrial[ni];
    ni++; X[2]=Xtrial[ni];
    ni++; X[3]=Xtrial[ni];

    return(0);
  }

}

// get spc coordinates (t,r,th,ph) minus target (t,r,th,ph) for a given X
static void X2spc(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff)
{
  FTYPE X[NDIM],spc_curr[NDIM];
  int i;
  int ni;


  ni=0; // redo Xguess so in normal HARM format
  ni++; X[0]=Xguess[ni];
  ni++; X[1]=Xguess[ni];
  ni++; X[2]=Xguess[ni];
  ni++; X[3]=Xguess[ni];


  for(i=0;i<=1;i++){
    if(i==1){
      // this avoids out of bounds attempts by iterative method
      if(X[0]-startx[0] < 0) X[0]=startx[0];
      if(X[0]-Xmax[0] >= 0) X[0]=Xmax[0];

      if(X[1]-startx[1] < 0) X[1]=startx[1];
      if(X[1]-Xmax[1] >= 0) X[1]=Xmax[1];

      if(X[2]-startx[2] < 0) X[2]=startx[2];
      if(X[2]-Xmax[2] >= 0) X[2]=Xmax[2];

      if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
        if(X[3]-startx[3] < 0) X[3]=X[3]+(Xmax[3]-startx[3]);
        if(X[3]-Xmax[3]  >= 0) X[3]=X[3]-(Xmax[3]-startx[3]);
      }
      else{
        // then normal non-periodic
        if(X[3]-startx[3] < 0) X[3]=startx[3];
        if(X[3]-Xmax[3] >= 0) X[3]=Xmax[3];
      }
    }


    // find current r,th,ph from current X
    old_x1232rthphcoord(X,&spc_curr[0],&spc_curr[1],&spc_curr[2],&spc_curr[3]);


    // spc_diff should contain function to be zeroed
    ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
    ni++; spc_diff[ni]=spc_curr[0]-spc_target[ni];
    ni++; spc_diff[ni]=spc_curr[1]-spc_target[ni];
    ni++; spc_diff[ni]=spc_curr[2]-spc_target[ni];
    ni++; spc_diff[ni]=spc_curr[3]-spc_target[ni];

    if((spc_diff[1]>BIG)||(spc_diff[2]>BIG)||(spc_diff[3]>BIG)||(spc_diff[4]>BIG)){
      // catches nans and infs
    }
    else break;
  }
}




// input parms, get returnparm
// finds X for a given r,theta
static int root_find_parameter1(FTYPE *parms, FTYPE *returnvalue)
{
  int notfoundsolution;


  // 0 means use numerical Jacobian, input dummy usrfun function that won't be used
  // parms is list of parms of FTYPE starting with parms[1] through whatever parms expected for resid_parameter1()
  // &returnvalue[-1] is sent so that accessing [1] element is actually returnvalue[0]
  newt(0, parms, &returnvalue[-1], 1, &notfoundsolution,resid_parameter1,usrfun_parameter1);

  return(notfoundsolution);
}
 
 
// like X2spc() for coordinates, but for a single parameter
static void resid_parameter1(int n, FTYPE *parms, FTYPE *guess, FTYPE *diff)
{

  FTYPE r0=guess[1];

  FTYPE xout=parms[1];
  FTYPE xin=parms[2];
  FTYPE end=parms[3];
  FTYPE start=parms[4];
  // see logradial_interpolation_interp.nb
  diff[1]=-xout + xin*log10(1.0+end/r0)/log10(1.0+start/r0);
  
}


// for single parameter
static int usrfun_parameter1(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
{

  // not really used

  return(0);

}


/* auxiliary function required by newt */
// used if have ANALYTICJAC 1 in global.joninterp.h
static int usrfun_joninterp2(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
{
  int i = 0, j = 0, k = 0;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE X[NDIM],spc_curr[NDIM];
  FTYPE V[NDIM];
  int ni;



  // X and Xguess different rank
  ni=0; // redo Xguess so in normal HARM format
  ni++; X[0]=Xguess[ni];
  ni++; X[1]=Xguess[ni];
  ni++; X[2]=Xguess[ni];
  ni++; X[3]=Xguess[ni];


  if(
     (oN0!=1)&&((X[0]-startx[0] < 0) || (X[0]-Xmax[0] >= 0)) ||
     (oN1!=1)&&((X[1]-startx[1] < 0) || (X[1]-Xmax[1] >= 0)) ||
     (oN2!=1)&&((X[2]-startx[2] < 0) || (X[2]-Xmax[2] >= 0)) ||
     (oN3!=1)&&((X[3]-startx[3] < 0) || (X[3]-Xmax[3] >= 0))
     ){
    // then beyond bounds
  }




  // find current r,th from current X
  old_x1232rthphcoord(X,&spc_curr[0],&spc_curr[1],&spc_curr[2],&spc_curr[3]);

  // spc_diff should contain function to be zeroed
  ni=0; // newt() needs things to start at 1 not 0, so set ni=0 so next ni++ iterates ni to 1
  ni++; spc_diff[ni]=spc_curr[0]-spc_target[ni];
  ni++; spc_diff[ni]=spc_curr[1]-spc_target[ni];
  ni++; spc_diff[ni]=spc_curr[2]-spc_target[ni];
  ni++; spc_diff[ni]=spc_curr[3]-spc_target[ni];


  // DEBUG
  //  for(i=1;i<=4;i++) fprintf(stderr,"%d: Xguess=%21.15g spc_diff=%21.15g spc_target=%21.15g\n",i,Xguess[i],spc_diff[i],spc_target[i]);

  //////////////////
  //
  // calculate dxdxp
  V[0]=spc_curr[0];
  V[1]=spc_curr[1];
  V[2]=spc_curr[2];
  V[3]=spc_curr[3];

  // DEBUG:
  //  for(i=0;i<4;i++) fprintf(stderr,"V[%d]=%21.15g\n",i,V[i]);

  // Using dxdxprim() (which may or may not have analytical dxdxp[][]'s)
  // This requires X be defined as in HARM (and so only uniform coordinate that maps to V)
  // And requires V to be as in HARM (and so oldgridtype, often spherical polar coordinates)
  // These requirements are generally satisfied when using HARM-produced data.
  // This code has nothing directly to do with newgridtype, which was only used to produce an X that we are now trying to get a V from to know what input data to use to interpolate to the new grid.
  dxdxprim(X,V,dxdxp);
  // assign to alpha (didn't use alpha directly since rank of alpha is smaller than dxdxp)
  for (j = 0; j < n; j++){
    for (k = 0; k < n; k++){
      alpha[j+1][k+1]=dxdxp[j][k];
    }
  }

  
  return (0);
}
 
 
#define NORMMETHOD 1
// 0: no norm method
// 1: special norm method
 
/* auxiliary function required by mnewt */
static int usrfun_joninterp(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *beta, FTYPE **alpha, FTYPE*norm)
{
  int h = 0, i = 0, j = 0, k = 0;
  FTYPE t,r,th,ph,spc_curr[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE X[NDIM];
  FTYPE V[NDIM];
  int numnormterms;
  int ni;


  // X and Xguess different rank
  ni=0; // redo Xguess so in normal HARM format
  ni++; X[0]=Xguess[ni];
  ni++; X[1]=Xguess[ni];
  ni++; X[2]=Xguess[ni];
  ni++; X[3]=Xguess[ni];


  // find current r,th from current X
  old_x1232rthphcoord(X,&t,&r,&th,&ph);

  V[0]=spc_curr[0]=t;
  V[1]=spc_curr[1]=r;
  V[2]=spc_curr[2]=th;
  V[3]=spc_curr[3]=ph;


  //////////////////
  //
  // calculate dxdxp
  dxdxprim(X,V,dxdxp);
  // assign to alpha (didn't use alpha directly since rank of alpha is smaller than dxdxp)
  for (j = 0; j < n; j++){
    for (k = 0; k < n; k++){
      alpha[j+1][k+1]=dxdxp[j][k];
    }
  }
  

  /////////////////////
  // normalization
  //
#if(NORMMETHOD==1)
  *norm=0.0;
  numnormterms=0;
  for (j = 0; j < n; j++){
    for (k = 0; k < n; k++){
      if(fabs(alpha[j + 1][k + 1])>NUMEPSILON){
        *norm+=fabs(alpha[j + 1][k + 1]);
        numnormterms++;
      }
    }
  }
  *norm=(FTYPE)(numnormterms)/(*norm); // (i.e. inverse of average of absolute values)
#else
  *norm=X[1];
  //*norm=1.0;
#endif
  // apply normalization
  for (j = 0; j < n; j++)
    for (k = 0; k < n; k++)
      alpha[j + 1][k + 1] *= (*norm);
  // determine normalized error
  for (k = 0; k < n; k++)
    beta[k + 1] = (spc_curr[k] - spc_target[k+1]) *(*norm);
  
  
  return (0);
}






//////////////////////////////
//
//////////////////// NEXT SEVERAL FUNCTIONS RELATED TO INTERPOLATION, SUPERSAMPLING, OR FILTERING
//
///////////////////////////////

//#define NUMIJK (2*2) // 2D
//#define NUMIJK (2*2*2) // 3D
#define NUMIJK (2*2*2*2) // 4D

#define NUMIJKMEM (NUMIJK+1)

// rref and thref are true grid location's r and theta
// x1,x2,x3 are 3 points x position
// y1,y2,y3 are 3 points y position
// z1,z2,z3 are 3 points values
// quad_interp is fully 3D capable
static int quad_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref, int hold, int iold, int jold, int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE*dist, FTYPE *mytc, FTYPE*myxc, FTYPE*myyc, FTYPE*myzc, FTYPE(*myfun)[NUMIJKMEM],  FTYPE *tcref,FTYPE*xcref,FTYPE*ycref,FTYPE *zcref)
{
  FTYPE t1,t2,t3, x1,x2,x3,  y1,y2,y3,  z1,z2,z3;
  int p;
  FTYPE myt[NUMIJKMEM],myr[NUMIJKMEM],myth[NUMIJKMEM],myph[NUMIJKMEM];
  //
  //
  int holdp,ioldp,joldp,koldp;
  int hsel[NUMIJKMEM],isel[NUMIJKMEM],jsel[NUMIJKMEM],ksel[NUMIJKMEM];
  int which,whichc;
  FTYPE interpolate;
  int atboundary;
  FTYPE ftemp;



  holdp = hold+1 ;
  ioldp = iold+1 ;
  joldp = jold+1 ;
  koldp = kold+1 ;

  //  atboundary=0;
  // if refining, then oN1 is really old image size, not refined, which is correct

  

  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    if((hold>=oN0)||(hold<0)||(iold>=oN1)||(iold<0)||(jold>=oN2)||(jold<0)||(kold<0)){
      fprintf(stderr,"shouldn't reach here: icoord() should truncate hold/iold/jold/kold\n");
      myexit(1);
    }
  }
  else{
    if((hold>=oN0)||(hold<0)||(iold>=oN1)||(iold<0)||(jold>=oN2)||(jold<0)||(kold>=oN3)||(kold<0)){
      fprintf(stderr,"shouldn't reach here: icoord() should truncate hold/iold/jold/kold\n");
      myexit(1);
    }
  }

  if((holdp<0)||(ioldp<0)||(joldp<0)||(koldp<0)){
    fprintf(stderr,"shouldn't reach here: holdp/ioldp/joldp/koldp can't be less than truncated hold/iold/jold/kold which are at minimum=0\n");
    myexit(1);
  }


  if(BOUNDARYEXTRAP==0){
    // only case that has to be corrected here
    if(holdp>=oN0){ holdp=oN0-1; hold=holdp-1;}
    if(ioldp>=oN1){ ioldp=oN1-1; iold=ioldp-1;}
    if(joldp>=oN2){ joldp=oN2-1; jold=joldp-1;}
  }

  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    // koldp can be up to =oN3
    if(koldp>=oN3+1){
      koldp=koldp-oN3; kold=koldp-1;
      phref=phref-(Xmax[3]-startx[3]);
    }

  }
  else{
    if(BOUNDARYEXTRAP==0){
      if(koldp>=oN3){ koldp=oN3-1; kold=koldp-1;}
    }
  }


  // otherwise do good interpolation
  
  // get distances and function values
  if(NUMIJK>=4){ // these only used in 2D
    hsel[1]=hold; isel[1]=iold;  jsel[1]=jold;  ksel[1]=kold;
    hsel[2]=hold; isel[2]=ioldp; jsel[2]=jold;  ksel[2]=kold;
    hsel[3]=hold; isel[3]=iold;  jsel[3]=joldp; ksel[3]=kold;
    hsel[4]=hold; isel[4]=ioldp; jsel[4]=joldp; ksel[4]=kold;
  }

  if(NUMIJK>=8){ // these only used in 3D
    hsel[5]=hold; isel[5]=iold;  jsel[5]=jold;  ksel[5]=koldp;
    hsel[6]=hold; isel[6]=ioldp; jsel[6]=jold;  ksel[6]=koldp;
    hsel[7]=hold; isel[7]=iold;  jsel[7]=joldp; ksel[7]=koldp;
    hsel[8]=hold; isel[8]=ioldp; jsel[8]=joldp; ksel[8]=koldp;
  }

  if(NUMIJK>=16){ // these only used in 4D
    hsel[9]=holdp;  isel[9]=iold;   jsel[9]=jold;   ksel[9]=kold;
    hsel[10]=holdp; isel[10]=ioldp; jsel[10]=jold;  ksel[10]=kold;
    hsel[11]=holdp; isel[11]=iold;  jsel[11]=joldp; ksel[11]=kold;
    hsel[12]=holdp; isel[12]=ioldp; jsel[12]=joldp; ksel[12]=kold;

    hsel[13]=holdp; isel[13]=iold;  jsel[13]=jold;  ksel[13]=koldp;
    hsel[14]=holdp; isel[14]=ioldp; jsel[14]=jold;  ksel[14]=koldp;
    hsel[15]=holdp; isel[15]=iold;  jsel[15]=joldp; ksel[15]=koldp;
    hsel[16]=holdp; isel[16]=ioldp; jsel[16]=joldp; ksel[16]=koldp;
  }

  old_xyzcoord(tref,rref,thref,phref,tcref,xcref,ycref,zcref);

  if(DEBUGINTERP)   fprintf(stderr,"tref=%g rref=%g thref=%g phref=%g  tcref=%g xcref=%g ycref=%g zcref=%g\n",tref,rref,thref,phref,*tcref,*xcref,*ycref,*zcref);

  int coli;
  for(p=1;p<=NUMIJK;p++){
    old_ijk2rthph(hsel[p],isel[p],jsel[p],ksel[p], &myt[p],&myr[p],&myth[p],&myph[p]) ;
    old_xyzcoord(myt[p],myr[p],myth[p],myph[p], &mytc[p],&myxc[p],&myyc[p],&myzc[p]);
    dist[p]=old_distance(mytc[p],myxc[p],myyc[p],myzc[p],*tcref,*xcref,*ycref,*zcref);
    if(DATATYPE==0){
      for(coli=0;coli<numoutputcols;coli++) myfun[coli][p]=(FTYPE)((unsigned char)oldimage[coli][hsel[p]][isel[p]][jsel[p]][ksel[p]]);
      if(DEBUGINTERP)       fprintf(stderr,"p=%d hsel=%d isel=%d jsel=%d ksel=%d :: myt=%g myr=%g myth=%g myph=%g :: mytc=%g myxc=%g myyc=%g myzc=%g :: fun=%g\n",p,hsel[p],isel[p],jsel[p],ksel[p], myt[p],myr[p],myth[p],myph[p], mytc[p],myxc[p],myyc[p],myzc[p], myfun[0][p]); fflush(stderr); // COLIMARK: only showing coli=0 for now
    }
    else{
      for(coli=0;coli<numoutputcols;coli++) myfun[coli][p]=olddata[coli][hsel[p]][isel[p]][jsel[p]][ksel[p]];
    }
  }

  // feeds back dist, mytc, myxc, myyc, myzc,  myfun[][]
  return(0);

}



// interpolate on an (N-1)D plane
// tref, rref, thref, phref are true grid location's t, r, theta, phi
// t1,t2,t3 are 3 points t position
// x1,x2,x3 are 3 points x position
// y1,y2,y3 are 3 points y position
// z1,z2,z3 are 3 points z position
// f1,f2,f3 are 3 points function values
// plane_interp is NOT fully 3D capable because assumes only using x-y space and so essentially nearest neighbor in k-space
static int plane_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref, int hold, int iold, int jold, int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp)
{
  FTYPE t1,t2,t3, x1,x2,x3, y1,y2,y3, z1,z2,z3;
  FTYPE f1,f2,f3;
  FTYPE a3,b3,c3,aoc,boc,newvalue;
  int p;
  FTYPE myt[NUMIJKMEM],myr[NUMIJKMEM],myth[NUMIJKMEM],myph[NUMIJKMEM],myfun[MAXCOLS][NUMIJKMEM],dist[NUMIJKMEM],mytc[NUMIJKMEM],myxc[NUMIJKMEM],myyc[NUMIJKMEM],myzc[NUMIJKMEM];
  FTYPE tcref,xcref,ycref,zcref;
  int quad_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref, int hold, int iold, int jold, int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE*dist, FTYPE *mytc, FTYPE*myxc, FTYPE*myyc, FTYPE*myzc, FTYPE(*myfun)[NUMIJKMEM], FTYPE *tcref, FTYPE*xcref,FTYPE*ycref,FTYPE *zcref);

  int nearest_interp(int hold, int iold,int jold,int kold,unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
  int which,whichc;
  int a,b,c;
  FTYPE interpolate;
  FTYPE ftemp;


  if(quad_interp(tref, rref, thref, phref,  hold, iold, jold, kold, oldimage,olddata, dist,  mytc, myxc, myyc, myzc, myfun, &tcref, &xcref,&ycref,&zcref)>=1){
    return(nearest_interp(hold,iold,jold,kold,oldimage,olddata, newdatatemp));
  }

  //////////////////////
  //
  // rest is specific to plane interpolation
  //
  // rest doesn't care about hold/iold/jold/kold, just actual positions and function value
  //


  which=1;
  for(p=1;p<=4;p++){ // GODMARK3D
    if(dist[p]>dist[which]) which=p;
  }
  // which now contains the point furthest away from our i,j location
  // as opposed to choosing any 3 of the 4
  whichc=1;
  for(p=1;p<=4;p++){ // GODMARK3D
    if(dist[p]<dist[whichc]) whichc=p;
  }
  // whichc now contains the closest point
  
  if(DEBUGINTERP)   fprintf(stderr,"hold=%d iold=%d jold=%d kold=%d :: which=%d whichc=%d\n",hold,iold,jold,kold, which,whichc);


  // GODMARK3D
  if(which==1){         a=2; b=3; c=4;}
  else if(which==2){    a=1; b=3; c=4;}
  else if(which==3){    a=1; b=2; c=4;}
  else if(which==4){    a=1; b=2; c=3;}

  //a=whichc;b=2;c=3;
  //  a=1;b=2;c=3;


  int coli;
  for(coli=0;coli<numoutputcols;coli++){ // over columns of data

    // GODMARK3D
    x1=myxc[a]; y1=myyc[a]; f1=myfun[coli][a];
    x2=myxc[b]; y2=myyc[b]; f2=myfun[coli][b];
    x3=myxc[c]; y3=myyc[c]; f3=myfun[coli][c];
    
    if(1){
      // GODMARK3D
      a3=y2* f1 - y3* f1 - y1* f2 + y3* f2 + y1* f3 - y2* f3;
      b3=-(x2* f1 - x3* f1 - x1* f2 + x3* f2 + x1* f3 - x2* f3);
      c3=x2* y1 - x3* y1 - x1* y2 + x3* y2 + x1* y3 - x2* y3;
      aoc=a3/c3;
      boc=b3/c3;
      
      if(DEBUGINTERP)     fprintf(stderr,"%d %d %d %d :: %g %g %g %g %g\n",hold,iold,jold,kold,a3,b3,c3,aoc,boc); fflush(stderr);
      
      newvalue=f1;
      interpolate=-aoc*(xcref-x1)-boc*(ycref-y1); // references off point (x1,y1) with value f1
      newvalue+=interpolate;
    }
    else{
      newvalue=((y3 - ycref)*(x2*f1 - x1*f2) + x3*(-(y2*f1) + ycref*f1 + y1*f2 - ycref*f2) + 
                (-(x2*y1) + x1*y2 - x1*ycref + x2*ycref)*f3 + xcref*(y2*f1 - y3*f1 - y1*f2 + y3*f2 + y1*f3 - y2*f3))/(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
    }
    if(DEBUGINTERP)   fprintf(stderr,"%d %d %d %d :: %g %g %g\n",hold,iold,jold,kold, newvalue,interpolate,f1); fflush(stderr);
    
    if(DATATYPE==0){
      // note that for images the function value may go one image value up
      // or down on a whim.  For a given pallette file (say john.pal) the
      // bottom then can have sharp colors changes to black.  Let's up
      // everything by 0.5 to avoid this.
      newvalue+=0.5;
      if(newvalue>255.0) newvalue=255;
      else if(newvalue<0.0) newvalue=0;

    }

    // assign for return out of this function
    newdatatemp[coli]=newvalue;
    
  }// over columns

  return(0);
}




// simple bi-linear interpolation wrapper
int bilinear_interp_simple(int hold, int iold, int jold, int kold, FTYPE *X, FTYPE *startx, FTYPE *dx, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp)
{
  FTYPE ftemp;
  int atboundary[NDIM],nearboundary[NDIM];
  FTYPE dhold,diold,djold,dkold;
  int nearest_interp_ij(int hold, int iold,int jold,int kold,unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
  int bilinear_interp_ij(int hold, int iold, int jold,int kold,  FTYPE dhold, FTYPE diold,FTYPE djold,FTYPE dkold,unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);


  is_atboundary(&hold,&iold,&jold,&kold,atboundary,nearboundary);



  if(BILINEARREDUCE2NEARESTATBOUNDARY && (nearboundary[0]||nearboundary[1]||nearboundary[2]||nearboundary[3]) ){
    // then just use nearest neighbor
    return(nearest_interp_ij(hold,iold,jold,kold,oldimage,olddata,newdatatemp));
  }
  else{

    // then need to shift down since going to use upper point inside interp
    // don't need to shift up since already shifted inside atboundary
    if(nearboundary[0]==1) hold=hold -1*(oN0>1);
    if(nearboundary[1]==1) iold=iold -1*(oN1>1);
    if(nearboundary[2]==1) jold=jold -1*(oN2>1);
    if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
      // then no need to fix since boundary condtion data is at upper k
    }
    else{
      if(nearboundary[3]==1) kold=kold -1*(oN3>1);
    }

    // below valid at boundary
    dhold = ((X[0]-startx[0])/dx[0] - 0.5) - hold;
    diold = ((X[1]-startx[1])/dx[1] - 0.5) - iold;
    djold = ((X[2]-startx[2])/dx[2] - 0.5) - jold;
    dkold = ((X[3]-startx[3])/dx[3] - 0.5) - kold;



    return(bilinear_interp_ij(hold,iold,jold,kold, dhold,diold,djold,dkold, oldimage,olddata, newdatatemp));
  }
  
  return(0);

}


// check if at the boundary of the grid of data
int is_atboundary(int *holdvar, int *ioldvar, int *joldvar, int *koldvar, int *atboundary, int *nearboundary)
{
  int hold, iold, jold, kold;
  int holdp, ioldp,joldp,koldp;
  int jj;



  hold=*holdvar;
  iold=*ioldvar;
  jold=*joldvar;
  kold=*koldvar;


  DLOOPA(jj){
    nearboundary[jj]=0;
    atboundary[jj]=0;
  }


  if(BOUNDARYEXTRAP==1){
    // don't modify old grid positions even if near boundary
    // Only modify to ensure not completely off
    // h,i,j,k vary from 0 to oN-1 but for bilinear can go from -1 to oN
    if(hold<-1) hold=-1;
    if(hold>oN0) hold=oN0;

    if(iold<-1) iold=-1;
    if(iold>oN1) iold=oN1;

    if(jold<-1) jold=-1;
    if(jold>oN2) jold=oN2;

    if(kold<-1) kold=-1;
    if(kold>oN3) kold=oN3;

  }
  else{
    // ensure bi-linear doesn't use bad values by shifting stencil effectively.  Can lead to poor interpolation (since really extrapolation)
    // if refining, then oN1 is really old image size, not refined, which is correct
    // take care of boundary effects
    if(hold>=oN0){ hold=oN0-1; atboundary[0]=1;}
    if(hold>=oN0-1){ hold=oN0-1; nearboundary[0]=1;}
    if(hold<0){ hold=0; atboundary[0]=-1;}

    if(iold>=oN1){ iold=oN1-1; atboundary[1]=1;}
    if(iold>=oN1-1){ iold=oN1-1; nearboundary[1]=1;}
    if(iold<0){ iold=0; atboundary[1]=-1;}

    if(jold>=oN2){ jold=oN2-1; atboundary[2]=1;}
    if(jold>=oN2-1){ jold=oN2-1; nearboundary[2]=1;}
    if(jold<0){ jold=0; atboundary[2]=-1;}

    if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
      // then ok for kold to be oN3-1 since oN3 data exists
      if(kold>=oN3){ kold=kold-oN3;}
      if(kold<0){ kold=kold+oN3;}
    }
    else{
      if(kold>=oN3-1){ kold=oN3-1; nearboundary[3]=1;}
      if(kold>=oN3){ kold=oN3-1; atboundary[3]=1;}
      if(kold<0){ kold=0; atboundary[3]=-1;}
    }

  }





  // assign final values
  *holdvar=hold;
  *ioldvar=iold;
  *joldvar=jold;
  *koldvar=kold;


  return(0);

}



// bi-linear interpolation (used by bi-linear wrapper)
static int bilinear_interp_ij(int hold, int iold, int jold,int kold,  FTYPE dhold, FTYPE diold,FTYPE djold,FTYPE dkold,  unsigned char*****oldimage,FTYPE*****olddata,FTYPE *newdatatemp)
{
  int holdp,ioldp,joldp,koldp;
  FTYPE newvalue;
  FTYPE myfun[MAXCOLS][NUMIJKMEM],dist[NUMIJKMEM];
  int hsel[NUMIJKMEM],isel[NUMIJKMEM],jsel[NUMIJKMEM],ksel[NUMIJKMEM];
  int p;
  FTYPE ftemp;
  FTYPE totaldist;


  holdp = hold+1*(oN0>1) ;
  ioldp = iold+1*(oN1>1) ;
  joldp = jold+1*(oN2>1) ;
  koldp = kold+1*(oN3>1) ;

  // did have boundary check here

  // below ?sel[] assignments are same as in quad_interp()
  // get distances and function values
  if(NUMIJK>=4){ // these only used in 2D
    hsel[1]=hold; isel[1]=iold;  jsel[1]=jold;  ksel[1]=kold;
    hsel[2]=hold; isel[2]=ioldp; jsel[2]=jold;  ksel[2]=kold;
    hsel[3]=hold; isel[3]=iold;  jsel[3]=joldp; ksel[3]=kold;
    hsel[4]=hold; isel[4]=ioldp; jsel[4]=joldp; ksel[4]=kold;
  }

  if(NUMIJK>=8){ // these only used in 3D
    hsel[5]=hold; isel[5]=iold;  jsel[5]=jold;  ksel[5]=koldp;
    hsel[6]=hold; isel[6]=ioldp; jsel[6]=jold;  ksel[6]=koldp;
    hsel[7]=hold; isel[7]=iold;  jsel[7]=joldp; ksel[7]=koldp;
    hsel[8]=hold; isel[8]=ioldp; jsel[8]=joldp; ksel[8]=koldp;
  }

  if(NUMIJK>=16){ // these only used in 4D
    hsel[9]=holdp;  isel[9]=iold;   jsel[9]=jold;   ksel[9]=kold;
    hsel[10]=holdp; isel[10]=ioldp; jsel[10]=jold;  ksel[10]=kold;
    hsel[11]=holdp; isel[11]=iold;  jsel[11]=joldp; ksel[11]=kold;
    hsel[12]=holdp; isel[12]=ioldp; jsel[12]=joldp; ksel[12]=kold;

    hsel[13]=holdp; isel[13]=iold;  jsel[13]=jold;  ksel[13]=koldp;
    hsel[14]=holdp; isel[14]=ioldp; jsel[14]=jold;  ksel[14]=koldp;
    hsel[15]=holdp; isel[15]=iold;  jsel[15]=joldp; ksel[15]=koldp;
    hsel[16]=holdp; isel[16]=ioldp; jsel[16]=joldp; ksel[16]=koldp;
  }

  // get distances
  if(NUMIJK>=4){ // these only used in 2D
    dist[1]=(1.-dhold)*(1. - diold)*(1.-djold)*(1.-dkold);
    dist[2]=(1.-dhold)*(diold)*(1. - djold)*(1.-dkold);
    dist[3]=(1.-dhold)*(1. - diold)*(djold)*(1.-dkold);
    dist[4]=(1.-dhold)*(diold)*(djold)*(1.-dkold);
  }

  if(NUMIJK>=8){
    dist[5]=(1.-dhold)*(1. - diold)*(1.-djold)*(dkold);
    dist[6]=(1.-dhold)*(diold)*(1. - djold)*(dkold);
    dist[7]=(1.-dhold)*(1. - diold)*(djold)*(dkold);
    dist[8]=(1.-dhold)*(diold)*(djold)*(dkold);
  }

  if(NUMIJK>=16){
    dist[9]=(dhold)*(1. - diold)*(1.-djold)*(1.-dkold);
    dist[10]=(dhold)*(diold)*(1. - djold)*(1.-dkold);
    dist[11]=(dhold)*(1. - diold)*(djold)*(1.-dkold);
    dist[12]=(dhold)*(diold)*(djold)*(1.-dkold);

    dist[13]=(dhold)*(1. - diold)*(1.-djold)*(dkold);
    dist[14]=(dhold)*(diold)*(1. - djold)*(dkold);
    dist[15]=(dhold)*(1. - diold)*(djold)*(dkold);
    dist[16]=(dhold)*(diold)*(djold)*(dkold);
  }


  int coli;
  for(coli=0;coli<numoutputcols;coli++){ // over all independent columsn of data
  
    for(p=1;p<=NUMIJK;p++){
      if(DATATYPE==0){
        myfun[coli][p]=(double)((unsigned char)oldimage[hsel[p]][isel[p]][jsel[p]][ksel[p]]);
      }
      else{
        myfun[coli][p]=olddata[coli][hsel[p]][isel[p]][jsel[p]][ksel[p]];
      }
      //    fprintf(stderr,"bi: p=%d dist=%10.5g myfun=%10.5g\n",p,dist[p],myfun[0][p]);
    }
    
    // get total distance and cumulative value
    newvalue=0.0;
    totaldist=0.0;
    for(p=1;p<=NUMIJK;p++){
      newvalue+=dist[p]*myfun[coli][p];
      totaldist+=dist[p];
    }
    
    // normalize value
    newvalue = newvalue/totaldist;
    
    if(DATATYPE==0){
      // note that for images the function value may go one image value up
      // or down on a whim.  For a given pallette file (say john.pal) the
      // bottom then can have sharp colors changes to black.  Let's up
      // everything by 0.5 to avoid this.
      newvalue+=0.5;
      if(newvalue>255.0) newvalue=255;
      else if(newvalue<0.0) newvalue=0;
    }

    newdatatemp[coli]=newvalue;
  }

  return(0);
}



// nearest neighbor
static int nearest_interp_ij(int hold, int iold,int jold,int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp)
{
  FTYPE newvalue;
  int atboundary;
  int coli;


  if(BOUNDARYEXTRAP==0){
    atboundary=0;
    // if refining, then oN1 is really old image size, not refined, which is correct
    /* take care of boundary effects */
    if(hold>=oN0){ hold=oN0-1;atboundary=1;}
    if(hold<0){ hold=0;atboundary=1;}

    if(iold>=oN1){ iold=oN1-1;atboundary=1;}
    if(iold<0){ iold=0;atboundary=1;}

    if(jold>=oN2){ jold=oN2-1;atboundary=1;}
    if(jold<0){ jold=0;atboundary=1;}
  }

  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    if(kold<0){ kold=kold+oN3;}
    if(kold>=oN3){ kold=kold-oN3;}
  }
  else{
    if(BOUNDARYEXTRAP==0){
      if(kold<0){ kold=0;atboundary=1;}
      if(kold>=oN3){ kold=oN3-1;atboundary=1;}
    }
  }

  if(DEBUGINTERP){
    for(coli=0;coli<numoutputcols;coli++) fprintf(stderr,"Problem?: hold=%d iold=%d jold=%d kold=%d   oN0=%d oN1=%d oN2=%d oN3=%d :: coli=%d %21.15g\n",hold,iold,jold,kold,oN0,oN1,oN2,oN3,coli,olddata[coli][hold][iold][jold][kold]);
  }

  
  for(coli=0;coli<numoutputcols;coli++){ // over all independent columsn of data
    if(DATATYPE==0)   newdatatemp[coli]=(FTYPE)oldimage[coli][hold][iold][jold][kold] ;
    else    newdatatemp[coli]=(FTYPE)olddata[coli][hold][iold][jold][kold] ;
  }


  return(0);
}




//
// don't use, not working
// also not setup for oN0>1 or oN3>1
//
// bi-linear interpolation using real distances instead of grid-based distances
static int bilinear_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref, int hold, int iold, int jold, int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp)
{
  FTYPE f1,f2;
  int holdp,ioldp,joldp,koldp;
  FTYPE newvalue;
  FTYPE dist[NUMIJKMEM],mytc[NUMIJKMEM],myxc[NUMIJKMEM],myyc[NUMIJKMEM],myzc[NUMIJKMEM],myfun[MAXCOLS][NUMIJKMEM];
  FTYPE tcref,xcref,ycref,zcref;
  int hsel[NUMIJKMEM],isel[NUMIJKMEM],jsel[NUMIJKMEM],ksel[NUMIJKMEM];
  int p;
  FTYPE ftemp;
  int nearest_interp(int hold, int iold,int jold,int kold,unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp);
  int quad_interp(FTYPE tref, FTYPE rref, FTYPE thref, FTYPE phref,  int hold, int iold, int jold, int kold,  unsigned char*****oldimage,FTYPE*****olddata,  FTYPE*dist,  FTYPE *mytc, FTYPE*myxc, FTYPE*myyc, FTYPE*myzc, FTYPE(*myfun)[NUMIJKMEM],  FTYPE *tcref,FTYPE*xcref,FTYPE*ycref,FTYPE *zcref);




  if(quad_interp(tref, rref, thref, phref,  hold, iold, jold, kold,  oldimage,olddata,  dist, mytc, myxc, myyc, myzc,  myfun,  &tcref,&xcref,&ycref,&zcref)>=1){
    nearest_interp(hold,iold,jold,kold,oldimage,olddata,newdatatemp);
    return(0);
  }


  int coli;
  for(coli=0;coli<numoutputcols;coli++){ // over all independent columsn of data
      

    // now assign actual weighted value, bilinear filtering
    newvalue=0.0;
    // first 3 are weights by true distance
    if(0){
      // =
      // half of _ and the other half of -
      // 1-2 and 3-4
      p=1; newvalue =0.5*(1.0-dist[p]/(dist[1]+dist[2]))*myfun[coli][p];
      p=2; newvalue+=0.5*(1.0-dist[p]/(dist[1]+dist[2]))*myfun[coli][p];
      p=3; newvalue+=0.5*(1.0-dist[p]/(dist[3]+dist[4]))*myfun[coli][p];
      p=4; newvalue+=0.5*(1.0-dist[p]/(dist[3]+dist[4]))*myfun[coli][p];
    }
    else if(0){
      // || // 1-3 and 2-4
      p=1; newvalue =0.5*(1.0-dist[p]/(dist[1]+dist[3]))*myfun[coli][p];
      p=3; newvalue+=0.5*(1.0-dist[p]/(dist[1]+dist[3]))*myfun[coli][p];
      p=2; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[4]))*myfun[coli][p];
      p=4; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[4]))*myfun[coli][p];
    }
    else if(0){
      // X // 1-4 and 2-3
      p=1; newvalue =0.5*(1.0-dist[p]/(dist[1]+dist[4]))*myfun[coli][p];
      p=4; newvalue+=0.5*(1.0-dist[p]/(dist[1]+dist[4]))*myfun[coli][p];
      p=2; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[3]))*myfun[coli][p];
      p=3; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[3]))*myfun[coli][p];
    }
    else if(1){
      // weight by x and y
      p=1; newvalue =(1.0-fabs(xcref-myxc[1])/fabs(myxc[2]-myxc[1]))*(1.0-fabs(ycref-myyc[1])/fabs(myyc[2]-myyc[1]))*myfun[coli][p];
      p=2; newvalue +=(1.0-fabs(myxc[2]-xcref)/fabs(myxc[2]-myxc[1]))*(1.0-fabs(myyc[2]-ycref)/fabs(myyc[2]-myyc[1]))*myfun[coli][p];
      p=3; newvalue +=(1.0-fabs(xcref-myxc[3])/fabs(myxc[3]-myxc[4]))*(1.0-fabs(ycref-myyc[3])/fabs(myyc[3]-myyc[4]))*myfun[coli][p];
      p=4; newvalue +=(1.0-fabs(myxc[4]-xcref)/fabs(myxc[3]-myxc[4]))*(1.0-fabs(myyc[4]-ycref)/fabs(myyc[3]-myyc[4]))*myfun[coli][p];
    }
    else if(0){
      if(myxc[2]>myxc[1]){
        f1=(xcref-myxc[1])*(myfun[coli][2]-myfun[coli][1])/(myxc[2]-myxc[1])+myfun[coli][1];
      }
      else{
        f1=(xcref-myxc[2])*(myfun[coli][1]-myfun[coli][2])/(myxc[1]-myxc[2])+myfun[coli][2];
      }
      if(myxc[4]>myxc[3]){
        f2=(xcref-myxc[3])*(myfun[coli][4]-myfun[coli][3])/(myxc[4]-myxc[3])+myfun[coli][3];
      }
      else{
        f2=(xcref-myxc[4])*(myfun[coli][3]-myfun[coli][4])/(myxc[3]-myxc[4])+myfun[coli][4];
      }
      if(myyc[3]>myyc[1]){ // 1 associated with f1, 3 associated with f2
        newvalue=(ycref-myyc[1])*(f2-f1)/(myyc[3]-myyc[1])+f1;
      }
      else{
        newvalue=(ycref-myyc[3])*(f1-f2)/(myyc[1]-myyc[3])+f2;
      }
    }
    
    if(DATATYPE==0){
      // note that for images the function value may go one image value up
      // or down on a whim.  For a given pallette file (say john.pal) the
      // bottom then can have sharp colors changes to black.  Let's up
      // everything by 0.5 to avoid this.
      newvalue+=0.5;
      if(newvalue>255.0) newvalue=255;
      else if(newvalue<0.0) newvalue=0;
    }

    newdatatemp[coli]=newvalue;
  }

  return(0);
}





static int nearest_interp(int hold,int iold,int jold,int kold, unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp)
{
  return(nearest_interp_ij(hold,iold,jold,kold,oldimage,olddata,newdatatemp));
}



// bicubic interpolation wrapper
// treats t and \phi directions as nearest neighbor
static int bicubic_interp_wrap(int nt, int nx, int ny, int nz,  int hold, int iold, int jold, int kold,  FTYPE x0, FTYPE x1, FTYPE x2,FTYPE x3,  unsigned char*****oldimage,FTYPE*****olddata, FTYPE *newdatatemp)
{
  int j,k;
  static int firsttime=1;
  static int lasth,lastk;
  int recomputederivs;
  FTYPE Xget[NDIM];

  static FTYPE **y1a,**y2a,**y12a;
  static FTYPE **ya;
  void dervs_for_bicubic(int nx, int ny, FTYPE **ya, FTYPE **y1a, FTYPE **y2a, FTYPE **y12a);
  FTYPE bicubic_interp(int nx,int ny,int hl,int il,int jl,int kl, FTYPE *Xget,FTYPE **ya,FTYPE **y1a, FTYPE **y2a, FTYPE**y12a,FTYPE dx1,FTYPE dx2);


  if(firsttime){
    lasth=hold;
    lastk=kold;


    //  fprintf(stderr,"got here1\n"); fflush(stderr);

    y1a = fmatrix(0,nx-1,0,ny-1) ;
    if((y1a==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }
    y2a = fmatrix(0,nx-1,0,ny-1) ;
    if((y2a==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }
    y12a = fmatrix(0,nx-1,0,ny-1) ;
    if((y12a==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }

    ya = fmatrix(0,nx-1,0,ny-1) ;
    if((ya==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }

  }

  
  int priorcoli=-100;
  int coli;
  for(coli=0;coli<numoutputcols;coli++){ // over all independent columsn of data

    // always copy values
    if(DATATYPE==0){
      for(j=0;j<nx;j++) for(k=0;k<ny;k++){
          ya[j][k]=(FTYPE)oldimage[coli][hold][j][k][kold]; // GODMARK3D -- nearest neighbor in h-direction and k-direction
        }
    }
    // else if(DATATYPE==1) ya=olddata;
    else{
      for(j=0;j<nx;j++) for(k=0;k<ny;k++){
          ya[j][k]=olddata[coli][hold][j][k][kold]; // GODMARK3D -- nearest neighbor in h-direction and k-direction
        }
    }
    
    
    if(firsttime || lastk!=kold|| lasth!=hold || priorcoli!=coli){
      // required to recompute if coli changed
      // GODMARK3D -- for now redo derivs for each new kold and each new hold
      recomputederivs=1;
    }
    else recomputederivs=0;


    if(recomputederivs){
      lasth=hold; // store hold into lasth for which recomputed derivatives
      lastk=kold; // store kold into lastk for which recomputed derivatives
      priorcoli=coli;
      
      
      // compute derivatives of low resolution data
      dervs_for_bicubic(nx,ny,ya,y1a,y2a,y12a);
    }
    
    
    
    Xget[0]=x0;
    Xget[1]=x1;
    Xget[2]=x2;
    Xget[3]=x3;
    
    
      
    newdatatemp[coli]=bicubic_interp(nx,ny,hold,iold,jold,kold,Xget,ya,y1a, y2a, y12a,dx[1],dx[2]);
  }


  // no longer using firsttime after below line
  if(firsttime) firsttime=0;


  return(0);
}




// compute derivatives of data
// COLIMARK: Done per data or column point outside this function
void dervs_for_bicubic(int nx, int ny, FTYPE **ya, FTYPE **y1a, FTYPE **y2a, FTYPE **y12a)
{
  int j,k;


  for(k=0;k<ny;k++) for(j=0;j<nx;j++){


      if((j>0)&&(j<nx-1)&&(k>0)&&(k<ny-1)){
        y1a[j][k]=(ya[j+1][k]-ya[j-1][k])/(2.0*dx[1]);
        y2a[j][k]=(ya[j][k+1]-ya[j][k-1])/(2.0*dx[2]);
        y12a[j][k]=(ya[j+1][k+1]-ya[j+1][k-1]-ya[j-1][k+1]+ya[j-1][k-1])/(4.0*dx[1]*dx[2]);
      }
      else if((j==0)&&(k==0)){
        y1a[j][k]=(ya[j+1][k]-ya[j][k])/dx[1];
        y2a[j][k]=(ya[j][k+1]-ya[j][k])/dx[2];
        y12a[j][k]=(ya[j+1][k+1]-ya[j+1][k]-ya[j][k+1]+ya[j][k])/(dx[1]*dx[2]);      
      }
      else if((j==nx-1)&&(k==0)){
        y1a[j][k]=(ya[j][k]-ya[j-1][k])/(dx[1]);
        y2a[j][k]=(ya[j][k+1]-ya[j][k])/(dx[2]);
        y12a[j][k]=(ya[j][k+1]-ya[j][k]-ya[j-1][k+1]+ya[j-1][k])/(dx[1]*dx[2]);
      }
      else if((j==0)&&(k==ny-1)){
        y1a[j][k]=(ya[j+1][k]-ya[j][k])/(dx[1]);
        y2a[j][k]=(ya[j][k]-ya[j][k-1])/(dx[2]);
        y12a[j][k]=(ya[j+1][k]-ya[j+1][k-1]-ya[j][k+1]+ya[j][k-1])/(dx[1]*dx[2]);
      }
      else if(j==0){
        y1a[j][k]=(ya[j+1][k]-ya[j][k])/(dx[1]);
        y2a[j][k]=(ya[j][k+1]-ya[j][k-1])/(2.0*dx[2]);
        y12a[j][k]=(ya[j+1][k+1]-ya[j+1][k-1]-ya[j][k+1]+ya[j][k-1])/(2.0*dx[1]*dx[2]);
      }
      else if(j==nx-1){
        y1a[j][k]=(ya[j][k]-ya[j-1][k])/(dx[1]);
        y2a[j][k]=(ya[j][k+1]-ya[j][k-1])/(2.0*dx[2]);
        y12a[j][k]=(ya[j][k+1]-ya[j][k-1]-ya[j-1][k+1]+ya[j-1][k-1])/(2.0*dx[1]*dx[2]);
      }
      else if(k==0){
        y1a[j][k]=(ya[j+1][k]-ya[j-1][k])/(2.0*dx[1]);
        y2a[j][k]=(ya[j][k+1]-ya[j][k])/(dx[2]);
        y12a[j][k]=(ya[j+1][k+1]-ya[j+1][k]-ya[j-1][k+1]+ya[j-1][k])/(2.0*dx[1]*dx[2]);
      }
      else if(k==ny-1){
        y1a[j][k]=(ya[j+1][k]-ya[j-1][k])/(2.0*dx[1]);
        y2a[j][k]=(ya[j][k]-ya[j][k-1])/(dx[2]);
        y12a[j][k]=(ya[j+1][k]-ya[j+1][k-1]-ya[j-1][k]+ya[j-1][k-1])/(2.0*dx[1]*dx[2]);
      }
      else{
        fprintf(stderr,"No such j=%d k=%d condition\n",j,k);
        exit(1);
      }

#if(0) //test
      y1a[j][k]=0;
      y2a[j][k]=0;
      y12a[j][k]=0;
    
#endif
      // fprintf(stdout,"%g %g %g %g\n",ya[j][k],y1a[j][k],y2a[j][k],y12a[j][k]);

    }
  //  exit(0);


}



// see Numerical recipies figure 3.6.1 and section 3.6 on page 123
//      ftemp=bicubic_interp(nxlow,nylow,il,jl,ftempi,ftempj,ya,y1a,y2a,y12a,dx[1],dx[2]);
// COLIMARK: Done per column of data, not over all data
FTYPE bicubic_interp(int nx,int ny, int hl,int il,int jl,int kl, FTYPE *Xget,FTYPE **ya,FTYPE **y1a, FTYPE **y2a, FTYPE**y12a,FTYPE dx1,FTYPE dx2)
{
  FTYPE Xp[4+1][NDIM];
  int pointsh[4+1],pointsi[4+1],pointsj[4+1],pointsk[4+1];
  FTYPE yap[4+1],y1ap[4+1],y2ap[4+1],y12ap[4+1];
  int i;
  extern void bcuint(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE x1l,
                     FTYPE x1u, FTYPE x2l, FTYPE x2u, FTYPE x1, FTYPE x2, FTYPE *ansy,
                     FTYPE *ansy1, FTYPE *ansy2);
  FTYPE answer,answerd1,answerd2;
  int atboundary;

  


  atboundary=0;
  // should just shift stencil near outer boundary

  pointsi[1]=il;   pointsj[1]=jl;
  if(pointsi[1]>nx-1){
    atboundary=1;
    pointsi[1]=nx-1;
  }
  if(pointsj[1]>ny-1){
    atboundary=1;
    pointsj[1]=ny-1;
  }

  pointsi[2]=il+1; pointsj[2]=jl;
  if(pointsi[2]>nx-1) atboundary=1;
  if(pointsj[2]>ny-1) atboundary=1;

  pointsi[3]=il+1; pointsj[3]=jl+1;
  if(pointsi[3]>nx-1) atboundary=1;
  if(pointsj[3]>ny-1) atboundary=1;

  pointsi[4]=il;   pointsj[4]=jl+1;
  if(pointsi[4]>nx-1) atboundary=1;
  if(pointsj[4]>ny-1) atboundary=1;

  if(atboundary){
    // then just use nearest neighbor
    return(ya[pointsi[1]][pointsj[1]]);
  }
  // otherwise do good interpolation

  //  for(k=0;k<ny;k++) for(j=0;j<nx;j++) {
  //    fprintf(stdout,"%g %g %g %g\n",ya[j][k],y1a[j][k],y2a[j][k],y12a[j][k]);
  //  }
  //  exit(0);


  // GODMARK3D
  pointsh[1]=pointsh[2]=pointsh[3]=pointsh[4]=hl;
  pointsk[1]=pointsk[2]=pointsk[3]=pointsk[4]=kl;

  // get surrounding points positions.
  for(i=1;i<=4;i++){  // GODMARK3D
    //    fprintf(stdout,"%d\n",ya);
    oldf_ijk2x123(pointsh[i], pointsi[i], pointsj[i], pointsk[i], Xp[i]);
    //    for(j=0;j<4;j++) Xp[i][j]=0;
    yap[i]=ya[pointsi[i]][pointsj[i]];
    y1ap[i]=y1a[pointsi[i]][pointsj[i]];
    y2ap[i]=y2a[pointsi[i]][pointsj[i]];
    y12ap[i]=y12a[pointsi[i]][pointsj[i]];
    //    fprintf(stdout,"i=%d j=%d Xp[i][1]=%g Xp[i][2]=%g yap=%g y1ap=%g y2ap=%g y12ap=%g\n",i,j,Xp[i][1],Xp[i][2],yap[i],y1ap[i],y2ap[i],y12ap[i]);
    //    fprintf(stdout,"%d %d %d %g %g %g %g %g %g %g\n",i,pointsi[i],pointsj[i],Xp[i][1],Xp[i][2],yap[i],y1ap[i],y2ap[i],y12ap[i],ya[pointsi[i]][pointsj[i]]);
  }

  // too much information in Xp, only using relevant information.  Assumes grid is rectangular, which it is.
  bcuint(yap,y1ap,y2ap,y12ap,Xp[1][1],Xp[2][1],Xp[1][2],Xp[4][2],Xget[1],Xget[2],&answer,&answerd1,&answerd2);

  return(answer);

}



static void interpicoord(FTYPE *X,int loc, int *h, int *i, int *j, int *k)
{

  
  lowericoord(X,loc,h,i,j,k);

  // restrict to old data size
  if(*h<0) *h=0;
  if(*h>oN0-1) *h=oN0-1;

  if(*i<0) *i=0;
  if(*i>oN1-1) *i=oN1-1;

  if(*j<0) *j=0;
  if(*j>oN2-1) *j=oN2-1;

  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    if(*k<0){
      *k=*k+oN3;
      // must shift X too so consistent when obtaining, e.g., dkold for bilinear
      X[3]=X[3]+(Xmax[3]-startx[3]);
    }
    if(*k>oN3-1){
      *k=*k-oN3;
      // must shift X too so consistent when obtaining, e.g., dkold for bilinear
      X[3]=X[3]-(Xmax[3]-startx[3]);
    }
  }
  else{
    if(*k<0) *k=0;
    if(*k>oN3-1) *k=oN3-1;
  }



}




// unlike icoord() in coord.c, this finds lower integer so that old and oldp bound ifloat(X)
static int lowericoord(FTYPE *X,int loc, int *h, int *i, int *j, int *k)
{

  if(loc == CENT){
    *h = ROUND2INT((X[0]-startx[0])/dx[0] + 0.0) - 1 ;
    *i = ROUND2INT((X[1]-startx[1])/dx[1] + 0.0) - 1 ;
    *j = ROUND2INT((X[2]-startx[2])/dx[2] + 0.0) - 1 ;
    *k = ROUND2INT((X[3]-startx[3])/dx[3] + 0.0) - 1 ;
  }


  return(0);

}









// refine data
void refine_data(void)
{
  int coli,h,i,j,k;

  //////////////////////////
  //
  // check for refinement of old data
  //
  didrefine=0;
  if(fabs(refinefactor-1.0)>0.1){
    didrefine=1;
    fprintf(stderr,"refine\n"); fflush(stderr);

    // then refine or derefine
    if(oN0!=1) roN0=(int)(refinefactor*(FTYPE)oN0);
    else roN0=oN0;
    if(oN1!=1) roN1=(int)(refinefactor*(FTYPE)oN1);
    else roN1=oN1;
    if(oN2!=1) roN2=(int)(refinefactor*(FTYPE)oN2);
    else roN2=oN2;
    if(oN3!=1) roN3=(int)(refinefactor*(FTYPE)oN3);
    else roN3=oN3;

    fprintf(stderr,"refine from %dX%dX%dX%d to %dX%dX%dX%d\n",oN0,oN1,oN2,oN3, roN0,roN1,roN2,roN3); fflush(stderr);

    if(refinefactor>1.0){
      if(DATATYPE==0){
        oldimage = c5matrix(0,numoutputcols-1,0,roN0-1,0,roN1-1,0,roN2-1,0,roN3-1) ;
        for(coli=0;coli<numoutputcols;coli++) for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)      for(i=0;i<oN1;i++) for(h=0;h<oN0;h++) oldimage[coli][h][i][j][k]=oldimage0[coli][h][i][j][k];
      }
      else{
        olddata = f5matrix(0,numoutputcols-1,0,roN0-1,0,roN1-1,0,roN2-1,0,roN3-1) ;
        for(coli=0;coli<numoutputcols;coli++) for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)      for(i=0;i<oN1;i++) for(h=0;h<oN0;h++) olddata[coli][h][i][j][k]=olddata0[coli][h][i][j][k];
      }
      //writeimage("jontest1.r8",oldimage,roN0,roN1,roN2,roN3);
      //      exit(0);
      low2high(oN0, oN1, oN2, oN3, roN0, roN1, roN2, roN3, oldimage,olddata);
    }
    else{
      high2low(oN0, oN1, oN2, oN3, roN0, roN1, roN2, roN3, oldimage0,olddata0);
      if(DATATYPE==0){
        oldimage = c5matrix(0,numoutputcols-1,0,roN0-1,0,roN1-1,0,roN2-1,0,roN3-1) ;
        for(coli=0;coli<numoutputcols;coli++) for(k=0;k<roN3;k++) for(j=0;j<roN2;j++)      for(i=0;i<roN1;i++)  for(h=0;h<roN0;h++) oldimage[coli][h][i][j][k]=oldimage0[coli][h][i][j][k];
      }
      else{
        olddata = f5matrix(0,numoutputcols-1,0,roN0-1,0,roN1-1,0,roN2-1,0,roN3-1) ;
        for(coli=0;coli<numoutputcols;coli++) for(k=0;k<roN3;k++)  for(j=0;j<roN2;j++)      for(i=0;i<roN1;i++) for(h=0;h<roN0;h++) olddata[coli][h][i][j][k]=olddata0[coli][h][i][j][k];
      }
    }
    // reset size
    oN0=roN0;
    oN1=roN1;
    oN2=roN2;
    oN3=roN3;
    fprintf(stderr,"done refining\n"); fflush(stderr);
  }
  else{
    // no refinement
    if(DATATYPE==0) oldimage=oldimage0;
    else{
      olddata=olddata0;
    }
  }
}






// assumes data is of size nxhigh*nyhigh, but only the smallest portion is filled with nxlow*nylow data
// assumes square grid data and going from integer sizes to larger integer size grid.
void low2high(int ntlow, int nxlow, int nylow, int nzlow,  int nthigh, int nxhigh, int nyhigh, int nzhigh,  unsigned char*****oldimage,FTYPE*****olddata)
{
  FTYPE **Ilowf;
  unsigned char **Ilowc;
  int hh, ih,jh,kh;
  int hl, il,jl,kl;
  FTYPE dhl, dil, djl,dkl,ftemp,ftemph,ftempi,ftempj,ftempk;
  FTYPE t,r,th,ph;
  unsigned char uctemp;

  int j,k;

  FTYPE **y1a,**y2a,**y12a;
  FTYPE **ya;
  FTYPE Xget[NDIM];


  if(numoutputcols>1){
    fprintf(stderr,"low2high not setup for numoutputcols=%d>1\n",numoutputcols);
    exit(1);
  }


  // GODMARK3D
  hl=0;
  kl=0;

  hh=0;
  kh=0;

  t=0;
  ph=0;

  dhl=0;
  dkl=0;

  ftemph=0;
  ftempk=0;
  // GODMARK3D


  //  fprintf(stderr,"got here0\n"); fflush(stderr);

  if(DATATYPE==0){
    /* make arrays for images */    
    Ilowc = cmatrix(0,nxlow-1,0,nylow-1) ;
    if((Ilowc==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }
  }
  else{
    /* make arrays for dumps */
    Ilowf = fmatrix(0,nxlow-1,0,nylow-1) ;
    if((Ilowf==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }
  }

  // setup Ilow
  for(il=0;il<nxlow;il++){
    for(jl=0;jl<nylow;jl++){
      if(DATATYPE==0) Ilowc[il][jl]=oldimage[0][0][il][jl][0];// GODMARK3D // COLIMARK
      else Ilowf[il][jl]=olddata[0][0][il][jl][0]; // GODMARK3D // COLIMARK
    }
  }

  //  fprintf(stderr,"got here1\n"); fflush(stderr);

  y1a = fmatrix(0,nxlow-1,0,nylow-1) ;
  if((y1a==NULL)){
    fprintf(stderr,"Cannot allocate memory\n");
    exit(1);
  }
  y2a = fmatrix(0,nxlow-1,0,nylow-1) ;
  if((y2a==NULL)){
    fprintf(stderr,"Cannot allocate memory\n");
    exit(1);
  }
  y12a = fmatrix(0,nxlow-1,0,nylow-1) ;
  if((y12a==NULL)){
    fprintf(stderr,"Cannot allocate memory\n");
    exit(1);
  }

  if(DATATYPE==0){
    ya = fmatrix(0,nxlow-1,0,nylow-1) ;
    if((ya==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }
    for(j=0;j<nxlow;j++) for(k=0;k<nylow;k++){
        ya[j][k]=(FTYPE)oldimage[0][0][j][k][0]; // GODMARK3D // COLIMARK
      }
  }
  // else if(DATATYPE==1) ya=olddata;
  else{
    ya = fmatrix(0,nxlow-1,0,nylow-1) ;
    if((ya==NULL)){
      fprintf(stderr,"Cannot allocate memory\n");
      exit(1);
    }
    for(j=0;j<nxlow;j++) for(k=0;k<nylow;k++){
        ya[j][k]=olddata[0][0][j][k][0]; // GODMARK3D // COLIMARK
      }

  }

  //  fprintf(stderr,"got here2\n"); fflush(stderr);


  // compute derivatives of low resolution data
  dervs_for_bicubic(nxlow,nylow,ya,y1a,y2a,y12a);


  //  for(j=0;j<nxlow;j++) for(k=0;k<nylow;k++){
  //    fprintf(stdout,"%g %g %g %g\n",ya[j][k],y1a[j][k],y2a[j][k],y12a[j][k]);
  //  }
  //  exit(0);

  //  fprintf(stderr,"got here3\n"); fflush(stderr);

  
  /*
    if(DATATYPE==1){
    if(1){ // needs to be set
    for(jl=0;jl<nxlow;jl++)      for(il=0;il<nylow;il++) {
    ftemp=olddata[0][0][il][jl][0]; // GODMARK3D // COLIMARK
    if(ftemp<0.0) ftemp=0.0;
    if(ftemp>255.0) ftemp=255.0;
    uctemp=(unsigned char)ftemp;
    oldimage[0][0][ih][jh][0]=uctemp; // GODMARK3D // COLIMARK
    }
    }
    }
  */

  //  writeimage("jontest00.r8",Ilowc,nxlow,nylow,nzlow);
  //  writeimage("jontest01.r8",oldimage,nxhigh,nyhigh,nzhigh);
  //  writeimage("jontest02.r8",oldimage,nxlow,nylow,nzlow);
  //  exit(0);
  // determine high resolution version
  for(jh=0;jh<nyhigh;jh++){
    for(ih=0;ih<nxhigh;ih++) {
      // determine this pixels bilinear filtered value
      // determine nearest neighbor in low space

      // 0: NN // generates symmetric result
      // 1: bilinear // generates symmetric result, but uses different "stencil" and choices
      // 2: planar 
      // 3: bicubic 
#define HIGHLOWINTERPTYPE 3


#if(HIGHLOWINTERPTYPE==0)
      ///////////////////////////////////////////
      //
      ///////////////// Nearest neighbor interpolation
      //
      ///////////////////////////////////////////
      // determine nearest neighbor in low space
#define SHIFTNN (0.0)
      ftempi=((FTYPE)(ih+(FTYPE)SHIFTNN)*(FTYPE)nxlow/((FTYPE)nxhigh));
      il=(int)(ftempi-SHIFTNN);
      ftempj=((FTYPE)(jh+(FTYPE)SHIFTNN)*(FTYPE)nylow/((FTYPE)nyhigh));
      jl=(int)(ftempj-SHIFTNN);
      dil=(FTYPE)ftempi-(FTYPE)(il+SHIFTNN);
      djl=(FTYPE)ftempj-(FTYPE)(jl+SHIFTNN);

      FTYPE blob[MAXCOLS]
        nearest_interp_ij(hl, il, jl, kl ,Ilowc,Ilowf,blob);
      ftemp=blob[0]; // GODCOLIMARK
#elif(HIGHLOWINTERPTYPE==1)
      ///////////////////////////////////////////
      //
      ///////////////// bilinear interpolation
      //
      ///////////////////////////////////////////
      // determine nearest neighbor in low space
#define SHIFTBL (0.0)
      ftempi=((FTYPE)(ih+(FTYPE)SHIFTBL)*(FTYPE)nxlow/((FTYPE)nxhigh));
      il=(int)(ftempi-SHIFTBL);
      ftempj=((FTYPE)(jh+(FTYPE)SHIFTBL)*(FTYPE)nylow/((FTYPE)nyhigh));
      jl=(int)(ftempj-SHIFTBL);
      dil=(FTYPE)ftempi-(FTYPE)(il+SHIFTBL);
      djl=(FTYPE)ftempj-(FTYPE)(jl+SHIFTBL);

      //      fprintf(stderr,"%d %d %d %d %g %g\n",ih,jh,il,jl,dil,djl);

      ftemp=bilinear_interp_ij(hl,il,jl,kl, dhl,dil,djl,dkl, Ilowc,Ilowf);
#elif(HIGHLOWINTERPTYPE==2)
      ///////////////////////////////////////////
      //
      ///////////////// Planar interpolation
      //
      ///////////////////////////////////////////
      // determine nearest neighbor in low space
#define SHIFTPLANE (0.5) //symmetric for this choice
      ftempi=((FTYPE)(ih+(FTYPE)SHIFTPLANE)*(FTYPE)nxlow/((FTYPE)nxhigh))-(FTYPE)SHIFTPLANE;
      il=(int)(ftempi);
      ftempj=((FTYPE)(jh+(FTYPE)SHIFTPLANE)*(FTYPE)nylow/((FTYPE)nyhigh))-(FTYPE)SHIFTPLANE;
      jl=(int)(ftempj);
      dil=(FTYPE)(ftempi+SHIFTPLANE)-(FTYPE)il;
      djl=(FTYPE)(ftempj+SHIFTPLANE)-(FTYPE)jl;

      oldf_ijk2rthph(ftemph,ftempi,ftempj,ftempk, &t,&r,&th,&ph);
      // fprintf(stderr,"%d %d : : %g %g : %g %g\n",jh,ih,ftempi,ftempj,r,th); fflush(stderr);
      ftemp=plane_interp(t,r,th,ph, hl,il,jl,kl, Ilowc,Ilowf);
#elif(HIGHLOWINTERPTYPE==3)
      ///////////////////////////////////////////
      //
      ///////////////// bicubic interpolation
      //
      ///////////////////////////////////////////
      //#define SHIFTBC (-0.5) // symmetric with this choice
#define SHIFTBC (0.5) // symmetric with this choice
      ftempi=((FTYPE)(ih+(FTYPE)SHIFTBC)*(FTYPE)nxlow/((FTYPE)nxhigh))-(FTYPE)SHIFTBC;
      il=(int)(ftempi);
      ftempj=((FTYPE)(jh+(FTYPE)SHIFTBC)*(FTYPE)nylow/((FTYPE)nyhigh))-(FTYPE)SHIFTBC;
      jl=(int)(ftempj);


      //      fprintf(stderr,"got here4 %d %d:: %d %d\n",ih,jh,il,jl); fflush(stderr);

      // get the desired points position
      oldf_ijk2x123(ftemph,ftempi,ftempj,ftempk,Xget);
      ftemp=bicubic_interp(nxlow,nylow,hl,il,jl,kl,Xget,ya,y1a,y2a,y12a,dx[1],dx[2]);
#endif


      //      fprintf(stderr,"ih=%d jh=%d il=%d jl=%d dil=%g djl=%g ftemp=%10.5g\n",ih,jh,il,jl,dil,djl,ftemp); fflush(stderr);


      // just write back into the given multi-pointers
      if(DATATYPE==0) oldimage[0][0][ih][jh][0]=(unsigned char)ftemp;// GODMARK3D // COLIMARK
      else olddata[0][0][ih][jh][0]=ftemp; // GODMARK3D // COLIMARK

      //fprintf(stderr,"olddata[%d][%d][%d][%d]=%g\n",0,ih,jh,0,olddata[0][0][ih][jh][0]); // GODMARK3D // COLIMARK
    }
  }
  if(0){// debug (need to setup this memory stuff -- segfaults for data currently
    // in principle could output in different order if wanted
    if(DATATYPE!=0){
      if(1){ // needs to be set
        for(jh=0;jh<nyhigh;jh++)      for(ih=0;ih<nxhigh;ih++) {
            ftemp=olddata[0][0][ih][jh][0];// GODMARK3D // COLIMARK
            if(ftemp<0.0) ftemp=0.0;
            if(ftemp>255.0) ftemp=255.0;
            uctemp=(unsigned char)ftemp;
            oldimage[0][0][ih][jh][0]=uctemp; // GODMARK3D // COLIMARK
          }
      }
    }
    writeimage("jontest.r8",oldimage,nthigh,nxhigh,nyhigh,nzhigh);
    //exit(0);
  }
 

  if(DATATYPE==0){
    free_cmatrix(Ilowc,0,nxlow-1,0,nylow-1) ;
  }
  else{
    free_fmatrix(Ilowf,0,nxlow-1,0,nylow-1) ;
  }
  //  exit(0);

}


// area weighted based interpolation from high resolution to lower resolution
// assumes feed in high resolution and put low resolution version into same buffer at end
void high2low(int nthigh, int nxhigh, int nyhigh, int nzhigh,  int ntlow, int nxlow, int nylow, int nzlow,  unsigned char *****oldimage,FTYPE*****olddata)
{
  FTYPE *Ihigh;
  FTYPE *Ilow;
  FTYPE *IlowW;
  int ih,jh;
  int il,jl;
  FTYPE ilfrac,jlfrac;
  FTYPE ioldbcl,joldbcl;
  FTYPE deltaih,deltajh;
  FTYPE W;


  if(numoutputcols>1){
    fprintf(stderr,"high2low not setup for numoutputcols=%d>1\n",numoutputcols);
    exit(1);
  }

  // +1's are just so I[1] is first and I[N] exists
  Ihigh=(FTYPE*)malloc(sizeof(FTYPE)*(nxhigh+1)*(nyhigh+1));
  Ilow=(FTYPE*)malloc(sizeof(FTYPE)*(nxlow+1)*(nylow+1));
  IlowW=(FTYPE*)malloc(sizeof(FTYPE)*(nxlow+1)*(nylow+1));
  //if((Ihigh==NULL)||(Ilow==NULL)||(IlowW==NULL)){
  if((Ilow==NULL)||(IlowW==NULL)||(Ihigh==NULL)){
    fprintf(stderr,"Cannot allocate memory\n");
    exit(1);
  }

  for(jh=1;jh<=nyhigh;jh++){
    for(ih=1;ih<=nxhigh;ih++){
      Ihigh[jh*nxhigh+ih]=olddata[0][0][ih][jh][0]; // GODMARK3D // COLIMARK
    }
  }

  for(jl=1;jl<=nylow;jl++){
    for(il=1;il<=nxlow;il++){
      Ilow[jl*nxlow+il]=0;
      IlowW[jl*nxlow+il]=0;
    }
  }

  for(jh=1;jh<=nyhigh;jh++){
    for(ih=1;ih<=nxhigh;ih++) {
      ilfrac=ih*(FTYPE)nxlow/(FTYPE)nxhigh;
      jlfrac=jh*(FTYPE)nylow/(FTYPE)nyhigh;
      ioldbcl=(ceil(ilfrac)-1)*nxhigh/nxlow;
      joldbcl=(ceil(jlfrac)-1)*nyhigh/nylow;
      deltaih=fabs(ioldbcl-(ih-0.5));
      deltajh=fabs(joldbcl-(jh-0.5));
      if(deltaih>1.0) deltaih=1.0;
      if(deltajh>1.0) deltajh=1.0;
      
      W=deltaih*deltajh;
      
      il=(int)ceil(ilfrac);
      jl=(int)ceil(jlfrac);
      
      Ilow[jl*nxlow+il]+=Ihigh[jh*nxhigh+ih]*W;
      IlowW[jl*nxlow+il]+=W;
      /*
        if((il==1)&&(jl==1)){
        printf("%g %g %g %g %g %g %g %d %d %g %g %g\n",ilfrac,jlfrac,ioldbcl,joldbcl,deltaih,deltajh,W,il,jl,Ilow[jl*nxlow+il],IlowW[jl*nxlow+il],Ihigh[jh*nxhigh+ih]*W);
        }
      */
    }
  }
  for(jl=1;jl<=nylow;jl++){
    for(il=1;il<=nxlow;il++){
      // go back to where we came from, just filling the partial array
      olddata[0][0][il][jl][0]=Ilow[jl*nxlow+il]/IlowW[jl*nxlow+il]; // GODMARK3D // COLIMARK
    }
  }

  free(Ihigh);
  free(Ilow);
  free(IlowW);



}
