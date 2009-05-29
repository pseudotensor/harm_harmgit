#include "decs.h"

 

// local declarations of local functions
static void new_xycoord(int i, int j, int k, FTYPE *xc, FTYPE *yc, FTYPE *zc);
static void new_coord(int i,int j, int k, FTYPE *r, FTYPE *th, FTYPE *ph) ;
static int old_coord(FTYPE r,FTYPE th,FTYPE ph,FTYPE *X) ;
static FTYPE old_distance(FTYPE x1, FTYPE y1, FTYPE z1,FTYPE x2, FTYPE y2, FTYPE z2);
static FTYPE new_distances(FTYPE r1, FTYPE th1, FTYPE ph1, FTYPE r2, FTYPE th2, FTYPE ph2);
static void old_xyzcoord(FTYPE r, FTYPE th, FTYPE ph, FTYPE *xc, FTYPE*yc, FTYPE *zc);
static FTYPE plane_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage,FTYPE***olddata);
static FTYPE bilinear_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage, FTYPE***olddata);
static FTYPE bilinear_interp_ij(int iold, int jold,int kold,FTYPE diold,FTYPE djold,FTYPE dkold,unsigned char***oldimage,FTYPE***olddata);
static FTYPE bilinear_interp_simple(int iold, int jold, int kold, FTYPE *X, FTYPE *startx, FTYPE *dx, unsigned char***oldimage,FTYPE***olddata);
static void interpicoord(FTYPE *X,int loc, int *i, int *j, int *k);
static int is_atboundary(int *ioldvar, int *joldvar, int *koldvar, int *atboundary, int *nearboundary);
static FTYPE nearest_interp(int iold,int jold, int kold, unsigned char***oldimage,FTYPE***olddata);
static FTYPE bicubic_interp_wrap(int nx, int ny, int nz, int iold, int jold, int kold, FTYPE x1, FTYPE x2, FTYPE x3, unsigned char***oldimage,FTYPE***olddata);
static void X2spc(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_curr);

static int root_find1(FTYPE r,FTYPE th, FTYPE ph, FTYPE*X) ;
static int root_find2(FTYPE r, FTYPE th, FTYPE ph, FTYPE *X);  

static FTYPE nearest_interp_ij(int iold,int jold,int kold,unsigned char***oldimage,FTYPE***olddata);
static void old_ijk2rthph(int iold, int jold, int kold, FTYPE*r, FTYPE*th, FTYPE *ph);
static void oldf_ijk2rthph(FTYPE iold, FTYPE jold, FTYPE kold, FTYPE*r, FTYPE*th, FTYPE *ph);
static void dervs_for_bicubic(int nx, int ny, FTYPE **ya, FTYPE **fy1a, FTYPE **fy2a, FTYPE **fy12a);
static void oldf_ijk2x123(FTYPE iold, FTYPE jold, FTYPE kold, FTYPE*X);
static void old_ijk2x123(int iold, int jold, int kold, FTYPE*X);
static void old_x1232rthphcoord(FTYPE *X,FTYPE *r, FTYPE *th, FTYPE *ph);

static int root_find_parameter1(FTYPE *parms, FTYPE *returnvalue);
static void resid_parameter1(int n, FTYPE *parms, FTYPE *guess, FTYPE *diff);
static int usrfun_parameter1(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha);
static int usrfun_joninterp(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *beta, FTYPE **alpha, FTYPE*norm);
static int usrfun_joninterp2(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha);


static void low2high(int nxlow, int nylow, int nzlow, int nxhigh, int nyhigh, int nzhigh, unsigned char***oldimage,FTYPE***olddata);
static void high2low(int nxhigh, int nyhigh, int nzhigh, int nxlow, int nylow, int nzlow, unsigned char***oldimage,FTYPE***olddata);

static int lowericoord(FTYPE *X,int loc, int *i, int *j, int *k);







// just copy over data from old -> new
//unsigned char***oldimage,FTYPE***olddata,unsigned char***newimage,FTYPE***newdata
void copy_old2new(void)
{
  int i,j,k;

  if(DATATYPE==0) LOOPINTERP newimage[i][j][k]=oldimage[i][j][k];
  else if(DATATYPE==1) LOOPINTERP newdata[i][j][k]=olddata[i][j][k];


}



// function that actually does interpolation
//unsigned char***oldimage,FTYPE***olddata,unsigned char***newimage,FTYPE***newdata
int compute_spatial_interpolation(void)
{
  //
  // used for image specific details
  int interptypetodo,getinterp;
  //
  int iold,jold,kold ;
  /////////////////////////////
  // coordinate stuff
  int i,j,k;
  int ii,jj,kk;
  //  FTYPE Xmax[NDIM];
  FTYPE r,th,r2,th2,ph,ph2;

  /////////////////////////////
  // special interpolation stuff
  int nosolution;
  FTYPE diold,djold,dkold;
  FTYPE ftemp;
  unsigned char uctemp;

  
  if(DATATYPE==0)   fprintf(stderr,"start to interpolate image\n");
  else if(DATATYPE==1)   fprintf(stderr,"start to interpolate data\n");
  fflush(stderr);
  /* interpolate to new image */

  if(DEBUGINTERP){
    fprintf(stderr,"%g %g :: %g %g :: %g %g\n",startx[1],Xmax[1],startx[2],Xmax[2],startx[3],Xmax[3]) ; fflush(stderr);
  }


#if(0)
  kold=0;
  if(oN3!=1){
    // pick kold (nearest neighbor)
    for(jj=0;jj<oN2;jj++)           for(ii=0;ii<oN1;ii++){
	olddata[ii][jj]=olddata3d[ii][jj][kold];
      }
  }
#endif
  

  LOOPINTERP{

    if(DEBUGINTERP){
      if(j==0) fprintf(stderr,"%d %d %d\n",i,j,k); fflush(stderr);
    }

    ////////////////////////////
    // coordinate solver
    ///////////////////////////

    // find new real coordinates (e.g. real r,th,ph) for i,j,k
    new_coord(i,j,k,&r,&th,&ph) ;

    //HACK
    //    if(ph>0 && ph<M_PI/2.0){
    //      r=100;
    //    }


    // find old grid's X from that new grid's r,th,ph
    // input must be those coordinates (here r,th,ph) corresponding to metric's original coordinates after transforming back using dxdxp from PRIMECOORDS.  Otherwise, not seeking correct solution in root finder that obtains X
    nosolution=old_coord(r,th,ph,X) ;


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
	newimage[i][j][k] = defaultvalue ;
	//fprintf(stderr,"badgot here i=%d j=%d\n",i,j); fflush(stderr);
	getinterp=0;
      }
      else if(
	      (oN1!=1)&&((X[1]-startx[1] < 0) || (X[1]-Xmax[1] >= 0)) ||
	      (oN2!=1)&&((X[2]-startx[2] < 0) || (X[2]-Xmax[2] >= 0)) ||
	      (oN3!=1)&&((X[3]-startx[3] < 0) || (X[3]-Xmax[3] >= 0))
	      ) {

	if(EXTRAPOLATE){
	  getinterp=1;
	  interptypetodo=0;
	}
	else{
	  newimage[i][j][k] = defaultvalue ;
	  //fprintf(stderr,"badgot here i=%d j=%d\n",i,j); fflush(stderr);
	  getinterp=0;
	  interptypetodo=0;
	}
      }
      else{
	getinterp=1;
	interptypetodo=INTERPTYPE;
      }
    }
    else if(DATATYPE==1){
      // make defaultvalue only if no solution
      if(nosolution) {
	newdata[i][j][k]=defaultvalue;
	getinterp=0;
      }
      else if(
	      (oN1!=1)&&((X[1]-startx[1] < 0) || (X[1]-Xmax[1] >= 0)) ||
	      (oN2!=1)&&((X[2]-startx[2] < 0) || (X[2]-Xmax[2] >= 0)) ||
	      (oN3!=1)&&((X[3]-startx[3] < 0) || (X[3]-Xmax[3] >= 0))
	      ){
	if(EXTRAPOLATE){
	  getinterp=1;
	  interptypetodo=0;
	}
	else{
	  newdata[i][j][k]=defaultvalue;
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

    if(DEBUGINTERP){
      fprintf(stderr,"get=%d nosol=%d :: %d %d %d :: %g %g %g :: %g %g %g\n",getinterp,nosolution,i,j,k,r,th,ph,X[1],X[2],X[3]) ; fflush(stderr);
    }


    if(getinterp){
      //      fprintf(stderr,"got here i=%d j=%d\n",i,j); fflush(stderr);

      // get old grid's i,j,k from determined X
      interpicoord(X,CENT,&iold,&jold,&kold);
      // GODMARK: for now only use nearest neighbor in k -> kold

      //      if(DEBUGINTERP){
      //	if(kold==15 || kold==0){
      //	  fprintf(stderr,"interp: %d %d %d\n",iold,jold,kold);
      //	  fflush(stderr);
      //	}
      //      }

      // check if symmetric -- yes, it is
      //      if(j==nN2/3 && (i==nN1/4 || i==3*nN1/4-1)){
      //      	fprintf(stderr,"i=%d j=%d r=%21.15g th=%21.15g X[1]=%21.15g X[2]=%21.15g iold=%d jold=%d\n",i,j,r,th,X[1],X[2],iold,jold);
      //      }

      if(DEBUGINTERP){
	fprintf(stderr,"doing iold=%d jold=%d kold=%d i=%d j=%d k=%d\n",iold,jold,kold,i,j,k); fflush(stderr);
      }

      //////////////////////////////
      // TRANSPLANT
      //////////////////////////////

      if(DEBUGINTERP){
	fprintf(stderr,"Done copy\n"); fflush(stderr);
      }

      if(interptypetodo==0){
	ftemp=nearest_interp(iold,jold,kold,oldimage,olddata);
      }
      else if(interptypetodo==1){
	ftemp=bilinear_interp_simple(iold,jold,kold,X,startx,dx,oldimage,olddata);
	//	ftemp=bilinear_interp(r,th,ph,iold, jold, kold, oldimage,olddata); // not working
      }
      else if(interptypetodo==2){
	ftemp=plane_interp(r, th, ph, iold, jold, kold, oldimage,olddata);
      }
      else if(interptypetodo==3){
	ftemp=bicubic_interp_wrap(oN1,oN2,oN3,iold,jold,kold,X[1],X[2],X[3],oldimage,olddata);
      }
      else{
	fprintf(stderr,"no such interptypetodo=%d\n",interptypetodo);
	exit(1);
      }

      if(DEBUGINTERP){
	fprintf(stderr,"Done interp i=%d j=%d k=%d\n",i,j,k); fflush(stderr);
      }

      // assign (interpolated) old value to new grid position
      //      fprintf(stderr,"%d %d %d\n",i,j,(unsigned char)(ftemp));
      if(DATATYPE==0) newimage[i][j][k] = FLOAT2IMAGE(ftemp);
      else newdata[i][j][k] = (FTYPE)(ftemp) ;

      if(DEBUGINTERP){
	fprintf(stderr,"Done set newdata\n"); fflush(stderr);
      }

      // test symmetry -- yes!
      //      if(j==nN2/3 && (i==nN1/4 || i==3*nN1/4-1)){
      //      	fprintf(stderr,"i=%d j=%d data=%21.15g\n",i,j,newdata[i][j][k]);
      //      }


      // DEBUG
#if(0)
      if((r>95)&&(r<105)&&(fabs(th-M_PI*0.5)<0.1)){
	fprintf(stderr,"%d %d : %g %g : %g %g : %d %d : %g\n",i,j,r,th,X[1],X[2],iold,jold,ftemp);
      }
#endif

      if(DEBUGINTERP){
	if(r>40){
	  fprintf(stderr,"%d %d %d: %g %g %g: %g %g %g: %d %d %d : %g\n",i,j,k,r,th,ph,X[1],X[2],X[3],iold,jold,kold,ftemp);
	}
      }


      /*
	fprintf(stderr,"newim:%d %d %d %d %c\n",i,j,iold,jold,newimage[i][j][k]) ;
      */
    }
    else{
      if(DATATYPE==0)      newimage[i][j][k] = defaultvalue ;
      else newdata[i][j][k]=defaultvalue;
    }// end else

  }// end loop over new zones

  // report that done
  if(DATATYPE==0)   fprintf(stderr,"done interpolating image\n");
  else if(DATATYPE==1)   fprintf(stderr,"done interpolating data\n");
  fflush(stderr);


  return(0);

}





// set dxc,dyc,dyz only used by new_xycoord() and outputted to file as $dx1,2,3
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
  if(newgridtype==GRIDTYPECART){
    dxc = (endxc-startxc)/(FTYPE)nN1 ;
    dyc = (endyc-startyc)/(FTYPE)nN2 ;
    dzc = (endzc-startzc)/(FTYPE)nN3 ;
    //    fprintf(stderr,"%g %g %g: Zin=%g Zout=%g\n",dxc,dyc,dzc,Zin,Zout);
    //    exit(0);
    fakedxc=dxc;
    fakedyc=dyc;
    fakedzc=dzc;
  }
  else if(newgridtype==GRIDTYPELOGLOGCYL){
    fakedxc= (endxc-startxc)/(FTYPE)nN1 ;
    fakedyc= (endyc-startyc)/(FTYPE)nN2 ;
    fakedzc= (endzc-startzc)/(FTYPE)nN3 ;

    fakeRin=Rin;
    // log or split log for each direction
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

    fakedxc= (endxc-startxc)/(FTYPE)nN1 ;
    fakedyc= (endyc-startyc)/(FTYPE)nN2 ;
    fakedzc= (endzc-startzc)/(FTYPE)nN3 ;

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
  else if(newgridtype==GRIDTYPENOCHANGE){
    // here d?c, start?c, faked?c define internal (not real) coordinate system that is used (with bl_coord()) to get real coordinates

    //    dxc = (Rout-Rin)/(FTYPE)nN1 ;
    //    dyc = (M_PI-0)/(FTYPE)nN2 ;
    //    dzc = (2.0*M_PI-0)/(FTYPE)nN3 ;
    dxc=dX[1]*totalsize[1]/nN1;
    dyc=dX[2]*totalsize[2]/nN2;
    dzc=dX[3]*totalsize[3]/nN3;
    startxc=startx[1];
    startyc=startx[2];
    startzc=startx[3]; // was 0
    //	     &t,&totalsize[1],&totalsize[2],&startx[1],&startx[2],&dX[1],&dX[2],&readnstep,&gam,&spin,&R0,&Rin,&Rout,&hslope,&dt,&defcoord,&MBH,&QBH);
    fakedxc=dxc;
    fakedyc=dyc;
    fakedzc=dzc;
    //fprintf(stderr,"%g %g %g\n",dxc,dyc,dzc);
    //exit(0);
  }
  else{
    dualfprintf(fail_file,"No such newgridtype=%d\n",newgridtype);
    myexit(24326264);
  }

  //  N1=nN1;
  //N2=nN2;
  //N3=nN3;
}


// Get new real coordinates (often real Cartesian x,y,z) from new grid's i,j,k
static void new_xycoord(int i, int j, int k, FTYPE *xc, FTYPE *yc, FTYPE *zc)
{
  FTYPE xnew,ynew,znew;
  FTYPE newX1,newX2,newX3;
  FTYPE newX1u,newX1d;
  FTYPE newX2u,newX2d;
  FTYPE newX3u,newX3d;
  FTYPE xt,yt,zt;
  FTYPE X[NDIM], V[NDIM];

  /*
    xnew = (i+0.5)*dxc ;
    // *y = (j-nN2/2+0.5)*dyc ;
    ynew = (nN2/2-(j+0.5))*dyc ;

    *x= sqrt(xnew*xnew + ynew*ynew) ;
    *y = atan2(xnew,ynew) ;	// deliberately interchanged args
    */

  if(newgridtype==GRIDTYPELOGSPC){
    // for this grid, only outer radius corresponds to startxc/endxc, etc., while in between xc,yc,zc are arbitrary
    *xc = startxc+(i+0.5)*dxc ; // varies from startxc to endxc
    *yc = startyc+(j+0.5)*dyc ;
    *zc = startzc+(k+0.5)*dzc ;

    //    *xc = startxc+(endxc-startxc)*(i+0.5)/((FTYPE)nN1); // varies from startxc to endxc
    //    *yc = startyc+(endyc-startyc)*(j+0.5)/((FTYPE)nN2); // varies from startyc to endyc
    //    *zc = startzc+(endzc-startzc)*(k+0.5)/((FTYPE)nN3); // varies from startzc to endzc
  }
  else if(newgridtype==GRIDTYPECART){
    // GODMARK: Assumes input data at CENT
    *xc = startxc+(i+0.5)*dxc ;
    *yc = startyc+(j+0.5)*dyc ;
    *zc = startzc+(k+0.5)*dzc ;
  }
  else if(newgridtype==GRIDTYPELOGLOGCYL){

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
  else if(newgridtype==GRIDTYPENOCHANGE){
    // GODMARK: Assumes input data at CENT

    X[0] = 0.0;
    X[1] = startxc+(i+0.5)*dxc;
    X[2] = startyc+(j+0.5)*dyc ;
    X[3] = startzc+(k+0.5)*dzc ;
    
    // use bl_coord()
    bl_coord(X,V);

    *xc = V[1];
    *yc = V[2];
    *zc = V[3];
  }
  else{
    dualfprintf(fail_file,"No such newgridtype=%d\n",newgridtype);
    myexit(2467346);
  }
  //  fprintf(stderr,"%d %d %g %g : %g %g %g : %g %g %g\n",i,j,*xc,*yc,newX1,log10(fakeRin),(i+0.5)*dxc, pow(newX1,10.0),pow(newX2u,10.0),-pow(newX2d,10.0)); fflush(stderr);
 
}


// Get original-type real coordinates (e.g. for spherical polar : new values of r,th,ph) from new grid's i,j,k (implicitly from xc,yc,zc in Cartesian)
// This must really give those coordinates (here labelled r,th,ph) that correspond to original metric's coordinates that are obtained after transforming from PRIMECOORDS to normal coords using dxdxp matrix
// of course, input i,j,k or result from new_xycoord() can be xc,yc,zc that are manipulated so that new i,j,k doesn't have to correspond to x,y,z but some other x',y',z'
static void new_coord(int i,int j,int k, FTYPE *r,FTYPE *th,FTYPE *ph)
{
  FTYPE xc,yc,zc,Rc ;
  FTYPE rspc;



  // get new real coordiantes for new i,j,k
  new_xycoord(i,j,k,&xc,&yc,&zc);



  if(newgridtype==GRIDTYPELOGSPC){
    // assume original grid was SPC
    // here xc,yc,zc are related to r,th,ph via
    //
    // xc like R-cyl = x'
    // yc like Z-cyl = z'
    // zc like y-direction = y'
    //
    // i.e. xc = x' = log_{10}(1+r) sin\theta cos\phi
    // i.e. zc = y' = log_{10}(1+r) sin\theta sin\phi
    // i.e. yc = z' = log_{10}(1+r) cos\theta
    //
    // tan(\theta) =  x'/z' = xc/yc
    // r = -1 + 10^{\sqrt{x'^2 + y'^2 + z'^2}}
    // tan(\phi) = x'/y' = xc/zc
    //
    // assume startxc,etc. correspond to outermost radial or x,y,z coordinate in Cart coords

    //    *r = -1.0 + pow(sqrt( xc*xc + yc*yc + zc*zc),10.0);
    FTYPE rfake = sqrt(xc*xc + yc*yc + zc*zc);
    *r = gridr0global*(pow(rfake/gridAAglobal,10.0)-1.0); // wrong unless range in startxc=startyc, etc.

    Rc = sqrt(xc*xc+zc*zc); // R' = \sqrt(x'^2 + y'^2)

    *th = atan2(Rc,yc) ;	/* deliberately interchanged args */
    if(*th<0) *th+=M_PI;
    if(*th<0) *th+=M_PI;
    if(*th<0) *th+=M_PI;
    if(*th>M_PI) *th-=M_PI;
    if(*th>M_PI) *th-=M_PI;
    if(*th>M_PI) *th-=M_PI;

    *ph = atan2(xc,zc) ; // Note that since \hat{i}\times\hat{j} = -\hat{k}, then \hat{k}=-\hat{y} so we flip xc,zc = x',y' to get back normal r,\theta,\phi coordinates (or is this just consistent with how we setup vis5d+ files?)
    if(*ph<0) *ph+=2.0*M_PI;
    if(*ph<0) *ph+=2.0*M_PI;
    if(*ph<0) *ph+=2.0*M_PI;
    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
	  

  }
  else{


    if(newgridtype==GRIDTYPENOCHANGE){
      // then doesn't matter what input grid type is
      *r = xc;
      *th = yc;
      *ph = zc;
    }
    else{
      // for Cart:
      // xc like R-cyl
      // yc like Z-cyl
      // zc like y-direction
	
      if(oldgridtype==GRIDTYPECART){ // Cart input
	*r=xc;
	*th=yc;
	*ph=zc;
      }
      else if(oldgridtype==GRIDTYPESPC){ // Spherical polar input (normal)
	if(defcoord!=666){
	  *r = sqrt(xc*xc + yc*yc + zc*zc) ;
	  Rc=sqrt(xc*xc+zc*zc);
	  *th = atan2(Rc,yc) ;	/* deliberately interchanged args */
	  if(*th<0) *th+=M_PI;
	  if(*th<0) *th+=M_PI;
	  if(*th<0) *th+=M_PI;
	  if(*th>M_PI) *th-=M_PI;
	  if(*th>M_PI) *th-=M_PI;
	  if(*th>M_PI) *th-=M_PI;
	  //	    *th = atan(Rc/yc) ;
	  //dualfprintf(fail_file," got here %d %d\n",newgridtype,oN3);
	  if(newgridtype==GRIDTYPECART && oN3==1){
	    *ph = atan2(fabs(xc),fabs(zc));
	    if(*ph<0) *ph+=M_PI;
	    if(*ph>M_PI) *ph-=M_PI;
	    //		dualfprintf(fail_file,"newcoord here ph=%g\n",*ph);
	  }
	  else{
	    *ph = atan2(xc,zc) ;
	    if(*ph<0) *ph+=2.0*M_PI;
	    if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
	  }
	}
	else{
	  *r=xc;
	  *th=yc;
	  *ph=zc;
	}
      }// end if oldgridtype==GRIDTYPESPC
      else if(oldgridtype==GRIDTYPECYL){ // Cylindrical coordinates input

	//	    rspc = sqrt(xc*xc + yc*yc + zc*zc) ;
	Rc=sqrt(xc*xc+zc*zc);
	    
	*r = Rc; // Cyl radius
	*th = yc; // Cyl height
	*ph = atan2(xc,zc) ;	// Cyl angle
	if(*ph<0) *ph+=2.0*M_PI;
	if(*ph>2.0*M_PI) *ph-=2.0*M_PI;
      }
    }
  }



  return ;
}


// This function only used for distances for final interpolation, so not crucial how input and output related, but most correct gives best quality interpolation result
// notice that y and z are switched compared to in Griffiths book final page
// Cartesian x,y,z from spherical polar r,th,ph
// always the same for any coordinates unless setup x,y,z to not be Cartesian elsewhere also
// used for Cartesian distance, so really want Cartesian coordinates!
static void old_xyzcoord(FTYPE r, FTYPE th, FTYPE ph, FTYPE *xc, FTYPE*yc, FTYPE *zc)
{

  if(oldgridtype==GRIDTYPECART){ // Cart input
    *xc = r;
    *yc = th;
    *zc = ph;
  }
  else if(oldgridtype==GRIDTYPESPC){ // spc input
    //  if(defcoord!=666){
    
    if(newgridtype==GRIDTYPENOCHANGE){
      // test for fake distance
      *xc = r;
      *yc = th;
      *zc = ph;
    }
    else{

      if(newgridtype==GRIDTYPELOGSPC){
	// then fake Cartesian generated
	*xc = log(r+1.0)*sin(th)*cos(ph);
	*yc = log(r+1.0)*cos(th); // my z
	*zc = log(r+1.0)*sin(th)*sin(ph); // my y
      }
      else{

	if(newgridtype==GRIDTYPECART && oN3==1){
	}
	else{
	  *xc = fabs(r*sin(th)*cos(ph));
	}
	*yc = r*cos(th); // my z
	*zc = r*sin(th)*sin(ph); // my y
      }
    }
  }
  else if(oldgridtype==GRIDTYPECYL){ // Cyl input

    *xc = r*sin(ph);
    *zc = r*cos(ph);
    *yc = th;
  }
  //  }
  //  else{
  //    *xc = r;
  //    *yc = th;
  //  }

}

// Cartesian distance between 2 spherical polar points
static FTYPE new_distance(FTYPE r1, FTYPE th1, FTYPE ph1, FTYPE r2, FTYPE th2, FTYPE ph2)
{
  FTYPE xc1,yc1,zc1,xc2,yc2,zc2;

  old_xyzcoord(r1,th1,ph1,&xc1,&yc1,&zc1);
  old_xyzcoord(r2,th2,ph2,&xc2,&yc2,&zc2);

  return(old_distance(xc1,yc1,zc1,xc2,yc2,zc2)); // can use same Cartesian function
  //return(sqrt((xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2)+(zc1-zc2)*(zc1-zc2)));
} 

// Cartesian distance between 2 Cartesian points
static FTYPE old_distance(FTYPE xc1, FTYPE yc1, FTYPE zc1, FTYPE xc2, FTYPE yc2, FTYPE zc2)
{
  return(sqrt((xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2)+(zc1-zc2)*(zc1-zc2)));
} 

// Get r,th,ph from old grid's i,j,k
static void old_ijk2rthph(int iold, int jold, int kold, FTYPE*r, FTYPE*th, FTYPE *ph)
{
  FTYPE X[NDIM];

  old_ijk2x123(iold,jold,kold,X);
  old_x1232rthphcoord(X,r,th,ph);

}

// Get r,th,ph from old grid's i,j,k
static void oldf_ijk2rthph(FTYPE iold, FTYPE jold, FTYPE kold, FTYPE*r, FTYPE*th, FTYPE *ph)
{
  FTYPE X[NDIM];


  oldf_ijk2x123(iold,jold,kold,X);
  old_x1232rthphcoord(X,r,th,ph);

}

// get r,th,ph for given x1,x2,x3 for old grid
static void old_x1232rthphcoord(FTYPE *X,FTYPE *r, FTYPE *th, FTYPE *ph)
{
  FTYPE V[NDIM];

  bl_coord(X,V);
  *r=V[1];
  *th=V[2];
  *ph=V[3];

}


// get X1, X2 from i,j for old uniform X-grid 
static void old_ijk2x123(int iold, int jold, int kold, FTYPE*X)
{
  coord(iold,jold,kold,CENT,X);
}

// get X1, X2 from i,j for old uniform X-grid 
static void oldf_ijk2x123(FTYPE iold, FTYPE jold, FTYPE kold, FTYPE*X)
{
  coordf(iold,jold,kold,CENT,X);
}

// get i,j from X1, X2 for old uniform X-grid
static void old_x1232ijk(FTYPE *X,int *iold, int *jold, int *kold)
{
  lowericoord(X,CENT,iold,jold,kold);
}

// find old grid's X for given r,th,ph
static int old_coord(FTYPE r, FTYPE th, FTYPE ph, FTYPE *X)
{

#if(ROOTMETHOD<=1)
  return(root_find1(r,th,ph,X)) ;
#else
  return(root_find2(r,th,ph,X)) ;
#endif

  fprintf(stderr,"shouldn't reach end of old_coord()\n");
  myexit(1);
  return(-1) ;
}










// finds X for a given r,theta
static int root_find1(FTYPE r, FTYPE th, FTYPE ph, FTYPE *X)
{
  int jj;
  int ntrial,mintrial;
  FTYPE tolx,tolf,tolxallowed,tolfallowed,tolxreport,tolfreport;
  FTYPE Xtrial[NDIM];
  int notfoundsolution;

  // target r,th,ph
  spc_target[1]=r;
  spc_target[2]=th;
  spc_target[3]=ph;


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
  Xtrial[1]=X[1];
  Xtrial[2]=X[2];
  Xtrial[3]=X[3];
    
  // iterate to find solution
  nrerrorflag=0;
  notfoundsolution=mnewt(ntrial, mintrial, Xtrial, 3, tolx, tolf, tolxallowed, tolfallowed, tolxreport, tolfreport, NULL, usrfun_joninterp);
  if(notfoundsolution) SLOOPA(jj) fprintf(stderr,"didn't find solution X[%d]=%21.15g\n",jj,X[jj]);
  if(nrerrorflag) SLOOPA(jj) fprintf(stderr,"nrerror: didn't find solution X[%d]=%21.15g\n",jj,X[jj]);
  if(notfoundsolution||nrerrorflag){
    return(1);
  }
  else{
    // assign answer
    X[1]=Xtrial[1];
    X[2]=Xtrial[2];
    X[3]=Xtrial[3];
    return(0);
  }


}


#define NUMATTEMPTS 100

// finds X for a given r,theta
static int root_find2(FTYPE r, FTYPE th, FTYPE ph, FTYPE *X)
{
  int ntrial,mintrial;
  FTYPE tolx,tolf,tolxallowed,tolfallowed,tolxreport,tolfreport;
  FTYPE Xtrial[NDIM];
  int notfoundsolution;
  int attempt;

  // target r,th
  spc_target[1]=r;
  spc_target[2]=th;
  spc_target[3]=ph;

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
    Xtrial[1]=X[1];
    Xtrial[2]=X[2];
    Xtrial[3]=X[3];
    
    // iterate to find solution (only need X2spc, newt computes Jacobian)
    //          newt(ANALYTICJAC,root guess, num dimension, &check,vecfun,usrfun)
    nrerrorflag=0;
#if(ROOTMETHOD==2)
    newt  (ANALYTICJAC,NULL,Xtrial, 3, &notfoundsolution,X2spc,usrfun_joninterp2);
#elif(ROOTMETHOD==3)
    broydn(ANALYTICJAC,NULL,Xtrial, 3, &notfoundsolution,X2spc,usrfun_joninterp2);
#endif

    // quick message, and try again
    if(SIMPLEDEBUGINTERP){
      if(notfoundsolution) fprintf(stderr,"ns  X[%d]=%21.15g\n",jj,X[jj]);
      if(nrerrorflag) fprintf(stderr,"nr  X[%d]=%21.15g\n",jj,X[jj]);
    }

    if((notfoundsolution==0)&&(nrerrorflag==0)) break; // otherwise try until NUMATTEMPTS attempts have been tried
    else{
      // try a new starting position, which is ultimately why failed typically
      coord((int)((oN1/2-1)*ranc(0,0)),(int)((oN2/2-1)*ranc(0,0)),(int)((oN3/2-1)*ranc(0,0)),CENT,X);
    }
  }

  if(notfoundsolution||nrerrorflag){
    if(SIMPLEDEBUGINTERP){
      fprintf(stderr,"Couldn't find solution after %d attempts\n",attempt); fflush(stderr);
    }
    return(1);
  }
  else{
    // assign answer
    X[1]=Xtrial[1];
    X[2]=Xtrial[2];
    X[3]=Xtrial[3];
    return(0);
  }

}

// get spc coordinates (r,th) minus target (r,th) for a given X
static void X2spc(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff)
{
  FTYPE X[NDIM],spc_curr[NDIM];
  int i;

  for(i=0;i<=1;i++){
    if(i==1){
      // this avoids out of bounds attempts by iterative method
      if(Xguess[1]-startx[1] < 0) Xguess[1]=startx[1];
      if(Xguess[1]-Xmax[1] >= 0) Xguess[1]=Xmax[1];
      if(Xguess[2]-startx[2] < 0) Xguess[2]=startx[2];
      if(Xguess[2]-Xmax[2] >= 0) Xguess[2]=Xmax[2];

      if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
	if(Xguess[3]-startx[3] < 0) Xguess[3]=Xguess[3]+(Xmax[3]-startx[3]);
	if(Xguess[3]-Xmax[3]  >= 0) Xguess[3]=Xguess[3]-(Xmax[3]-startx[3]);
      }
      else{
	// then normal non-periodic
	if(Xguess[3]-startx[3] < 0) Xguess[3]=startx[3];
	if(Xguess[3]-Xmax[3] >= 0) Xguess[3]=Xmax[3];
      }
    }


    // X and Xguess different rank
    X[0]=0;
    X[1]=Xguess[1];
    X[2]=Xguess[2];
    X[3]=Xguess[3];
    // find current r,th,ph from current X
    old_x1232rthphcoord(X,&spc_curr[1],&spc_curr[2],&spc_curr[3]);

    // spc_diff should contain function to be zeroed
    spc_diff[1]=spc_curr[1]-spc_target[1];
    spc_diff[2]=spc_curr[2]-spc_target[2];
    spc_diff[3]=spc_curr[3]-spc_target[3];

    if((spc_diff[1]>BIG)||(spc_diff[2]>BIG)||(spc_diff[3]>BIG)){
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
 
 
// get spc coordinates (r,th) minus target (r,th) for a given X
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

  // X and Xguess different rank
  X[0]=0;
  X[1]=Xguess[1];
  X[2]=Xguess[2];
  X[3]=Xguess[3];

  // find current r,th from current X
  old_x1232rthphcoord(X,&spc_curr[1],&spc_curr[2],&spc_curr[3]);

  // spc_diff should contain function to be zeroed
  spc_diff[1]=spc_curr[1]-spc_target[1];
  spc_diff[2]=spc_curr[2]-spc_target[2];
  spc_diff[3]=spc_curr[3]-spc_target[3];


  //////////////////
  //
  // calculate dxdxp
  //  dxdxprim(X,spc_curr[1],spc_curr[2],dxdxp);
  V[0]=0;
  V[1]=spc_curr[1];
  V[2]=spc_curr[2];
  V[3]=spc_curr[3];

  dxdxprim(X,V,dxdxp);
  // assign to alpha (didn't use alpha directly since rank of alpha is smaller than dxdxp)
  for (j = 0; j < n; j++){
    for (k = 0; k < n; k++){
      alpha[j+1][k+1]=dxdxp[j+1][k+1];
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
  int i = 0, j = 0, k = 0;
  FTYPE r,th,ph,spc_curr[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE X[NDIM];
  FTYPE V[NDIM];
  int numnormterms;

  // X and Xguess different rank
  X[0]=0;
  X[1]=Xguess[1];
  X[2]=Xguess[2];
  X[3]=Xguess[3];
  // find current r,th from current X
  old_x1232rthphcoord(X,&r,&th,&ph);

  V[1]=spc_curr[1]=r;
  V[2]=spc_curr[2]=th;
  V[3]=spc_curr[3]=ph;

  //////////////////
  //
  // calculate dxdxp
  //  dxdxprim(X,r,th,dxdxp);
  dxdxprim(X,V,dxdxp);
  // assign to alpha (didn't use alpha directly since rank of alpha is smaller than dxdxp)
  for (j = 0; j < n; j++){
    for (k = 0; k < n; k++){
      alpha[j+1][k+1]=dxdxp[j+1][k+1];
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
    beta[k + 1] = (spc_curr[k+1] - spc_target[k+1]) *(*norm);
  
  
  return (0);
}






//////////////////////////////
//
//////////////////// NEXT SEVERAL FUNCTIONS RELATED TO INTERPOLATION, SUPERSAMPLING, OR FILTERING
//
///////////////////////////////

//#define NUMIJK (4) // 2D
#define NUMIJK (8) // 3D

#define NUMIJKMEM (NUMIJK+1)

// rref and thref are true grid location's r and theta
// x1,x2,x3 are 3 points x position
// y1,y2,y3 are 3 points y position
// z1,z2,z3 are 3 points values
// quad_interp is fully 3D capable
static int quad_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage,FTYPE***olddata, FTYPE*dist, FTYPE*myx, FTYPE*myy, FTYPE*myz, FTYPE*myfun,FTYPE*xref,FTYPE*yref,FTYPE *zref)
{
  FTYPE x1,  x2,  x3,  y1,  y2,  y3,  z1,  z2,  z3;
  int p;
  FTYPE myr[NUMIJKMEM],myth[NUMIJKMEM],myph[NUMIJKMEM];
  //
  //
  int ioldp,joldp,koldp,isel[NUMIJKMEM],jsel[NUMIJKMEM],ksel[NUMIJKMEM];
  int which,whichc;
  FTYPE interpolate;
  int atboundary;
  FTYPE ftemp;



  ioldp = iold+1 ;
  joldp = jold+1 ;
  koldp = kold+1 ;

  //  atboundary=0;
  // if refining, then oN1 is really old image size, not refined, which is correct

  
  /*
  //take care of boundary effects
  if(iold>=oN1){ iold=oN1-1;atboundary=1;}
  if(iold<0){ iold=0;atboundary=1;}
  if(jold>=oN2){ jold=oN2-1;atboundary=1;}
  if(jold<0){ jold=0;atboundary=1;}
  
  if(ioldp>=oN1){ ioldp=oN1-1;atboundary=1;}
  if(ioldp<0){ ioldp=0;atboundary=1;}
  if(joldp>=oN2){ joldp=oN2-1;atboundary=1;}
  if(joldp<0){ joldp=0;    atboundary=1;}

  if(atboundary){
  return(1);
  }
  */

  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    if((iold>=oN1)||(iold<0)||(jold>=oN2)||(jold<0)||(kold<0)){
      fprintf(stderr,"shouldn't reach here: icoord() should truncate iold/jold/kold\n");
      myexit(1);
    }
  }
  else{
    if((iold>=oN1)||(iold<0)||(jold>=oN2)||(jold<0)||(kold>=oN3)||(kold<0)){
      fprintf(stderr,"shouldn't reach here: icoord() should truncate iold/jold/kold\n");
      myexit(1);
    }
  }

  if((ioldp<0)||(joldp<0)||(koldp<0)){
    fprintf(stderr,"shouldn't reach here: ioldp/joldp/koldp can't be less than truncated iold/jold/kold which are at minimum=0\n");
    myexit(1);
  }

  // only case that has to be corrected
  if(ioldp>=oN1){ ioldp=oN1-1; iold=ioldp-1;}
  if(joldp>=oN2){ joldp=oN2-1; jold=joldp-1;}
  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    // koldp can be up to =oN3
    if(koldp>=oN3+1){
      koldp=koldp-oN3; kold=koldp-1;
      phref=phref-(Xmax[3]-startx[3]);
    }

  }
  else{
    if(koldp>=oN3){ koldp=oN3-1; kold=koldp-1;}
  }


  // otherwise do good interpolation
  
  // get distances and function values
  isel[1]=iold;  jsel[1]=jold;  ksel[1]=kold;
  isel[2]=ioldp; jsel[2]=jold;  ksel[2]=kold;
  isel[3]=iold;  jsel[3]=joldp; ksel[3]=kold;
  isel[4]=ioldp; jsel[4]=joldp; ksel[4]=kold;

  if(NUMIJK==8){
    // these only used in 3D
    isel[5]=iold;  jsel[5]=jold;  ksel[5]=koldp;
    isel[6]=ioldp; jsel[6]=jold;  ksel[6]=koldp;
    isel[7]=iold;  jsel[7]=joldp; ksel[7]=koldp;
    isel[8]=ioldp; jsel[8]=joldp; ksel[8]=koldp;
  }

  old_xyzcoord(rref,thref,phref,xref,yref,zref);

  if(DEBUGINTERP)   fprintf(stderr,"rref=%g thref=%g phref=%g xref=%g yref=%g zref=%g\n",rref,thref,phref,*xref,*yref,*zref);

  for(p=1;p<=NUMIJK;p++){
    old_ijk2rthph(isel[p],jsel[p],ksel[p],&myr[p],&myth[p],&myph[p]) ;
    //    new_coord(isel[p],jsel[p],&myr[p],&myth[p]) ;
    old_xyzcoord(myr[p],myth[p],myph[p],&myx[p],&myy[p],&myz[p]);
    //dist[p]=new_distance(myr[p],myth[p],r,th);
    dist[p]=old_distance(myx[p],myy[p],myz[p],*xref,*yref,*zref);
    if(DATATYPE==0){
      myfun[p]=(FTYPE)((unsigned char)oldimage[isel[p]][jsel[p]][ksel[p]]);
      if(DEBUGINTERP)       fprintf(stderr,"p=%d isel%d jsel=%d myr=%g myth=%g myx=%g myy=%g myz=%g\n",p,isel[p],jsel[p],myr[p],myth[p],myx[p],myy[p],myfun[p]); fflush(stderr);
    }
    else{
      myfun[p]=olddata[isel[p]][jsel[p]][ksel[p]];
    }
  }

  // feeds back dist, myx, myy, myz, myfun
  return(0);

}


// rref and thref are true grid location's r and theta
// x1,x2,x3 are 3 points x position
// y1,y2,y3 are 3 points y position
// z1,z2,z3 are 3 points values
// plane_interp is NOT fully 3D capable because assumes only using x-y space and so essentially nearest neighbor in k-space
static FTYPE plane_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage,FTYPE***olddata)
{
  FTYPE x1,  x2,  x3,  y1,  y2,  y3,  z1,  z2,  z3;
  FTYPE f1,f2,f3;
  FTYPE a3,b3,c3,aoc,boc,newvalue;
  int p;
  FTYPE myr[NUMIJKMEM],myth[NUMIJKMEM],myph[NUMIJKMEM],myfun[NUMIJKMEM],dist[NUMIJKMEM],myx[NUMIJKMEM],myy[NUMIJKMEM],myz[NUMIJKMEM];
  FTYPE xref,yref,zref;
  int quad_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage,FTYPE***olddata, FTYPE*dist, FTYPE*myx, FTYPE*myy, FTYPE*myz, FTYPE*myfun,FTYPE*xref,FTYPE*yref,FTYPE *zref);

  FTYPE nearest_interp(int iold,int jold,int kold,unsigned char***oldimage,FTYPE***olddata);
  int which,whichc;
  int a,b,c;
  FTYPE interpolate;
  FTYPE ftemp;


  if(quad_interp(rref, thref, phref, iold, jold, kold, oldimage,olddata, dist, myx, myy, myz, myfun, &xref,&yref,&zref)>=1){
    ftemp=nearest_interp(iold,jold,kold,oldimage,olddata);
    return(ftemp);
  }

  //////////////////////
  //
  // rest is specific to plane interpolation
  //
  // rest doesn't care about iold/jold/kold, just actual positions and function value
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
  
  if(DEBUGINTERP)   fprintf(stderr,"iold=%d jold=%d kold=%d which=%d whichc=%d\n",iold,jold,kold,which,whichc);


  // GODMARK3D
  if(which==1){         a=2; b=3; c=4;}
  else if(which==2){    a=1; b=3; c=4;}
  else if(which==3){    a=1; b=2; c=4;}
  else if(which==4){    a=1; b=2; c=3;}

  //a=whichc;b=2;c=3;
  //  a=1;b=2;c=3;

  // GODMARK3D
  x1=myx[a]; y1=myy[a]; f1=myfun[a];
  x2=myx[b]; y2=myy[b]; f2=myfun[b];
  x3=myx[c]; y3=myy[c]; f3=myfun[c];

  if(1){
    // GODMARK3D
    a3=y2* f1 - y3* f1 - y1* f2 + y3* f2 + y1* f3 - y2* f3;
    b3=-(x2* f1 - x3* f1 - x1* f2 + x3* f2 + x1* f3 - x2* f3);
    c3=x2* y1 - x3* y1 - x1* y2 + x3* y2 + x1* y3 - x2* y3;
    aoc=a3/c3;
    boc=b3/c3;
    
    if(DEBUGINTERP)     fprintf(stderr,"%d %d %d :: %g %g %g %g %g\n",iold,jold,kold,a3,b3,c3,aoc,boc); fflush(stderr);
    
    newvalue=f1;
    interpolate=-aoc*(xref-x1)-boc*(yref-y1); // references off point (x1,y1) with value f1
    newvalue+=interpolate;
  }
  else{
    newvalue=((y3 - yref)*(x2*f1 - x1*f2) + x3*(-(y2*f1) + yref*f1 + y1*f2 - yref*f2) + 
	      (-(x2*y1) + x1*y2 - x1*yref + x2*yref)*f3 + xref*(y2*f1 - y3*f1 - y1*f2 + y3*f2 + y1*f3 - y2*f3))/(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
  }
  if(DEBUGINTERP)   fprintf(stderr,"%d %d %g %g %g\n",iold,jold,newvalue,interpolate,f1); fflush(stderr);

  if(DATATYPE==0){
    // note that for images the function value may go one image value up
    // or down on a whim.  For a given pallette file (say john.pal) the
    // bottom then can have sharp colors changes to black.  Let's up
    // everything by 0.5 to avoid this.
    newvalue+=0.5;
    if(newvalue>255.0) newvalue=255;
    else if(newvalue<0.0) newvalue=0;
  }
 
  return(newvalue);
}



static FTYPE bilinear_interp_simple(int iold, int jold, int kold, FTYPE *X, FTYPE *startx, FTYPE *dx, unsigned char***oldimage,FTYPE***olddata)
{
  FTYPE ftemp;
  int atboundary[NDIM],nearboundary[NDIM];
  FTYPE diold,djold,dkold;
  FTYPE nearest_interp_ij(int iold,int jold,int kold,unsigned char***oldimage,FTYPE***olddata);
  FTYPE bilinear_interp_ij(int iold, int jold,int kold,FTYPE diold,FTYPE djold,FTYPE dkold,unsigned char***oldimage,FTYPE***olddata);


  is_atboundary(&iold,&jold,&kold,atboundary,nearboundary);

  if(BILINEARREDUCE2NEARESTATBOUNDARY && (nearboundary[1]||nearboundary[2]||nearboundary[3]) ){
    // then just use nearest neighbor
    ftemp=nearest_interp_ij(iold,jold,kold,oldimage,olddata);
    return(ftemp);
  }
  else{

    // then need to shift down since going to use upper point inside interp
    // don't need to shift up since already shifted inside atboundary
    if(nearboundary[1]==1) iold--;
    if(nearboundary[2]==1) jold--;
    if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
      // then no need to fix since boundary condtion data is at upper k
    }
    else{
      if(nearboundary[3]==1) kold--;
    }

    // below valid at boundary
    diold = ((X[1]-startx[1])/dx[1] - 0.5) - iold;
    djold = ((X[2]-startx[2])/dx[2] - 0.5) - jold;
    dkold = ((X[3]-startx[3])/dx[3] - 0.5) - kold;

    ftemp=bilinear_interp_ij(iold,jold,kold,diold,djold,dkold,oldimage,olddata);
    return(ftemp);

  }
  

}


int is_atboundary(int *ioldvar, int *joldvar, int *koldvar, int *atboundary, int *nearboundary)
{
  int iold, jold, kold;
  int ioldp,joldp,koldp;
  int jj;

  iold=*ioldvar;
  jold=*joldvar;
  kold=*koldvar;

  SLOOPA(jj){
    nearboundary[jj]=0;
    atboundary[jj]=0;
  }
  // if refining, then oN1 is really old image size, not refined, which is correct
  /* take care of boundary effects */
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


#if(0)
  // this shouldn't be necessary  
  ioldp = iold+1 ;
  joldp = jold+1 ;
  koldp = kold+1 ;
  
  if(ioldp>=oN1){ ioldp=oN1-1;atboundary[1]=1;}
  if(ioldp<0){ ioldp=0;atboundary[1]=1;}
  if(joldp>=oN2){ joldp=oN2-1;atboundary[2]=1;}
  if(joldp<0){ joldp=0;    atboundary[2]=1;}
  if(koldp>=oN3){ koldp=oN3-1;atboundary[3]=1;}
  if(koldp<0){ koldp=0;    atboundary[3]=1;}
#endif


  *ioldvar=iold;
  *joldvar=jold;
  *koldvar=kold;


  return(0);

}


static FTYPE bilinear_interp_ij(int iold, int jold,int kold,FTYPE diold,FTYPE djold,FTYPE dkold,unsigned char***oldimage,FTYPE***olddata)
{
  int ioldp,joldp,koldp;
  FTYPE newvalue;
  FTYPE myfun[NUMIJKMEM],dist[NUMIJKMEM];
  int isel[NUMIJKMEM],jsel[NUMIJKMEM],ksel[NUMIJKMEM];
  int p;
  FTYPE ftemp;
  FTYPE totaldist;

  ioldp = iold+1 ;
  joldp = jold+1 ;
  koldp = kold+1 ;

  // did have boundary check here

  // get locations
  isel[1]=iold;  jsel[1]=jold;  ksel[1]=kold;
  isel[2]=ioldp; jsel[2]=jold;  ksel[2]=kold;
  isel[3]=iold;  jsel[3]=joldp; ksel[3]=kold;
  isel[4]=ioldp; jsel[4]=joldp; ksel[4]=kold;

  if(NUMIJK==8){
    // these only used in 3D
    isel[5]=iold;  jsel[5]=jold;  ksel[5]=koldp;
    isel[6]=ioldp; jsel[6]=jold;  ksel[6]=koldp;
    isel[7]=iold;  jsel[7]=joldp; ksel[7]=koldp;
    isel[8]=ioldp; jsel[8]=joldp; ksel[8]=koldp;
  }

  // get distances
  dist[1]=(1. - diold)*(1.-djold)*(1.-dkold);
  dist[2]=(diold)*(1. - djold)*(1.-dkold);
  dist[3]=(1. - diold)*(djold)*(1.-dkold);
  dist[4]=(diold)*(djold)*(1.-dkold);

  if(NUMIJK==8){
    dist[5]=(1. - diold)*(1.-djold)*(dkold);
    dist[6]=(diold)*(1. - djold)*(dkold);
    dist[7]=(1. - diold)*(djold)*(dkold);
    dist[8]=(diold)*(djold)*(dkold);
  }


  //  fprintf(stderr,"bi1: %g %g\n",diold,djold);

  for(p=1;p<=NUMIJK;p++){
    if(DATATYPE==0){
      myfun[p]=(double)((unsigned char)oldimage[isel[p]][jsel[p]][ksel[p]]);
    }
    else{
      myfun[p]=olddata[isel[p]][jsel[p]][ksel[p]];
    }
    //    fprintf(stderr,"bi: p=%d dist=%10.5g myfun=%10.5g\n",p,dist[p],myfun[p]);
  }
  newvalue=0.0;
  totaldist=0.0;
  for(p=1;p<=NUMIJK;p++){
    newvalue+=dist[p]*myfun[p];
    totaldist+=dist[p];
  }

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

  return(newvalue);
}


static FTYPE nearest_interp_ij(int iold,int jold,int kold,unsigned char***oldimage,FTYPE***olddata)
{
  FTYPE newvalue;
  int atboundary;

  atboundary=0;
  // if refining, then oN1 is really old image size, not refined, which is correct
  /* take care of boundary effects */
  if(iold>=oN1){ iold=oN1-1;atboundary=1;}
  if(iold<0){ iold=0;atboundary=1;}
  if(jold>=oN2){ jold=oN2-1;atboundary=1;}
  if(jold<0){ jold=0;atboundary=1;}
  if(PERIODICINPHI && oN3>1 && oldgridtype==GRIDTYPESPC){
    if(kold<0){ kold=kold+oN3;}
    if(kold>=oN3){ kold=kold-oN3;}
  }
  else{
    if(kold<0){ kold=0;atboundary=1;}
    if(kold>=oN3){ kold=oN3-1;atboundary=1;}
  }

#if(DEBUGINTERP)
  fprintf(stderr,"Problem?: iold=%d jold=%d kold=%d oN1=%d oN2=%d oN3=%d :: %21.15g\n",iold,jold,kold,oN1,oN2,oN3,olddata[iold][jold][kold]);
#endif


  if(DATATYPE==0)			newvalue=(FTYPE)oldimage[iold][jold][kold] ;
  else 			newvalue=(FTYPE)olddata[iold][jold][kold] ;
  return(newvalue);
}



// don't use, not working
static FTYPE bilinear_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage,FTYPE***olddata)
{
  FTYPE f1,f2;
  int ioldp,joldp,koldp;
  FTYPE newvalue;
  FTYPE dist[NUMIJKMEM],myx[NUMIJKMEM],myy[NUMIJKMEM],myz[NUMIJKMEM],myfun[NUMIJKMEM];
  FTYPE xref,yref,zref;
  int isel[NUMIJKMEM],jsel[NUMIJKMEM],ksel[NUMIJKMEM];
  int p;
  FTYPE ftemp;
  FTYPE nearest_interp(int iold,int jold,int kold,unsigned char***oldimage,FTYPE***olddata);
  int quad_interp(FTYPE rref, FTYPE thref, FTYPE phref, int iold, int jold, int kold, unsigned char***oldimage,FTYPE***olddata, FTYPE*dist, FTYPE*myx, FTYPE*myy, FTYPE*myz, FTYPE*myfun,FTYPE*xref,FTYPE*yref,FTYPE *zref);





  if(quad_interp(rref, thref, phref, iold, jold, kold, oldimage,olddata, dist, myx, myy, myz,myfun,&xref,&yref,&zref)>=1){
    ftemp=nearest_interp(iold,jold,kold,oldimage,olddata);
    return(ftemp);
  }


  // now assign actual weighted value, bilinear filtering
  newvalue=0.0;
  // first 3 are weights by true distance
  if(0){
    // =
    // half of _ and the other half of -
    // 1-2 and 3-4
    p=1; newvalue =0.5*(1.0-dist[p]/(dist[1]+dist[2]))*myfun[p];
    p=2; newvalue+=0.5*(1.0-dist[p]/(dist[1]+dist[2]))*myfun[p];
    p=3; newvalue+=0.5*(1.0-dist[p]/(dist[3]+dist[4]))*myfun[p];
    p=4; newvalue+=0.5*(1.0-dist[p]/(dist[3]+dist[4]))*myfun[p];
  }
  else if(0){
    // || // 1-3 and 2-4
    p=1; newvalue =0.5*(1.0-dist[p]/(dist[1]+dist[3]))*myfun[p];
    p=3; newvalue+=0.5*(1.0-dist[p]/(dist[1]+dist[3]))*myfun[p];
    p=2; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[4]))*myfun[p];
    p=4; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[4]))*myfun[p];
  }
  else if(0){
    // X // 1-4 and 2-3
    p=1; newvalue =0.5*(1.0-dist[p]/(dist[1]+dist[4]))*myfun[p];
    p=4; newvalue+=0.5*(1.0-dist[p]/(dist[1]+dist[4]))*myfun[p];
    p=2; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[3]))*myfun[p];
    p=3; newvalue+=0.5*(1.0-dist[p]/(dist[2]+dist[3]))*myfun[p];
  }
  else if(1){
    // weight by x and y
    p=1; newvalue =(1.0-fabs(xref-myx[1])/fabs(myx[2]-myx[1]))*(1.0-fabs(yref-myy[1])/fabs(myy[2]-myy[1]))*myfun[p];
    p=2; newvalue +=(1.0-fabs(myx[2]-xref)/fabs(myx[2]-myx[1]))*(1.0-fabs(myy[2]-yref)/fabs(myy[2]-myy[1]))*myfun[p];
    p=3; newvalue +=(1.0-fabs(xref-myx[3])/fabs(myx[3]-myx[4]))*(1.0-fabs(yref-myy[3])/fabs(myy[3]-myy[4]))*myfun[p];
    p=4; newvalue +=(1.0-fabs(myx[4]-xref)/fabs(myx[3]-myx[4]))*(1.0-fabs(myy[4]-yref)/fabs(myy[3]-myy[4]))*myfun[p];
  }
  else if(0){
    if(myx[2]>myx[1]){
      f1=(xref-myx[1])*(myfun[2]-myfun[1])/(myx[2]-myx[1])+myfun[1];
    }
    else{
      f1=(xref-myx[2])*(myfun[1]-myfun[2])/(myx[1]-myx[2])+myfun[2];
    }
    if(myx[4]>myx[3]){
      f2=(xref-myx[3])*(myfun[4]-myfun[3])/(myx[4]-myx[3])+myfun[3];
    }
    else{
      f2=(xref-myx[4])*(myfun[3]-myfun[4])/(myx[3]-myx[4])+myfun[4];
    }
    if(myy[3]>myy[1]){ // 1 associated with f1, 3 associated with f2
      newvalue=(yref-myy[1])*(f2-f1)/(myy[3]-myy[1])+f1;
    }
    else{
      newvalue=(yref-myy[3])*(f1-f2)/(myy[1]-myy[3])+f2;
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

  return(newvalue);
}





static FTYPE nearest_interp(int iold,int jold,int kold,unsigned char***oldimage,FTYPE***olddata)
{
  return(nearest_interp_ij(iold,jold,kold,oldimage,olddata));
}



//	ftemp=bicubic_interp_wrap(oN1,oN2,iold,jold,X[1],X[2],oldimage,olddata);
static FTYPE bicubic_interp_wrap(int nx, int ny, int nz, int iold, int jold, int kold,FTYPE x1, FTYPE x2,FTYPE x3,unsigned char***oldimage,FTYPE***olddata)
{
  int j,k;
  static int firsttime=1;
  static int lastk;
  int recomputederivs;
  FTYPE Xget[NDIM];

  static FTYPE **y1a,**y2a,**y12a;
  static FTYPE **ya;
  void dervs_for_bicubic(int nx, int ny, FTYPE **ya, FTYPE **y1a, FTYPE **y2a, FTYPE **y12a);
  FTYPE bicubic_interp(int nx,int ny,int il,int jl,int kl, FTYPE *Xget,FTYPE **ya,FTYPE **y1a, FTYPE **y2a, FTYPE**y12a,FTYPE dx1,FTYPE dx2);


  if(firsttime){
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

  
  // always copy values
  if(DATATYPE==0){
    for(j=0;j<nx;j++) for(k=0;k<ny;k++){
	ya[j][k]=(FTYPE)oldimage[j][k][kold]; // GODMARK3D -- nearest neighbor in k-direction
      }
  }
  // else if(DATATYPE==1) ya=olddata;
  else{
    for(j=0;j<nx;j++) for(k=0;k<ny;k++){
	ya[j][k]=olddata[j][k][kold]; // GODMARK3D -- nearest neighbor in k-direction
      }
  }


  if(firsttime || lastk!=kold){
    // GODMARK3D -- for now redo derivs for each new kold
    recomputederivs=1;
  }
  else recomputederivs=0;


  if(recomputederivs){
    lastk=kold; // store kold into lastk for which recomputed derivatives


    //  fprintf(stderr,"got here2\n"); fflush(stderr);


    // compute derivatives of low resolution data
    dervs_for_bicubic(nx,ny,ya,y1a,y2a,y12a);
  }
  


  Xget[0]=0;
  Xget[1]=x1;
  Xget[2]=x2;
  Xget[3]=x3;
 

  // no longer using firsttime after below line
  if(firsttime) firsttime=0;


 
  return(bicubic_interp(nx,ny,iold,jold,kold,Xget,ya,y1a, y2a, y12a,dx[1],dx[2]));



}




// compute derivatives of data
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
FTYPE bicubic_interp(int nx,int ny,int il,int jl,int kl, FTYPE *Xget,FTYPE **ya,FTYPE **y1a, FTYPE **y2a, FTYPE**y12a,FTYPE dx1,FTYPE dx2)
{
  FTYPE Xp[4+1][NDIM];
  int pointsi[4+1],pointsj[4+1],pointsk[4+1];
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
  pointsk[1]=pointsk[2]=pointsk[3]=pointsk[4]=kl;

  // get surrounding points positions.
  for(i=1;i<=4;i++){  // GODMARK3D
    //    fprintf(stdout,"%d\n",ya);
    oldf_ijk2x123(pointsi[i], pointsj[i], pointsk[i], Xp[i]);
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



static void interpicoord(FTYPE *X,int loc, int *i, int *j, int *k)
{

  
  lowericoord(X,loc,i,j,k);

  // restrict to old data size
  if(*i<0) *i=0;
  if(*j<0) *j=0;
  if(*i>oN1-1) *i=oN1-1;
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
static int lowericoord(FTYPE *X,int loc, int *i, int *j, int *k)
{

  if(loc == CENT){
    *i = ROUND2INT((X[1]-startx[1])/dx[1] + 0.0) - 1 ;
    *j = ROUND2INT((X[2]-startx[2])/dx[2] + 0.0) - 1 ;
    *k = ROUND2INT((X[3]-startx[3])/dx[3] + 0.0) - 1 ;
  }


  return(0);

}









// refine data
void refine_data(void)
{
  int i,j,k;

  //////////////////////////
  //
  // check for refinement of old data
  //
  didrefine=0;
  if(fabs(refinefactor-1.0)>0.1){
    didrefine=1;
    fprintf(stderr,"refine\n"); fflush(stderr);
    // then refine or derefine
    if(oN1!=1) roN1=(int)(refinefactor*(FTYPE)oN1);
    else roN1=oN1;
    if(oN2!=1) roN2=(int)(refinefactor*(FTYPE)oN2);
    else roN2=oN2;
    if(oN3!=1) roN3=(int)(refinefactor*(FTYPE)oN3);
    else roN3=oN3;
    fprintf(stderr,"refine from %dX%dX%d to %dX%dX%d\n",oN1,oN2,oN3,roN1,roN2,roN3); fflush(stderr);
    if(refinefactor>1.0){
      if(DATATYPE==0){
	oldimage = c3matrix(0,roN1-1,0,roN2-1,0,roN3-1) ;
	for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)      for(i=0;i<oN1;i++) oldimage[i][j][k]=oldimage0[i][j][k];
      }
      else{
	// debug below 1 line
	//	if(DOUBLEWORK)	oldimage = cmatrix(0,roN1-1,0,roN2-1) ; 
	olddata = f3matrix(0,roN1-1,0,roN2-1,0,roN3-1) ;
	for(k=0;k<oN3;k++) for(j=0;j<oN2;j++)      for(i=0;i<oN1;i++) olddata[i][j][k]=olddata0[i][j][k];
	//	for(j=0;j<oN2;j++)      for(i=0;i<oN1;i++)	fprintf(stderr,"%d %d %g\n",j,i,olddata[i][j]);
	//	fprintf(stderr,"done refine copy\n"); fflush(stderr);
	//if(DOUBLEWORK) for(j=0;j<oN2;j++)      for(i=0;i<oN1;i++) oldimage[i][j]=(unsigned char)olddata0[i][j];
      }
      //writeimage("jontest1.r8",oldimage,roN1,roN2,roN3);
      //      exit(0);
      low2high(oN1, oN2, oN3, roN1, roN2, roN3, oldimage,olddata);
      //      for(j=0;j<roN2;j++)      for(i=0;i<roN1;i++)	fprintf(stderr,"%d %d %g\n",j,i,olddata[i][j]);
    }
    else{
      high2low(oN1, oN2, oN3, roN1, roN2, roN3, oldimage0,olddata0);
      if(DATATYPE==0){
	oldimage = c3matrix(0,roN1-1,0,roN2-1,0,roN3-1) ;
	for(k=0;k<roN3;k++) for(j=0;j<roN2;j++)      for(i=0;i<roN1;i++) oldimage[i][j][k]=oldimage0[i][j][k];
      }
      else{
	olddata = f3matrix(0,roN1-1,0,roN2-1,0,roN3-1) ;
	for(k=0;k<roN3;k++)  for(j=0;j<roN2;j++)      for(i=0;i<roN1;i++) olddata[i][j][k]=olddata0[i][j][k];
      }
    }
    // reset size
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
void low2high(int nxlow, int nylow, int nzlow, int nxhigh, int nyhigh, int nzhigh, unsigned char***oldimage,FTYPE***olddata)
{
  FTYPE **Ilowf;
  unsigned char **Ilowc;
  int ih,jh,kl;
  int il,jl;
  FTYPE dil, djl,dkl,ftemp,ftempi,ftempj,ftempk;
  FTYPE r,th,ph;
  unsigned char uctemp;

  int j,k;

  FTYPE **y1a,**y2a,**y12a;
  FTYPE **ya;
  FTYPE Xget[NDIM];


  // GODMARK3D
  kl=0;
  ph=0;
  dkl=0;
  ftempk=0; // GODMARK3D

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
      if(DATATYPE==0) Ilowc[il][jl]=oldimage[il][jl][0];// GODMARK3D
      else if(DATATYPE==1) Ilowf[il][jl]=olddata[il][jl][0]; // GODMARK3D
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
	ya[j][k]=(FTYPE)oldimage[j][k][0]; // GODMARK3D
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
	ya[j][k]=olddata[j][k][0]; // GODMARK3D
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
    ftemp=olddata[il][jl][0]; // GODMARK3D
    if(ftemp<0.0) ftemp=0.0;
    if(ftemp>255.0) ftemp=255.0;
    uctemp=(unsigned char)ftemp;
    oldimage[ih][jh][0]=uctemp; // GODMARK3D
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

      ftemp=nearest_interp_ij(il, jl, kl ,Ilowc,Ilowf);
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

      ftemp=bilinear_interp_ij(il, jl , kl, dil,djl,dkl, Ilowc,Ilowf);
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

      oldf_ijk2rthph(ftempi, ftempj, ftempk, &r, &th, &ph);
      // fprintf(stderr,"%d %d : : %g %g : %g %g\n",jh,ih,ftempi,ftempj,r,th); fflush(stderr);
      ftemp=plane_interp(r, th, ph, il, jl,kl, Ilowc,Ilowf);
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
      oldf_ijk2x123(ftempi,ftempj,ftempk,Xget);
      ftemp=bicubic_interp(nxlow,nylow,il,jl,kl,Xget,ya,y1a,y2a,y12a,dx[1],dx[2]);
#endif


      //      fprintf(stderr,"ih=%d jh=%d il=%d jl=%d dil=%g djl=%g ftemp=%10.5g\n",ih,jh,il,jl,dil,djl,ftemp); fflush(stderr);


      // just write back into the given multi-pointers
      if(DATATYPE==0) oldimage[ih][jh][0]=(unsigned char)ftemp;// GODMARK3D
      else olddata[ih][jh][0]=ftemp; // GODMARK3D

      //fprintf(stderr,"olddata[%d][%d][%d]=%g\n",ih,jh,0,olddata[ih][jh][0]); // GODMARK3D
    }
  }
  if(0){// debug (need to setup this memory stuff -- segfaults for data currently
    // in principle could output in different order if wanted
    if(DATATYPE==1){
      if(1){ // needs to be set
	for(jh=0;jh<nyhigh;jh++)      for(ih=0;ih<nxhigh;ih++) {
	    ftemp=olddata[ih][jh][0];// GODMARK3D
	    if(ftemp<0.0) ftemp=0.0;
	    if(ftemp>255.0) ftemp=255.0;
	    uctemp=(unsigned char)ftemp;
	    oldimage[ih][jh][0]=uctemp; // GODMARK3D
	  }
      }
    }
    writeimage("jontest.r8",oldimage,nxhigh,nyhigh,nzhigh);
    //exit(0);
  }
 

  if(DATATYPE==0){
    free_cmatrix(Ilowc,0,nxlow-1,0,nylow-1) ;
  }
  else if(DATATYPE==1){
    free_fmatrix(Ilowf,0,nxlow-1,0,nylow-1) ;
  }
  //  exit(0);

}


// area weighted based interpolation from high resolution to lower resolution
// assumes feed in high resolution and put low resolution version into same buffer at end
void high2low(int nxhigh, int nyhigh, int nzhigh, int nxlow, int nylow, int nzlow, unsigned char ***oldimage,FTYPE***olddata)
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
      Ihigh[jh*nxhigh+ih]=olddata[ih][jh][0]; // GODMARK3D
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
      olddata[il][jl][0]=Ilow[jl*nxlow+il]/IlowW[jl*nxlow+il]; // GODMARK3D
    }
  }

  free(Ihigh);
  free(Ilow);
  free(IlowW);



}
