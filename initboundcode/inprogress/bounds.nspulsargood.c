
#include "decs.h"

/* bound array containing entire set of primitive variables */

#define MAXORDEREXTRAP 4

#define EXTRAP 1
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

int bound_prim(int boundstage, FTYPE prim[][N2M][NPR])
{
  static FTYPE *myBsq;
  static int firsttime=1;

  extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *ptrgeom,FTYPE *p2interp);
  SFTYPE Dt=-1.0; // no boundary zone should ever be in accounting
  int i, j, k;
  struct of_geom geom,rgeom;
  FTYPE vcon[NDIM]; // coordinate basis vcon
  FTYPE ucon[NDIM];
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j
  FTYPE prescale[NPR];
  FTYPE rescaled[NBIGBND][NPR];
  FTYPE Bcon[NDIM],Bcov[NDIM];
  FTYPE copied[NBIGBND][NPR];
  FTYPE det,sqrtdet;
  FTYPE utmin;
  FTYPE signsolution;

  FTYPE mytop,mybottom;
  FTYPE uu0,vu3,vu1o,vu2o;
  FTYPE vold[NDIM];
  FTYPE avu1free1,avu1free2,avu2free1,avu2free2;

  FTYPE uu0a,uu0b,ucov[NDIM];

  FTYPE Bsq;

  FTYPE r,th,X[NDIM];
  FTYPE xc[MAXORDEREXTRAP+2],yc[MAXORDEREXTRAP+2],AA,BB,CC,DD,EE;

  FTYPE newB3;
  FTYPE plusminus;

  int extraps;

  // 0: normal other field components not relaxed
  // 1: fix B^2
  // 2: fix B^r,  relax B_\phi and B^\theta
  // 3: fix B^r,  relax B^\phi and B^\theta
  int WHICHFIXED;

  // 0: B^\phi
  // 1: B_\phi
  // 2: B^\theta
  // 3: B^r
  int WHICHTOEXTRAP;


  // 0: 0th
  // 2: parabolic
  // 3: cubic
  int ORDEREXTRAP;



  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    if (mycpupos[1] == 0) {
      if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
	/* inner r boundary condition: u, just copy */
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--)	  PBOUNDLOOP	  prim[i][j][k] = prim[ri][rj][k];
#elif(EXTRAP==1)

	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    for(k=RHO;k<=UU;k++){
	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    }
	    k=U1;	    // treat U1 as special
	    prim[i][j][k] = prim[ri][rj][k] * (1. - (i-ri)*dx[1]) ;
	    for(k=U2;k<=U3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. + (i-ri)*dx[1]) ;
	    }
	    k=B1; // treat B1 special
	    prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    for(k=B2;k<=B3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. + (i-ri)*dx[1]) ;
	    }
	  }
#elif(EXTRAP==2)
	  ri=0;
	  rj=j;
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  for(i=-1;i>=-N1BND;i--){
	    // set guess
	    PBOUNDLOOP prim[i][j][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);	    
	    rescale(-1,1,prim[i][j],&geom,prescale);
	  }
#endif
	  for(i=-1;i>=-N1BND;i--){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom, Dt) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom, Dt) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#endif
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom,Dt) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom,Dt)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	  }
	}
      }
      else if((BCtype[X1DN]==NSSURFACE)){
	// same as above but don't modify magnetic field
	// fix field
	// outflow v and densities 

	if(firsttime){
	  myBsq=malloc(sizeof(FTYPE)*N1BND*N2);
	  if(myBsq==NULL){
	    dualfprintf(fail_file,"Couldn't open myBsq\n");
	    myexit(1);
	  }
	  // shift it so same as boundary region
	  myBsq+=N1BND*N2;
	  // fill with correct values
	  for(i=-1;i>=-N1BND;i--) for (j = 0; j < N2; j++){

	    // use initial/fixed field
	    Bcon[0]=0;
	    Bcon[1]=p[i][j][B1];
	    Bcon[2]=p[i][j][B2];
	    Bcon[3]=p[i][j][B3];

	    get_geometry(i, j, CENT, &geom);
	    lower(Bcon,&geom,Bcov);
	    Bsq=0.0;
	    for(k=0;k<=3;k++) Bsq+=Bcon[k]*Bcov[k];

	    myBsq[i*N2+j]=Bsq; // safe and sound storage of B^2
	    //	    dualfprintf(fail_file,"i=%d j=%d Bsq=%g\n",i,j,myBsq[i*N2+j]);
	  }// end loop
	}// end if firsttime




	/* inner r boundary condition: u, just copy */
	for (j = 0; j < N2; j++) {


#if(EXTRAP==1)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k]=prim[ri][rj][k];
	  }
	  // have B^\phi copied from above
	  // Force B^r and B^\theta to be analytic solution
	  for(i=-1;i>=-N1BND;i--){

	    get_geometry(i, j, CENT, &geom);
	    coord(i,j,CENT,X);
	    bl_coord(X,&r,&th);

#if(0)
	    //	    if( ((int)(t/20.0))%2==0){// relax B^\theta first
	    if(t<100.0){// relax B^\theta first
	      // this will overwrite fixed B^\theta field.
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=2;
	      ORDEREXTRAP=2;
	    }
	    else if(t<150.0){
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=2;
	      ORDEREXTRAP=4;
	    }
	    else if(t<200.0){
	      // this uses new adjusted (as above) fixed B^\theta and adjusts B_\phi
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=0;
	      ORDEREXTRAP=2;
	    }
	    else{
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=0;
	      ORDEREXTRAP=4;
	    }

#elif(0)
	    //	    if( ((int)(t/20.0))%2==0){// relax B^\theta first
	    if(t<200.0){// relax B^\theta first
	      // this will overwrite fixed B^\theta field.
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=2;
	      ORDEREXTRAP=2;
	    }
	    else if(t<400.0){
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=0;
	      ORDEREXTRAP=2;
	    }
	    else{
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=0;
	      ORDEREXTRAP=2;
	    }

#elif(0)
	    // this doesn't work because field grows instead of getting smaller, so B^\phi is forced to be 0 or nonexistent
	    WHICHFIXED=1;  // fix B^2
	    WHICHTOEXTRAP=2; // B^\theta (can't extrap B_\phi and fix B^2 (yet))
	    ORDEREXTRAP=2; // parabolic
#elif(0)
	      WHICHFIXED=2;  // fix ONLY B^r
	      // WHICHTOEXTRAP==2,0
	      ORDEREXTRAP=2; 
	    }
#elif(1)
	    if(t<100.0){
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=2;
	      ORDEREXTRAP=2;	      
	    }
	    else if(t<200.0){
	      WHICHFIXED=0;
	      WHICHTOEXTRAP=0;
	      ORDEREXTRAP=2;
	    }
	    else{
	      WHICHFIXED=2;  // fix ONLY B^r
	      // WHICHTOEXTRAP==2,0
	      ORDEREXTRAP=2; 
	    }
#endif

	    //	    if((th<2.89377)&&(th>M_PI-2.89377)){ // then fix
	    //	    if(1||(th<M_PI-.774609)&&(th>0.774609)){ // then fix
	    if(0&&(th<M_PI-.702388)&&(th>0.702388)){ // then fix
	      prim[i][j][B3]=0.0;
	    }
	    else{

	      for(extraps=0;extraps<=1;extraps++){
		if(WHICHFIXED<2) extraps=3; // get out after once
		else if(WHICHFIXED==2){
		  if(extraps==0) WHICHTOEXTRAP=2;
		  else WHICHTOEXTRAP=0;
		}
		else if(WHICHFIXED==3){
		  if(extraps==0) WHICHTOEXTRAP=2;
		  else WHICHTOEXTRAP=1;
		}

	      if(ORDEREXTRAP>=2){
		if(WHICHTOEXTRAP==0){
		  if(WHICHFIXED==0) for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
		  else if(WHICHFIXED==1){
		    // relaxing B^\phi, assume fixing B^r, and solve for B^\theta
		    prim[i][j][B1]=p[i][j][B1];
		  }
		  else if(WHICHFIXED==2) for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];

		  // extrapolation of B^\phi
		  for(k=1;k<=ORDEREXTRAP+1;k++) yc[k]=prim[ri+(k-1)][rj][B3];
		}
		else if(WHICHTOEXTRAP==1){
		  if(WHICHFIXED==0) for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
		  else if(WHICHFIXED==2) for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		  // extrapolation of B_\phi
		  for(k=1;k<=ORDEREXTRAP+1;k++){

		    Bcon[0]=0;
		    Bcon[1]=prim[ri+(k-1)][rj][B1];
		    Bcon[2]=prim[ri+(k-1)][rj][B2];
		    Bcon[3]=prim[ri+(k-1)][rj][B3];
		    get_geometry(ri+(k-1), rj, CENT, &geom);
		    lower(Bcon,&geom,Bcov);
		    yc[k]=Bcov[3]*geom.gcov[0][0];
		  }
		}
		else if(WHICHTOEXTRAP==2){
		  if(WHICHFIXED==0){
		    for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		    for(k=B3;k<=B3;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  else if(WHICHFIXED==1){
		    // relaxing B^\theta, assume fixing B^r, and solve for B^\phi
		    prim[i][j][B1]=p[i][j][B1];
		  }
		  else if(WHICHFIXED==2){
		    for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  for(k=1;k<=ORDEREXTRAP+1;k++) yc[k]=prim[ri+(k-1)][rj][B2];
		}
		else if(WHICHTOEXTRAP==3){
		  if(WHICHFIXED==0){
		    for(k=B2;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
		    for(k=B3;k<=B3;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  else if(WHICHFIXED==1){
		    // relaxing B^r, assume fixing B^theta, and solve for B^\phi
		    prim[i][j][B2]=p[i][j][B2];
		  }
		  else if(WHICHFIXED==2){ // never used...
		    for(k=B2;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  for(k=1;k<=ORDEREXTRAP+1;k++) yc[k]=prim[ri+(k-1)][rj][B1];
		}

		// location of boundary zone
		get_geometry(i, j, CENT, &geom);
		coord(i,j,CENT,X);
		bl_coord(X,&r,&th);
		xc[0]=r;

		for(k=1;k<=ORDEREXTRAP+1;k++){
		  // zones to extrapolate from
		  get_geometry(ri+(k-1), rj, CENT, &geom);
		  coord(ri+(k-1),rj,CENT,X);
		  bl_coord(X,&r,&th);
		  xc[k]=r;
		}

		if(ORDEREXTRAP==2){

		  AA=( yc[1]*xc[2]*xc[3]*(xc[2]-xc[3]) +  yc[2]*xc[1]*xc[3]*(xc[3]-xc[1]) + yc[3]*xc[1]*xc[2]*(xc[1]-xc[2]) )/((xc[1]-xc[2])*(xc[1]-xc[3])*(xc[2]-xc[3]));
		  
		  BB=( yc[1]*(xc[3]*xc[3] - xc[2]*xc[2]) +  yc[2]*(xc[1]*xc[1] - xc[3]*xc[3]) + yc[3]*(xc[2]*xc[2] - xc[1]*xc[1]) )/((xc[1]-xc[2])*(xc[1]-xc[3])*(xc[2]-xc[3]));
		  
		  //CC=((yc[1]-yc[3])/(xc[1]-xc[3]) - (yc[2]-yc[3])/(xc[2]-xc[3]))/(xc[1]-xc[2]);
		  CC=( yc[1]*(xc[2] - xc[3]) +  yc[2]*(xc[3] - xc[1]) + yc[3]*(xc[1] - xc[2]) )/((xc[1]-xc[2])*(xc[1]-xc[3])*(xc[2]-xc[3]));

		  yc[0]=AA+BB*xc[0]+CC*xc[0]*xc[0];

		}
		else if(ORDEREXTRAP==3){
		  AA=(xc[1]*(xc[1] - xc[3])*xc[3]*(xc[1] - xc[4])*(xc[3] - xc[4])*xc[4]*yc[2] + xc[2]*pow(xc[4],2)*(-(pow(xc[3],3)*yc[1]) + pow(xc[3],2)*xc[4]*yc[1] + pow(xc[1],2)*(xc[1] - xc[4])*yc[3]) + 
		      pow(xc[1],2)*xc[2]*pow(xc[3],2)*(-xc[1] + xc[3])*yc[4] + pow(xc[2],3)*(xc[4]*(-(pow(xc[3],2)*yc[1]) + xc[3]*xc[4]*yc[1] + xc[1]*(xc[1] - xc[4])*yc[3]) + xc[1]*xc[3]*(-xc[1] + xc[3])*yc[4]) + 
		      pow(xc[2],2)*(xc[1]*xc[4]*(-pow(xc[1],2) + pow(xc[4],2))*yc[3] + pow(xc[3],3)*(xc[4]*yc[1] - xc[1]*yc[4]) + xc[3]*(-(pow(xc[4],3)*yc[1]) + pow(xc[1],3)*yc[4])))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4]));


		  BB=(pow(xc[1],2)*(xc[1] - xc[4])*pow(xc[4],2)*(yc[2] - yc[3]) + pow(xc[3],3)*(pow(xc[4],2)*(yc[1] - yc[2]) + pow(xc[1],2)*(yc[2] - yc[4])) + 
		      pow(xc[2],2)*(pow(xc[4],3)*(yc[1] - yc[3]) + pow(xc[1],3)*(yc[3] - yc[4]) + pow(xc[3],3)*(-yc[1] + yc[4])) + pow(xc[3],2)*(pow(xc[4],3)*(-yc[1] + yc[2]) + pow(xc[1],3)*(-yc[2] + yc[4])) + 
		      pow(xc[2],3)*(pow(xc[4],2)*(-yc[1] + yc[3]) + pow(xc[3],2)*(yc[1] - yc[4]) + pow(xc[1],2)*(-yc[3] + yc[4])))/((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4]));

		  CC=(-(xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] + xc[4])*(yc[2] - yc[3])) + xc[3]*(pow(xc[4],3)*(yc[1] - yc[2]) + pow(xc[1],3)*(yc[2] - yc[4])) + pow(xc[3],3)*(-(xc[4]*yc[1]) - xc[1]*yc[2] + xc[4]*yc[2] + xc[1]*yc[4]) + 
		      pow(xc[2],3)*(xc[4]*yc[1] + xc[1]*yc[3] - xc[4]*yc[3] - xc[1]*yc[4] + xc[3]*(-yc[1] + yc[4])) + xc[2]*(pow(xc[4],3)*(-yc[1] + yc[3]) + pow(xc[3],3)*(yc[1] - yc[4]) + pow(xc[1],3)*(-yc[3] + yc[4])))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4]));

		  DD=(xc[1]*(xc[1] - xc[4])*xc[4]*(yc[2] - yc[3]) + pow(xc[3],2)*(xc[4]*yc[1] + xc[1]*yc[2] - xc[4]*yc[2] - xc[1]*yc[4]) + pow(xc[2],2)*(-(xc[4]*yc[1]) - xc[1]*yc[3] + xc[4]*yc[3] + xc[3]*(yc[1] - yc[4]) + xc[1]*yc[4]) + 
		      xc[2]*(pow(xc[4],2)*(yc[1] - yc[3]) + pow(xc[1],2)*(yc[3] - yc[4]) + pow(xc[3],2)*(-yc[1] + yc[4])) + xc[3]*(pow(xc[4],2)*(-yc[1] + yc[2]) + pow(xc[1],2)*(-yc[2] + yc[4])))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4]));


		  yc[0]=AA+BB*xc[0]+CC*xc[0]*xc[0]+DD*xc[0]*xc[0]*xc[0];


		}
		else if(ORDEREXTRAP==4){
		  AA=(-(xc[1]*(xc[1] - xc[3])*xc[3]*(xc[1] - xc[4])*(xc[3] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[3] - xc[5])*(xc[4] - xc[5])*xc[5]*yc[2]) + 
		      pow(xc[2],2)*(-(xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[4] - xc[5])*xc[5]*(xc[1]*xc[4] + (xc[1] + xc[4])*xc[5])*yc[3]) + xc[3]*pow(xc[5],3)*(pow(xc[4],3)*(xc[4] - xc[5])*yc[1] + pow(xc[1],3)*(-xc[1] + xc[5])*yc[4]) + 
				   pow(xc[1],3)*xc[3]*(xc[1] - xc[4])*pow(xc[4],3)*yc[5] + pow(xc[3],4)*(pow(xc[4],3)*xc[5]*yc[1] - xc[4]*pow(xc[5],3)*yc[1] - pow(xc[1],3)*xc[5]*yc[4] + xc[1]*pow(xc[5],3)*yc[4] + pow(xc[1],3)*xc[4]*yc[5] - 
											  xc[1]*pow(xc[4],3)*yc[5]) + pow(xc[3],3)*(-(pow(xc[4],4)*xc[5]*yc[1]) + xc[4]*pow(xc[5],4)*yc[1] + pow(xc[1],4)*xc[5]*yc[4] - xc[1]*pow(xc[5],4)*yc[4] - pow(xc[1],4)*xc[4]*yc[5] + 
															    xc[1]*pow(xc[4],4)*yc[5])) + pow(xc[2],4)*(-(xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[4] - xc[5])*xc[5]*yc[3]) + xc[3]*pow(xc[5],2)*(pow(xc[4],2)*(xc[4] - xc[5])*yc[1] + pow(xc[1],2)*(-xc[1] + xc[5])*yc[4]) + 
																			       pow(xc[1],2)*xc[3]*(xc[1] - xc[4])*pow(xc[4],2)*yc[5] + pow(xc[3],3)*(xc[5]*(xc[4]*(xc[4] - xc[5])*yc[1] + xc[1]*(-xc[1] + xc[5])*yc[4]) + xc[1]*(xc[1] - xc[4])*xc[4]*yc[5]) + 
																			       pow(xc[3],2)*(xc[1]*(xc[1] - xc[5])*xc[5]*(xc[1] + xc[5])*yc[4] + pow(xc[4],3)*(-(xc[5]*yc[1]) + xc[1]*yc[5]) + xc[4]*(pow(xc[5],3)*yc[1] - pow(xc[1],3)*yc[5]))) + 
		      pow(xc[2],3)*(xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[4] - xc[5])*xc[5]*(xc[1] + xc[4] + xc[5])*yc[3] + 
				   pow(xc[3],4)*(xc[5]*(-(pow(xc[4],2)*yc[1]) + xc[4]*xc[5]*yc[1] + xc[1]*(xc[1] - xc[5])*yc[4]) + xc[1]*xc[4]*(-xc[1] + xc[4])*yc[5]) + 
				   xc[3]*(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(xc[1] + xc[5])*yc[4] + pow(xc[4],4)*(-(pow(xc[5],2)*yc[1]) + pow(xc[1],2)*yc[5]) + pow(xc[4],2)*(pow(xc[5],4)*yc[1] - pow(xc[1],4)*yc[5])) + 
				   pow(xc[3],2)*(xc[1]*xc[5]*(-pow(xc[1],3) + pow(xc[5],3))*yc[4] + pow(xc[4],4)*(xc[5]*yc[1] - xc[1]*yc[5]) + xc[4]*(-(pow(xc[5],4)*yc[1]) + pow(xc[1],4)*yc[5]))) + 
		      xc[2]*(pow(xc[1],2)*(xc[1] - xc[4])*pow(xc[4],2)*(xc[1] - xc[5])*(xc[4] - xc[5])*pow(xc[5],2)*yc[3] + 
			  pow(xc[3],4)*(pow(xc[5],2)*(-(pow(xc[4],3)*yc[1]) + pow(xc[4],2)*xc[5]*yc[1] + pow(xc[1],2)*(xc[1] - xc[5])*yc[4]) + pow(xc[1],2)*pow(xc[4],2)*(-xc[1] + xc[4])*yc[5]) + 
			  pow(xc[3],2)*(pow(xc[5],3)*(-(pow(xc[4],4)*yc[1]) + pow(xc[4],3)*xc[5]*yc[1] + pow(xc[1],3)*(xc[1] - xc[5])*yc[4]) + pow(xc[1],3)*pow(xc[4],3)*(-xc[1] + xc[4])*yc[5]) + 
			  pow(xc[3],3)*(pow(xc[1],2)*pow(xc[5],2)*(-pow(xc[1],2) + pow(xc[5],2))*yc[4] + pow(xc[4],4)*(pow(xc[5],2)*yc[1] - pow(xc[1],2)*yc[5]) + 
				       pow(xc[4],2)*(-(pow(xc[5],4)*yc[1]) + pow(xc[1],4)*yc[5]))))/((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4])*(xc[1] - xc[5])*(xc[2] - xc[5])*(xc[3] - xc[5])*(xc[4] - xc[5]));



		  BB=(pow(xc[1],2)*(xc[1] - xc[4])*pow(xc[4],2)*(xc[1] - xc[5])*(xc[4] - xc[5])*pow(xc[5],2)*(yc[2] - yc[3]) + 
		      pow(xc[3],3)*(-(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(xc[1] + xc[5])*(yc[2] - yc[4])) + pow(xc[4],2)*(pow(xc[5],4)*(yc[1] - yc[2]) + pow(xc[1],4)*(yc[2] - yc[5])) + 
				   pow(xc[4],4)*(pow(xc[5],2)*(-yc[1] + yc[2]) + pow(xc[1],2)*(-yc[2] + yc[5]))) + 
		      pow(xc[3],4)*(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(yc[2] - yc[4]) + pow(xc[4],3)*(pow(xc[5],2)*(yc[1] - yc[2]) + pow(xc[1],2)*(yc[2] - yc[5])) + 
				   pow(xc[4],2)*(pow(xc[5],3)*(-yc[1] + yc[2]) + pow(xc[1],3)*(-yc[2] + yc[5]))) + 
		      pow(xc[3],2)*(pow(xc[1],3)*(xc[1] - xc[5])*pow(xc[5],3)*(yc[2] - yc[4]) + pow(xc[4],4)*(pow(xc[5],3)*(yc[1] - yc[2]) + pow(xc[1],3)*(yc[2] - yc[5])) + 
				   pow(xc[4],3)*(pow(xc[5],4)*(-yc[1] + yc[2]) + pow(xc[1],4)*(-yc[2] + yc[5]))) + 
		      pow(xc[2],3)*(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(xc[1] + xc[5])*(yc[3] - yc[4]) + pow(xc[4],4)*(pow(xc[5],2)*(yc[1] - yc[3]) + pow(xc[1],2)*(yc[3] - yc[5])) + 
				   pow(xc[3],2)*(pow(xc[5],4)*(yc[1] - yc[4]) + pow(xc[1],4)*(yc[4] - yc[5]) + pow(xc[4],4)*(-yc[1] + yc[5])) + pow(xc[4],2)*(pow(xc[5],4)*(-yc[1] + yc[3]) + pow(xc[1],4)*(-yc[3] + yc[5])) + 
				   pow(xc[3],4)*(pow(xc[5],2)*(-yc[1] + yc[4]) + pow(xc[4],2)*(yc[1] - yc[5]) + pow(xc[1],2)*(-yc[4] + yc[5]))) + 
		      pow(xc[2],4)*(-(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(yc[3] - yc[4])) + pow(xc[4],2)*(pow(xc[5],3)*(yc[1] - yc[3]) + pow(xc[1],3)*(yc[3] - yc[5])) + 
				   pow(xc[3],3)*(pow(xc[5],2)*(yc[1] - yc[4]) + pow(xc[1],2)*(yc[4] - yc[5]) + pow(xc[4],2)*(-yc[1] + yc[5])) + pow(xc[4],3)*(pow(xc[5],2)*(-yc[1] + yc[3]) + pow(xc[1],2)*(-yc[3] + yc[5])) + 
				   pow(xc[3],2)*(pow(xc[5],3)*(-yc[1] + yc[4]) + pow(xc[4],3)*(yc[1] - yc[5]) + pow(xc[1],3)*(-yc[4] + yc[5]))) + 
		      pow(xc[2],2)*(-(pow(xc[1],3)*(xc[1] - xc[5])*pow(xc[5],3)*(yc[3] - yc[4])) + pow(xc[4],3)*(pow(xc[5],4)*(yc[1] - yc[3]) + pow(xc[1],4)*(yc[3] - yc[5])) + 
				   pow(xc[3],4)*(pow(xc[5],3)*(yc[1] - yc[4]) + pow(xc[1],3)*(yc[4] - yc[5]) + pow(xc[4],3)*(-yc[1] + yc[5])) + pow(xc[4],4)*(pow(xc[5],3)*(-yc[1] + yc[3]) + pow(xc[1],3)*(-yc[3] + yc[5])) + 
				   pow(xc[3],3)*(pow(xc[5],4)*(-yc[1] + yc[4]) + pow(xc[4],4)*(yc[1] - yc[5]) + pow(xc[1],4)*(-yc[4] + yc[5]))))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4])*(xc[1] - xc[5])*(xc[2] - xc[5])*(xc[3] - xc[5])*(xc[4] - xc[5]));

		  CC=(-(xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[4] - xc[5])*xc[5]*(xc[4]*xc[5] + xc[1]*(xc[4] + xc[5]))*(yc[2] - yc[3])) + 
		      pow(xc[3],4)*(-(xc[1]*(xc[1] - xc[5])*xc[5]*(xc[1] + xc[5])*(yc[2] - yc[4])) + xc[4]*(pow(xc[5],3)*(yc[1] - yc[2]) + pow(xc[1],3)*(yc[2] - yc[5])) + pow(xc[4],3)*(xc[5]*(-yc[1] + yc[2]) + xc[1]*(-yc[2] + yc[5]))) + 
		      xc[3]*(-(pow(xc[1],3)*(xc[1] - xc[5])*pow(xc[5],3)*(yc[2] - yc[4])) + pow(xc[4],3)*(pow(xc[5],4)*(yc[1] - yc[2]) + pow(xc[1],4)*(yc[2] - yc[5])) + 
			  pow(xc[4],4)*(pow(xc[5],3)*(-yc[1] + yc[2]) + pow(xc[1],3)*(-yc[2] + yc[5]))) + 
		      pow(xc[3],3)*(xc[1]*xc[5]*(pow(xc[1],3) - pow(xc[5],3))*(yc[2] - yc[4]) + pow(xc[4],4)*(xc[5]*(yc[1] - yc[2]) + xc[1]*(yc[2] - yc[5])) + xc[4]*(pow(xc[5],4)*(-yc[1] + yc[2]) + pow(xc[1],4)*(-yc[2] + yc[5]))) + 
		      pow(xc[2],4)*(xc[1]*(xc[1] - xc[5])*xc[5]*(xc[1] + xc[5])*(yc[3] - yc[4]) + pow(xc[4],3)*(xc[5]*yc[1] + xc[1]*yc[3] - xc[5]*yc[3] - xc[1]*yc[5]) + 
				   pow(xc[3],3)*(-(xc[5]*yc[1]) - xc[1]*yc[4] + xc[5]*yc[4] + xc[4]*(yc[1] - yc[5]) + xc[1]*yc[5]) + xc[3]*(pow(xc[5],3)*(yc[1] - yc[4]) + pow(xc[1],3)*(yc[4] - yc[5]) + pow(xc[4],3)*(-yc[1] + yc[5])) + 
				   xc[4]*(pow(xc[5],3)*(-yc[1] + yc[3]) + pow(xc[1],3)*(-yc[3] + yc[5]))) + xc[2]*
		      (pow(xc[1],3)*(xc[1] - xc[5])*pow(xc[5],3)*(yc[3] - yc[4]) + pow(xc[4],4)*(pow(xc[5],3)*(yc[1] - yc[3]) + pow(xc[1],3)*(yc[3] - yc[5])) + 
		       pow(xc[3],3)*(pow(xc[5],4)*(yc[1] - yc[4]) + pow(xc[1],4)*(yc[4] - yc[5]) + pow(xc[4],4)*(-yc[1] + yc[5])) + pow(xc[4],3)*(pow(xc[5],4)*(-yc[1] + yc[3]) + pow(xc[1],4)*(-yc[3] + yc[5])) + 
		       pow(xc[3],4)*(pow(xc[5],3)*(-yc[1] + yc[4]) + pow(xc[4],3)*(yc[1] - yc[5]) + pow(xc[1],3)*(-yc[4] + yc[5]))) + 
		      pow(xc[2],3)*(-(xc[1]*xc[5]*(pow(xc[1],3) - pow(xc[5],3))*(yc[3] - yc[4])) + xc[4]*(pow(xc[5],4)*(yc[1] - yc[3]) + pow(xc[1],4)*(yc[3] - yc[5])) + 
				   pow(xc[3],4)*(xc[5]*(yc[1] - yc[4]) + xc[1]*(yc[4] - yc[5]) + xc[4]*(-yc[1] + yc[5])) + pow(xc[4],4)*(xc[5]*(-yc[1] + yc[3]) + xc[1]*(-yc[3] + yc[5])) + 
				   xc[3]*(pow(xc[5],4)*(-yc[1] + yc[4]) + pow(xc[4],4)*(yc[1] - yc[5]) + pow(xc[1],4)*(-yc[4] + yc[5]))))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4])*(xc[1] - xc[5])*(xc[2] - xc[5])*(xc[3] - xc[5])*(xc[4] - xc[5]));


		  DD=(xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[4] - xc[5])*xc[5]*(xc[1] + xc[4] + xc[5])*(yc[2] - yc[3]) + 
		      pow(xc[3],2)*(-(xc[1]*xc[5]*(pow(xc[1],3) - pow(xc[5],3))*(yc[2] - yc[4])) + xc[4]*(pow(xc[5],4)*(yc[1] - yc[2]) + pow(xc[1],4)*(yc[2] - yc[5])) + 
				   pow(xc[4],4)*(xc[5]*(-yc[1] + yc[2]) + xc[1]*(-yc[2] + yc[5]))) + pow(xc[3],4)*
		      (xc[1]*(xc[1] - xc[5])*xc[5]*(yc[2] - yc[4]) + pow(xc[4],2)*(xc[5]*(yc[1] - yc[2]) + xc[1]*(yc[2] - yc[5])) + xc[4]*(pow(xc[5],2)*(-yc[1] + yc[2]) + pow(xc[1],2)*(-yc[2] + yc[5]))) + 
		      xc[3]*(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(xc[1] + xc[5])*(yc[2] - yc[4]) + pow(xc[4],4)*(pow(xc[5],2)*(yc[1] - yc[2]) + pow(xc[1],2)*(yc[2] - yc[5])) + 
			  pow(xc[4],2)*(pow(xc[5],4)*(-yc[1] + yc[2]) + pow(xc[1],4)*(-yc[2] + yc[5]))) + 
		      pow(xc[2],2)*(xc[1]*xc[5]*(pow(xc[1],3) - pow(xc[5],3))*(yc[3] - yc[4]) + pow(xc[4],4)*(xc[5]*(yc[1] - yc[3]) + xc[1]*(yc[3] - yc[5])) + 
				   xc[3]*(pow(xc[5],4)*(yc[1] - yc[4]) + pow(xc[1],4)*(yc[4] - yc[5]) + pow(xc[4],4)*(-yc[1] + yc[5])) + xc[4]*(pow(xc[5],4)*(-yc[1] + yc[3]) + pow(xc[1],4)*(-yc[3] + yc[5])) + 
				   pow(xc[3],4)*(xc[5]*(-yc[1] + yc[4]) + xc[4]*(yc[1] - yc[5]) + xc[1]*(-yc[4] + yc[5]))) + 
		      pow(xc[2],4)*(-(xc[1]*(xc[1] - xc[5])*xc[5]*(yc[3] - yc[4])) + xc[4]*(pow(xc[5],2)*(yc[1] - yc[3]) + pow(xc[1],2)*(yc[3] - yc[5])) + pow(xc[4],2)*(-(xc[5]*yc[1]) - xc[1]*yc[3] + xc[5]*yc[3] + xc[1]*yc[5]) + 
				   pow(xc[3],2)*(xc[5]*yc[1] + xc[1]*yc[4] - xc[5]*yc[4] - xc[1]*yc[5] + xc[4]*(-yc[1] + yc[5])) + xc[3]*(pow(xc[5],2)*(-yc[1] + yc[4]) + pow(xc[4],2)*(yc[1] - yc[5]) + pow(xc[1],2)*(-yc[4] + yc[5]))) + 
		      xc[2]*(-(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(xc[1] + xc[5])*(yc[3] - yc[4])) + pow(xc[4],2)*(pow(xc[5],4)*(yc[1] - yc[3]) + pow(xc[1],4)*(yc[3] - yc[5])) + 
			  pow(xc[3],4)*(pow(xc[5],2)*(yc[1] - yc[4]) + pow(xc[1],2)*(yc[4] - yc[5]) + pow(xc[4],2)*(-yc[1] + yc[5])) + pow(xc[4],4)*(pow(xc[5],2)*(-yc[1] + yc[3]) + pow(xc[1],2)*(-yc[3] + yc[5])) + 
			  pow(xc[3],2)*(pow(xc[5],4)*(-yc[1] + yc[4]) + pow(xc[4],4)*(yc[1] - yc[5]) + pow(xc[1],4)*(-yc[4] + yc[5]))))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4])*(xc[1] - xc[5])*(xc[2] - xc[5])*(xc[3] - xc[5])*(xc[4] - xc[5]));



		  EE= (xc[1]*(xc[1] - xc[4])*xc[4]*(xc[1] - xc[5])*(xc[4] - xc[5])*xc[5]*(-yc[2] + yc[3]) + pow(xc[3],3)*(-(xc[1]*(xc[1] - xc[5])*xc[5]*(yc[2] - yc[4])) + xc[4]*(pow(xc[5],2)*(yc[1] - yc[2]) + pow(xc[1],2)*(yc[2] - yc[5])) + 
											pow(xc[4],2)*(xc[5]*(-yc[1] + yc[2]) + xc[1]*(-yc[2] + yc[5]))) + xc[3]*(-(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(yc[2] - yc[4])) + 
																	   pow(xc[4],2)*(pow(xc[5],3)*(yc[1] - yc[2]) + pow(xc[1],3)*(yc[2] - yc[5])) + pow(xc[4],3)*(pow(xc[5],2)*(-yc[1] + yc[2]) + pow(xc[1],2)*(-yc[2] + yc[5]))) + 
		       pow(xc[3],2)*(xc[1]*(xc[1] - xc[5])*xc[5]*(xc[1] + xc[5])*(yc[2] - yc[4]) + pow(xc[4],3)*(xc[5]*(yc[1] - yc[2]) + xc[1]*(yc[2] - yc[5])) + xc[4]*(pow(xc[5],3)*(-yc[1] + yc[2]) + pow(xc[1],3)*(-yc[2] + yc[5]))) + 
		       pow(xc[2],3)*(xc[1]*(xc[1] - xc[5])*xc[5]*(yc[3] - yc[4]) + pow(xc[4],2)*(xc[5]*yc[1] + xc[1]*yc[3] - xc[5]*yc[3] - xc[1]*yc[5]) + pow(xc[3],2)*(-(xc[5]*yc[1]) - xc[1]*yc[4] + xc[5]*yc[4] + xc[4]*(yc[1] - yc[5]) + xc[1]*yc[5]) + 
				    xc[3]*(pow(xc[5],2)*(yc[1] - yc[4]) + pow(xc[1],2)*(yc[4] - yc[5]) + pow(xc[4],2)*(-yc[1] + yc[5])) + xc[4]*(pow(xc[5],2)*(-yc[1] + yc[3]) + pow(xc[1],2)*(-yc[3] + yc[5]))) + 
		       xc[2]*(pow(xc[1],2)*(xc[1] - xc[5])*pow(xc[5],2)*(yc[3] - yc[4]) + pow(xc[4],3)*(pow(xc[5],2)*(yc[1] - yc[3]) + pow(xc[1],2)*(yc[3] - yc[5])) + 
			   pow(xc[3],2)*(pow(xc[5],3)*(yc[1] - yc[4]) + pow(xc[1],3)*(yc[4] - yc[5]) + pow(xc[4],3)*(-yc[1] + yc[5])) + pow(xc[4],2)*(pow(xc[5],3)*(-yc[1] + yc[3]) + pow(xc[1],3)*(-yc[3] + yc[5])) + 
			   pow(xc[3],3)*(pow(xc[5],2)*(-yc[1] + yc[4]) + pow(xc[4],2)*(yc[1] - yc[5]) + pow(xc[1],2)*(-yc[4] + yc[5]))) + 
		       pow(xc[2],2)*(-(xc[1]*(xc[1] - xc[5])*xc[5]*(xc[1] + xc[5])*(yc[3] - yc[4])) + xc[4]*(pow(xc[5],3)*(yc[1] - yc[3]) + pow(xc[1],3)*(yc[3] - yc[5])) + pow(xc[4],3)*(-(xc[5]*yc[1]) - xc[1]*yc[3] + xc[5]*yc[3] + xc[1]*yc[5]) + 
				    pow(xc[3],3)*(xc[5]*yc[1] + xc[1]*yc[4] - xc[5]*yc[4] - xc[1]*yc[5] + xc[4]*(-yc[1] + yc[5])) + xc[3]*(pow(xc[5],3)*(-yc[1] + yc[4]) + pow(xc[4],3)*(yc[1] - yc[5]) + pow(xc[1],3)*(-yc[4] + yc[5]))))/
		    ((xc[1] - xc[2])*(xc[1] - xc[3])*(xc[2] - xc[3])*(xc[1] - xc[4])*(xc[2] - xc[4])*(xc[3] - xc[4])*(xc[1] - xc[5])*(xc[2] - xc[5])*(xc[3] - xc[5])*(xc[4] - xc[5]));


		  yc[0]=AA+BB*xc[0]+CC*xc[0]*xc[0]+DD*xc[0]*xc[0]*xc[0]+EE*xc[0]*xc[0]*xc[0]*xc[0];

		}


	    

		if(WHICHTOEXTRAP==1){
		  // convert back to B^\phi
		  prim[i][j][B3]=(yc[0]/geom.gcov[0][0]-(prim[i][j][B1]*geom.gcov[3][1]+prim[i][j][B2]*geom.gcov[3][2]))/geom.gcov[3][3];
		}
		else if(WHICHTOEXTRAP==0){
		  prim[i][j][B3]=yc[0];
		}
		else if(WHICHTOEXTRAP==2){
		  prim[i][j][B2]=yc[0];
		}
		else if(WHICHTOEXTRAP==3){
		  prim[i][j][B1]=yc[0];
		}
	      }
	      else{// 0th order extrap
		if((WHICHTOEXTRAP==1)||(WHICHTOEXTRAP==0)){
		  if(WHICHFIXED==0){
		    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  else if(WHICHFIXED==1){
		    // relaxing B^\phi , fixing B^r, and solving for B^\theta
		    prim[i][j][B1]=p[i][j][B1];
		  }
		  else if(WHICHFIXED==2){
		    for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		  }

		  prim[i][j][B3]=prim[ri][rj][B3];
		}
		else if(WHICHTOEXTRAP==2){
		  if(WHICHFIXED==0){
		    for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		    for(k=B3;k<=B3;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  else if(WHICHFIXED==1){
		    // relaxing B^\theta , fixing B^r, and solving for B^\phi
		    prim[i][j][B1]=p[i][j][B1];
		  }
		  else if(WHICHFIXED==2){
		    for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		  }

		  prim[i][j][B2]=prim[ri][rj][B2];
		}
		else if(WHICHTOEXTRAP==3){
		  if(WHICHFIXED==0){
		    for(k=B2;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
		    for(k=B3;k<=B3;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  else if(WHICHFIXED==1){
		    // relaxing B^r , fixing B^\theta, and solving for B^\phi
		    prim[i][j][B2]=p[i][j][B2];
		  }
		  else if(WHICHFIXED==2){
		    for(k=B1;k<=B1;k++)  prim[i][j][k] = p[i][j][k];
		  }
		  prim[i][j][B1]=prim[ri][rj][B1];
		}
	      }

	      //	    dualfprintf(fail_file,"%g %g %g %g\n",prim[i][j][B3],yc[1],yc[2],yc[3]);

	      //  dualfprintf(fail_file,"yc[1]=%g yc[2]=%g yc[3]=%g \n xc[1]=%g xc[2]=%g xc[3]=%g \n AA=%g BB=%g CC=%g ans=%g\n",yc[1],yc[2],yc[3],xc[1],xc[2],xc[3],AA,BB,CC,prim[i][j][B3]);
	      } // end over extraps loop
	    } // end else doing extrapolation

	    // if WHICHFIXED==0, then done, else need to solve for other component
	    if(WHICHFIXED==1){
	      Bsq=myBsq[i*N2+j];
	      get_geometry(i, j, CENT, &geom);
	      coord(i,j,CENT,X);
	      bl_coord(X,&r,&th);

	      if(WHICHTOEXTRAP==0){
		// then solve for B^\theta
		// B^2 was fixed
		plusminus=1.0; // always positive for this case (except across poles, which is not assigned here)
		det=pow(2.*prim[i][j][B1]*geom.gcov[2-1][3-1] + 2.*prim[i][j][B3]*geom.gcov[3-1][4-1],2) - 4.*geom.gcov[3-1][3-1]*(-1.*Bsq + pow(prim[i][j][B1],2)*geom.gcov[2-1][2-1] + 2.*prim[i][j][B1]*prim[i][j][B3]*geom.gcov[2-1][4-1] + pow(prim[i][j][B3],2)*geom.gcov[4-1][4-1]);
		if(det<0.0) det=0.0;
		prim[i][j][B2]=(0.5*(-2.*prim[i][j][B1]*geom.gcov[2-1][3-1] - 2.*prim[i][j][B3]*geom.gcov[3-1][4-1] + plusminus*sqrt(det)))/geom.gcov[3-1][3-1];
	      }
	      else if((WHICHTOEXTRAP>=2)&&(WHICHTOEXTRAP<=3)){
		// then solve for B^\phi
		// B^2 was fixed
		// for positive rotation and dipolar field
		if(th<M_PI*0.5) plusminus=-1.0;
		else plusminus=1.0;
		det=pow(2.*prim[i][j][B1]*geom.gcov[2-1][4-1] + 2.*prim[i][j][B2]*geom.gcov[3-1][4-1],2) - 4.*(-1.*Bsq + pow(prim[i][j][B1],2)*geom.gcov[2-1][2-1] + 2.*prim[i][j][B1]*prim[i][j][B2]*geom.gcov[2-1][3-1] + pow(prim[i][j][B2],2)*geom.gcov[3-1][3-1])*geom.gcov[4-1][4-1];
		if(det<0.0) det=0.0;
		prim[i][j][B3] = (0.5*(-2.*prim[i][j][B1]*geom.gcov[2-1][4-1] - 2.*prim[i][j][B2]*geom.gcov[3-1][4-1] + plusminus*sqrt(det)))/geom.gcov[4-1][4-1];

		//		for(k=1;k<=3;k++) dualfprintf(fail_file,"i=%d j=%d : B[%d]=%g\n",i,j,k,prim[i][j][B1+k-1]);
		
	      }
	      else{
		dualfprintf(fail_file,"No defined WHICHFIXED=%d with WHICHTOEXTRAP=%d\n",WHICHFIXED,WHICHTOEXTRAP);
		myexit(1);
	      }
	    } // end if WHICHFIXED==1

	  }// end loop over i
	  





	  // generic -- determine velocity once field is determined

	  for(i=-1;i>=-N1BND;i--){
	    get_geometry(i, j, CENT, &geom);

	    Bcon[0]=0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    Bcon[3]=prim[i][j][B3];

	    lower(Bcon,&geom,Bcov);
	    Bsq=0.0;
	    for(k=0;k<=3;k++) Bsq+=Bcon[k]*Bcov[k];



	    vcon[1]=-Bcon[1]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
	    vcon[2]=-Bcon[2]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;
	    vcon[3]=Omegastar-Bcon[3]*(Bcov[0]+Omegastar*Bcov[3])/Bsq;

	    //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"Bcon[%d]=%21.15g Bcov[%d]=%21.15g\n",k,Bcon[k],k,Bcov[k]);

	    //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"vcon[%d]=%21.15g Bsq=%21.15g\n",k,vcon[k],Bsq);

	    //	    for(k=1;k<=3;k++) if(isnan(vcon[k])){
	    //	      dualfprintf(fail_file,"nan encountered: k=%d vcon=%21.15g\n",k,vcon[k]);
	    //	    }
	    // sticks primitive velocity into prim[i][j][U1->U3]
	    if(vcon2pr(WHICHVEL,vcon,&geom,prim[i][j])>=1){
	      dualfprintf(fail_file,"vcon2pr(rescaled): v^\phi[rescaled]=%g u^t=%g : v^\phi[fixed]=%g\n",ucon[PH]/ucon[TT],ucon[TT],vcon[PH]);
	      return(1);
	    }

	    for(k=U1;k<=B3;k++) if(isnan(prim[i][j][k])){
	      dualfprintf(fail_file,"nan encountered: k=%d prim=%21.15g\n",k,prim[i][j][k]);
	    }
	    
	    //	    dualfprintf(fail_file,"i=%d j=%d\n",geom.i,geom.j);


	  }



#elif(EXTRAP==0)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    for(k=RHO;k<=UU;k++)  prim[i][j][k] = prim[ri][rj][k];


	    for(k=U1;k<=U3;k++)  prim[i][j][k] = prim[ri][rj][k];

	    // fixed field
	    //	    for(k=B1;k<=B3;k++)  prim[i][j][k] = p[i][j][k];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];

	    for(k=B3;k<=B3;k++)  prim[i][j][k] = prim[ri][rj][k];


	    get_geometry(i, j, CENT, &geom);

	    Bcon[0]=0.0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    Bcon[3]=prim[i][j][B3];
	    lower(Bcon,&geom,Bcov);


	    pr2ucon(WHICHVEL,prim[i][j],&geom,ucon);

	    if(Bcov[TH]!=0.0){
	      vcon[RR]=ucon[RR]/ucon[TT];
	      vcon[TH]=-(vcon[RR]*Bcov[RR]+vcon[PH]*Bcov[PH]+Bcov[TT])/(Bcov[TH]);
	    }
	    else{
	      vcon[TH]=ucon[TH]/ucon[TT];
	      vcon[RR]=-(vcon[TH]*Bcov[TH]+vcon[PH]*Bcov[PH]+Bcov[TT])/(Bcov[RR]);
	    }
	    //vcon[TH]=ucon[TH]/ucon[TT]; // "" for this component
	    vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	    vcon2pr(WHICHVEL, vcon,&geom,prim[i][j]); // get in terms of primitive velocity

	    /*
	    // outflow field, NS field completely reconnects through surface
	    for(k=B3;k<=B3;k++){
	      //	      prim[i][j][k] = prim[ri][rj][k];
	      if(i==-1) prim[i][j][k] = prim[ri][rj][k] -(prim[ri+1][rj][k]- prim[ri][rj][k]) ;
	      else if(i==-2) prim[i][j][k] = prim[ri][rj][k]-2.0(prim[ri+2][rj][k]- prim[ri+1][rj][k]) ;
	      else if(i==-3) prim[i][j][k] = prim[ri][rj][k]-3.0(prim[ri+3][rj][k]- prim[ri+2][rj][k]) ;

	    }
	    */
	  }
#elif(EXTRAP==111)

	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    for(k=RHO;k<=UU;k++){
	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    }
	    // set  	    ucon[TT,etc.]
	    vcon[RR]=0; // surface that completely dissipates normal direction momentum
	    vcon[TH]=0; // "" for this component
	    vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	    get_geometry(i, j, CENT, &geom);
	    vcon2pr(WHICHVEL,vcon,&geom,prim[i][j]); // get in terms of primitive velocity

	    /*
	    // outflow field, NS field completely reconnects through surface
	    for(k=B3;k<=B3;k++){
	      //	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	      if(i==-1) prim[i][j][k] = prim[ri][rj][k] -(prim[ri+1][rj][k] * gdet[ri+1][rj][CENT]- prim[ri][rj][k] * gdet[ri][rj][CENT])/gdet[i][j][CENT] ;
	      else if(i==-2) prim[i][j][k] = prim[ri][rj][k] -2.0*(prim[ri+2][rj][k] * gdet[ri+2][rj][CENT]- prim[ri+1][rj][k] * gdet[ri+1][rj][CENT])/gdet[i][j][CENT] ;
	      else if(i==-3) prim[i][j][k] = prim[ri][rj][k] -3.0*(prim[ri+3][rj][k] * gdet[ri+3][rj][CENT]- prim[ri+2][rj][k] * gdet[ri+2][rj][CENT])/gdet[i][j][CENT] ;
	    }
	    */
	    // fixed field
	    for(k=B1;k<=B3;k++)  prim[i][j][k] = p[i][j][k];

	  }
#elif(EXTRAP==122)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--){
	    PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k]=prim[ri][rj][k];
	  }
	  //	  dualfprintf(fail_file,"got here2\n");
	  // now rescaled[][NPR] contains interpolated quantities as primitive-type
	  // these are guarenteed to preserve v.B=0 and E.B=0, can't screw that up

	  // now work with them

	  // have B^\phi copied from above
	  // Force B^r and B^\theta to be analytic solution
	  for(i=-1;i>=-N1BND;i--){
	    //	    prim[i][j][B3]=rescaled[i+N1BND][B3];
	    //test
	    //	    for(k=B1;k<=B3;k++) prim[i][j][k]=rescaled[i+N1BND][k];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
	  }

	  // forced
	  vu3=vold[3]=vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	  
	  for(i=-1;i>=-N1BND;i--){
	    get_geometry(i, j, CENT, &geom);
	    
	    //	    Bcon[0]=0.0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    //	    Bcon[3]=prim[i][j][B3];
	    //	    lower(Bcon,&geom,Bcov);

	    // get u^t from interpolation
	    if(ucon_calc(rescaled[i+N1BND],&geom,ucon)>=1){
	      dualfprintf(fail_file,"Shouldn't fail: ucon_calc(rescaled)1\n");
	      return(1);
	    }

	    // fix uu2
	    ucon[2]=0.0;
	    ucon[1]=0.0;

	    // 2 possible u^t's for new v^\phi
	    
	    uu0a=(0.5*(2.*geom.gcov[1-1][2-1]*ucon[1] + 2.*geom.gcov[1-1][3-1]*ucon[2] + 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 + 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3 - 
			1.*sqrt(pow(-2.*geom.gcov[1-1][2-1]*ucon[1] - 2.*geom.gcov[1-1][3-1]*ucon[2] - 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 - 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3,2) - 
			       4.*(-1. - 1.*geom.gcov[2-1][2-1]*pow(ucon[1],2) - 2.*geom.gcov[2-1][3-1]*ucon[1]*ucon[2] - 1.*geom.gcov[3-1][3-1]*pow(ucon[2],2))*
			       (-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2)))))/(-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2));



	    uu0b=(0.5*(2.*geom.gcov[1-1][2-1]*ucon[1] + 2.*geom.gcov[1-1][3-1]*ucon[2] + 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 + 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3 + 
		       sqrt(pow(-2.*geom.gcov[1-1][2-1]*ucon[1] - 2.*geom.gcov[1-1][3-1]*ucon[2] - 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 - 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3,2) - 
			    4.*(-1. - 1.*geom.gcov[2-1][2-1]*pow(ucon[1],2) - 2.*geom.gcov[2-1][3-1]*ucon[1]*ucon[2] - 1.*geom.gcov[3-1][3-1]*pow(ucon[2],2))*
			    (-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2)))))/(-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2));


	    // pick closest u^t
	    if(fabs(ucon[0]-uu0a)<fabs(ucon[0]-uu0b)) ucon[0]=uu0a;
	    else ucon[0]=uu0b;

	    // fix new u^\phi
	    ucon[PH]=vcon[PH]*ucon[0];
	    
	    lower(ucon,&geom,ucov);

	    if(ucov[3]!=0.0){
	      prim[i][j][B3]=-(ucov[1]*Bcon[1]+ucov[2]*Bcon[2])/ucov[3];
	    }
	    else prim[i][j][B3]=0.0;


	    // sticks primitive velocity into prim[i][j][U1->U3]
	    ucon2pr(WHICHVEL,ucon,&geom,prim[i][j]);


	  }



#elif(EXTRAP==20022)
	  ri=0;
	  rj=j;
	  //	  dualfprintf(fail_file,"got here1\n");
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  //	  dualfprintf(fail_file,"got here1.5\n");
	  //	  PBOUNDLOOP dualfprintf(fail_file,"list: k=%d prim=%g prescale=%g\n",k,prim[ri][rj][k],prescale[k]);
	  for(i=-1;i>=-N1BND;i--){
	    // set guess and get copied version
	    PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    if(rescale(-1,1,rescaled[i+N1BND],&geom,prescale)>=1){
	      dualfprintf(fail_file,"failed to interpolate conserved quantity, using copy of primitives\n");
	      failed=0; // don't let fail
	      // then interpolate of conserved quantities didn't work out, so revert to copied
	      PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k];
	    }
	  }
	  //	  dualfprintf(fail_file,"got here2\n");
	  // now rescaled[][NPR] contains interpolated quantities as primitive-type
	  // these are guarenteed to preserve v.B=0 and E.B=0, can't screw that up

	  // now work with them

	  // have B^\phi copied from above
	  // Force B^r and B^\theta to be analytic solution
	  for(i=-1;i>=-N1BND;i--){
	    prim[i][j][B3]=rescaled[i+N1BND][B3];
	    //test
	    //	    for(k=B1;k<=B3;k++) prim[i][j][k]=rescaled[i+N1BND][k];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
	  }

	  // forced
	  vu3=vold[3]=vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	  
	  for(i=-1;i>=-N1BND;i--){
	    get_geometry(i, j, CENT, &geom);
	    
	    //	    Bcon[0]=0.0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    //	    Bcon[3]=prim[i][j][B3];
	    //	    lower(Bcon,&geom,Bcov);

	    // get u^t from interpolation
	    if(ucon_calc(rescaled[i+N1BND],&geom,ucon)>=1){
	      dualfprintf(fail_file,"Shouldn't fail: ucon_calc(rescaled)1\n");
	      return(1);
	    }

	    // 2 possible u^t's for new v^\phi
	    
	    uu0a=(0.5*(2.*geom.gcov[1-1][2-1]*ucon[1] + 2.*geom.gcov[1-1][3-1]*ucon[2] + 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 + 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3 - 
			1.*sqrt(pow(-2.*geom.gcov[1-1][2-1]*ucon[1] - 2.*geom.gcov[1-1][3-1]*ucon[2] - 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 - 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3,2) - 
			       4.*(-1. - 1.*geom.gcov[2-1][2-1]*pow(ucon[1],2) - 2.*geom.gcov[2-1][3-1]*ucon[1]*ucon[2] - 1.*geom.gcov[3-1][3-1]*pow(ucon[2],2))*
			       (-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2)))))/(-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2));



	    uu0b=(0.5*(2.*geom.gcov[1-1][2-1]*ucon[1] + 2.*geom.gcov[1-1][3-1]*ucon[2] + 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 + 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3 + 
		       sqrt(pow(-2.*geom.gcov[1-1][2-1]*ucon[1] - 2.*geom.gcov[1-1][3-1]*ucon[2] - 2.*geom.gcov[2-1][4-1]*ucon[1]*vu3 - 2.*geom.gcov[3-1][4-1]*ucon[2]*vu3,2) - 
			    4.*(-1. - 1.*geom.gcov[2-1][2-1]*pow(ucon[1],2) - 2.*geom.gcov[2-1][3-1]*ucon[1]*ucon[2] - 1.*geom.gcov[3-1][3-1]*pow(ucon[2],2))*
			    (-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2)))))/(-1.*geom.gcov[1-1][1-1] - 2.*geom.gcov[1-1][4-1]*vu3 - 1.*geom.gcov[4-1][4-1]*pow(vu3,2));


	    // pick closest u^t
	    if(fabs(ucon[0]-uu0a)<fabs(ucon[0]-uu0b)) ucon[0]=uu0a;
	    else ucon[0]=uu0b;

	    // fix new u^\phi
	    ucon[PH]=vcon[PH]*ucon[0];
	    
	    lower(ucon,&geom,ucov);

	    if(ucov[3]!=0.0){
	      prim[i][j][B3]=-(ucov[1]*Bcon[1]+ucov[2]*Bcon[2])/ucov[3];
	    }
	    else prim[i][j][B3]=0.0;


	    // sticks primitive velocity into prim[i][j][U1->U3]
	    ucon2pr(WHICHVEL,ucon,&geom,prim[i][j]);


	  }



#elif(EXTRAP==222222)
	  ri=0;
	  rj=j;
	  //	  dualfprintf(fail_file,"got here1\n");
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  //	  dualfprintf(fail_file,"got here1.5\n");
	  //	  PBOUNDLOOP dualfprintf(fail_file,"list: k=%d prim=%g prescale=%g\n",k,prim[ri][rj][k],prescale[k]);
	  for(i=-1;i>=-N1BND;i--){
	    // set guess and get copied version
	    PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    if(rescale(-1,1,rescaled[i+N1BND],&geom,prescale)>=1){
	      dualfprintf(fail_file,"failed to interpolate conserved quantity, using copy of primitives\n");
	      failed=0; // don't let fail
	      // then interpolate of conserved quantities didn't work out, so revert to copied
	      PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k];
	    }
	  }
	  //	  dualfprintf(fail_file,"got here2\n");
	  // now rescaled[][NPR] contains interpolated quantities as primitive-type
	  // these are guarenteed to preserve v.B=0 and E.B=0, can't screw that up

	  // now work with them

	  // have B^\phi copied from above
	  // Force B^r and B^\theta to be analytic solution
	  for(i=-1;i>=-N1BND;i--){
	    prim[i][j][B3]=rescaled[i+N1BND][B3];
	    //test
	    //	    for(k=B1;k<=B3;k++) prim[i][j][k]=rescaled[i+N1BND][k];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
	  }

	  // forced
	  vu3=vold[3]=vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	  
	  for(i=-1;i>=-N1BND;i--){
	    get_geometry(i, j, CENT, &geom);
	    
	    Bcon[0]=0.0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    Bcon[3]=prim[i][j][B3];
	    lower(Bcon,&geom,Bcov);

	    // get u^t from interpolation
	    if(ucon_calc(rescaled[i+N1BND],&geom,ucon)>=1){
	      dualfprintf(fail_file,"Shouldn't fail: ucon_calc(rescaled)1\n");
	      return(1);
	    }

	    
	    // minimum u^t allowed


	    utmin=sqrt((-1.*pow(Bcov[2],2)*geom.gcov[2-1][2-1] + 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] - 1.*pow(Bcov[1],2)*geom.gcov[3-1][3-1])/
	      ((-1.*pow(geom.gcov[2-1][3-1],2) + geom.gcov[2-1][2-1]*geom.gcov[3-1][3-1])*pow(Bcov[0] + Bcov[3]*vu3,2) + 
		   2.*Bcov[1]*(Bcov[0] + Bcov[3]*vu3)*(-1.*geom.gcov[3-1][3-1]*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) + geom.gcov[2-1][3-1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3)) + 
		   pow(Bcov[2],2)*(-1.*pow(geom.gcov[1-1][2-1],2) + geom.gcov[1-1][1-1]*geom.gcov[2-1][2-1] - 2.*geom.gcov[1-1][2-1]*geom.gcov[2-1][4-1]*vu3 + 
				   vu3*(2.*geom.gcov[1-1][4-1]*geom.gcov[2-1][2-1] - 1.*pow(geom.gcov[2-1][4-1],2)*vu3 + geom.gcov[2-1][2-1]*geom.gcov[4-1][4-1]*vu3)) + 
		   pow(Bcov[1],2)*(-1.*pow(geom.gcov[1-1][3-1],2) + geom.gcov[1-1][1-1]*geom.gcov[3-1][3-1] - 2.*geom.gcov[1-1][3-1]*geom.gcov[3-1][4-1]*vu3 + 
				   vu3*(2.*geom.gcov[1-1][4-1]*geom.gcov[3-1][3-1] - 1.*pow(geom.gcov[3-1][4-1],2)*vu3 + geom.gcov[3-1][3-1]*geom.gcov[4-1][4-1]*vu3)) + 
		   2.*Bcov[2]*((Bcov[0] + Bcov[3]*vu3)*(geom.gcov[2-1][3-1]*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) - 1.*geom.gcov[2-1][2-1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3)) + 
			       Bcov[1]*(-1.*geom.gcov[2-1][3-1]*(geom.gcov[1-1][1-1] + 2.*geom.gcov[1-1][4-1]*vu3) + geom.gcov[1-1][2-1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) + 
					vu3*(geom.gcov[1-1][3-1]*geom.gcov[2-1][4-1] + geom.gcov[2-1][4-1]*geom.gcov[3-1][4-1]*vu3 - 1.*geom.gcov[2-1][3-1]*geom.gcov[4-1][4-1]*vu3)))));



	    

	    if(utmin>ucon[0]) ucon[0]=utmin+1E-10;


	    //	    dualfprintf(fail_file,"i=%d j=%d\n",geom.i,geom.j);
	    //	    dualfprintf(fail_file,"ut=%21.15g utmin=%21.15g mytop=%21.15g mybottom=%21.15g\n",ucon[0],utmin,mytop,mybottom);


	    uu0=ucon[0];
	    vu1o=vold[1]=ucon[1]/ucon[0];
	    vu2o=vold[2]=ucon[2]/ucon[0];


	    // determine new v's



	    avu1free1=(-1.*Bcov[1]*geom.gcov[3-1][3-1]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3) - 1.*pow(Bcov[2],2)*pow(uu0,2)*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) + 
		       Bcov[2]*pow(uu0,2)*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3)) + 
		       0.5*sqrt(4.*pow(uu0,4)*pow(Bcov[1]*geom.gcov[3-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + pow(Bcov[2],2)*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) - 
						  1.*Bcov[2]*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3)),2) - 
				4.*(pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2)*
				(geom.gcov[3-1][3-1]*pow(uu0,2)*pow(Bcov[0] + Bcov[3]*vu3,2) - 2.*Bcov[2]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3)*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) + 
				 pow(Bcov[2],2)*(1. + pow(uu0,2)*(geom.gcov[1-1][1-1] + vu3*(2.*geom.gcov[1-1][4-1] + geom.gcov[4-1][4-1]*vu3))))))/
	      ((pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2));

	    avu1free2=(-1.*(Bcov[1]*geom.gcov[3-1][3-1]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3) + pow(Bcov[2],2)*pow(uu0,2)*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) - 
			    1.*Bcov[2]*pow(uu0,2)*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3)) + 
			    0.5*sqrt(4.*pow(uu0,4)*pow(Bcov[1]*geom.gcov[3-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + pow(Bcov[2],2)*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) - 
						       1.*Bcov[2]*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[1]*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3)),2) - 
				     4.*(pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2)*
				     (geom.gcov[3-1][3-1]*pow(uu0,2)*pow(Bcov[0] + Bcov[3]*vu3,2) - 2.*Bcov[2]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3)*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) + 
				      pow(Bcov[2],2)*(1. + pow(uu0,2)*(geom.gcov[1-1][1-1] + vu3*(2.*geom.gcov[1-1][4-1] + geom.gcov[4-1][4-1]*vu3)))))))/
	      ((pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2));


	    avu2free1=(-1.*Bcov[2]*geom.gcov[2-1][2-1]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3) - 1.*pow(Bcov[1],2)*pow(uu0,2)*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) + 
		       Bcov[1]*pow(uu0,2)*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[2]*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3)) + 
		       0.5*sqrt(4.*pow(uu0,4)*pow(Bcov[2]*geom.gcov[2-1][2-1]*(Bcov[0] + Bcov[3]*vu3) + pow(Bcov[1],2)*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) - 
						  1.*Bcov[1]*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[2]*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3)),2) - 
				4.*(pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2)*
				(geom.gcov[2-1][2-1]*pow(uu0,2)*pow(Bcov[0] + Bcov[3]*vu3,2) - 2.*Bcov[1]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3)*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) + 
				 pow(Bcov[1],2)*(1. + pow(uu0,2)*(geom.gcov[1-1][1-1] + vu3*(2.*geom.gcov[1-1][4-1] + geom.gcov[4-1][4-1]*vu3))))))/
	      ((pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2));



	    avu2free2=(-1.*(Bcov[2]*geom.gcov[2-1][2-1]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3) + pow(Bcov[1],2)*pow(uu0,2)*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) - 
			    1.*Bcov[1]*pow(uu0,2)*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[2]*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3)) + 
			    0.5*sqrt(4.*pow(uu0,4)*pow(Bcov[2]*geom.gcov[2-1][2-1]*(Bcov[0] + Bcov[3]*vu3) + pow(Bcov[1],2)*(geom.gcov[1-1][3-1] + geom.gcov[3-1][4-1]*vu3) - 
						       1.*Bcov[1]*(geom.gcov[2-1][3-1]*(Bcov[0] + Bcov[3]*vu3) + Bcov[2]*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3)),2) - 
				     4.*(pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2)*
				     (geom.gcov[2-1][2-1]*pow(uu0,2)*pow(Bcov[0] + Bcov[3]*vu3,2) - 2.*Bcov[1]*pow(uu0,2)*(Bcov[0] + Bcov[3]*vu3)*(geom.gcov[1-1][2-1] + geom.gcov[2-1][4-1]*vu3) + 
				      pow(Bcov[1],2)*(1. + pow(uu0,2)*(geom.gcov[1-1][1-1] + vu3*(2.*geom.gcov[1-1][4-1] + geom.gcov[4-1][4-1]*vu3)))))))/
	      ((pow(Bcov[2],2)*geom.gcov[2-1][2-1] - 2.*Bcov[1]*Bcov[2]*geom.gcov[2-1][3-1] + pow(Bcov[1],2)*geom.gcov[3-1][3-1])*pow(uu0,2));


	    //	    dualfprintf(fail_file,"%g %g %g %g\n",avu1free1,avu1free2,avu2free1,avu2free2);


	    if(Bcov[RR]!=0.0){
	      // then copy v^\theta and fix v^r

	      // use avu2free's
	      
	      // find free quantity v^\theta

	      // get version closest to old version
	      if(fabs(avu2free1-vu2o)<fabs(avu2free2-vu2o)) vcon[TH]=avu2free1;
	      else vcon[TH]=avu2free2;

	      // now fix v^r, which should be equal to vcon[RR]=avu2free*vu1o (for 1 or 2, whatever used above)
	      vcon[RR]=-(vcon[TH]*Bcov[TH]+vcon[PH]*Bcov[PH]+Bcov[TT])/(Bcov[RR]);

	    }
	    else if(Bcov[TH]!=0.0){
	      // then copy v^r and fix v^\theta

	      // get version closest to old version
	      if(fabs(avu1free1-vu1o)<fabs(avu1free2-vu1o)) vcon[RR]=avu1free1;
	      else vcon[RR]=avu1free2;

	      
	      vcon[TH]=-(vcon[RR]*Bcov[RR]+vcon[PH]*Bcov[PH]+Bcov[TT])/(Bcov[TH]);


	    }
	    else{
	      dualfprintf(fail_file,"So Bcov[TH,RR]=0, but v^\phi!=0, so B_\phi=0?\n");
	      myexit(1);
	    }

	    //	    for(k=1;k<=3;k++) dualfprintf(fail_file,"vold[%d]=%21.15g vnew[%d]=%21.15g\n",k,vold[k],k,vcon[k]);

	    //	    for(k=0;k<=3;k++) dualfprintf(fail_file,"Bcov[%d]=%21.15g\n",k,Bcov[k]);


	    ///////////////
	    //
	    // some checks and convert to primitive
	    //
	    ///////////////

	    for(k=1;k<=3;k++) if(isnan(vcon[k])){
	      dualfprintf(fail_file,"nan encountered: k=%d vcon=%21.15g\n",k,vcon[k]);
	    }
	    // sticks primitive velocity into prim[i][j][U1->U3]
	    if(vcon2pr(WHICHVEL,vcon,&geom,prim[i][j])>=1){
	      dualfprintf(fail_file,"vcon2pr(rescaled): v^\phi[rescaled]=%g u^t=%g : v^\phi[fixed]=%g\n",ucon[PH]/ucon[TT],ucon[TT],vcon[PH]);
	      return(1);
	    }

	    for(k=U1;k<=B3;k++) if(isnan(prim[i][j][k])){
	      dualfprintf(fail_file,"nan encountered: k=%d prim=%21.15g\n",k,prim[i][j][k]);
	    }


	  }



#elif(EXTRAP==11111111111)
	  ri=0;
	  rj=j;
	  //	  dualfprintf(fail_file,"got here1\n");
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  //	  dualfprintf(fail_file,"got here1.5\n");
	  //	  PBOUNDLOOP dualfprintf(fail_file,"list: k=%d prim=%g prescale=%g\n",k,prim[ri][rj][k],prescale[k]);
	  for(i=-1;i>=-N1BND;i--){
	    // set guess and get copied version
	    PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    if(rescale(-1,1,rescaled[i+N1BND],&geom,prescale)>=1){
	      dualfprintf(fail_file,"failed to interpolate conserved quantity, using copy of primitives\n");
	      failed=0; // don't let fail
	      // then interpolate of conserved quantities didn't work out, so revert to copied
	      PBOUNDLOOP rescaled[i+N1BND][k]=copied[i+N1BND][k];
	    }
	  }
	  //	  dualfprintf(fail_file,"got here2\n");
	  // now rescaled[][NPR] contains interpolated quantities as primitive-type
	  // these are guarenteed to preserve v.B=0 and E.B=0, can't screw that up

	  // now work with them

	  // have B^\phi copied from above
	  // Force B^r and B^\theta to be analytic solution
	  for(i=-1;i>=-N1BND;i--){
	    prim[i][j][B3]=rescaled[i+N1BND][B3];
	    //test
	    //	    for(k=B1;k<=B3;k++) prim[i][j][k]=rescaled[i+N1BND][k];
	    for(k=B1;k<=B2;k++)  prim[i][j][k] = p[i][j][k];
	  }

	  // forced
	  vcon[PH]=Omegastar; // surface rotates with angular frequency Omegastar to observer at infinity
	  
	  for(i=-1;i>=-N1BND;i--){
	    get_geometry(i, j, CENT, &geom);
	    
	    Bcon[0]=0.0;
	    Bcon[1]=prim[i][j][B1];
	    Bcon[2]=prim[i][j][B2];
	    Bcon[3]=prim[i][j][B3];
	    lower(Bcon,&geom,Bcov);

	    // get u^t from interpolation
	    if(ucon_calc(rescaled[i+N1BND],&geom,ucon)>=1){
	      dualfprintf(fail_file,"Shouldn't fail: ucon_calc(rescaled)1\n");
	      return(1);
	    }
	    // make sure u^t not any smaller than needed by v^\phi
	    
	    utmin=sqrt(-(-(pow(Bcov[TH],2)*(geom.gcov[2-1][2-1])) + 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][3-1]) - pow(Bcov[RR],2)*(geom.gcov[3-1][3-1])))/
	      sqrt(-(-((pow((geom.gcov[2-1][3-1]),2) - (geom.gcov[2-1][2-1])*(geom.gcov[3-1][3-1]))*pow(Bcov[TT] + Bcov[PH]*vcon[PH],2)) + 
		     2.*Bcov[RR]*(Bcov[TT] + Bcov[PH]*vcon[PH])*(-((geom.gcov[3-1][3-1])*((geom.gcov[1-1][2-1]) + (geom.gcov[2-1][4-1])*vcon[PH])) + (geom.gcov[2-1][3-1])*((geom.gcov[1-1][3-1]) + (geom.gcov[3-1][4-1])*vcon[PH])) + 
		     pow(Bcov[TH],2)*(-pow((geom.gcov[1-1][2-1]),2) + (geom.gcov[1-1][1-1])*(geom.gcov[2-1][2-1]) - 2.*(geom.gcov[1-1][2-1])*(geom.gcov[2-1][4-1])*vcon[PH] + vcon[PH]*(2.*(geom.gcov[1-1][4-1])*(geom.gcov[2-1][2-1]) - pow((geom.gcov[2-1][4-1]),2)*vcon[PH] + (geom.gcov[2-1][2-1])*(geom.gcov[4-1][4-1])*vcon[PH])) + 
		     pow(Bcov[RR],2)*(-pow((geom.gcov[1-1][3-1]),2) + (geom.gcov[1-1][1-1])*(geom.gcov[3-1][3-1]) - 2.*(geom.gcov[1-1][3-1])*(geom.gcov[3-1][4-1])*vcon[PH] + vcon[PH]*(2.*(geom.gcov[1-1][4-1])*(geom.gcov[3-1][3-1]) - pow((geom.gcov[3-1][4-1]),2)*vcon[PH] + (geom.gcov[3-1][3-1])*(geom.gcov[4-1][4-1])*vcon[PH])) + 
		     2.*Bcov[TH]*((Bcov[TT] + Bcov[PH]*vcon[PH])*((geom.gcov[2-1][3-1])*((geom.gcov[1-1][2-1]) + (geom.gcov[2-1][4-1])*vcon[PH]) - (geom.gcov[2-1][2-1])*((geom.gcov[1-1][3-1]) + (geom.gcov[3-1][4-1])*vcon[PH])) + 
				  Bcov[RR]*(-((geom.gcov[2-1][3-1])*((geom.gcov[1-1][1-1]) + 2.*(geom.gcov[1-1][4-1])*vcon[PH])) + (geom.gcov[1-1][2-1])*((geom.gcov[1-1][3-1]) + (geom.gcov[3-1][4-1])*vcon[PH]) + vcon[PH]*((geom.gcov[1-1][3-1])*(geom.gcov[2-1][4-1]) + (geom.gcov[2-1][4-1])*(geom.gcov[3-1][4-1])*vcon[PH] - (geom.gcov[2-1][3-1])*(geom.gcov[4-1][4-1])*vcon[PH])))));
	    
//	    utmin=1.0/sqrt(-geom.gcov[TT][TT]-vcon[PH]*vcon[PH]*geom.gcov[PH][PH]-2.0*vcon[PH]*geom.gcov[TT][PH]);
	    if(ucon[TT]<utmin) ucon[TT]=utmin*(1.0+1E-13); // 1E-13 to avoid roundoff problems

	    //	    dualfprintf(fail_file,"got here3\n");

	    // test
	    //	    vcon[PH]=ucon[PH]/ucon[TT];

	    //	    ucon[PH]=ucon[TT]*vcon[PH]; // redo u^\phi

	    if(Bcov[RR]!=0.0){
	      // then copy v^\theta and fix v^r
	      //vcon[TH]=ucon[TH]/ucon[TT];
	      
	      if(geom.j>totalsize[2]) signsolution=-1.0;
	      else signsolution=1.0;
	      //	      signsolution=-1.0;

	      det=(pow(-2.*Bcov[RR]*Bcov[TH]*(geom.gcov[1-1][2-1])*pow(ucon[TT],2.0) + 2.*pow(Bcov[RR],2.0)*(geom.gcov[1-1][3-1])*pow(ucon[TT],2.0) + 
         2.*Bcov[TT]*Bcov[TH]*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) - 2.*Bcov[TT]*Bcov[RR]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) + 2.*Bcov[TH]*Bcov[PH]*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0)*vcon[PH] - 
         2.*Bcov[RR]*Bcov[PH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] - 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + 
         2.*pow(Bcov[RR],2.0)*(geom.gcov[3-1][4-1])*pow(ucon[TT],2.0)*vcon[PH],2.0) - 
       4*(pow(Bcov[TH],2.0)*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) - 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) + pow(Bcov[RR],2.0)*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0))*
        (pow(Bcov[RR],2.0) + pow(Bcov[RR],2.0)*(geom.gcov[1-1][1-1])*pow(ucon[TT],2.0) - 2.*Bcov[TT]*Bcov[RR]*(geom.gcov[1-1][2-1])*pow(ucon[TT],2.0) + 
          pow(Bcov[TT],2.0)*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) - 2.*Bcov[RR]*Bcov[PH]*(geom.gcov[1-1][2-1])*pow(ucon[TT],2.0)*vcon[PH] + 
          2.*pow(Bcov[RR],2.0)*(geom.gcov[1-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + 2.*Bcov[TT]*Bcov[PH]*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0)*vcon[PH] - 
          2.*Bcov[TT]*Bcov[RR]*(geom.gcov[2-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + pow(Bcov[PH],2.0)*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0)*pow(vcon[PH],2.0) - 
	 2.*Bcov[RR]*Bcov[PH]*(geom.gcov[2-1][4-1])*pow(ucon[TT],2.0)*pow(vcon[PH],2.0) + pow(Bcov[RR],2.0)*(geom.gcov[4-1][4-1])*pow(ucon[TT],2.0)*pow(vcon[PH],2.0)));

	      sqrtdet=sqrt(det);

	      //	      dualfprintf(fail_file,"det=%21.15g Bd0=%g Bd1=%g Bd2=%g Bd3=%g uu0=%g\n",det,Bcov[TT],Bcov[RR],Bcov[TH],Bcov[PH],ucon[TT]);

	      vcon[TH]=(2.*Bcov[RR]*Bcov[TH]*(geom.gcov[1-1][2-1])*pow(ucon[TT],2.0) - 2.*pow(Bcov[RR],2.0)*(geom.gcov[1-1][3-1])*pow(ucon[TT],2.0) - 2.*Bcov[TT]*Bcov[TH]*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) + 
     2.*Bcov[TT]*Bcov[RR]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) - 2.*Bcov[TH]*Bcov[PH]*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0)*vcon[PH] + 2.*Bcov[RR]*Bcov[PH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] + 
     2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] - 2.*pow(Bcov[RR],2.0)*(geom.gcov[3-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] +
			
			signsolution*sqrtdet)/
   (2.*(pow(Bcov[TH],2.0)*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) - 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) + pow(Bcov[RR],2.0)*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0)));
	      
	      vcon[RR]=-(vcon[TH]*Bcov[TH]+vcon[PH]*Bcov[PH]+Bcov[TT])/(Bcov[RR]);

	    }
	    else if(Bcov[TH]!=0.0){
	      // then copy v^r and fix v^\theta
	      //vcon[RR]=ucon[RR]/ucon[TT];

	      if(geom.j>totalsize[2]) signsolution=-1.0;
	      else signsolution=1.0;

	      det=pow(2.*pow(Bcov[TH],2.0)*(geom.gcov[1-1][2-1])*pow(ucon[TT],2.0) - 
         2.*Bcov[RR]*Bcov[TH]*(geom.gcov[1-1][3-1])*pow(ucon[TT],2.0) - 2.*Bcov[TT]*Bcov[TH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) + 2.*Bcov[TT]*Bcov[RR]*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0) - 
         2.*Bcov[TH]*Bcov[PH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] + 2.*pow(Bcov[TH],2.0)*(geom.gcov[2-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + 
         2.*Bcov[RR]*Bcov[PH]*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] - 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[3-1][4-1])*pow(ucon[TT],2.0)*vcon[PH],2.0) - 
       4*(pow(Bcov[TH],2.0)*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) - 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) + pow(Bcov[RR],2.0)*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0))*
        (pow(Bcov[TH],2.0) + pow(Bcov[TH],2.0)*(geom.gcov[1-1][1-1])*pow(ucon[TT],2.0) - 2.*Bcov[TT]*Bcov[TH]*(geom.gcov[1-1][3-1])*pow(ucon[TT],2.0) + 
          pow(Bcov[TT],2.0)*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0) - 2.*Bcov[TH]*Bcov[PH]*(geom.gcov[1-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] + 
          2.*pow(Bcov[TH],2.0)*(geom.gcov[1-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + 2.*Bcov[TT]*Bcov[PH]*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] - 
          2.*Bcov[TT]*Bcov[TH]*(geom.gcov[3-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + pow(Bcov[PH],2.0)*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0)*pow(vcon[PH],2.0) - 
	 2.*Bcov[TH]*Bcov[PH]*(geom.gcov[3-1][4-1])*pow(ucon[TT],2.0)*pow(vcon[PH],2.0) + pow(Bcov[TH],2.0)*(geom.gcov[4-1][4-1])*pow(ucon[TT],2.0)*pow(vcon[PH],2.0));


	      sqrtdet=sqrt(det);

	      vcon[RR]=(-2.*pow(Bcov[TH],2.0)*(geom.gcov[1-1][2-1])*pow(ucon[TT],2.0) + 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[1-1][3-1])*pow(ucon[TT],2.0) + 2.*Bcov[TT]*Bcov[TH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) - 
     2.*Bcov[TT]*Bcov[RR]*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0) + 2.*Bcov[TH]*Bcov[PH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] - 
     2.*pow(Bcov[TH],2.0)*(geom.gcov[2-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] - 2.*Bcov[RR]*Bcov[PH]*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0)*vcon[PH] + 
     2.*Bcov[RR]*Bcov[TH]*(geom.gcov[3-1][4-1])*pow(ucon[TT],2.0)*vcon[PH] + signsolution*sqrtdet)/
   (2.*(pow(Bcov[TH],2.0)*(geom.gcov[2-1][2-1])*pow(ucon[TT],2.0) - 2.*Bcov[RR]*Bcov[TH]*(geom.gcov[2-1][3-1])*pow(ucon[TT],2.0) + pow(Bcov[RR],2.0)*(geom.gcov[3-1][3-1])*pow(ucon[TT],2.0)));
	      
	      vcon[TH]=-(vcon[RR]*Bcov[RR]+vcon[PH]*Bcov[PH]+Bcov[TT])/(Bcov[TH]);


	    }
	    else{
	      dualfprintf(fail_file,"So Bcov[TH,RR]=0, but v^\phi!=0, so B_\phi=0?\n");
	      myexit(1);
	    }

	    // ucon[RR]=vcon[RR]*ucon[TT];
	    //ucon[TH]=vcon[TH]*ucon[TT];

#if(1)
	    // sticks primitive velocity into prim[i][j][U1->U3]
	    if(vcon2pr(WHICHVEL,vcon,&geom,prim[i][j])>=1){
	      dualfprintf(fail_file,"vcon2pr(rescaled): v^\phi[rescaled]=%g u^t=%g : v^\phi[fixed]=%g\n",ucon[PH]/ucon[TT],ucon[TT],vcon[PH]);
	      return(1);
	    }
#else
	    // sticks primitive velocity into prim[i][j][U1->U3]
	    ucon2pr(WHICHVEL,ucon,&geom,prim[i][j]);
#endif


	  }






#endif







	  for(i=-1;i>=-N1BND;i--){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom, Dt) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom, Dt) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#endif
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom,Dt) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom,Dt)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	  }
	}
      }

    }

    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){
	/* outer r BC: outflow */
      
	for (j = 0; j < N2; j++) {
#if(EXTRAP==0)
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1-1+N1BND;i++)	  PBOUNDLOOP prim[i][j][k] = prim[ri][rj][k];
#elif(EXTRAP==1)
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1-1+N1BND;i++){
	    for(k=RHO;k<=UU;k++){
	      prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    }
	    k=U1; // treat U1 as special
	    prim[i][j][k] = prim[ri][rj][k] * (1. - 2*(i-ri)*dx[1]) ;
	    for(k=U2;k<=U3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. - (i-ri)*dx[1]) ;
	    }
	    k=B1; // treat B1 special
	    prim[i][j][k] = prim[ri][rj][k] * gdet[ri][rj][CENT]/gdet[i][j][CENT] ;
	    for(k=B2;k<=B3;k++){
	      prim[i][j][k] = prim[ri][rj][k] * (1. - (i-ri)*dx[1]) ;
	    }
	  }
#elif(EXTRAP==2)
	  ri=N1-1;
	  rj=j;
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prescale);
	  for(i=N1;i<=N1-1+N1BND;i++){
	    // set guess
	    PBOUNDLOOP prim[i][j][k]=prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    rescale(-1,1,prim[i][j],&geom,prescale);
	  }
#endif

	  for(i=N1;i<=N1-1+N1BND;i++){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom,Dt) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom,Dt) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#endif
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom,Dt) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom, Dt)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
	  }
	}
      }
      /* if fixed BC: do nothing */
    }

    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      for (i = -N1BND; i <=N1-1+N1BND; i++){
	ri=i;
	rj=0;
	for(j=-N2BND;j<=-1;j++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj+(rj-j-1)][k];
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
       for (i = -N1BND; i <=N1-1+N1BND; i++){
	ri=i;
	rj=N2-1;
	for(j=N2;j<=N2-1+N2BND;j++) PBOUNDLOOP  prim[i][j][k] = prim[ri][rj+(rj-j+1)][k];
      }
    }
    /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
    /* inner pole */
    if (mycpupos[2] == 0) {
      for (i = -N1BND; i <= N1-1+N1BND; i++) {
	for (j = -N2BND; j < 0; j++) {
	  if(POSDEFMETRIC==0){
	    // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	    // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
#if(1)
	    prim[i][j][U2] *= 1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= 1.;
	    prim[i][j][B3] *= 1.;
#else
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= -1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= -1.;
#endif
	  }
	  else{
	    prim[i][j][U2] *= 1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= 1.;
	    prim[i][j][B3] *= 1.;
	  }
	}
      }
    }
    /* outer pole */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -N1BND; i <= N1-1+N1BND; i++) {
	for (j = N2; j <= N2-1+N2BND; j++) {
	  // normal axis
#if(1)
	  prim[i][j][U1] *= 1.;
	  prim[i][j][U2] *= 1.;
	  prim[i][j][U3] *= 1.;
	  prim[i][j][B1] *= 1.;
	  prim[i][j][B2] *= 1.;
	  prim[i][j][B3] *= 1.;
#else
	  prim[i][j][U2] *= -1.;
	  prim[i][j][U3] *= -1.;
	  prim[i][j][B2] *= -1.;
	  prim[i][j][B3] *= -1.;

#endif


	  // monopole
#if(0)
	  prim[i][j][U1] *= 1.;
	  prim[i][j][U2] *= -1.;
	  prim[i][j][U3] *= 1.;
	  prim[i][j][B1] *= -1.;
	  prim[i][j][B2] *= 1.;
	  prim[i][j][B3] *= -1.;
#endif
#if(0)
	  // split monopole
	  prim[i][j][U1] *= 1.;
	  prim[i][j][U2] *= -1.;
	  prim[i][j][U3] *= 1.;
	  prim[i][j][B1] *= 1.;
	  prim[i][j][B2] *= -1.;
	  prim[i][j][B3] *= 1.;
#endif
	}
      }
    }

  }// end if stage0 or stagem1

  if (USEMPI) bound_mpi(boundstage, prim);


  firsttime=0;
  return (0);
}
