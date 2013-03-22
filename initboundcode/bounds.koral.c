
#include "decs.h"

int bound_x1dn_radbeamflatinflow(
								 int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
								 int *inboundloop,
								 int *outboundloop,
								 int *innormalloop,
								 int *outnormalloop,
								 int (*inoutlohi)[NUMUPDOWN][NDIM],
								 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
								 int *dosetbc,
								 int enerregion,
								 int *localenerpos
								 );

int bound_radbeam2dbeaminflow(int dir,
						  int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
						  int *inboundloop,
						  int *outboundloop,
						  int *innormalloop,
						  int *outnormalloop,
						  int (*inoutlohi)[NUMUPDOWN][NDIM],
						  int riin, int riout, int rjin, int rjout, int rkin, int rkout,
						  int *dosetbc,
						  int enerregion,
						  int *localenerpos
						  );

int bound_radbeam2dflowinflow(int dir,
							   int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
							   int *inboundloop,
							   int *outboundloop,
							   int *innormalloop,
							   int *outnormalloop,
							   int (*inoutlohi)[NUMUPDOWN][NDIM],
							   int riin, int riout, int rjin, int rjout, int rkin, int rkout,
							   int *dosetbc,
							   int enerregion,
							   int *localenerpos
							   );


int bound_radshadowinflow(int dir,
						  int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
						  int *inboundloop,
						  int *outboundloop,
						  int *innormalloop,
						  int *outnormalloop,
						  int (*inoutlohi)[NUMUPDOWN][NDIM],
						  int riin, int riout, int rjin, int rjout, int rkin, int rkout,
						  int *dosetbc,
						  int enerregion,
						  int *localenerpos
						  );

int bound_radatmbeaminflow(int dir,
						  int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
						  int *inboundloop,
						  int *outboundloop,
						  int *innormalloop,
						  int *outnormalloop,
						  int (*inoutlohi)[NUMUPDOWN][NDIM],
						  int riin, int riout, int rjin, int rjout, int rkin, int rkout,
						  int *dosetbc,
						  int enerregion,
						  int *localenerpos
						  );

int bound_radwallinflow(int dir,
								 int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
								 int *inboundloop,
								 int *outboundloop,
								 int *innormalloop,
								 int *outnormalloop,
								 int (*inoutlohi)[NUMUPDOWN][NDIM],
								 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
								 int *dosetbc,
								 int enerregion,
								 int *localenerpos
								 );


int bound_radbondiinflow(int dir,
							   int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
							   int *inboundloop,
							   int *outboundloop,
							   int *innormalloop,
							   int *outnormalloop,
							   int (*inoutlohi)[NUMUPDOWN][NDIM],
							   int riin, int riout, int rjin, int rjout, int rkin, int rkout,
							   int *dosetbc,
							   int enerregion,
							   int *localenerpos
							   );

int bound_raddot(
                 int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                 int *inboundloop,
                 int *outboundloop,
                 int *innormalloop,
                 int *outnormalloop,
                 int (*inoutlohi)[NUMUPDOWN][NDIM],
                 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                 int *dosetbc,
                 int enerregion,
                 int *localenerpos
                 );

int bound_radnt(int dir,
                int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                int *inboundloop,
                int *outboundloop,
                int *innormalloop,
                int *outnormalloop,
                int (*inoutlohi)[NUMUPDOWN][NDIM],
                int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                int *dosetbc,
                int enerregion,
                int *localenerpos
                );


int bound_x1dn_cylaxis(
                       int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                       int *inboundloop,
                       int *outboundloop,
                       int *innormalloop,
                       int *outnormalloop,
                       int (*inoutlohi)[NUMUPDOWN][NDIM],
                       int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                       int *dosetbc,
                       int enerregion,
                       int *localenerpos
                       );


int bound_x1up_radcylbeam(
                          int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                          int *inboundloop,
                          int *outboundloop,
                          int *innormalloop,
                          int *outnormalloop,
                          int (*inoutlohi)[NUMUPDOWN][NDIM],
                          int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                          int *dosetbc,
                          int enerregion,
                          int *localenerpos
                          );

/* bound array containing entire set of primitive variables */



static int bound_prim_user_general(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int ispstag, int* dirprim, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);



int bound_prim_user_dir(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int dirprim[NPR];
  int pl,pliter;


  // specify location of primitives
  PALLLOOP(pl) dirprim[pl]=CENT;
  //  dualfprintf(fail_file,"start bound_prim\n"); // CHANGINGMARK
  bound_prim_user_general(boundstage, finalstep, boundtime, whichdir, boundvartype, BOUNDPRIMLOC, dirprim, prim);
  //  dualfprintf(fail_file,"end bound_prim\n"); // CHANGINGMARK

  return(0);
}


// assume single user function takes care of primitive locations
int bound_pstag_user_dir(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{

  int dirprim[NPR];
  int pl,pliter;




  if(FLUXB!=FLUXCTSTAG) return(0); // nothing to do


  // specify location of primitives
  PALLLOOP(pl) dirprim[pl]=CENT;
  dirprim[B1]=FACE1;
  dirprim[B2]=FACE2;
  dirprim[B3]=FACE3;

#if(0)
  // Assume JCM setup bound_prim_user_general() correctly, so turned this off
  // GODMARK: assume non-field velocity not set and have to have something reasonable
  // use global values of non-field parts at present time
  // note that bound_pstag() setup loops to be over only B1..B3, but user may violate this and just stick in something so no failures even if not using data
  FULLLOOP PLOOPNOB1(pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(p,i,j,k,pl);
  FULLLOOP PLOOPNOB2(pl) MACP0A1(prim,i,j,k,pl)=MACP0A1(p,i,j,k,pl);
#endif


  // assume before calling this that bound_pstag() setup PLOOPINTERP so only doing B1,B2,B3 (even though user may not respect this in bound_prim_user_general() -- which is ok since non-field quantities in pstag aren't needed -- may be problem if user_general() assumes primitive is reasonable)
  //  dualfprintf(fail_file,"start bound_pstag\n"); // CHANGINGMARK
  bound_prim_user_general(boundstage, finalstep, boundtime, whichdir, boundvartype, BOUNDPSTAGLOC, dirprim, prim);
  //  dualfprintf(fail_file,"end bound_pstag\n"); // CHANGINGMARK



  return (0);
}



// user boundary routine
int bound_prim_user_general(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int ispstag, int* dirprim, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int inboundloop[NDIM];
  int outboundloop[NDIM];
  int innormalloop[NDIM];
  int outnormalloop[NDIM];
  int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
  int riin,riout,rjin,rjout,rkin,rkout;
  int dosetbc[COMPDIM*2];
  int enerregion;
  int *localenerpos;
  int dir;
  int donebc[COMPDIM*2];



  DIRLOOP(dir) donebc[dir]=0; // assume BC not done yet


  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);


  // for CZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  //  enerregion=TRUEGLOBALENERREGION;
  enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
  localenerpos=enerposreg[enerregion];
  // if(WITHINENERREGION(localenerpos,i,j,k)){
  // note that localenerpos[X1DN] sets first evolved cell




  // first do non-directional internal grid assignments overwritting any evolution (ok to do for each whichdir)
  if(WHICHPROBLEM==RADDOT){
    bound_raddot(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
  }
  


  if(whichdir==1){

    if((dosetbc[X1DN] || dosetbc[X1UP]) && (donebc[X1DN]==0 && donebc[X1UP]==0)){
      if( (BCtype[X1DN]==PERIODIC)&&(BCtype[X1UP]==PERIODIC) ){
        bound_x1_periodic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[X1DN]=donebc[X1UP]=1;
      }
    }

    dir=X1DN;
    if(dosetbc[dir] && donebc[dir]==0){
      if(BCtype[dir]==HORIZONOUTFLOW){
        bound_x1dn_outflow(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if((BCtype[dir]==OUTFLOW)||(BCtype[dir]==FIXEDOUTFLOW)||(BCtype[dir]==FREEOUTFLOW)){
        bound_x1dn_outflow_simple(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if((BCtype[dir]==R0SING)||(BCtype[dir]==SYMM)||(BCtype[dir]==ASYMM) ){
        bound_x1dn_sym(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==FIXEDUSEPANALYTIC){
        bound_x1dn_analytic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADBEAMFLATINFLOW){
        bound_x1dn_radbeamflatinflow(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADSHADOWINFLOW){
        bound_radshadowinflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADATMBEAMINFLOW){
        bound_radatmbeaminflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADWALLINFLOW){
        bound_radwallinflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==CYLAXIS){
        bound_x1dn_cylaxis(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else{
        dualfprintf(fail_file,"No x1dn boundary condition specified: %d\n",BCtype[dir]);
        myexit(7598730);
      }
    }


    dir=X1UP;
    if(dosetbc[dir] && donebc[dir]==0){
      if((BCtype[dir]==OUTFLOW)||(BCtype[dir]==FIXEDOUTFLOW)||(BCtype[dir]==FREEOUTFLOW)){
        bound_x1up_outflow_simple(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADBEAM2DFLOWINFLOW){
        bound_radbeam2dflowinflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADATMBEAMINFLOW){
        bound_radatmbeaminflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==FIXEDUSEPANALYTIC){
        bound_x1up_analytic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADBONDIINFLOW){
        bound_radbondiinflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADNTBC){
        bound_radnt(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADCYLBEAMBC){
        bound_x1up_radcylbeam(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else{
        dualfprintf(fail_file,"No x1up boundary condition specified: %d\n",BCtype[dir]);
        myexit(7598731);
      }
    }


  }
  else if(whichdir==2){


    if((dosetbc[X2DN] || dosetbc[X2UP]) && (donebc[X2DN]==0 && donebc[X2UP]==0)){
      if( (BCtype[X2DN]==PERIODIC)&&(BCtype[X2UP]==PERIODIC) ){
        bound_x2_periodic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[X2DN]=donebc[X2UP]=1;
      }
    }


    dir=X2DN;
    if(dosetbc[dir] && donebc[dir]==0){
      if((BCtype[dir]==OUTFLOW)||(BCtype[dir]==FIXEDOUTFLOW)||(BCtype[dir]==FREEOUTFLOW)){
        bound_x2dn_outflow_simple(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==POLARAXIS && special3dspc){
        int whichcall=1;
        bound_x2dn_polaraxis_full3d(whichcall,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if((BCtype[dir]==POLARAXIS)||(BCtype[dir]==SYMM)||(BCtype[dir]==ASYMM) ){
        bound_x2dn_polaraxis(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==FIXEDUSEPANALYTIC){
        bound_x2dn_analytic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim);
        donebc[dir]=1;
      }
      else{
        dualfprintf(fail_file,"No x2dn boundary condition specified: %d\n",BCtype[dir]);
        myexit(7598732);
      }
    }


    dir=X2UP;
    if(dosetbc[dir] && donebc[dir]==0){
      if((BCtype[dir]==OUTFLOW)||(BCtype[dir]==FIXEDOUTFLOW)||(BCtype[dir]==FREEOUTFLOW)){
        bound_x2up_outflow_simple(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==POLARAXIS && special3dspc){
        int whichcall=1;
        bound_x2up_polaraxis_full3d(whichcall,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if((BCtype[dir]==POLARAXIS)||(BCtype[dir]==SYMM)||(BCtype[dir]==ASYMM) ){
        bound_x2up_polaraxis(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==FIXEDUSEPANALYTIC){
        bound_x2up_analytic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADSHADOWINFLOWX2UP){
        bound_radshadowinflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADWALLINFLOW){
        bound_radwallinflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADNTBC){
        // do ASYMM first
        BCtype[dir]=ASYMM;
        bound_x2up_polaraxis(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        // now do anything else that needs to be done
        BCtype[dir]=RADNTBC;       
        bound_radnt(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else{
        dualfprintf(fail_file,"No x2dn boundary condition specified: %d\n",BCtype[dir]);
        myexit(7598733);
      }
    }




  }
  else if(whichdir==3){


    
    if((dosetbc[X3DN] || dosetbc[X3UP]) && (donebc[X3DN]==0 && donebc[X3UP]==0)){
      if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
        bound_x3_periodic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[X3DN]=donebc[X3UP]=1;
      }
    }


    dir=X3DN;
    if(dosetbc[dir] && donebc[dir]==0){
      if((BCtype[dir]==OUTFLOW)||(BCtype[dir]==FIXEDOUTFLOW)||(BCtype[dir]==FREEOUTFLOW)){
		bound_x3dn_outflow_simple(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
		donebc[dir]=1;
      }
      else if(BCtype[dir]==FIXEDUSEPANALYTIC){
        bound_x3dn_analytic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADBEAM2DBEAMINFLOW){
        bound_radbeam2dbeaminflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);	
        donebc[dir]=1;
      }
      else{
        dualfprintf(fail_file,"No x3dn boundary condition specified: %d\n",BCtype[dir]);
        myexit(34672546);
      }
    }


    dir=X3UP;
    if(dosetbc[dir] && donebc[dir]==0){
      if((BCtype[dir]==OUTFLOW)||(BCtype[dir]==FIXEDOUTFLOW)||(BCtype[dir]==FREEOUTFLOW)){
		bound_x3up_outflow_simple(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
		donebc[dir]=1;
      }
      else if(BCtype[dir]==FIXEDUSEPANALYTIC){
        bound_x3up_analytic(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim);
        donebc[dir]=1;
      }
      else{
        dualfprintf(fail_file,"No x3up boundary condition specified: %d\n",BCtype[dir]);
        myexit(34672547);
      }
    }





  }
  else{
    dualfprintf(fail_file,"No such whichdir=%d\n",whichdir);
    myexit(2436262);
  }
  






  return(0);
}





// see interpline.c
int apply_bc_line(int nprlocalstart, int nprlocalend, int*nprlocallist, int doinverse, int iterdir, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])
{
  int flip_y(int nprlocalstart, int nprlocalend, int*nprlocallist, int iterdir, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM]);

  if(doinverse==0){
    flip_y(nprlocalstart,nprlocalend,nprlocallist,iterdir, recontype, bs, be, yin);
  }
  else{
    flip_y(nprlocalstart,nprlocalend,nprlocallist,iterdir, recontype, bs, be, yin);
    flip_y(nprlocalstart,nprlocalend,nprlocallist,iterdir, recontype, bs, be, yout);
  }

  return(0);

}


#include "reconstructeno.h"

int flip_y(int nprlocalstart, int nprlocalend, int*nprlocallist, int iterdir, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM])
{
  int pllocal,pl,myi;


#if( WENO_DIR_FLIP_CONS_SIGN_DN )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterdir == WENO_DIR_FLIP_CONS_SIGN_DN && (recontype == CVT_C2A || recontype == CVT_A2C) && mycpupos[iterdir] == 0 ) { 
    NUMPRIMLOOP(pllocal,pl) 
      for( myi = bs; myi < 0; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif
	
#if( WENO_DIR_FLIP_CONS_SIGN_UP )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterdir == WENO_DIR_FLIP_CONS_SIGN_UP && (recontype == CVT_C2A || recontype == CVT_A2C)  && mycpupos[iterdir] == numbercpu[iterdir] - 1 ) { 
    NUMPRIMLOOP(pllocal,pl) 
      for( myi = N1*(iterdir==1) + N2*(iterdir==2) + N3*(iterdir==3); myi <= be; myi++ ) {
        y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif


  return(0);

}

void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], 
              FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
              FTYPE *p2interp_l, FTYPE *p2interp_r )
{
}

void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, 
                FTYPE (*p2interp)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*dq)[NSTORE2][NSTORE3][NPR2INTERP], 
                FTYPE (*pleft)[NSTORE2][NSTORE3][NPR2INTERP], FTYPE (*pright)[NSTORE2][NSTORE3][NPR2INTERP], 
                FTYPE *p2interp_l, FTYPE *p2interp_r )
{
}


///Called after the MPI boundary routines
// many things here are copied from above
int bound_prim_user_after_mpi_dir(int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR])
{
  int dirprim[NPR];
  int pliter,pl;

  int inboundloop[NDIM];
  int outboundloop[NDIM];
  int innormalloop[NDIM];
  int outnormalloop[NDIM];
  int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
  int riin,riout,rjin,rjout,rkin,rkout;
  int dosetbc[COMPDIM*2];
  int enerregion;
  int *localenerpos;
  int dir;




  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout, dosetbc);


  // for CZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  //  enerregion=TRUEGLOBALENERREGION;
  enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
  localenerpos=enerposreg[enerregion];
  // if(WITHINENERREGION(localenerpos,i,j,k)){
  // note that localenerpos[X1DN] sets first evolved cell


  if(ispstag==BOUNDPRIMLOC){
    // specify location of primitives
    PALLLOOP(pl) dirprim[pl]=CENT;
  }
  
  if(ispstag==BOUNDPSTAGLOC){
    // specify location of primitives
    PALLLOOP(pl) dirprim[pl]=CENT;
    dirprim[B1]=FACE1;
    dirprim[B2]=FACE2;
    dirprim[B3]=FACE3;
  }


  // use post-MPI values to operate poledeath() and/or polesmooth() since otherwise boundary values poledeath uses are not set yet
  if(whichdir==2){

    dir=X2DN;
    if(dosetbc[dir]){
      if(BCtype[dir]==POLARAXIS && special3dspc){
        bound_x2dn_polaraxis_full3d(2,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
      }
    }


    dir=X2UP;
    if(dosetbc[dir]){
      if(BCtype[dir]==POLARAXIS && special3dspc){
        bound_x2up_polaraxis_full3d(2,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
      }
    }

  }


  // can only really check boundaries after MPI is done
  // e.g. periodicx3 and ncpux3==2 and then MPi required to set k<0 and k>=ncpux2*N3 cells
  if(whichdir==1 && N2==1 && N3==1 || N3==1 && whichdir==2 || N3>1 && whichdir==3){ // not completely general conditional
    bound_checks1(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
  }


  return(0);
}















////////////////////////////////////
//
//  Koral specific physical boundary conditions
//
////////////////////////////////////



// X1 lower for radiation beam injection
int bound_x1dn_radbeamflatinflow(
								 int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
								 int *inboundloop,
								 int *outboundloop,
								 int *innormalloop,
								 int *outnormalloop,
								 int (*inoutlohi)[NUMUPDOWN][NDIM],
								 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
								 int *dosetbc,
								 int enerregion,
								 int *localenerpos
								 )

{


#pragma omp parallel  // assume don't require EOS
  {

    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if(BCtype[X1DN]==RADBEAMFLATINFLOW && totalsize[1]>1 && mycpupos[1] == 0 ) {



      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



        ri=riin;
        rj=j;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
        LOOPBOUND1IN{

    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);


          if(ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);

            // set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates
            int whichvel;
            whichvel=VEL4;
            int whichcoord;
            whichcoord=CARTMINKMETRIC2;

            extern FTYPE RADBEAMFLAT_FRATIO,RADBEAMFLAT_ERAD, RADBEAMFLAT_RHO, RADBEAMFLAT_UU;

            FTYPE ERADAMB=(RADBEAMFLAT_ERAD/RHOBAR);
            FTYPE uradx=1.0/sqrt(1.0 - RADBEAMFLAT_FRATIO*RADBEAMFLAT_FRATIO); // radiation 4-velocity
            FTYPE ERADINJ;
            ERADINJ=1000.0*(RADBEAMFLAT_ERAD/RHOBAR);


            //primitives in whichvel,whichcoord
            if(V[2]>.4 && V[2]<.6){//beam to be imposed

              MACP0A1(prim,i,j,k,URAD0) = ERADINJ;
              MACP0A1(prim,i,j,k,URAD1) = uradx;
              MACP0A1(prim,i,j,k,URAD2) = 0.;
              MACP0A1(prim,i,j,k,URAD3) = 0.;
            }
            else{ //no beam
              MACP0A1(prim,i,j,k,URAD0) = ERADAMB;
              MACP0A1(prim,i,j,k,URAD1) = 0.;
              MACP0A1(prim,i,j,k,URAD2) = 0.;
              MACP0A1(prim,i,j,k,URAD3) = 0.;
            }


            // get all primitives in WHICHVEL/PRIMECOORDS value
            primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).


          }// end if not staggered field


        }// end loop over inner i's
      }
    }


  }// end parallel region

  return(0);
} 




// X1 lower for beam injection to create shadow
int bound_radshadowinflow(int dir,
						  int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
						  int *inboundloop,
						  int *outboundloop,
						  int *innormalloop,
						  int *outnormalloop,
						  int (*inoutlohi)[NUMUPDOWN][NDIM],
						  int riin, int riout, int rjin, int rjout, int rkin, int rkout,
						  int *dosetbc,
						  int enerregion,
						  int *localenerpos
						  )
  
{

  if(dir==X2DN){
	dualfprintf(fail_file,"Shouldn't be in bound_radshadowinflow() with dir=%d -- should be using ASYMM for RADDBLSHADOW\n");
	myexit(87243534);
  }


  FTYPE RHOAMB=1.e-4; // matches init.koral.c
  FTYPE TAMB=1.e7/TEMPBAR; // matches init.koral.c
  FTYPE BLOBW=0.22;
  FTYPE RHOBLOB=1.e3;
  //
  FTYPE NLEFT,angle;
  if(WHICHPROBLEM==RADSHADOW){
	NLEFT=0.99999;
	angle=0.0;
  }
  else if(WHICHPROBLEM==RADDBLSHADOW){
	//NLEFT=0.99999; // very hard on code -- only MINM with jon choice for CASES works.
	NLEFT=0.99;
	//  NLEFT=0.7;
	//  NLEFT=0.93;
	angle=0.4;
  }
  //
  FTYPE gammax=1.0/sqrt(1.0-NLEFT*NLEFT);
  FTYPE uradx=NLEFT*gammax/sqrt(1.0+angle*angle);
  FTYPE urady=-NLEFT*gammax*angle/sqrt(1.0+angle*angle);
  FTYPE TLEFT=TAMB*100.0;



#pragma omp parallel  // assume don't require EOS
  {
	// make rho,u consistent with on-domain values
	FTYPE rho;
	FTYPE uint;
	FTYPE Trad;

	int i,j,k,pl,pliter;
	FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
	int failreturn;
#endif
	int ri, rj, rk; // reference i,j,k
	FTYPE prescale[NPR];
	int jj,kk;
	struct of_geom geomdontuse[NPR];
	struct of_geom *ptrgeom[NPR];
	struct of_geom rgeomdontuse[NPR];
	struct of_geom *ptrrgeom[NPR];

	// assign memory
	PALLLOOP(pl){
	  ptrgeom[pl]=&(geomdontuse[pl]);
	  ptrrgeom[pl]=&(rgeomdontuse[pl]);
	}



	//////////////////////
	// dir==X1DN
	//////////////////////
	if(dir==X1DN && BCtype[X1DN]==RADSHADOWINFLOW && (totalsize[1]>1) && (mycpupos[1] == 0) ){


	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


		ri=riin;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
		LOOPBOUND1IN{
    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		  if(ispstag==0){
			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	  
			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);

			// set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates
			int whichvel;
			whichvel=VEL4;
			int whichcoord;
			whichcoord=CARTMINKMETRIC2;



			FTYPE xx,yy,zz,rsq;
			coord(i, j, k, CENT, X);
			bl_coord(X, V);
			xx=V[1];
			yy=V[2];
			zz=V[3];

			  
			if(yy>0.3 && WHICHPROBLEM==RADDBLSHADOW || WHICHPROBLEM==RADSHADOW ){

			  rsq=xx*xx+yy*yy+zz*zz;
			  rho=(RHOBLOB-RHOAMB)*exp(-sqrt(rsq)/(BLOBW*BLOBW))+RHOAMB;      
			  //			rho=RHOAMB;
			  Trad=TAMB*RHOAMB/rho;
			
			  uint=calc_PEQ_ufromTrho(Trad,rho);
			  FTYPE ERAD;
			  ERAD=calc_LTE_EfromT(TLEFT);
			  //			  FTYPE FxRAD;
			  //			  FxRAD=NLEFT*ERAD;
			  FTYPE ux=0.0; // orthonormal 4-velocity.  Matches init.koral.c

			  //			dualfprintf(fail_file,"i=%d j=%d k=%d TLEFT=%g ERAD=%g FxRAD=%g : ARAD_CODE=%g\n",i,j,k,TLEFT,ERAD,FxRAD,ARAD_CODE);

			  MACP0A1(prim,i,j,k,RHO) = rho;
			  MACP0A1(prim,i,j,k,UU) = uint;
			  MACP0A1(prim,i,j,k,U1) = ux/sqrt(ptrgeom[U1]->gcov[GIND(1,1)]); // assumed no spatial mixing


			  //E, F^i
			  MACP0A1(prim,i,j,k,URAD0) = ERAD;
			  //	      MACP0A1(prim,i,j,k,URAD1) = 0.;
			  //			  MACP0A1(prim,i,j,k,URAD1) = FxRAD;
			  MACP0A1(prim,i,j,k,URAD1) = uradx; //FxRAD;
			  //	      MACP0A1(prim,i,j,k,URAD2) = RADBEAMFLAT_FRATIO*MACP0A1(prim,i,j,k,URAD0);
			  MACP0A1(prim,i,j,k,URAD2) = urady;
			  MACP0A1(prim,i,j,k,URAD3) = 0.;

			  //			  dualfprintf(fail_file,"BC: i=%d j=%d rho=%g Trad=%g uint=%g ERAD=%g\n",i,j,rho,Trad,uint,ERAD);

              // KORALTODO GODMARK: ERAD is really fluid frame value, not radiation frame!  Need the below to account for that.  Currently just adjusted ERAD so injection is similar to expected.

              // get all primitives in WHICHVEL/PRIMECOORDS value
              primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
			} // over spatially relevant region

		  }// end if not staggered fields

		}// end loop over inner i's

	  }// end over loop
	}// end if correct boundary condition and core
	  
	  
	////////////////////////////////
	// dir==X2DN handled via AYMM BC
	////////////////////////////////
	  
	  
	//////////////////////
	// dir==X2UP
	//////////////////////
	if(dir==X2UP && BCtype[X2UP]==RADSHADOWINFLOWX2UP && (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ){
		

	  OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
	  ////////	LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


		ri=i;
		rj=rjout;
		rk=k;
		  
		  
		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
		LOOPBOUND2OUT{
    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		  if(ispstag==0){
			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	  
			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);

			// set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates
			int whichvel;
			whichvel=VEL4;
			int whichcoord;
			whichcoord=CARTMINKMETRIC2;



			FTYPE xx,yy,zz,rsq;
			coord(i, j, k, CENT, X);
			bl_coord(X, V);
			xx=V[1];
			yy=V[2];
			zz=V[3];


			if( WHICHPROBLEM==RADDBLSHADOW || WHICHPROBLEM==RADSHADOW ){

			  rsq=xx*xx+yy*yy+zz*zz;
			  rho=(RHOBLOB-RHOAMB)*exp(-sqrt(rsq)/(BLOBW*BLOBW))+RHOAMB;      
			  //			rho=RHOAMB;
			  Trad=TAMB*RHOAMB/rho;
			
			  uint=calc_PEQ_ufromTrho(Trad,rho);
			  FTYPE ERAD;
			  ERAD=calc_LTE_EfromT(TLEFT);
			  //			  FTYPE FxRAD;
			  //			  FxRAD=NLEFT*ERAD;
			  FTYPE ux=0.0; // orthonormal 4-velocity.  Matches init.koral.c

			  //			dualfprintf(fail_file,"i=%d j=%d k=%d TLEFT=%g ERAD=%g FxRAD=%g : ARAD_CODE=%g\n",i,j,k,TLEFT,ERAD,FxRAD,ARAD_CODE);

			  MACP0A1(prim,i,j,k,RHO) = rho;
			  MACP0A1(prim,i,j,k,UU) = uint;
			  MACP0A1(prim,i,j,k,U1) = ux/sqrt(ptrgeom[U1]->gcov[GIND(1,1)]); // assumed no spatial mixing


			  //E, F^i
			  MACP0A1(prim,i,j,k,URAD0) = ERAD;
			  //	      MACP0A1(prim,i,j,k,URAD1) = 0.;
			  //			  MACP0A1(prim,i,j,k,URAD1) = FxRAD;
			  MACP0A1(prim,i,j,k,URAD1) = uradx; //FxRAD;
			  //	      MACP0A1(prim,i,j,k,URAD2) = RADBEAMFLAT_FRATIO*MACP0A1(prim,i,j,k,URAD0);
			  MACP0A1(prim,i,j,k,URAD2) = urady;
			  MACP0A1(prim,i,j,k,URAD3) = 0.;

			  //			  dualfprintf(fail_file,"BC: i=%d j=%d rho=%g Trad=%g uint=%g ERAD=%g\n",i,j,rho,Trad,uint,ERAD);


              // get all primitives in WHICHVEL/PRIMECOORDS value
              primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
			}// if spatially relevant region

		  }// end if not staggered fields
			
		}// end loop over inner j's
	  }// end over loop
	}// end if correct boundary condition and core
	  
	  
  }// end parallel region
  


 



  return(0);
} 





// X3 lower for radiation beam injection
int bound_radbeam2dbeaminflow(int dir,
							int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
							int *inboundloop,
							int *outboundloop,
							int *innormalloop,
							int *outnormalloop,
							int (*inoutlohi)[NUMUPDOWN][NDIM],
							int riin, int riout, int rjin, int rjout, int rkin, int rkout,
							int *dosetbc,
							int enerregion,
							int *localenerpos
							)

{


#pragma omp parallel  // assume don't require EOS
  {

    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if(dir==X3DN && BCtype[X3DN]==RADBEAM2DBEAMINFLOW && totalsize[3]>1 && mycpupos[3] == 0 ){

      extern int BEAMNO,FLATBACKGROUND;
      FTYPE RHOAMB=1.e0/RHOBAR;
	  FTYPE TAMB=1e7/TEMPBAR;
	  FTYPE PAR_D=1./RHOBAR;
	  FTYPE PAR_E=1e-4/RHOBAR;

      // BEAM PROPERTIES
	  int IFBEAM=1; // whether to have a beam
	  FTYPE TLEFT=1e9/TEMPBAR;
      FTYPE NLEFT=0.99;
      //	  FTYPE NLEFT=0.999;
      //	  FTYPE NLEFT=0.999999; // paper says this, while koral code says 0.999

	  FTYPE BEAML,BEAMR;
	  if (BEAMNO==1){
		BEAML=2.9;
		BEAMR=3.1;
	  }
	  else if (BEAMNO==2){
		BEAML=5.8;
		BEAMR=6.2;
	  }
	  else if (BEAMNO==3){
		BEAML=15.5;
		BEAMR=16.5;
	  }
	  else if (BEAMNO==4){
		BEAML=37;
		BEAMR=43;
	  }

	  

	  OPENMPBCLOOPVARSDEFINELOOPX3DIR; OPENMPBCLOOPSETUPLOOPX3DIR;
	  ////////	LOOPX3dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX3DIR(i,j);

		ri=i;
		rj=j;
		rk=rkin;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
		LOOPBOUND3IN{

    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
		  
          // NOTEMARK: only really makes sense near the hole if in KSCOORDS
		  if(pl==U3) if(MACP0A1(prim,i,j,k,U3)>0.0) MACP0A1(prim,i,j,k,U3)=0.0; // limit so no arbitrary fluid inflow
		  if(pl==URAD3) if(MACP0A1(prim,i,j,k,URAD3)>0.0) MACP0A1(prim,i,j,k,URAD3)=0.0; // limit so no arbitrary radiative inflow


		  // only overwrite copy if not inside i<0 since want to keep consistent with outflow BCs used for other \phi at those i
		  if(1      &&     startpos[1]+i>=0 && ispstag==0){ // only do something special with non-field primitives

			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);

			// set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates
			int whichvel;
			whichvel=VEL4;
			int whichcoord;
			whichcoord=MCOORD;

			// get metric grid geometry for these ICs
			int getprim=0;
			struct of_geom geomrealdontuse;
			struct of_geom *ptrgeomreal=&geomrealdontuse;
			gset(getprim,whichcoord,i,j,k,ptrgeomreal);


			FTYPE ERADAMB;
			FTYPE rho,uint,Vr;
			FTYPE uradx,urady,uradz;
            FTYPE uradcon[NDIM],othersrad[NUMOTHERSTATERESULTS];


            // get coordinate basis in VEL4 format
            ucon_calc(&MACP0A1(prim,i,j,k,URAD1-U1),ptrgeom[URAD1],uradcon,othersrad);
            // get coordinate basis in MCOORD basis
            uradx=uradcon[1]*sqrt(fabs(ptrgeom[URAD1]->gcov[GIND(1,1)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(1,1)]));
            urady=uradcon[2]*sqrt(fabs(ptrgeom[URAD2]->gcov[GIND(2,2)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(2,2)]));
            uradz=uradcon[3]*sqrt(fabs(ptrgeom[URAD3]->gcov[GIND(3,3)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(3,3)]));
            if(uradz>0.0) uradz=0.0; // limit so no arbitrary radiative inflow


			if(FLATBACKGROUND){
			  Vr=0.0;
			  rho=RHOAMB;
			  uint=calc_PEQ_ufromTrho(TAMB,rho);
			  ERADAMB=calc_LTE_EfromT(TAMB);

			  // override so like outflow conditions to avoid shear at boundary
			  ERADAMB=MACP0A1(prim,i,j,k,URAD0);

			}
			else{
			  //zaczynam jednak od profilu analitycznego:   
			  FTYPE r=V[1];
			  FTYPE mD=PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
			  FTYPE mE=PAR_E/(pow(r*r*sqrt(2./r),gamideal)*pow(1.-2./r,(gamideal+1.)/4.));
			  Vr=sqrt(2./r)*(1.-2./r);


			  FTYPE W=1./sqrt(1.-Vr*Vr*ptrgeomreal->gcov[GIND(1,1)]); // assumes RHO location is good for all these quantities
			  rho=PAR_D/(r*r*sqrt(2./r));
			  FTYPE T=TAMB;
			  //			FTYPE ERAD=calc_LTE_EfromT(T);
			  uint=mE/W;
			  ERADAMB=calc_LTE_Efromurho(uint,rho);
			  
			  // override so like outflow conditions to avoid shear at boundary
			  ERADAMB=MACP0A1(prim,i,j,k,URAD0);
			}

            PBOUNDLOOP(pliter,pl){

              //primitives in whichvel,whichcoord
              if(V[1]>BEAML && V[1]<BEAMR && IFBEAM){//beam to be imposed
                
                // beam injection
                // override URAD0
                FTYPE ERADINJ;
                ERADINJ=calc_LTE_EfromT(TLEFT);
                // override uradz
                uradz=1.0/sqrt(1.0 - NLEFT*NLEFT);
                uradx=urady=0.0;
                
                //                dualfprintf(fail_file,"i=%d k=%d ERADAMB=%g ERADINJ=%g uradx=%g urady=%g uradz=%g\n",i,k,ERADAMB,ERADINJ,uradx,urady,uradz);

                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD0) = ERADINJ;
                if(pl==URAD1) MACP0A1(prim,i,j,k,URAD1) = uradx;
                if(pl==URAD2) MACP0A1(prim,i,j,k,URAD2) = urady;
                if(pl==URAD3) MACP0A1(prim,i,j,k,URAD3) = uradz;
              }
              else{ //no beam
                
                //              dualfprintf(fail_file,"i=%d k=%d ERADAMB=%g\n",i,k,ERADAMB);
                
                //			  MACP0A1(prim,i,j,k,URAD0) = ERADAMB;
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD0) = ERADAMB; // so matches outer radial boundary when no beam
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD1) = uradx;
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD2) = urady;
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD3) = uradz;
              }
            } // over allowed pl's to bound


            //            PLOOP(pliter,pl) dualfprintf(fail_file,"BEFOREBC: pl=%d prim=%g\n",pl,MACP0A1(prim,i,j,k,pl));

            //            if(i==10 && k==0){
            //              PLOOP(pliter,pl) dualfprintf(fail_file,"BEFOREBC: pl=%d prim=%g\n",pl,MACP0A1(prim,i,j,k,pl));
            //              if(MACP0A1(prim,i,j,k,UU)>0.1) dualfprintf(fail_file,"BEFORE BC DEATHUU: ijk=%d %d %d u=%g\n",i,j,k,MACP0A1(prim,i,j,k,UU));
            //              if(fabs(MACP0A1(prim,i,j,k,U1))>1.0) dualfprintf(fail_file,"BEFORE BC DEATHU1: ijk=%d %d %d u=%g\n",i,j,k,MACP0A1(prim,i,j,k,U1));
            //            }

            // KORALTODO: ERADINJ is in fluid frame, need to convert, but probably ok.

			// get all primitives in WHICHVEL/PRIMECOORDS value
			primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).

            //            MACP0A1(prim,i,j,k,URAD1)=0.0;

            //            PLOOP(pliter,pl) dualfprintf(fail_file,"AFTERBC: pl=%d prim=%g\n",pl,MACP0A1(prim,i,j,k,pl));


            //            if(i==10 && k==0){
            //              PLOOP(pliter,pl) dualfprintf(fail_file,"AFTERBC: pl=%d prim=%g\n",pl,MACP0A1(prim,i,j,k,pl));
            //              if(MACP0A1(prim,i,j,k,UU)>0.1) dualfprintf(fail_file,"AFTER BC DEATHUU: ijk=%d %d %d u=%g\n",i,j,k,MACP0A1(prim,i,j,k,UU));
            //              if(fabs(MACP0A1(prim,i,j,k,U1))>1.0) dualfprintf(fail_file,"AFTER BC DEATHU1: ijk=%d %d %d u=%g\n",i,j,k,MACP0A1(prim,i,j,k,U1));
            //            }

		  }// end if not staggered field


		}// end loop over inner i's
	  } // over block
	}// end if correct BC and should be doind BC for this core
  }// end parallel region

  return(0);
} 





// X1 upper for inflow
int bound_radbeam2dflowinflow(int dir,
							   int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
							   int *inboundloop,
							   int *outboundloop,
							   int *innormalloop,
							   int *outnormalloop,
							   int (*inoutlohi)[NUMUPDOWN][NDIM],
							   int riin, int riout, int rjin, int rjout, int rkin, int rkout,
							   int *dosetbc,
							   int enerregion,
							   int *localenerpos
							   )

{


#pragma omp parallel  // assume don't require EOS
  {

    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if(dir==X1UP && BCtype[X1UP]==RADBEAM2DFLOWINFLOW && totalsize[1]>1 && mycpupos[1] == ncpux1-1 ){

      extern int FLATBACKGROUND;
	  FTYPE RHOAMB=1.e0/RHOBAR;
	  FTYPE TAMB=1e7/TEMPBAR;
	  FTYPE PAR_D=1./RHOBAR;
	  FTYPE PAR_E=1e-4/RHOBAR;


	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



		ri=riout;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
		LOOPBOUND1OUT{

    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
		  

		  if(ispstag==0){ // only do something special with non-field primitives

			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);

			// set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates
			int whichvel;
			whichvel=VEL4;
			int whichcoord;
			whichcoord=MCOORD;


			FTYPE ERADAMB;
			FTYPE rho,uint,Vr;
			if(FLATBACKGROUND){
			  Vr=0.0;
			  rho=RHOAMB;
			  uint=calc_PEQ_ufromTrho(TAMB,rho);
			  ERADAMB=calc_LTE_EfromT(TAMB);
			}
			else{
			  //zaczynam jednak od profilu analitycznego:   
			  FTYPE r=V[1];
			  FTYPE mD=PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
			  FTYPE mE=PAR_E/(pow(r*r*sqrt(2./r),gamideal)*pow(1.-2./r,(gamideal+1.)/4.));
			  Vr=sqrt(2./r)*(1.-2./r);

			  // get metric grid geometry for these ICs
			  int getprim=0;
			  struct of_geom geomrealdontuse;
			  struct of_geom *ptrgeomreal=&geomrealdontuse;
			  gset(getprim,whichcoord,i,j,k,ptrgeomreal);

			  FTYPE W=1./sqrt(1.-Vr*Vr*ptrgeomreal->gcov[GIND(1,1)]); // assumes RHO location is good for all these quantities
			  rho=PAR_D/(r*r*sqrt(2./r));
			  FTYPE T=TAMB;
			  //			FTYPE ERAD=calc_LTE_EfromT(T);
			  uint=mE/W;
			  ERADAMB=calc_LTE_Efromurho(uint,rho);
			}
			FTYPE uradx,urady,uradz;
			uradx=urady=uradz=0.0;

			// set quantities at outer radial edge
			MACP0A1(prim,i,j,k,RHO) = rho;
			MACP0A1(prim,i,j,k,UU)  = uint;
			MACP0A1(prim,i,j,k,U1)  = -Vr;
			MACP0A1(prim,i,j,k,U2)  = 0.;
			MACP0A1(prim,i,j,k,U3)  = 0.;
			MACP0A1(prim,i,j,k,PRAD0) = ERADAMB;
			MACP0A1(prim,i,j,k,PRAD1) = uradx;
			MACP0A1(prim,i,j,k,PRAD2) = urady;
			MACP0A1(prim,i,j,k,PRAD3) = uradz;

            // KORALTODO: ERADAMB is in fluid frame, need to convert, but probably ok.

			
			//			dualfprintf(fail_file,"IC: ijk=%d %d %d : rho=%g u=%g Vr=%g erad=%g\n",i,j,k,rho,uint,-Vr,ERAD);

			// get all primitives in WHICHVEL/PRIMECOORDS value
			if (bl2met2metp2v(whichvel, whichcoord,MAC(prim,i,j,k), i,j,k) >= 1){
			  FAILSTATEMENT("bounds.koral.c:bound_radbeam2dflowinflow()", "bl2ks2ksp2v()", 1);
			}
			
			//			dualfprintf(fail_file,"POSTIC: ijk=%d %d %d : rho=%g u=%g Vr=%g erad=%g\n",i,j,k,rho,uint,-Vr,ERAD);


		  }// end if not staggered field


		}// end loop over inner i's
	  }// over block
	}// if correct BC and core
  }// end parallel region

  return(0);
} 





// X1 lower and upper for RADATM
int bound_radatmbeaminflow(int dir,
						  int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
						  int *inboundloop,
						  int *outboundloop,
						  int *innormalloop,
						  int *outnormalloop,
						  int (*inoutlohi)[NUMUPDOWN][NDIM],
						  int riin, int riout, int rjin, int rjout, int rkin, int rkout,
						  int *dosetbc,
						  int enerregion,
						  int *localenerpos
						  )
  
{





#pragma omp parallel  // assume don't require EOS
  {

	int i,j,k,pl,pliter;
	FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
	int failreturn;
#endif
	int ri, rj, rk; // reference i,j,k
	FTYPE prescale[NPR];
	int jj,kk;
	struct of_geom geomdontuse[NPR];
	struct of_geom *ptrgeom[NPR];
	struct of_geom rgeomdontuse[NPR];
	struct of_geom *ptrrgeom[NPR];

	// assign memory
	PALLLOOP(pl){
	  ptrgeom[pl]=&(geomdontuse[pl]);
	  ptrrgeom[pl]=&(rgeomdontuse[pl]);
	}


    extern FTYPE RADATM_MDOTEDD;
    extern FTYPE RADATM_LUMEDD;
    extern int RADATM_THINRADATM;
    extern FTYPE RADATM_FERATIO;
    extern FTYPE RADATM_FRATIO;
    extern FTYPE RADATM_RHOAMB;
    extern FTYPE RADATM_TAMB;


    FTYPE MINX=Rin_array[1];
    FTYPE kappaesperrho=calc_kappaes_user(1,0, 0,0,0);
    FTYPE FLUXLEFT=RADATM_FRATIO/kappaesperrho/pow(MINX,2.0);

    //at boundary
    FTYPE f = (FTYPE)kappaesperrho*FLUXLEFT*MINX*MINX;

    FTYPE p0=RADATM_RHOAMB*RADATM_TAMB;
    FTYPE KKK=p0/pow(RADATM_RHOAMB,gamideal);
    FTYPE C3=gamideal*KKK/(gamideal-1.)*pow(RADATM_RHOAMB,gamideal-1.)-(1.-f)*(1./MINX+0.*1./MINX/MINX+0.*4./3./MINX/MINX/MINX);

    //    dualfprintf(fail_file,"IT: %g %g %g : %g : %g %g %g\n",MINX,kappaesperrho,FLUXLEFT,f,p0,KKK,C3);
	  
	if(dir==X1DN && BCtype[X1DN]==RADATMBEAMINFLOW && (totalsize[1]>1) && (mycpupos[1] == 0) ){


	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


		ri=riin;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

        FTYPE *pr;
		LOOPBOUND1IN{

          pr = &MACP0A1(prim,i,j,k,0);
    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		  if(ispstag==0){
			// local PRIMECOORDS geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	  
			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);

            FTYPE xx,yy,zz,rsq;
            coord(i, j, k, CENT, X);
            bl_coord(X, V);
            xx=V[1];
            yy=V[2];
            zz=V[3];


            FTYPE rho=pow((gamideal-1.0)/gamideal/KKK*(C3+(1.-f)*(1./xx+0.*1./xx/xx+0.*4./3./xx/xx/xx)),1./(gamideal-1.0));

            FTYPE pre=KKK*pow(rho,gamideal);

            FTYPE uint=pre/(gamideal-1.0);

            FTYPE Fz=0;
            FTYPE Fy=0.;
            FTYPE Fx=FLUXLEFT*(MINX/xx)*(MINX/xx);

            FTYPE ERAD;
            if(RADATM_THINRADATM){
              ERAD=Fx/RADATM_FERATIO;
            }
            else{
              ERAD=calc_LTE_EfromT(calc_PEQ_Tfromurho(uint,rho));
            }

            //            dualfprintf(fail_file,"BC: i=%d j=%d rho=%g Trad=%g uint=%g ERAD=%g\n",i,j,rho,uint,ERAD);


            pr[RHO] = rho;
            pr[UU] = uint;
            pr[U1] = 0.0; // static
            pr[U2] = 0.0;
            pr[U3] = 0.0;

            //E, F^i in orthonormal fluid frame
            FTYPE pradffortho[NPR];
            pradffortho[PRAD0] = ERAD;
            pradffortho[PRAD1] = Fx;
            pradffortho[PRAD2] = Fy;
            pradffortho[PRAD3] = Fz;


            int whichvel=VEL4; // in which vel U1-U3 set
            int whichcoordfluid=MCOORD; // in which coordinates U1-U3 set
            int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
            whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

            //            PLOOP(pliter,pl) dualfprintf(fail_file,"ijk=%d %d %d pl=%d pr=%g\n",i,j,k,pl,pr[pl]);

		  }// end if not staggered fields

		}// end loop over inner i's

	  }// end over loop
	}// end if correct boundary condition and core




	if(dir==X1UP && BCtype[X1UP]==RADATMBEAMINFLOW && (totalsize[1]>1) && (mycpupos[1] == ncpux1-1) ){


	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


		ri=riout;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        FTYPE *pr;
		LOOPBOUND1OUT{
          
          pr = &MACP0A1(prim,i,j,k,0);

    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		  if(ispstag==0){
			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	  
			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);


            FTYPE xx,yy,zz,rsq;
            coord(i, j, k, CENT, X);
            bl_coord(X, V);
            xx=V[1];
            yy=V[2];
            zz=V[3];



            FTYPE rho=pow((gamideal-1.0)/gamideal/KKK*(C3+(1.-f)*(1./xx+0.*1./xx/xx+0.*4./3./xx/xx/xx)),1./(gamideal-1.0));

            FTYPE pre=KKK*pow(rho,gamideal);

            FTYPE uint=pre/(gamideal-1.0);

            FTYPE Fz=0;
            FTYPE Fy=0.;
            FTYPE Fx=FLUXLEFT*(MINX/xx)*(MINX/xx);

            FTYPE ERAD;
            if(RADATM_THINRADATM){
              ERAD=Fx/RADATM_FERATIO;
            }
            else{
              ERAD=calc_LTE_EfromT(calc_PEQ_Tfromurho(uint,rho));
            }


            pr[RHO] = rho;
            pr[UU] = uint;
            pr[U1] = 0.0; // static
            pr[U2] = 0.0;
            pr[U3] = 0.0;

            //E, F^i in orthonormal fluid frame
            FTYPE pradffortho[NPR];
            pradffortho[PRAD0] = ERAD;
            pradffortho[PRAD1] = Fx;
            pradffortho[PRAD2] = Fy;
            pradffortho[PRAD3] = Fz;

            
            int whichvel=VEL4;
            int whichcoordfluid=MCOORD;
            int whichcoordrad=whichcoordfluid;
            whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

		  }// end if not staggered fields

		}// end loop over outer i's

	  }// end over loop
	}// end if correct boundary condition and core
	  
	  
	  
  }// end parallel region
  


 



  return(0);
} 





// X1 lower and upper X2
int bound_radwallinflow(int dir,
                        int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                        int *inboundloop,
                        int *outboundloop,
                        int *innormalloop,
                        int *outnormalloop,
                        int (*inoutlohi)[NUMUPDOWN][NDIM],
                        int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                        int *dosetbc,
                        int enerregion,
                        int *localenerpos
                        )

{


#pragma omp parallel  // assume don't require EOS
  {

    FTYPE LTEFACTOR=1.;
    FTYPE URFX=100.;
    FTYPE URFY=-30.;


    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if(dir==X1DN && BCtype[X1DN]==RADWALLINFLOW && (totalsize[1]>1) && (mycpupos[1] == 0) ){



      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



        ri=riin;
        rj=j;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
        LOOPBOUND1IN{

    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);


          if(ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);

            int whichvel;
            whichvel=VEL4;
            int whichcoord;
            whichcoord=CARTMINKMETRIC2;

            PBOUNDLOOP(pliter,pl){
              if(pl==RHO) MACP0A1(prim,i,j,k,RHO) = 1.0;
              if(pl==UU) MACP0A1(prim,i,j,k,UU) = 1.0;
              if(pl==U1) MACP0A1(prim,i,j,k,U1) = 0.0;


              if(V[2]>0.3){
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD0) = pow(100,4.0);
                if(pl==URAD1) MACP0A1(prim,i,j,k,URAD1) = URFX;
                if(pl==URAD2) MACP0A1(prim,i,j,k,URAD2) = URFY;
                if(pl==URAD3) MACP0A1(prim,i,j,k,URAD3) = 0.;
              }
              else{
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD0) = 1.0;
                if(pl==URAD1) MACP0A1(prim,i,j,k,URAD1) = 0.0;
                if(pl==URAD2) MACP0A1(prim,i,j,k,URAD2) = 0.;
                if(pl==URAD3) MACP0A1(prim,i,j,k,URAD3) = 0.;
              }
            }// end over pl's allowed for bounding

            // get all primitives in WHICHVEL/PRIMECOORDS value
            if(1) primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).

          }// end if not staggered field


        }// end loop over inner i's
      }
    }




    if(dir==X2UP && BCtype[X2UP]==RADWALLINFLOW && (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ){



      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      ////////	LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);



        ri=i;
        rj=rjout;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
        LOOPBOUND2OUT{

    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);


          if(ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);

            int whichvel;
            whichvel=VEL4;
            int whichcoord;
            whichcoord=CARTMINKMETRIC2;

            PBOUNDLOOP(pliter,pl){
              if(pl==RHO) MACP0A1(prim,i,j,k,RHO) = 1.0;
              if(pl==UU) MACP0A1(prim,i,j,k,UU) = 1.0;
              if(pl==U1) MACP0A1(prim,i,j,k,U1) = 0.0;

              if(1){
                if(pl==URAD0) MACP0A1(prim,i,j,k,URAD0) = pow(100,4.0);
                if(pl==URAD1) MACP0A1(prim,i,j,k,URAD1) = URFX;
                if(pl==URAD2) MACP0A1(prim,i,j,k,URAD2) = URFY;
                if(pl==URAD3) MACP0A1(prim,i,j,k,URAD3) = 0.;
              }

            }// over allowed pl's


            // get all primitives in WHICHVEL/PRIMECOORDS value
            if(1) primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).


          }// end if not staggered field


        }// end loop over outer j's
      }
    }








  }// end parallel region

  return(0);
} 









// X1 upper for inflow
int bound_radbondiinflow(int dir,
							   int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
							   int *inboundloop,
							   int *outboundloop,
							   int *innormalloop,
							   int *outnormalloop,
							   int (*inoutlohi)[NUMUPDOWN][NDIM],
							   int riin, int riout, int rjin, int rjout, int rkin, int rkout,
							   int *dosetbc,
							   int enerregion,
							   int *localenerpos
							   )

{


#pragma omp parallel  // assume don't require EOS
  {

    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if(dir==X1UP && BCtype[X1UP]==RADBONDIINFLOW && totalsize[1]>1 && mycpupos[1] == ncpux1-1 ){


      extern FTYPE RADBONDI_TESTNO;
      extern FTYPE RADBONDI_PRADGAS;
      extern FTYPE RADBONDI_TGAS0;
      extern FTYPE RADBONDI_MDOTPEREDD;
      extern FTYPE RADBONDI_MDOTEDD;
      extern FTYPE RADBONDI_RHOAMB;
      extern FTYPE RADBONDI_TAMB;
      extern FTYPE RADBONDI_MUGAS;
      extern FTYPE RADBONDI_MINX;
      extern FTYPE RADBONDI_MAXX;


	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



		ri=riout;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
        FTYPE *pr;
		LOOPBOUND1OUT{

          pr = &MACP0A1(prim,i,j,k,0);
          
    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
		  

		  if(ispstag==0){ // only do something special with non-field primitives

			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);

            // Identical to init, so could use analytic, but don't for now.

            FTYPE xx,yy,zz,rsq;
            xx=V[1];
            yy=V[2];
            zz=V[3];

            int whichvel=VEL3; // in which vel U1-U3 set
            int whichcoordfluid=MCOORD; // in which coordinates U1-U3 set
            // get metric grid geometry for these ICs
            int getprim=0;
            struct of_geom geomrealdontuse;
            struct of_geom *ptrgeomreal=&geomrealdontuse;
            gset(getprim,whichcoordfluid,i,j,k,ptrgeomreal);

    
            FTYPE rho,ERAD,uint;
            FTYPE rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;

            FTYPE Fx,Fy,Fz;
            Fx=Fy=Fz=0;

            //at outern boundary
            r=RADBONDI_MAXX;
            ur=-sqrt(2./r);
            rho0=-RADBONDI_MDOTPEREDD*RADBONDI_MDOTEDD/(4.*Pi*r*r*ur);
            Tgas0=RADBONDI_TGAS0;
            
            //at given cell
            r=xx;
            ur=-sqrt(2./r);    
            ut=sqrt((-1.-ur*ur*ptrgeomreal->gcov[GIND(1,1)])/ptrgeomreal->gcov[GIND(0,0)]);
            vx=ur/ut;  
            rho=-RADBONDI_MDOTPEREDD*RADBONDI_MDOTEDD/(4.*Pi*r*r*ur);
            Tgas=Tgas0*pow(rho/rho0,gam-1.);      

            uint=calc_PEQ_ufromTrho(Tgas,rho);

            pgas=rho*Tgas;
            prad=RADBONDI_PRADGAS*pgas;
            ERAD=prad*3.;
    

            pr[RHO] = rho ;
            pr[UU]  = uint;
            pr[U1]  = vx;
            pr[U2]  = 0 ;    
            pr[U3]  = 0 ;

            //E, F^i in orthonormal fluid frame
            FTYPE pradffortho[NPR];
            pradffortho[PRAD0] = ERAD;
            pradffortho[PRAD1] = Fx;
            pradffortho[PRAD2] = Fy;
            pradffortho[PRAD3] = Fz;

            dualfprintf(fail_file,"rho=%g uint=%g vx=%g ERAD=%g\n",rho,uint,vx,ERAD);

            int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
            whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

            PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pr=%g\n",pl,pr[pl]);

		  }// end if not staggered field


		}// end loop over inner i's
	  }// over block
	}// if correct BC and core
  }// end parallel region

  return(0);
} 






// on-grid bounding
int bound_raddot(
                 int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                 int *inboundloop,
                 int *outboundloop,
                 int *innormalloop,
                 int *outnormalloop,
                 int (*inoutlohi)[NUMUPDOWN][NDIM],
                 int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                 int *dosetbc,
                 int enerregion,
                 int *localenerpos
                 )

{


#pragma omp parallel  // assume don't require EOS
  {

    extern FTYPE RADDOT_XDOT;
    extern FTYPE RADDOT_YDOT;
    extern FTYPE RADDOT_ZDOT;
    extern int RADDOT_IDOT;
    extern int RADDOT_JDOT;
    extern int RADDOT_KDOT;
    extern FTYPE RADDOT_FYDOT;
    extern FTYPE RADDOT_LTEFACTOR;
    extern FTYPE RADDOT_URFX;
    extern FTYPE RADDOT_F1;
    extern FTYPE RADDOT_F2;


    int i,j,k,pl,pliter;
    FTYPE X[NDIM],V[NDIM]; 
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }

    OPENMP3DLOOPVARSDEFINE; OPENMP3DLOOPSETUPFULL;
  
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize)) nowait // can nowait since each fluxvec[dir] is set separately
    OPENMP3DLOOPBLOCK{
      OPENMP3DLOOPBLOCK2IJK(i,j,k);
      ////COMPFULLLOOP

      // KORALTODO: Koral has different IC and BC for the dot
      if(startpos[1]+i==RADDOT_IDOT && startpos[2]+j==RADDOT_JDOT && startpos[3]+k==RADDOT_KDOT){
        // NOTEMARK: It's normal for inversion to fail right on the dot, but doesn't affect evolution since the dot's primitives are overwritten and never will the updated conserved quantities for the dot be used.
        
        if(ispstag==0){ // only do something special with non-field primitives
        
          FTYPE *pr=&MACP0A1(prim,i,j,k,0);
        
          // get current location
          bl_coord_ijk_2(i,j,k,CENT,X, V);
        
          // local geom
          PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
        
          pr[RHO] = 1.0 ;
          pr[UU]  = 1.0;
          pr[U1]  = 0.0;
          pr[U2]  = 0 ;    
          pr[U3]  = 0 ;
        
          //E, F^i in orthonormal fluid frame
          FTYPE pradffortho[NPR];
          pradffortho[PRAD0] = RADDOT_LTEFACTOR*calc_LTE_Efromurho(pr[RHO],pr[UU]);
          pradffortho[PRAD1] = 0;
          pradffortho[PRAD2] = 0;
          pradffortho[PRAD3] = 0;
        
          //          dualfprintf(fail_file,"GOT BC DOT: nstep=%ld steppart=%d\n",nstep,steppart);
          if(N1==1) pradffortho[PRAD0] *= RADDOT_F1;
          else{
            pradffortho[PRAD0]*=RADDOT_F2;
            pradffortho[PRAD2]=RADDOT_FYDOT*pradffortho[PRAD0];
          }

          int whichvel=VEL4; // in which vel U1-U3 set
          int whichcoordfluid=MCOORD; // in which coordinates U1-U3 set
          int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
          whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);
        
        }// end if DOT
      }// end if not staggered field
    }// end loop over zones
 


  }// end parallel region

  return(0);
} 






// RADNT
int bound_radnt(int dir,
                int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
                int *inboundloop,
                int *outboundloop,
                int *innormalloop,
                int *outnormalloop,
                int (*inoutlohi)[NUMUPDOWN][NDIM],
                int riin, int riout, int rjin, int rjout, int rkin, int rkout,
                int *dosetbc,
                int enerregion,
                int *localenerpos
                )

{


#pragma omp parallel  // assume don't require EOS
  {

    extern FTYPE RADNT_MINX;
    extern FTYPE RADNT_MAXX;
    extern FTYPE RADNT_KKK;
    extern FTYPE RADNT_ELL;
    extern FTYPE RADNT_UTPOT;
    extern FTYPE RADNT_RHOATMMIN;
    extern FTYPE RADNT_UINTATMMIN;
    extern FTYPE RADNT_ERADATMMIN;
    extern FTYPE RADNT_NODONUT;
    extern FTYPE RADNT_INFLOWING;
    extern FTYPE RADNT_TGASATMMIN;
    extern FTYPE RADNT_TRADATMMIN;
    extern FTYPE RADNT_ROUT;
    extern FTYPE RADNT_OMSCALE;
    extern FTYPE RADNT_FULLPHI;


    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];

    // assign memory
    PALLLOOP(pl){
      ptrgeom[pl]=&(geomdontuse[pl]);
      ptrrgeom[pl]=&(rgeomdontuse[pl]);
    }


  
    if(dir==X1UP && BCtype[X1UP]==RADNTBC && totalsize[1]>1 && mycpupos[1] == ncpux1-1 ){

	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



		ri=riout;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

	  
        FTYPE *pr;
		LOOPBOUND1OUT{
          pr = &MACP0A1(prim,i,j,k,0);
    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
		  

		  if(ispstag==0){ // only do something special with non-field primitives

			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            FTYPE r,th,ph;
            coord(i, j, k, CENT, X);
            bl_coord(X, V);
            r=V[1];
            th=V[2];
            ph=V[3];

            // Identical to IC, except more involved check for inflow vs. outflow

            int whichcoord;
            if(WHICHPROBLEM==RADFLATDISK) whichcoord=MCOORD; // whatever else
            else whichcoord=BLCOORDS; // want to setup things in BLCOORDS
            int whichvel=VEL4;

            // get metric grid geometry for these ICs
            int getprim=0;
            struct of_geom geomrealdontuse;
            struct of_geom *ptrgeomreal=&geomrealdontuse;
            gset(getprim,whichcoord,i,j,k,ptrgeomreal);

            FTYPE uconlab[NDIM];
            FTYPE others[NUMOTHERSTATERESULTS];

            // see if fluid flow wants to go in or out
            // get PRIMECOORDS ucon
            ucon_calc(pr,ptrgeom[U1],uconlab,others);
            // get MCOORD
            metptomet(i,j,k,uconlab);
            // get whichcoord
            coordtrans(MCOORD,whichcoord,i,j,k,CENT,uconlab);
            if(uconlab[RR]<=0.0){ // check in whichvel whichcoord
              pr[RHO]=RADNT_RHOATMMIN*pow(r/RADNT_ROUT,-1.5);
              pr[UU]=RADNT_UINTATMMIN*pow(r/RADNT_ROUT,-2.5);
              set_zamo_velocity(whichvel,ptrgeomreal,pr); // only sets U1-U3 to zamo
            }
            else{
              uconlab[RR]=0.0;
              // overwrite pr[U1-U3] with this non-radially moving flow in whichvel whichcoord version
              ucon2pr(whichvel,uconlab,ptrgeomreal,pr);
            }


            // assume radiation in fluid frame with zero flux
            FTYPE pradffortho[NPR];

            pradffortho[PRAD0] = RADNT_ERADATMMIN; // fluid frame ortho!
            pradffortho[PRAD1] = 0;
            pradffortho[PRAD2] = 0;
            pradffortho[PRAD3] = 0;
            
            int whichcoordfluid;
            whichcoordfluid=whichcoord;
            int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
            whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

		  }// end if not staggered field


		}// end loop over inner i's
	  }// over block
	}// if correct BC and core







    if(dir==X2UP && BCtype[X2UP]==RADNTBC && (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ){

      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      ////////	LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);

        ri=i;
        rj=rjout;
        rk=k;

        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

        FTYPE *pr;
        LOOPBOUND2OUT{
          pr = &MACP0A1(prim,i,j,k,0);

          // ASYMM already done before got here, so only change what's necessary to change

          if(ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);
            FTYPE r;
            r=V[1];

            FTYPE rin;
            if(WHICHPROBLEM==RADFLATDISK) rin=15.;
            else rin=6.;

            //hot boundary KORALTODO: JCM confused why not much output in URAD0 out of disk
            if(r>rin){

              //E, F^i in orthonormal fluid frame
              FTYPE pradffortho[NPR];
              if(WHICHPROBLEM==RADFLATDISK) pradffortho[PRAD0] = calc_LTE_EfromT(1.e11/TEMPBAR);
              else pradffortho[PRAD0] = calc_LTE_EfromT(1.e11/TEMPBAR)*(1.-sqrt(rin/r))/pow(r,3.);
              // KORALTODO: in reality, can only constrain magnitude, not direction, of outflow away from plane.
              pradffortho[PRAD1] = 0;
              pradffortho[PRAD2] = -0.5*pradffortho[PRAD0];
              pradffortho[PRAD3] = 0;


              //pr[RHO] and pr[UU] remain same as from ASYMM condition as well as any field
              //Keplerian gas with no inflow or outflow
              // KORALTODO: in reality, can only constrain magnitude, not direction, of outflow away from plane.  But since magnitude is chosen to be zero, then no issue here.
              pr[U1]=pr[U2]=0.0; // have to be careful with this for VEL3 (must have rin>>rergo).
              if(WHICHPROBLEM==RADFLATDISK) pr[U3]=0.0;
              else pr[U3]=1./(a + pow(r,1.5));
	
              int whichvel;
              whichvel=VEL3; // VEL3 so can set Keplerian rotation rate

              int whichcoordfluid;
              if(WHICHPROBLEM==RADFLATDISK) whichcoordfluid=MCOORD; // whatever else
              else whichcoordfluid=BLCOORDS; // want to setup things in BLCOORDS

              int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
              whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

            } // end if actually doing something to boundary cells in "hot" boundary

          }// end if not staggered field


        }// end loop over outer j's
      }
    }







  }// end parallel region

  return(0);
} 







// X1 inner CYLAXIS
int bound_x1dn_cylaxis(
				   int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
				   int *inboundloop,
				   int *outboundloop,
				   int *innormalloop,
				   int *outnormalloop,
				   int (*inoutlohi)[NUMUPDOWN][NDIM],
				   int riin, int riout, int rjin, int rjout, int rkin, int rkout,
				   int *dosetbc,
				   int enerregion,
				   int *localenerpos
				   )
{


#pragma omp parallel  // assume don't require EOS
  {
    int i,j,k,pl,pliter;
    FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
    int failreturn;
#endif
    int ri, rj, rk; // reference i,j,k
    FTYPE prescale[NPR];
    int jj,kk;

    extern FTYPE RADNT_OMSCALE;
    extern FTYPE RADNT_FULLPHI;

  
    if( (BCtype[X1DN]==CYLAXIS) ){



      /* inner radial BC (preserves u^t rho and u) */
      if ( (totalsize[1]>1) && (mycpupos[1] == 0) ) {
		////////	LOOPX1dir{

		{ // start block
		  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
		  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
		  OPENMPBCLOOPBLOCK{
			OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);

			rj=j;
			if(RADNT_FULLPHI==2.0*Pi){
              rk=k+N3/2;
              if(rk>=N3) rk-=N3;
            }
            else rk=k;
			LOOPBOUND1IN{
			  PBOUNDLOOP(pliter,pl){
				// SECTIONMARK: assume r=0 singularity can't move
				if(dirprim[pl]==FACE1 || dirprim[pl]==CORN3 || dirprim[pl]==CORN2 || dirprim[pl]==CORNT ) ri = -i; // FACE1 values
				else ri=-i-1; // "CENT" values for purposes of this BC
				MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
			  }// over pl
			}// over boundary zones
		  }
		}// end block



		if( (BCtype[X1DN]==CYLAXIS) ){

		  /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
		  ////	LOOPX1dir{

		  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
		  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
		  OPENMPBCLOOPBLOCK{
			OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


			// SECTIONMARK: assume r=0 singularity can't move
			i=0;
			PBOUNDLOOP(pliter,pl){
			  if(pl==U1 || pl==B1){
				if(dirprim[pl]==FACE1 || dirprim[pl]==CORN3 || dirprim[pl]==CORN2 || dirprim[pl]==CORNT ){
				  MACP0A1(prim,i,j,k,pl) = 0.0;
				}
			  }// else don't do this pl
			} // end over pl
	
			LOOPBOUND1IN {
			  PBOUNDLOOP(pliter,pl){
				if(pl==U1 || pl==B1){
				  MACP0A1(prim,i,j,k,pl) *= -1.;
				}// end if right pl
			  } // end over pl
			} // end over boundary zones
		  }// end loop 23
		}
      } //end if inner CPU wall
    }
    else{
      dualfprintf(fail_file,"Shouldn't be here in bounds\n");
      myexit(3946840);
    }

  } // end parallel region

  // if full 2pi, then for MPI would have to deal with separately.

  return(0);

}








// X1 upper for RADCYLBEAM
int bound_x1up_radcylbeam(
						  int boundstage, int finalstep, SFTYPE boundtime, int whichdir, int boundvartype, int *dirprim, int ispstag, FTYPE (*prim)[NSTORE2][NSTORE3][NPR],
						  int *inboundloop,
						  int *outboundloop,
						  int *innormalloop,
						  int *outnormalloop,
						  int (*inoutlohi)[NUMUPDOWN][NDIM],
						  int riin, int riout, int rjin, int rjout, int rkin, int rkout,
						  int *dosetbc,
						  int enerregion,
						  int *localenerpos
						  )
  
{





#pragma omp parallel  // assume don't require EOS
  {

	int i,j,k,pl,pliter;
	FTYPE vcon[NDIM],X[NDIM],V[NDIM]; 
#if(WHICHVEL==VEL3)
	int failreturn;
#endif
	int ri, rj, rk; // reference i,j,k
	FTYPE prescale[NPR];
	int jj,kk;
	struct of_geom geomdontuse[NPR];
	struct of_geom *ptrgeom[NPR];
	struct of_geom rgeomdontuse[NPR];
	struct of_geom *ptrrgeom[NPR];

	// assign memory
	PALLLOOP(pl){
	  ptrgeom[pl]=&(geomdontuse[pl]);
	  ptrrgeom[pl]=&(rgeomdontuse[pl]);
	}


    extern FTYPE RADNT_OMSCALE;
    extern FTYPE RADNT_FULLPHI;


	if(BCtype[X1UP]==RADCYLBEAMBC && (totalsize[1]>1) && (mycpupos[1] == ncpux1-1) ){


	  OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
	  ////////	LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
	  OPENMPBCLOOPBLOCK{
		OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


		ri=riout;
		rj=j;
		rk=k;


		// ptrrgeom : i.e. ref geom
		PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        FTYPE *pr;
		LOOPBOUND1OUT{
          
          pr = &MACP0A1(prim,i,j,k,0);

    
		  //initially copying everything
		  PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

		  if(ispstag==0){
			// local geom
			PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
	  
			//coordinates of the ghost cell
			bl_coord_ijk_2(i,j,k,CENT,X, V);


            FTYPE xx,yy,zz,rsq;
            coord(i, j, k, CENT, X);
            bl_coord(X, V);
            xx=V[1];
            yy=V[2];
            zz=V[3];


            pr[RHO] = 1;
            pr[UU] = 0.1;

            //Keplerian gas
            FTYPE rCYL=V[1];
            FTYPE Om=RADNT_OMSCALE/(a+pow(rCYL,1.5));
            
            pr[U1] = 0.0;
            pr[U2] = 0.0;
            pr[U3] = Om;

            //E, F^i in orthonormal fluid frame
            FTYPE pradffortho[NPR];
            pradffortho[PRAD0] = calc_LTE_EfromT(1.e10/TEMPBAR);
            pradffortho[PRAD0] = 1.0;
            pradffortho[PRAD1] = 0;
            pradffortho[PRAD2] = 0;
            pradffortho[PRAD3] = 0;

            int whichvel=VEL3;
            int whichcoordfluid=MCOORD;
            int whichcoordrad=whichcoordfluid;
            whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

		  }// end if not staggered fields

		}// end loop over outer i's

	  }// end over loop
	}// end if correct boundary condition and core
	  
	  
	  
  }// end parallel region
  


 



  return(0);
} 
