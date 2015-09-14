

/*! \file bounds.koral.c
  \brief User Boundary conditions

  Also calls general/frequently-used functions in bounds.tools.c

*/

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


int bound_radbeam2dksvertbeaminflow(int dir,
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

int bound_radcylbeamcart(int dir,
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

int bound_staticset(int dir,
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




int bound_waldmono(int dir,
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

  if(WHICHPROBLEM==RADDONUT && OUTERDEATH==1){
    if(whichdir==1 && boundvartype!=BOUNDPSTAGTYPE && boundvartype!=BOUNDPSTAGSIMPLETYPE){// assumes always calls whichdir==1, which is true if N1>1
      extern void debugfixupaltdeath_bc(FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
      debugfixupaltdeath_bc(prim);
    }
  }

  return(0);
}


/// assume single user function takes care of primitive locations
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



/// user boundary routine
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
      else if(BCtype[dir]==HORIZONOUTFLOWSTATIC){
        BCtype[dir]=HORIZONOUTFLOW;
        bound_x1dn_outflow(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        BCtype[dir]=HORIZONOUTFLOWSTATIC;
        bound_staticset(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
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
      else if(BCtype[dir]==RADCYLBEAMCARTBC){
        bound_radcylbeamcart(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
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
      else if(BCtype[dir]==HORIZONOUTFLOW){
        bound_x1up_outflow(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        donebc[dir]=1;
      }
      else if(BCtype[dir]==HORIZONOUTFLOWSTATIC){
        BCtype[dir]=HORIZONOUTFLOW;
        bound_x1up_outflow(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        BCtype[dir]=HORIZONOUTFLOWSTATIC;
        bound_staticset(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
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
      else if(BCtype[dir]==RADCYLBEAMCARTBC){
        bound_radcylbeamcart(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
        donebc[dir]=1;
      }
      else if(BCtype[dir]==WALDMONOBC){
        bound_waldmono(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADCYLJETBC){
        bound_x1up_radcyljet(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
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
      else if(BCtype[dir]==RADCYLBEAMCARTBC){
        bound_radcylbeamcart(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADCYLJETBC){
        bound_x2dn_radcyljet(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
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
      else if(BCtype[dir]==RADBEAM2DKSVERTBEAMINFLOW){
        bound_radbeam2dksvertbeaminflow(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
        donebc[dir]=1;
      }
      else if(BCtype[dir]==RADCYLBEAMCARTBC){
        bound_radcylbeamcart(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
        donebc[dir]=1;
      }
      else if(BCtype[dir]==WALDMONOBC){
        // do ASYMM first
        BCtype[dir]=ASYMM;
        bound_x2up_polaraxis(boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos);
        // now do anything else that needs to be done
        BCtype[dir]=WALDMONOBC;
        bound_waldmono(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
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
      else if(BCtype[dir]==RADCYLBEAMCARTBC){
        bound_radcylbeamcart(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
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
      else if(BCtype[dir]==RADCYLBEAMCARTBC){
        bound_radcylbeamcart(dir,boundstage,finalstep,boundtime,whichdir,boundvartype,dirprim,ispstag,prim,inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi,riin,riout,rjin,rjout,rkin,rkout,dosetbc,enerregion,localenerpos); 
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





/// see interpline.c
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


////Called after the MPI boundary routines
/// many things here are copied from above
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















/////////////////////////////////////
///
///  Koral specific physical boundary conditions
///
/////////////////////////////////////



/// X1 lower for radiation beam injection
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
      //////// LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



        ri=riin;
        rj=j;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

   
        
        LOOPBOUND1IN{
          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);


          if(ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);

            // set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates

            extern FTYPE RADBEAMFLAT_FRATIO,RADBEAMFLAT_ERAD, RADBEAMFLAT_RHO, RADBEAMFLAT_UU;

            //            FTYPE ERADAMB=RADBEAMFLAT_ERAD;
            //            FTYPE ERADINJ=1000.0*ERADAMB;

            FTYPE ERADAMB=RADBEAMFLAT_ERAD;
            //            FTYPE ERADINJ=1000.0*ERADAMB; // old harm setup
            FTYPE ERADINJ=100.0*ERADAMB; // koral's current setup



            if(0){
              // correct version in general

              pr[RHO] = RADBEAMFLAT_RHO ;
              pr[UU] = RADBEAMFLAT_UU;
              SLOOPA(jj) pr[U1+jj-1] = 0.0;

              //E, F^i in orthonormal fluid frame
              FTYPE pradffortho[NPR];
              FTYPE Fx=0,Fy=0,Fz=0;
              //primitives in whichvel,whichcoord
              if(V[2]>.4 && V[2]<.6){//beam to be imposed
                Fx=RADBEAMFLAT_FRATIO*ERADINJ;
                Fy=Fz=0.0;
                
                pradffortho[PRAD0] = ERADINJ;
                pradffortho[PRAD1] = Fx;
                pradffortho[PRAD2] = Fy;
                pradffortho[PRAD3] = Fz;
              }
              else{ //no beam
                Fx=Fy=Fz=0.0;
                pradffortho[PRAD0] = ERADAMB;
                pradffortho[PRAD1] = Fx;
                pradffortho[PRAD2] = Fy;
                pradffortho[PRAD3] = Fz;
              }
            
              int whichvel=VEL4;
              int whichcoordfluid=MCOORD;
              int whichcoordrad=whichcoordfluid;
              whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);
            }
            else if(0){
              // set vradx

              pr[RHO] = RADBEAMFLAT_RHO ;
              pr[UU] = RADBEAMFLAT_UU;
              SLOOPA(jj) pr[U1+jj-1] = 0.0;

              // assume RADBEAMFLAT_FRATIO is just vradx
              FTYPE uradx=1.0/sqrt(1.0 - RADBEAMFLAT_FRATIO*RADBEAMFLAT_FRATIO); // radiation 4-velocity

              //primitives in whichvel,whichcoord
              if(V[2]>.4 && V[2]<.6){//beam to be imposed
                pr[URAD0] = ERADINJ;
                pr[URAD1] = uradx;
                pr[URAD2] = 0.;
                pr[URAD3] = 0.;
              }
              else{ //no beam
                pr[URAD0] = ERADAMB;
                pr[URAD1] = 0.;
                pr[URAD2] = 0.;
                pr[URAD3] = 0.;
              }

              // get all primitives in WHICHVEL/PRIMECOORDS value
              int whichvel=VEL4;
              int whichcoord=MCOORD;
              primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
            }
            else if(1){
              // koral mixed way (must use ff ortho choice in init.koral.c (i.e. first if(1))

              pr[RHO] = RADBEAMFLAT_RHO ;
              pr[UU] = RADBEAMFLAT_UU;
              SLOOPA(jj) pr[U1+jj-1] = 0.0;


              //primitives in whichvel,whichcoord
              if(V[2]>.4 && V[2]<.6){//beam to be imposed
                FTYPE dxdxp[NDIM][NDIM];
                dxdxprim_ijk(i, j, k, CENT, dxdxp);

                // like koral
                //                FTYPE uradx=ERADINJ*RADBEAMFLAT_FRATIO;

                // fix:
                FTYPE uradx=1.0/sqrt(1.0 - RADBEAMFLAT_FRATIO*RADBEAMFLAT_FRATIO); // radiation 4-velocity


                pr[URAD0] = ERADINJ;
                pr[URAD1] = uradx/dxdxp[1][1];
                pr[URAD2] = 0.;
                pr[URAD3] = 0.;
              }
              else{ //no beam and ERADAMB is assumed as radiation frame -- as in koral.
                pr[URAD0] = ERADAMB;
                pr[URAD1] = 0.;
                pr[URAD2] = 0.;
                pr[URAD3] = 0.;
              }

            }


          }// end if not staggered field


        }// end loop over inner i's
      }
    }


  }// end parallel region

  return(0);
} 




/// X1 lower for beam injection to create shadow
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
  extern FTYPE RADSHADOW_NLEFT,RADSHADOW_ANGLE;
  extern FTYPE RADSHADOW_TLEFTOTAMB;
  extern FTYPE RADSHADOW_BEAMY;

  extern FTYPE RADDBLSHADOW_NLEFT,RADDBLSHADOW_ANGLE;
  extern FTYPE RADDBLSHADOW_TLEFTOTAMB;
  extern FTYPE RADDBLSHADOW_BEAMY;

  FTYPE NLEFT,angle,TLEFT,BEAMY;

  if(WHICHPROBLEM==RADSHADOW){
    NLEFT=RADSHADOW_NLEFT;
    angle=RADSHADOW_ANGLE;
    TLEFT=TAMB*RADSHADOW_TLEFTOTAMB;
    BEAMY=RADSHADOW_BEAMY;
  }
  else if(WHICHPROBLEM==RADDBLSHADOW){
    NLEFT=RADDBLSHADOW_NLEFT;
    angle=RADDBLSHADOW_ANGLE;
    TLEFT=TAMB*RADDBLSHADOW_TLEFTOTAMB;
    BEAMY=RADDBLSHADOW_BEAMY;
  }


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
      //////// LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


        ri=riin;
        rj=j;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

   
        LOOPBOUND1IN{
          FTYPE *pr=&MACP0A1(prim,i,j,k,0);
    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          if(ispstag==0){
            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
   
            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);
            FTYPE xx,yy,zz,rsq;
            xx=V[1];
            yy=V[2];
            zz=V[3];

            // default for beam or not
            rsq=xx*xx+yy*yy+zz*zz;
            rho=(RHOBLOB-RHOAMB)*exp(-sqrt(rsq)/(BLOBW*BLOBW))+RHOAMB;
            pr[RHO] = rho;
            Trad=TAMB*RHOAMB/rho;
            uint=calc_PEQ_ufromTrho(Trad,rho);
            pr[UU] = uint;

     
            if(yy>BEAMY && WHICHPROBLEM==RADDBLSHADOW || WHICHPROBLEM==RADSHADOW ){

              FTYPE ERAD;
              ERAD=calc_LTE_EfromT(TLEFT);

              FTYPE ux=0.0; // orthonormal 4-velocity.  Matches init.koral.c
              pr[U1] = ux/sqrt(ptrgeom[U1]->gcov[GIND(1,1)]); // assumed no spatial mixing

              //E, F^i
              if(1){
                // correct way of using ERAD(TLEFT) and NLEFT
                FTYPE Fx=NLEFT*ERAD/sqrt(1+angle*angle);
                FTYPE Fy=-NLEFT*ERAD*angle/sqrt(1+angle*angle);
                FTYPE Fz=0.0;

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
              }
              else{
                // old harm way
                FTYPE gammax=1.0/sqrt(1.0-NLEFT*NLEFT);
                FTYPE uradx=NLEFT*gammax/sqrt(1.0+angle*angle);
                FTYPE urady=-NLEFT*gammax*angle/sqrt(1.0+angle*angle);
                
                pr[URAD0] = ERAD;
                pr[URAD1] = uradx;
                pr[URAD2] = urady;
                pr[URAD3] = 0.;

                // get all primitives in WHICHVEL/PRIMECOORDS value
                int whichvel;
                whichvel=VEL4;
                int whichcoord;
                whichcoord=MCOORD;
                primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
              }

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
      //////// LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


        ri=i;
        rj=rjout;
        rk=k;
    
    
        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

   
        LOOPBOUND2OUT{
          FTYPE *pr=&MACP0A1(prim,i,j,k,0);
    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          if(ispstag==0){
            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
   
            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);
            FTYPE xx,yy,zz,rsq;
            xx=V[1];
            yy=V[2];
            zz=V[3];



            // default for beam or not
            rsq=xx*xx+yy*yy+zz*zz;
            rho=(RHOBLOB-RHOAMB)*exp(-sqrt(rsq)/(BLOBW*BLOBW))+RHOAMB;
            pr[RHO] = rho;
            Trad=TAMB*RHOAMB/rho;
            uint= calc_PEQ_ufromTrho(Trad,rho);
            pr[UU] = uint;

     
            if(WHICHPROBLEM==RADDBLSHADOW || WHICHPROBLEM==RADSHADOW){

              FTYPE ERAD;
              ERAD=calc_LTE_EfromT(TLEFT);

              FTYPE ux=0.0; // orthonormal 4-velocity.  Matches init.koral.c
              pr[U1] = ux/sqrt(ptrgeom[U1]->gcov[GIND(1,1)]); // assumed no spatial mixing

              //E, F^i
              if(1){
                // correct way of using ERAD(TLEFT) and NLEFT
                FTYPE Fx=NLEFT*ERAD/sqrt(1+angle*angle);
                FTYPE Fy=-NLEFT*ERAD*angle/sqrt(1+angle*angle);
                FTYPE Fz=0.0;

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
              }
              else{
                // old harm way
                FTYPE gammax=1.0/sqrt(1.0-NLEFT*NLEFT);
                FTYPE uradx=NLEFT*gammax/sqrt(1.0+angle*angle);
                FTYPE urady=-NLEFT*gammax*angle/sqrt(1.0+angle*angle);
                
                pr[URAD0] = ERAD;
                pr[URAD1] = uradx;
                pr[URAD2] = urady;
                pr[URAD3] = 0.;

                // get all primitives in WHICHVEL/PRIMECOORDS value
                int whichvel;
                whichvel=VEL4;
                int whichcoord;
                whichcoord=MCOORD;
                primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
              }



            }// if spatially relevant region

          }// end if not staggered fields
   
        }// end loop over inner j's
      }// end over loop
    }// end if correct boundary condition and core
   
   
  }// end parallel region
  


 



  return(0);
} 





/// X3 lower for radiation beam injection
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
    FTYPE X[NDIM],V[NDIM]; 
    int ri, rj, rk; // reference i,j,k
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

    extern int RADBEAM2D_BEAMNO;
    extern int RADBEAM2D_FLATBACKGROUND;
    extern FTYPE RADBEAM2D_RHOAMB;
    extern FTYPE RADBEAM2D_TAMB;
    extern int RADBEAM2D_BLOB;
    extern FTYPE RADBEAM2D_BLOBW;
    extern FTYPE RADBEAM2D_BLOBP;
    extern FTYPE RADBEAM2D_BLOBX;
    extern FTYPE RADBEAM2D_BLOBZ;
    extern FTYPE RADBEAM2D_PAR_D;
    extern FTYPE RADBEAM2D_PAR_E;
    extern int RADBEAM2D_IFBEAM;
    extern FTYPE RADBEAM2D_TLEFT;
    extern FTYPE RADBEAM2D_NLEFT;
    extern FTYPE RADBEAM2D_BEAML;
    extern FTYPE RADBEAM2D_BEAMR;

  
    if(dir==X3DN && BCtype[X3DN]==RADBEAM2DBEAMINFLOW && totalsize[3]>1 && mycpupos[3] == 0 ){


   

      OPENMPBCLOOPVARSDEFINELOOPX3DIR; OPENMPBCLOOPSETUPLOOPX3DIR;
      //////// LOOPX3dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX3DIR(i,j);

        ri=i;
        rj=j;
        rk=rkin;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

        
        LOOPBOUND3IN{
          FTYPE *pr=&MACP0A1(prim,i,j,k,0);
    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) pr[pl] = MACP0A1(prim,ri,rj,rk,pl);


          if(0){
            // NOTEMARK: only really makes sense near the hole if in KSCOORDS
            PBOUNDLOOP(pliter,pl) if(pl==U3) if(pr[U3]>0.0) pr[U3]=0.0; // limit so no arbitrary fluid inflow
            PBOUNDLOOP(pliter,pl) if(pl==URAD3) if(pr[URAD3]>0.0) pr[URAD3]=0.0; // limit so no arbitrary radiative inflow
          }

          // store this outflow result
          FTYPE pr0[NPR];
          PBOUNDLOOP(pliter,pl) pr0[pl]=pr[pl];

          // only overwrite copy if not inside i<0 since want to keep consistent with outflow BCs used for other \phi at those i
          if(1      &&     startpos[1]+i>=0 && ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            int whichcoord;
            whichcoord=MCOORD;

            // get metric grid geometry for these ICs
            int getprim=0;
            struct of_geom geomrealdontuse;
            struct of_geom *ptrgeomreal=&geomrealdontuse;
            gset(getprim,whichcoord,i,j,k,ptrgeomreal);


            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);


            FTYPE ERADAMB;
            FTYPE rho,uint,Vr;
            if(RADBEAM2D_FLATBACKGROUND){
              Vr=0.0;
              rho=RADBEAM2D_RHOAMB;
              uint=calc_PEQ_ufromTrho(RADBEAM2D_TAMB,rho);
              ERADAMB=calc_LTE_EfromT(RADBEAM2D_TAMB);

              // override so like outflow conditions to avoid shear at boundary
              if(0) ERADAMB=pr[URAD0]; // only makes sense with urad method

            }
            else{
              //zaczynam jednak od profilu analitycznego:   
              FTYPE r=V[1];
              FTYPE mD=RADBEAM2D_PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
              FTYPE mE=RADBEAM2D_PAR_E/(pow(r*r*sqrt(2./r),gamideal)*pow(1.-2./r,(gamideal+1.)/4.));
              Vr=sqrt(2./r)*(1.-2./r);


              FTYPE W=1./sqrt(1.-Vr*Vr*ptrgeomreal->gcov[GIND(1,1)]); // assumes RHO location is good for all these quantities
              rho=RADBEAM2D_PAR_D/(r*r*sqrt(2./r));
              FTYPE T=RADBEAM2D_TAMB;
              //   FTYPE ERAD=calc_LTE_EfromT(T);
              uint=mE/W;
              ERADAMB=calc_LTE_Efromurho(uint,rho);
     
              // override so like outflow conditions to avoid shear at boundary
              if(0) ERADAMB=pr[URAD0];
            }


            // beam parameters
            FTYPE ERADINJ=calc_LTE_EfromT(RADBEAM2D_TLEFT);
            // original top-hat beam
            //FTYPE beamshape=(FTYPE)(V[1]>RADBEAM2D_BEAML && V[1]<RADBEAM2D_BEAMR && RADBEAM2D_IFBEAM);
            // gaussian-like beam:
            FTYPE powbeam=6.0;
            FTYPE beamhalfwidth=0.5*(RADBEAM2D_BEAMR-RADBEAM2D_BEAML);
            FTYPE beamcenter=(RADBEAM2D_BEAML+RADBEAM2D_BEAMR)*0.5;
            FTYPE beamshape=exp(-pow(V[1]-beamcenter,powbeam)/(2.0*pow(beamhalfwidth,powbeam)))*RADBEAM2D_IFBEAM;


#define INJECTINFLUIDFRAME 0 // need beam to go out as designed, not with fluid, so usually 0

            if(INJECTINFLUIDFRAME){
              //E, F^i in orthonormal fluid frame
              FTYPE Fx,Fy,Fz;
              // default flux
              Fx=Fy=Fz=0;
              // beam flux
              Fz=RADBEAM2D_NLEFT*ERADINJ;
              
              FTYPE pradffortho[NPR];
              pradffortho[PRAD0] = ERADAMB + ERADINJ*beamshape;
              pradffortho[PRAD1] = Fx*beamshape;
              pradffortho[PRAD2] = Fy*beamshape;
              pradffortho[PRAD3] = Fz*beamshape;

              int whichvel=WHICHVEL; // in which vel U1-U3 set
              int whichcoordfluid=PRIMECOORDS; // in which coordinates U1-U3 set
              int whichcoordrad=MCOORD; // in which coordinates E,F are orthonormal
              whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

        
#if(0)
              // try using outflow outside of beam
              if(){
              }
              else{
                pr[PRAD0]=pr0[PRAD0];
                pr[PRAD1]=pr0[PRAD1];
                pr[PRAD2]=pr0[PRAD2];
                pr[PRAD3]=pr0[PRAD3];
                if(pr[PRAD3]>0.0) pr[PRAD3]=0.0; // but don't let radiative inflow
              }
#endif

            }
            else{

              int whichvel;
              whichvel=VEL4;

              // get coordinate basis in VEL4 format
              FTYPE uradcon[NDIM],othersrad[NUMOTHERSTATERESULTS];
              ucon_calc(&pr[URAD1-U1],ptrgeom[URAD1],uradcon,othersrad);
              // get coordinate basis in MCOORD basis
              FTYPE uradx,urady,uradz;
              uradx=urady=0.0;

              uradx=uradcon[1]*sqrt(fabs(ptrgeom[URAD1]->gcov[GIND(1,1)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(1,1)]));
              urady=uradcon[2]*sqrt(fabs(ptrgeom[URAD2]->gcov[GIND(2,2)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(2,2)]));
              uradz=uradcon[3]*sqrt(fabs(ptrgeom[URAD3]->gcov[GIND(3,3)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(3,3)]));
              if(uradz>0.0) uradz=0.0; // limit so no arbitrary radiative inflow

              // override uradz
              if(V[1]>RADBEAM2D_BEAML && V[1]<RADBEAM2D_BEAMR && RADBEAM2D_IFBEAM) uradz=1.0/sqrt(1.0 - RADBEAM2D_NLEFT*RADBEAM2D_NLEFT);
              //              else uradz=0.0;

              PBOUNDLOOP(pliter,pl) if(pl==URAD0) pr[URAD0] = ERADINJ;
              PBOUNDLOOP(pliter,pl) if(pl==URAD1) pr[URAD1] = uradx;
              PBOUNDLOOP(pliter,pl) if(pl==URAD2) pr[URAD2] = urady;
              PBOUNDLOOP(pliter,pl) if(pl==URAD3) pr[URAD3] = uradz;

              // get all primitives in WHICHVEL/PRIMECOORDS value
              primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],MAC(prim,i,j,k),MAC(prim,i,j,k)); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
            }


          }// end if not staggered field


        }// end loop over inner i's
      } // over block
    }// end if correct BC and should be doind BC for this core
  }// end parallel region

  return(0);
} 





/// X1 upper for inflow
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


    extern int RADBEAM2D_BEAMNO;
    extern int RADBEAM2D_FLATBACKGROUND;
    extern FTYPE RADBEAM2D_RHOAMB;
    extern FTYPE RADBEAM2D_TAMB;
    extern int RADBEAM2D_BLOB;
    extern FTYPE RADBEAM2D_BLOBW;
    extern FTYPE RADBEAM2D_BLOBP;
    extern FTYPE RADBEAM2D_BLOBX;
    extern FTYPE RADBEAM2D_BLOBZ;
    extern FTYPE RADBEAM2D_PAR_D;
    extern FTYPE RADBEAM2D_PAR_E;
    extern int RADBEAM2D_IFBEAM;
    extern FTYPE RADBEAM2D_TLEFT;
    extern FTYPE RADBEAM2D_NLEFT;
    extern FTYPE RADBEAM2D_BEAML;
    extern FTYPE RADBEAM2D_BEAMR;
  
    if(dir==X1UP && BCtype[X1UP]==RADBEAM2DFLOWINFLOW && totalsize[1]>1 && mycpupos[1] == ncpux1-1 ){

      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);



        ri=riout;
        rj=j;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

   
        LOOPBOUND1OUT{
          FTYPE *pr=&MACP0A1(prim,i,j,k,0);
    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);
    

          if(ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);

            int whichcoord;
            whichcoord=MCOORD;


            FTYPE ERADAMB;
            FTYPE rho,uint,Vr;
            if(RADBEAM2D_FLATBACKGROUND){
              Vr=0.0;
              rho=RADBEAM2D_RHOAMB;
              uint=calc_PEQ_ufromTrho(RADBEAM2D_TAMB,rho);
              ERADAMB=calc_LTE_EfromT(RADBEAM2D_TAMB);
            }
            else{
              //zaczynam jednak od profilu analitycznego:   
              FTYPE r=V[1];
              FTYPE mD=RADBEAM2D_PAR_D/(r*r*sqrt(2./r*(1.-2./r)));
              FTYPE mE=RADBEAM2D_PAR_E/(pow(r*r*sqrt(2./r),gamideal)*pow(1.-2./r,(gamideal+1.)/4.));
              Vr=sqrt(2./r)*(1.-2./r);

              // get metric grid geometry for these ICs
              int getprim=0;
              struct of_geom geomrealdontuse;
              struct of_geom *ptrgeomreal=&geomrealdontuse;
              gset(getprim,whichcoord,i,j,k,ptrgeomreal);

              FTYPE W=1./sqrt(1.-Vr*Vr*ptrgeomreal->gcov[GIND(1,1)]); // assumes RHO location is good for all these quantities
              rho=RADBEAM2D_PAR_D/(r*r*sqrt(2./r));
              FTYPE T=RADBEAM2D_TAMB;
              //   FTYPE ERAD=calc_LTE_EfromT(T);
              uint=mE/W;
              ERADAMB=calc_LTE_Efromurho(uint,rho);
            }

            // set quantities at outer radial edge
            pr[RHO] = rho;
            pr[UU]  = uint;
            pr[U1]  = -Vr;
            pr[U2]  = 0.;
            pr[U3]  = 0.;


            FTYPE ERAD=ERADAMB;


            if(1){
              //E, F^i in orthonormal fluid frame
              FTYPE Fx,Fy,Fz;
              // default flux
              Fx=Fy=Fz=0;

              FTYPE pradffortho[NPR];
              pradffortho[PRAD0] = ERAD;
              pradffortho[PRAD1] = Fx;
              pradffortho[PRAD2] = Fy;
              pradffortho[PRAD3] = Fz;

              int whichvel=VEL4;
              int whichcoordfluid=whichcoord; // in which coordinates U1-U3 set
              int whichcoordrad=whichcoord; // in which coordinates E,F are orthonormal
              whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

            }
            else if(0){

              FTYPE uradx,urady,uradz;
              uradx=urady=uradz=0.0;
              
              pr[PRAD0] = ERAD;
              pr[PRAD1] = uradx;
              pr[PRAD2] = urady;
              pr[PRAD3] = uradz;

              // get all primitives in WHICHVEL/PRIMECOORDS value
              int whichvel;
              whichvel=VEL4;
              if (bl2met2metp2v(whichvel, whichcoord,pr, i,j,k) >= 1){
                FAILSTATEMENT("bounds.koral.c:bound_radbeam2dflowinflow()", "bl2ks2ksp2v()", 1);
              }
            }   


          }// end if not staggered field


        }// end loop over inner i's
      }// over block
    }// if correct BC and core
  }// end parallel region

  return(0);
} 





/// X1 lower and upper for RADATM
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
      //////// LOOPX1dir{
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
      //////// LOOPX1dir{
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





/// X1 lower and upper X2
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
      //////// LOOPX1dir{
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
      //////// LOOPX2dir{
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









/// X1 upper for inflow
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
      extern FTYPE RADBONDI_MINX;
      extern FTYPE RADBONDI_MAXX;


      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
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

            //            dualfprintf(fail_file,"rho=%g uint=%g vx=%g ERAD=%g\n",rho,uint,vx,ERAD);

            int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
            whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

            //            PLOOP(pliter,pl) dualfprintf(fail_file,"pl=%d pr=%g\n",pl,pr[pl]);

          }// end if not staggered field


        }// end loop over inner i's
      }// over block
    }// if correct BC and core
  }// end parallel region

  return(0);
} 






/// on-grid bounding
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






/// RADNT
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
    extern FTYPE RADNT_RHODONUT;
    extern FTYPE RADNT_UINTATMMIN;
    extern FTYPE RADNT_ERADATMMIN;
    extern FTYPE RADNT_DONUTTYPE;
    extern FTYPE RADNT_INFLOWING;
    extern FTYPE RADNT_TGASATMMIN;
    extern FTYPE RADNT_TRADATMMIN;
    extern FTYPE RADNT_ROUT;
    extern FTYPE RADNT_OMSCALE;
    extern FTYPE RADNT_FULLPHI;
    extern FTYPE RADNT_DONUTRADPMAX;
    extern FTYPE RADNT_HOVERR;
    extern FTYPE RADNT_LPOW;


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
      //////// LOOPX1dir{
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


            if(pr[PRAD1]<0.0){ // if active cell (now ghost cell) inflows, set values, keep outflow value
              if(1){
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
              }
              else{
                FTYPE gammamax=10.;
                FTYPE rout=2.; //normalize at r_BL=2
                // RADNT_ERADATMMIN here, as in koral, assumes in ut frame, so not really consistent with initial conditions
                pr[PRAD0]=RADNT_ERADATMMIN*(rout/r)*(rout/r)*(rout/r)*(rout/r);

                FTYPE ut[NDIM]={0.,-gammamax*pow(r/rout,1.),0.,0.}; // assume in BLCOORDS 4-vel or SPCMINKMETRIC if no gravity
                SLOOPA(jj) pr[URAD1+jj-1]=ut[jj]; 

                // get all primitives in WHICHVEL/PRIMECOORDS value
                if (bl2met2metp2v(whichvel, whichcoord,pr, i,j,k) >= 1){
                  FAILSTATEMENT("bounds.koral.c:bound_radnt()", "bl2ks2ksp2v()", 1);
                }
              }
            }
            else{
              // keep outflow value already set
            }
            

          }// end if not staggered field


        }// end loop over inner i's
      }// over block
    }// if correct BC and core







    if(dir==X2UP && BCtype[X2UP]==RADNTBC && (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ){

      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      //////// LOOPX2dir{
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

            FTYPE rin,rout;
            int conddisk;
            if(WHICHPROBLEM==RADFLATDISK){
              rin=15.;
              rout=25.;
              conddisk=(r<rout); // as in new koral
            }
            else{
              rin=6.;
              rout=1E10;
              conddisk=(r>rin);
            }

            //hot boundary
            if(conddisk){

              //E, F^i in orthonormal fluid frame
              FTYPE pradffortho[NPR];
              if(WHICHPROBLEM==RADFLATDISK) pradffortho[PRAD0] = calc_LTE_EfromT(1.e11/TEMPBAR);
              else pradffortho[PRAD0] = calc_LTE_EfromT(1.e11/TEMPBAR)*(1.-sqrt(rin/r))/pow(r,3.);
              // KORALTODO: in reality, can only constrain magnitude, not direction, of outflow away from plane.
              // KORALTODO: Also, in reality, can't really set vrad=0, else solution dominated by numerical diffusion.
              pradffortho[PRAD1] = 0;
              if(WHICHPROBLEM==RADFLATDISK) pradffortho[PRAD2] = 0.0; // current koral value
              else pradffortho[PRAD2] = -0.5*pradffortho[PRAD0];
              pradffortho[PRAD3] = 0;


              //pr[RHO] and pr[UU] remain same as from ASYMM condition as well as any field
              //Keplerian gas with no inflow or outflow
              // KORALTODO: in reality, can only constrain magnitude, not direction, of outflow away from plane.  But since magnitude is chosen to be zero, then no issue here.
              pr[U1]=pr[U2]=0.0; // have to be careful with this for VEL3 (must have rin>>rergo).
              if(WHICHPROBLEM==RADFLATDISK) pr[U3]=0.0; // current koral value
              else pr[U3]=1./(a + pow(r,1.5));
 
              int whichvel;
              whichvel=VEL3; // VEL3 so can set Keplerian rotation rate

              int whichcoordfluid;
              if(WHICHPROBLEM==RADFLATDISK) whichcoordfluid=MCOORD; // whatever else
              else whichcoordfluid=BLCOORDS; // want to setup things in BLCOORDS

              //              dualfprintf(fail_file,"FLAT: rho=%g uint=%g prad0=%g Eff=%g\n",pr[RHO],pr[UU]/pr[RHO],pr[PRAD0]/pr[RHO],pradffortho[PRAD0]/pr[RHO]);

              int whichcoordrad=whichcoordfluid; // in which coordinates E,F are orthonormal
              whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

              //              dualfprintf(fail_file,"FLATPOST: rho=%g uint=%g Eff=%g\n",pr[RHO],pr[UU]/pr[RHO],pr[PRAD0]/pr[RHO]);
              //              dualfprintf(fail_file,"opacity : %g\n",pr[RHO]*KAPPA_ES_CODE(rho,T)/1E14*0.1);
              //              dualfprintf(fail_file,"LBAR: %g\n",LBAR);

            } // end if actually doing something to boundary cells in "hot" boundary

          }// end if not staggered field


        }// end loop over outer j's
      }
    }







  }// end parallel region

  return(0);
} 







/// X1 inner CYLAXIS
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
        //////// LOOPX1dir{

        { // start block
          OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
          //////// LOOPX1dir{
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
          //// LOOPX1dir{

          OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
          //////// LOOPX1dir{
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








/// X1 upper for RADCYLBEAM
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
      //////// LOOPX1dir{
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

/// X1 upper for RADCYLJET
int bound_x1up_radcyljet(
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



  extern FTYPE RADCYLJET_TYPE;


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



    if(BCtype[X1UP]==RADCYLJETBC && (totalsize[1]>1) && (mycpupos[1] == ncpux1-1) ){


      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
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


            FTYPE pradffortho[NPR];
            if(RADCYLJET_TYPE==2){
              //            pr[RHO] = 1;
              //            pr[UU] = 0.1;

              pr[U1] = 0.0;
              pr[U2] = 0.0;
              pr[U3] = 0.0;

              //E, F^i in orthonormal fluid frame
              pradffortho[PRAD0] = calc_LTE_EfromT(1.e10/TEMPBAR);
              pradffortho[PRAD0] = 1.0; // should be compared to Ehatjet in init.koral.c
              pradffortho[PRAD1] = 0;
              pradffortho[PRAD2] = 0;
              pradffortho[PRAD3] = 0;
            }
            if(RADCYLJET_TYPE==3){
              pr[RHO] = 0.1; // made same as rhojet to keep density flatish
              pr[UU] = 0.1; // ""

              pr[U1] = 0.0;
              pr[U2] = 0.0;
              pr[U3] = 0.0;

              //E, F^i in orthonormal fluid frame
              pradffortho[PRAD0] = calc_LTE_EfromT(1.e10/TEMPBAR);
              pradffortho[PRAD0] = 1.0; // should be compared to Ehatjet in init.koral.c
              pradffortho[PRAD1] = -0.1*pradffortho[PRAD0];
              pradffortho[PRAD2] = 0;
              pradffortho[PRAD3] = 0;
            }


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


/// X2 lower for RADCYLJET
int bound_x2dn_radcyljet(
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



  extern FTYPE RADCYLJET_TYPE;


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



    if(BCtype[X2DN]==RADCYLJETBC && (totalsize[2]>1) && (mycpupos[2] == 0) ){


      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      //////// LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


        ri=i;
        rj=rjin;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        FTYPE *pr;
        LOOPBOUND2IN{
          
          pr = &MACP0A1(prim,i,j,k,0);

    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          if(ispstag==0){
            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            extern int jetbound(int i, int j, int k, int loc, FTYPE *prin, FTYPE *prflux, FTYPE (*prim)[NSTORE2][NSTORE3][NPR]);
            int insidejet=jetbound(i,j,k,CENT,MAC(prim,i,j,k),MAC(prim,i,j,k),prim);


          }// end if not staggered fields

        }// end loop over outer i's

      }// end over loop
    }// end if correct boundary condition and core
   
   
   
  }// end parallel region
  


 



  return(0);
} 











/// X2 upper for radiation beam injection
int bound_radbeam2dksvertbeaminflow(int dir,
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

    if(dir==X2UP && BCtype[X2UP]==RADBEAM2DKSVERTBEAMINFLOW && totalsize[2]>1 && mycpupos[2] == ncpux2-1 ){


      extern int RADBEAM2DKSVERT_BEAMNO,RADBEAM2D_FLATBACKGROUND;
      FTYPE RHOAMB=1.e0/RHOBAR;
      FTYPE TAMB=1e7/TEMPBAR;
      FTYPE PAR_D=1./RHOBAR;
      FTYPE PAR_E=1e-4/RHOBAR;

      // BEAM PROPERTIES
      int IFBEAM=1; // whether to have a beam
      FTYPE TLEFT=1e9/TEMPBAR;
      //      FTYPE NLEFT=0.99;
      FTYPE NLEFT=0.995;
      //   FTYPE NLEFT=0.999;
      //   FTYPE NLEFT=0.999999; // paper says this, while koral code says 0.999

      FTYPE BEAML,BEAMR;
      if (RADBEAM2DKSVERT_BEAMNO==1){
        BEAML=2.9;
        BEAMR=3.1;
      }
      else if (RADBEAM2DKSVERT_BEAMNO==2){
        BEAML=5.8;
        BEAMR=6.2;
      }
      else if (RADBEAM2DKSVERT_BEAMNO==3){
        BEAML=15.5;
        BEAMR=16.5;
      }
      else if (RADBEAM2DKSVERT_BEAMNO==4){
        BEAML=37;
        BEAMR=43;
      }
      else if (RADBEAM2DKSVERT_BEAMNO==5){
        BEAML=7.;
        BEAMR=9.;
      }

   

      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      //////// LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);

        ri=i;
        rj=rjout;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);

   
        LOOPBOUND2OUT{
          FTYPE *pr=&MACP0A1(prim,i,j,k,0);
    
          //initially copying everything
          PBOUNDLOOP(pliter,pl) pr[pl] = MACP0A1(prim,ri,rj,rk,pl);
    
          // NOTEMARK: only really makes sense near the hole if in KSCOORDS
          if(pr[U2]<0.0) pr[U2]=0.0; // limit so no arbitrary fluid inflow
          if(pr[URAD2]<0.0) pr[URAD2]=0.0; // limit so no arbitrary radiative inflow


          // only overwrite copy if not inside i<0 since want to keep consistent with outflow BCs used for other \theta at those i
          FTYPE ERADAMB;
          FTYPE rho,uint,Vr;
          if(1      &&     startpos[1]+i>=0 && ispstag==0){ // only do something special with non-field primitives

            // local geom
            PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);

            //coordinates of the ghost cell
            bl_coord_ijk_2(i,j,k,CENT,X, V);

            int whichcoord;
            whichcoord=MCOORD;

            // get metric grid geometry for these ICs
            int getprim=0;
            struct of_geom geomrealdontuse;
            struct of_geom *ptrgeomreal=&geomrealdontuse;
            gset(getprim,whichcoord,i,j,k,ptrgeomreal);


            if(RADBEAM2D_FLATBACKGROUND){
              Vr=0.0;
              rho=RHOAMB;
              uint=calc_PEQ_ufromTrho(TAMB,rho);
              ERADAMB=calc_LTE_EfromT(TAMB);

              // override so like outflow conditions to avoid shear at boundary
              // ERADAMB=pr[URAD0];

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
              //   FTYPE ERAD=calc_LTE_EfromT(T);
              uint=mE/W;
              ERADAMB=calc_LTE_Efromurho(uint,rho);
     
              // override so like outflow conditions to avoid shear at boundary
              //ERADAMB=pr[URAD0];
            }


            // beam parameters
            FTYPE ERADINJ=calc_LTE_EfromT(TLEFT);
            // original top-hat beam
            //FTYPE beamshape=(FTYPE)(V[1]>BEAML && V[1]<BEAMR && IFBEAM);
            // gaussian-like beam:
            FTYPE powbeam=6.0;
            FTYPE beamhalfwidth=0.5*(BEAMR-BEAML);
            FTYPE beamcenter=(BEAML+BEAMR)*0.5;
            FTYPE beamshape=exp(-pow(V[1]-beamcenter,powbeam)/(2.0*pow(beamhalfwidth,powbeam)))*IFBEAM;


            if(1){
              //E, F^i in orthonormal fluid frame
              FTYPE Fx,Fy,Fz;
              // default flux
              Fx=Fy=Fz=0;
              // beam flux
              Fy=-NLEFT*ERADINJ;

              //              dualfprintf(fail_file,"ERADINJ=%g Fy=%g beamshape=%g\n",ERADINJ, Fy,beamshape);
              
              FTYPE pradffortho[NPR];
              pradffortho[PRAD0] = ERADAMB + ERADINJ*beamshape;
              pradffortho[PRAD1] = Fx*beamshape;
              pradffortho[PRAD2] = Fy*beamshape;
              pradffortho[PRAD3] = Fz*beamshape;

              int whichvel=WHICHVEL; // in which vel U1-U3 set
              int whichcoordfluid=PRIMECOORDS; // in which coordinates U1-U3 set
              int whichcoordrad=MCOORD; // in which coordinates E,F are orthonormal
              whichfluid_ffrad_to_primeall(&whichvel, &whichcoordfluid, &whichcoordrad, ptrgeom[RHO], pradffortho, pr, pr);

              //              dualfprintf(fail_file,"primrad: %g %g %g %g\n",pr[PRAD0],pr[PRAD1],pr[PRAD2],pr[PRAD3]);

            }
            else{
              // old harm way

              // set radiation quantities as R^t_\nu in orthonormal fluid frame using whichvel velocity and whichcoord coordinates
              int whichvel;
              whichvel=VEL4;


              FTYPE uradx,urady,uradz;
              FTYPE uradcon[NDIM],othersrad[NUMOTHERSTATERESULTS];


              // get coordinate basis in VEL4 format
              ucon_calc(&pr[URAD1-U1],ptrgeom[URAD1],uradcon,othersrad);
              // get coordinate basis in MCOORD basis
              uradx=uradcon[1]*sqrt(fabs(ptrgeom[URAD1]->gcov[GIND(1,1)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(1,1)]));
              urady=uradcon[2]*sqrt(fabs(ptrgeom[URAD2]->gcov[GIND(2,2)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(2,2)]));
              uradz=uradcon[3]*sqrt(fabs(ptrgeom[URAD3]->gcov[GIND(3,3)]))/sqrt(fabs(ptrgeomreal->gcov[GIND(3,3)]));
              if(urady<0.0) urady=0.0; // limit so no arbitrary radiative inflow


              PBOUNDLOOP(pliter,pl){

                //primitives in whichvel,whichcoord
                if(V[1]>BEAML && V[1]<BEAMR && IFBEAM){//beam to be imposed
                
                  // override uradz
                  urady=-1.0/sqrt(1.0 - NLEFT*NLEFT);
                  uradx=uradz=0.0;
                

                  if(pl==URAD0) pr[URAD0] = ERADINJ;
                  if(pl==URAD1) pr[URAD1] = uradx;
                  if(pl==URAD2) pr[URAD2] = urady;
                  if(pl==URAD3) pr[URAD3] = uradz;

                  //                dualfprintf(fail_file,"GOT BEAM: ijk=%d %d %d : %g %g\n",ERADINJ,urady);
                }
                else{ //no beam
                
                  //                dualfprintf(fail_file,"GOTNOBEAM: ijk=%d %d %d : %g %g\n",ERADAMB,urady);
                
                  //     pr[URAD0] = ERADAMB;
                  if(pl==URAD0) pr[URAD0] = ERADAMB; // so matches outer radial boundary when no beam
                  if(pl==URAD1) pr[URAD1] = uradx;
                  if(pl==URAD2) pr[URAD2] = urady;
                  if(pl==URAD3) pr[URAD3] = uradz;
                }
              } // over allowed pl's to bound


              // KORALTODO: ERADINJ is in fluid frame, need to convert, but probably ok.

              // get all primitives in WHICHVEL/PRIMECOORDS value
              primefluid_EVrad_to_primeall(&whichvel, &whichcoord, ptrgeom[RHO],pr,pr); // assumes ptrgeom[RHO] is same location as all other primitives (as is currently true).
            }

          }// end if not staggered field


        }// end loop over inner i's
      } // over block
    }// end if correct BC and should be doind BC for this core
  }// end parallel region

  return(0);
}














/// all boundaries for RADCYLBEAMCART
int bound_radcylbeamcart(int dir,
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
    int ri, rj, rk; // reference i,j,k
    struct of_geom rgeomdontuse[NPR];
    struct of_geom *ptrrgeom[NPR];
    // assign memory
    PALLLOOP(pl)  ptrrgeom[pl]=&(rgeomdontuse[pl]);



    int get_radcylbeamcart(int dir, int *dirprim, int ispstag, int ri, int rj, int rk, int i, int j, int k, FTYPE *prref, FTYPE *pr);


    ////////
    // X1DN
    ////////
    if(dir==X1DN && totalsize[1]>1 && mycpupos[1] == 0){


      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


        ri=riin;
        rj=j;
        rk=k;

        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        LOOPBOUND1IN{

          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
          FTYPE *prref = &MACP0A1(prim,ri,rj,rk,0);
          get_radcylbeamcart(dir, dirprim, ispstag,ri,rj,rk,i,j,k,prref,pr);

        }// end loop over outer i's

      }// end over loop
    }// end if correct boundary condition and core



    ////////
    // X1UP
    ////////
    if(dir==X1UP && totalsize[1]>1 && mycpupos[1] == ncpux1-1){


      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX1DIR(j,k);


        ri=riout;
        rj=j;
        rk=k;

        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        LOOPBOUND1OUT{

          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
          FTYPE *prref = &MACP0A1(prim,ri,rj,rk,0);
          get_radcylbeamcart(dir, dirprim, ispstag,ri,rj,rk,i,j,k,prref,pr);

        }// end loop over outer i's

      }// end over loop
    }// end if correct boundary condition and core





    ////////
    // X2DN
    ////////
    if(dir==X2DN && totalsize[2]>1 && mycpupos[2] == 0){


      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      //////// LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


        ri=i;
        rj=rjin;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        LOOPBOUND2IN{
          
          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
          FTYPE *prref = &MACP0A1(prim,ri,rj,rk,0);
          get_radcylbeamcart(dir, dirprim, ispstag,ri,rj,rk,i,j,k,prref,pr);

        }// end inner loop

      }// end over loop
    }// end if correct boundary condition and core


    ////////
    // X2UP
    ////////
    if(dir==X2UP && totalsize[2]>1 && mycpupos[2] == ncpux2-1){


      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      //////// LOOPX2dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX2DIR(i,k);


        ri=i;
        rj=rjout;
        rk=k;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        LOOPBOUND2OUT{
          
          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
          FTYPE *prref = &MACP0A1(prim,ri,rj,rk,0);
          get_radcylbeamcart(dir, dirprim, ispstag,ri,rj,rk,i,j,k,prref,pr);

        }// end inner loop

      }// end over loop
    }// end if correct boundary condition and core

    ////////
    // X3DN
    ////////
    if(dir==X3DN && totalsize[3]>1 && mycpupos[3] == 0){


      OPENMPBCLOOPVARSDEFINELOOPX3DIR; OPENMPBCLOOPSETUPLOOPX3DIR;
      //////// LOOPX3dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX3DIR(i,j);


        ri=i;
        rj=j;
        rk=rkin;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        LOOPBOUND3IN{
          
          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
          FTYPE *prref = &MACP0A1(prim,ri,rj,rk,0);
          get_radcylbeamcart(dir, dirprim, ispstag,ri,rj,rk,i,j,k,prref,pr);

        }// end inner loop

      }// end over loop
    }// end if correct boundary condition and core


    ////////
    // X3UP
    ////////
    if(dir==X3UP && totalsize[3]>1 && mycpupos[3] == ncpux3-1){


      OPENMPBCLOOPVARSDEFINELOOPX3DIR; OPENMPBCLOOPSETUPLOOPX3DIR;
      //////// LOOPX3dir{
#pragma omp for schedule(OPENMPSCHEDULE(),OPENMPCHUNKSIZE(blocksize))
      OPENMPBCLOOPBLOCK{
        OPENMPBCLOOPBLOCK2IJKLOOPX3DIR(i,j);


        ri=i;
        rj=j;
        rk=rkout;


        // ptrrgeom : i.e. ref geom
        PALLLOOP(pl) get_geometry(ri, rj, rk, dirprim[pl], ptrrgeom[pl]);
        
        LOOPBOUND3OUT{
          
          FTYPE *pr = &MACP0A1(prim,i,j,k,0);
          FTYPE *prref = &MACP0A1(prim,ri,rj,rk,0);
          get_radcylbeamcart(dir, dirprim, ispstag,ri,rj,rk,i,j,k,prref,pr);

        }// end inner loop

      }// end over loop
    }// end if correct boundary condition and core


   
   
   
  }// end parallel region
  


 



  return(0);
} 



int get_radcylbeamcart(int dir, int *dirprim, int ispstag, int ri, int rj, int rk, int i, int j, int k, FTYPE *prref, FTYPE *pr)
{
  int pliter,pl;
  

  //initially copying everything
  PBOUNDLOOP(pliter,pl) pr[pl] = prref[pl];


  if(ispstag==1) return(0); // nothing else to do for now


  int src;
  if(dir==X1DN || dir==X1UP) src=(j>0.9*totalsize[2]/2. && j<1.1*totalsize[2]/2.);
  else if(dir==X2DN || dir==X2UP) src=(i>0.9*totalsize[1]/2. && i<1.1*totalsize[1]/2.);
  else if(dir==X3DN || dir==X3UP){
    dualfprintf(fail_file,"Shouldn't be here for X3 dir\n");
    myexit(87253245);
  }

  // densities
  pr[RHO] = 1;
  pr[UU] = 0.1;


  if(src){

    FTYPE X[NDIM],V[NDIM]; 
    int jj,kk;
    struct of_geom geomdontuse[NPR];
    struct of_geom *ptrgeom[NPR];
    // assign memory
    PALLLOOP(pl) ptrgeom[pl]=&(geomdontuse[pl]);

    // local geom
    PALLLOOP(pl) get_geometry(i, j, k, dirprim[pl], ptrgeom[pl]);
   
    //coordinates of the ghost cell
    bl_coord_ijk_2(i,j,k,CENT,X, V);
    FTYPE xx=V[1];
    FTYPE yy=V[2];
    FTYPE zz=V[3];

    //Keplerian gas
    extern FTYPE RADNT_OMSCALE;
    FTYPE Omx=-RADNT_OMSCALE/(a+pow(fabs(yy),1.5))*yy;
    FTYPE Omy=RADNT_OMSCALE/(a+pow(fabs(xx),1.5))*xx;

    //    dualfprintf(fail_file,"dir=%d RADNT_OMSCALE=%g xx=%g a=%g Omx=%g Omy=%g\n",dir,RADNT_OMSCALE,xx,a,Omx,Omy);

    if(dir==X1DN || dir==X1UP){
      pr[U1] = 0.0;
      pr[U2] = Omy;
    }
    else if(dir==X2DN || dir==X2UP){
      pr[U1] = Omx;
      pr[U2] = 0.0;
    }

    pr[U3] = 0.0;

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
  }
  else{
    pr[U1] = 0.0;
    pr[U2] = 0.0;
    pr[U3] = 0.0;

    pr[URAD0] = 0.001;
    pr[URAD1] = 0.0;
    pr[URAD2] = 0.0;
    pr[URAD3] = 0.0;
  }

  return(0);
}






/// X1 upper and lower static
int bound_staticset(int dir,
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


    if((BCtype[X1DN]==OUTFLOWSTATIC || BCtype[X1DN]==HORIZONOUTFLOWSTATIC) && (totalsize[1]>1) && (mycpupos[1] == 0) ){


      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
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
          //          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          if(ispstag==0){

            //initially copying everything
            PBOUNDLOOP(pliter,pl) if(pl==U1) MACP0A1(prim,i,j,k,pl) = 0.0; // static

          }// end if not staggered fields

        }// end loop over outer i's

      }// end over loop
    }// end if correct boundary condition and core



    if((BCtype[X1UP]==OUTFLOWSTATIC || BCtype[X1UP]==HORIZONOUTFLOWSTATIC) && (totalsize[1]>1) && (mycpupos[1] == ncpux1-1) ){


      OPENMPBCLOOPVARSDEFINELOOPX1DIR; OPENMPBCLOOPSETUPLOOPX1DIR;
      //////// LOOPX1dir{
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
          //          PBOUNDLOOP(pliter,pl) MACP0A1(prim,i,j,k,pl) = MACP0A1(prim,ri,rj,rk,pl);

          if(ispstag==0){

            //initially copying everything
            PBOUNDLOOP(pliter,pl) if(pl==U1) MACP0A1(prim,i,j,k,pl) = 0.0; // static

          }// end if not staggered fields

        }// end loop over outer i's

      }// end over loop
    }// end if correct boundary condition and core



   
   
   
  }// end parallel region
  


 



  return(0);
} 





/// WALDMONO
int bound_waldmono(int dir,
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


  


    if(dir==X2UP && BCtype[X2UP]==WALDMONOBC && (totalsize[2]>1) && (mycpupos[2] == ncpux2-1) ){

      OPENMPBCLOOPVARSDEFINELOOPX2DIR; OPENMPBCLOOPSETUPLOOPX2DIR;
      //////// LOOPX2dir{
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

            FTYPE rin,rout;
            int conddisk;
            rin=15.;
            rout=25.;
            conddisk=1;//(r<rout); // as in new koral

            //hot boundary
            if(conddisk){

              //pr[RHO] and pr[UU] remain same as from ASYMM condition as well as any field
              //Keplerian gas with no inflow or outflow
              pr[U1]=pr[U2]=0.0; // have to be careful with this for VEL3 (must have rin>>rergo).
              pr[U3]=0.0; // current koral value
 
            } // end if actually doing something to boundary cells in "hot" boundary

          }// end if not staggered field


        }// end loop over outer j's
      }
    }







  }// end parallel region

  return(0);
} 

