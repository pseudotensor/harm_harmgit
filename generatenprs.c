#include <stdio.h>
#include "global.realdef.h"
#include "global.nondepmnemonics.h"
#include "definit.h"
#include "init.h"


/* \file generatenprs.c
   \brief Generates various code files needed for HARM code
 */


// upper limit of number of MAXNPR
#define SUPERMAXNPR 20

int main(void)
{
  int pl,pliter;

  // temp vars
  int mydoyfl,mydoyl,mydoynu,mydoentropy,mydoextrainterp;

  // things to define:
  int npr,npr2interp,npr2notinterp,nprbound,nfluxbound,nprdump,nprinvert;
  int maxnpr;
  int yfl1,yfl2,yfl3,yfl4,yfl5,yl,ynu,entropy,vsq;
  int rad0,rad1,rad2,rad3,nrad;
  int orignprstart,orignprend,orignprlist[SUPERMAXNPR];
  int orignpr2interpstart,orignpr2interpend,orignpr2interplist[SUPERMAXNPR];
  int orignpr2notinterpstart,orignpr2notinterpend,orignpr2notinterplist[SUPERMAXNPR];
  int orignprboundstart,orignprboundend,orignprboundlist[SUPERMAXNPR];
  int orignprfluxboundstart,orignprfluxboundend,orignprfluxboundlist[SUPERMAXNPR];
  int orignprdumpstart,orignprdumpend,orignprdumplist[SUPERMAXNPR];
  int orignprinvertstart,orignprinvertend,orignprinvertlist[SUPERMAXNPR];
  int i;

  // file output
  FILE * defout;




  // npr is for list of all possible variables while nprlist holds loop for general PLOOP loop

  /////////
  //
  // these lists are like those in rest of code that assign variable names to list
  //
  /////////
  npr=1; orignprstart=0; orignprend=0; orignprlist[orignprend]=RHO;
  npr++; orignprend++; orignprlist[orignprend]=UU;
  npr++; orignprend++; orignprlist[orignprend]=U1;
  npr++; orignprend++; orignprlist[orignprend]=U2;
  npr++; orignprend++; orignprlist[orignprend]=U3;


  // inversion init and final: inversion only for original 5 things
  orignprinvertstart=orignprstart; orignprinvertend=orignprend;
  for(i=orignprinvertstart;i<=orignprinvertend;i++){
    orignprinvertlist[i]=orignprlist[i];
  }

  //////////////
  //
  //
  //////////////
  if(EOMTYPE>=0){
    npr++; orignprend++; orignprlist[orignprend]=B1;
    npr++; orignprend++; orignprlist[orignprend]=B2;
    npr++; orignprend++; orignprlist[orignprend]=B3;
  }

  // 2interp init
  orignpr2interpstart=orignprstart; orignpr2interpend=orignprend;
  for(i=orignpr2interpstart;i<=orignpr2interpend;i++){
    orignpr2interplist[i]=orignprlist[i];
  }

  // 2notinterp init (default is nothing to not interpolate)
  // 2notinterp is actually those things that would normally interpolate but don't for some new reason (so has nothing to do with those quantities we never interpolate)
  orignpr2notinterpstart=0; orignpr2notinterpend=-1;

  // bound init
  orignprboundstart=orignprstart; orignprboundend=orignprend;
  for(i=orignprboundstart;i<=orignprboundend;i++){
    orignprboundlist[i]=orignprlist[i];
  }

  // fluxbound init
  orignprfluxboundstart=orignprstart; orignprfluxboundend=orignprend;
  for(i=orignprfluxboundstart;i<=orignprfluxboundend;i++){
    orignprfluxboundlist[i]=orignprlist[i];
  }

  // dump init
  orignprdumpstart=orignprstart; orignprdumpend=orignprend;
  for(i=orignprdumpstart;i<=orignprdumpend;i++){
    orignprdumplist[i]=orignprlist[i];
  }



  ///////////////////
  //
  // Radiation things
  //
  ////////////////////
  if(EOMRADTYPE!=EOMRADNONE){
    npr++; rad0 = npr-1;
    orignprend++; orignprlist[orignprend]=rad0;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=rad0;
    orignprboundend++; orignprboundlist[orignprboundend]=rad0;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=rad0;
    orignprdumpend++; orignprdumplist[orignprdumpend]=rad0;

    npr++; rad1 = npr-1;
    orignprend++; orignprlist[orignprend]=rad1;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=rad1;
    orignprboundend++; orignprboundlist[orignprboundend]=rad1;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=rad1;
    orignprdumpend++; orignprdumplist[orignprdumpend]=rad1;

    npr++; rad2 = npr-1;
    orignprend++; orignprlist[orignprend]=rad2;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=rad2;
    orignprboundend++; orignprboundlist[orignprboundend]=rad2;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=rad2;
    orignprdumpend++; orignprdumplist[orignprdumpend]=rad2;

    npr++; rad3 = npr-1;
    orignprend++; orignprlist[orignprend]=rad3;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=rad3;
    orignprboundend++; orignprboundlist[orignprboundend]=rad3;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=rad3;
    orignprdumpend++; orignprdumplist[orignprdumpend]=rad3;
  }
  else{
    rad0 = VARNOTDEFINED; // indicates not defined
    rad1 = VARNOTDEFINED; // indicates not defined
    rad2 = VARNOTDEFINED; // indicates not defined
    rad3 = VARNOTDEFINED; // indicates not defined
  }

  // number density of radiation
  if(EVOLVENRAD){
    npr++; nrad = npr-1;
    orignprend++; orignprlist[orignprend]=nrad;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=nrad;
    orignprboundend++; orignprboundlist[orignprboundend]=nrad;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=nrad;
    orignprdumpend++; orignprdumplist[orignprdumpend]=nrad;
  }
  else{
    nrad = VARNOTDEFINED; // indicates not defined
  }

  ///////////////////
  //
  // Rest of variables are most often on/off
  //
  // These variables: yflx, yl, ynu, entropy, all have variable number assignments unlike more stable names
  //
  ////////////////////

  if(DOYFL!=DONOYFL){
    npr++; yfl1 = npr-1;

    orignprend++; orignprlist[orignprend]=yfl1;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=yfl1;
    orignprboundend++; orignprboundlist[orignprboundend]=yfl1;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=yfl1;
    orignprdumpend++; orignprdumplist[orignprdumpend]=yfl1;
  }
  else{
    yfl1 = VARNOTDEFINED; // indicates not defined
  }

  if(DOYFL!=DONOYFL){
    npr++; yfl2 = npr-1;

    orignprend++; orignprlist[orignprend]=yfl2;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=yfl2;
    orignprboundend++; orignprboundlist[orignprboundend]=yfl2;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=yfl2;
    orignprdumpend++; orignprdumplist[orignprdumpend]=yfl2;
  }
  else{
    yfl2 = VARNOTDEFINED; // indicates not defined
  }

  if(DOYFL!=DONOYFL){
    npr++; yfl3 = npr-1;

    orignprend++; orignprlist[orignprend]=yfl3;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=yfl3;
    orignprboundend++; orignprboundlist[orignprboundend]=yfl3;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=yfl3;
    orignprdumpend++; orignprdumplist[orignprdumpend]=yfl3;
  }
  else{
    yfl3 = VARNOTDEFINED; // indicates not defined
  }

  if(DOYFL!=DONOYFL && (EOMRADTYPE!=EOMRADNONE)){
    npr++; yfl4 = npr-1;

    orignprend++; orignprlist[orignprend]=yfl4;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=yfl4;
    orignprboundend++; orignprboundlist[orignprboundend]=yfl4;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=yfl4;
    orignprdumpend++; orignprdumplist[orignprdumpend]=yfl4;
  }
  else{
    yfl4 = VARNOTDEFINED; // indicates not defined
  }

  if(DOYFL!=DONOYFL && (EOMRADTYPE!=EOMRADNONE)){
    npr++; yfl5 = npr-1;

    orignprend++; orignprlist[orignprend]=yfl5;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=yfl5;
    orignprboundend++; orignprboundlist[orignprboundend]=yfl5;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=yfl5;
    orignprdumpend++; orignprdumplist[orignprdumpend]=yfl5;
  }
  else{
    yfl5 = VARNOTDEFINED; // indicates not defined
  }




  if(DOYL!=DONOYL){
    npr++; yl = npr-1;

    orignprend++; orignprlist[orignprend]=yl;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=yl;
    orignprboundend++; orignprboundlist[orignprboundend]=yl;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=yl;
    orignprdumpend++; orignprdumplist[orignprdumpend]=yl;
  }
  else{
    yl = VARNOTDEFINED; // indicates not defined
  }


  if(DOYNU!=DONOYNU){
    npr++; ynu = npr-1;

    orignprend++; orignprlist[orignprend]=ynu;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=ynu;
    orignprboundend++; orignprboundlist[orignprboundend]=ynu;
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=ynu;
    orignprdumpend++; orignprdumplist[orignprdumpend]=ynu;
  }
  else{
    ynu = VARNOTDEFINED;
  }

  if(DOENTROPY!=DONOENTROPY){ // in reality primitive for entropy is not used -- don't have to interpolate -- how to handle unless have label for each type of object GODMARK
    // don't care about entropy for primitive bounding since entropy just simple function of primitives and don't use entropy in ghost zones

    npr++; entropy = npr-1;

    orignprend++; orignprlist[orignprend]=entropy; // do treat conserved version of entropy
    //    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=entropy; // do not interpolate primitive version of entropy

    //    orignprboundend++; orignprboundlist[orignprboundend]=entropy; // don't bound primitive version of entropy
    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=entropy; // do bound flux of entropy (related to conserved quantity)
    //    orignprdumpend++; orignprdumplist[orignprdumpend]=entropy; // don't dump primitive version of entropy
  }
  else{
    entropy = VARNOTDEFINED;
  }



  // v^2 (or \gamma of some type) is not related to any conserved quantity and does not need to be bounded
  // adds only to things to interpolate
  if(DOEXTRAINTERP){
    npr++; vsq = npr-1;

    //    orignprend++; orignprlist[orignprend]=vsq;
    orignpr2interpend++; orignpr2interplist[orignpr2interpend]=vsq; // only for interpolation
    //    orignprboundend++; orignprboundlist[orignprboundend]=vsq;
    //    orignprfluxboundend++; orignprfluxboundlist[orignprfluxboundend]=vsq;
    //    orignprdumpend++; orignprdumplist[orignprdumpend]=vsq;
  }
  else{
    vsq = VARNOTDEFINED;
  }


  // Now npr is the number of all possible quantities
  npr2interp=orignpr2interpend+1;    // number of interpolated quantities (independent from actual list of primitives)

  // for simplicity of pointer assignments in flux.c (p2interp=pr) and fluxctstag.c (pstag=pr),
  // set "npr2interp" to be no smaller than npr.  Loop control will still avoid unwanted values
  // GODMARK: alternatively could define "memory" version of macro definition for pointer memory so always know how many interpolating things total we are supposed to have
  if(npr2interp<npr) npr2interp=npr;

  nprbound = orignprboundend+1;      // default # of prim quants to bound
  nfluxbound=orignprfluxboundend+1;  // must be equal to NPR for BOUNDFLUXTYPE // default # of flux quants to bound
  nprdump = orignprdumpend+1;           // default # of quants to dump
  nprinvert=orignprinvertend+1;      // number of quantities to invert in 5D inversion (5 always for now)


  // determine maximum over all NPR types
  maxnpr=MAX(npr,nprbound);
  maxnpr=MAX(maxnpr,nprdump);
  maxnpr=MAX(maxnpr,nfluxbound);
  maxnpr=MAX(maxnpr,npr2interp);
  maxnpr=MAX(maxnpr,nprinvert);


  if(npr>SUPERMAXNPR){
    fprintf(stderr,"Must increase SUPERMAXNPR=%d to be greater than npr=%d\n",SUPERMAXNPR,npr);
    exit(1);
  }

  ///////////////////////////////////////////
  //
  // Now generate #define's for HARM
  //
  /////////////////////////////////////////////

  defout=fopen("global.defnprs.h","wt");
  if(defout==NULL){
    fprintf(stderr,"Can't open global.defnprs.h\n");
    return(1);
  }

  // define number of variables
  fprintf(defout,"// THIS IS AN AUTOGENERATED FILE, DO NOT MODIFY\n\n\n\n");
  fprintf(defout,"#define MAXNPR %d\n",maxnpr);
  fprintf(defout,"#define NPR %d\n",npr);
  fprintf(defout,"#define NPRREALSET %d\n",npr); // set of "real" quantities // GODMARK: excessive to do npr, but otherwise need to choose carefully what system of equations being solved.
  fprintf(defout,"#define NPRBOUND %d\n",nprbound);
  fprintf(defout,"#define NPRDUMP %d\n",nprdump);
  fprintf(defout,"#define NFLUXBOUND %d\n",nfluxbound);
  fprintf(defout,"#define NPR2INTERP %d\n",npr2interp);
  fprintf(defout,"#define NPRINVERT %d\n",nprinvert);

  //SUPERGODMARK
  // define name of radiation terms for general conserved quantities
  fprintf(defout,"#define PRAD0 %d\n",rad0);
  fprintf(defout,"#define PRAD1 %d\n",rad1);
  fprintf(defout,"#define PRAD2 %d\n",rad2);
  fprintf(defout,"#define PRAD3 %d\n",rad3);

  // define name of radiation terms for primitives
  fprintf(defout,"#define URAD0 %d\n",rad0);
  fprintf(defout,"#define URAD1 %d\n",rad1);
  fprintf(defout,"#define URAD2 %d\n",rad2);
  fprintf(defout,"#define URAD3 %d\n",rad3);

  fprintf(defout,"#define NRAD %d\n",nrad);

  // define name of extra variables
  fprintf(defout,"#define YFL1 %d\n",yfl1);
  fprintf(defout,"#define YFL2 %d\n",yfl2);
  fprintf(defout,"#define YFL3 %d\n",yfl3);
  fprintf(defout,"#define YFL4 %d\n",yfl4);
  fprintf(defout,"#define YFL5 %d\n",yfl5);
  fprintf(defout,"#define YL %d\n",yl);
  fprintf(defout,"#define YE %d\n",yl); // treating YE as primitive for YL conservation
  fprintf(defout,"#define YNU %d\n",ynu);
  fprintf(defout,"#define ENTROPY %d\n",entropy);
  fprintf(defout,"#define VSQ %d\n",vsq);

  fclose(defout);



  ///////////////////////////////////////////
  //
  // Now generate default assignments for npr type lists
  //
  // this generates a function added to initbase.c
  //
  /////////////////////////////////////////////

  defout=fopen("initbase.defaultnprlists.c","wt");
  if(defout==NULL){
    fprintf(stderr,"Can't open initbase.defaultnprlists.c\n");
    return(1);
  }

  // define number of variables
  fprintf(defout,"\n\n\n\n// THIS IS AN AUTOGENERATED FILE, DO NOT MODIFY\n\n\n\n");
  fprintf(defout,"void set_default_nprlists(void)\n{\n\nint pl,pliter;\n\n\n");
  
  // npr
  fprintf(defout,"\n%s=%d; %s=%d;\n","nprstart",orignprstart,"nprend",orignprend);
  for(pl=orignprstart;pl<=orignprend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","nprlist",pl,orignprlist[pl]);
  }
  
  // 2interp
  fprintf(defout,"\n#pragma omp parallel\n{// must set npr2interp stuff inside parallel region since threadprivate\n");
  fprintf(defout,"\n%s=%d; %s=%d;\n","npr2interpstart",orignpr2interpstart,"npr2interpend",orignpr2interpend);
  for(pl=orignpr2interpstart;pl<=orignpr2interpend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","npr2interplist",pl,orignpr2interplist[pl]);
  }
  fprintf(defout,"\n}\n");

  // 2notinterp (should always just give default -1 to 0 range)
  fprintf(defout,"\n#pragma omp parallel\n{// must set npr2interp stuff inside parallel region since threadprivate\n");
  fprintf(defout,"\n%s=%d; %s=%d;\n","npr2notinterpstart",orignpr2notinterpstart,"npr2notinterpend",orignpr2notinterpend);
  for(pl=orignpr2notinterpstart;pl<=orignpr2notinterpend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","npr2notinterplist",pl,orignpr2notinterplist[pl]);
  }
  fprintf(defout,"\n}\n");

  // bound
  fprintf(defout,"\n%s=%d; %s=%d;\n","nprboundstart",orignprboundstart,"nprboundend",orignprboundend);
  for(pl=orignprboundstart;pl<=orignprboundend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","nprboundlist",pl,orignprboundlist[pl]);
  }

  // fluxbound
  fprintf(defout,"\n%s=%d; %s=%d;\n","nprfluxboundstart",orignprfluxboundstart,"nprfluxboundend",orignprfluxboundend);
  for(pl=orignprfluxboundstart;pl<=orignprfluxboundend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","nprfluxboundlist",pl,orignprfluxboundlist[pl]);
  }

  // dump
  fprintf(defout,"\n%s=%d; %s=%d;\n","nprdumpstart",orignprdumpstart,"nprdumpend",orignprdumpend);
  for(pl=orignprdumpstart;pl<=orignprdumpend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","nprdumplist",pl,orignprdumplist[pl]);
  }

  // invert 
  fprintf(defout,"\n%s=%d; %s=%d;\n","nprinvertstart",orignprinvertstart,"nprinvertend",orignprinvertend);
  for(pl=orignprinvertstart;pl<=orignprinvertend;pl++){
    fprintf(defout,"%s[%d]=%d;\n","nprinvertlist",pl,orignprinvertlist[pl]);
  }

  fprintf(defout,"\n\n\n}\n\n\n");

  // old code from prepre_init():
  //
  // below used both if SPLITNPR=0,1
  // choice for range of PLOOP
  //  nprstart=0;
  //  nprend=NPR-1;
  //  for(pl=nprstart;pl<=nprend;pl++) nprlist[pl]=pl;

  // choice for range of PLOOPINTERP
  //  npr2interpstart=0;
  //  npr2interpend=NPR2INTERP-1;
  //  for(pl=npr2interpstart;pl<=npr2interpend;pl++) npr2interplist[pl]=pl;

  // choice for range of PLOOPNOTINTERP
  //  npr2notinterpstart=0;
  //  npr2notinterpend=-1;
  //  npr2notinterplist[0]=0;

  // default choice for range of PBOUNDLOOP and PLOOPMPI
  //  nprboundstart=0;
  //  nprboundend=NPRBOUND-1;
  //  for(pl=nprboundstart;pl<=nprboundend;pl++) nprboundlist[pl]=pl;



  fclose(defout);



  return(0);


}
