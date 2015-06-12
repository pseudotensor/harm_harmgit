
/*! \file image.c
     \brief image dumping

   produces an "r8" file.
*/

#include "decs.h"

#define GAMMIESTARTI (0)
#define GAMMIEENDI (0)

#define JONSTARTI (0)
#define JONENDI (NPRDUMP-1)

#define NORMAL (0)
#define ZOOM (1)
// revert to only normal for now
//#define NUMLIMITS (2)
#define NUMLIMITS (1)

#define LOG (0)
#define LINEAR (1)
#define NUMSCALE (2)

#define DOCONS (1)
// whether to do (0=no, 1=yes) conserved quantities in image

#define MINVECTOR (SMALL)


// OPENMPMARK: Assume imagedefs() image() only called by master thread
// global variables for this file
static FTYPE imageparms[NUMIMAGEPARMS];
static int imagescale,imagewhichpl,imagevartype,imagelimits;



/// setup image dump
int image_dump(long dump_cnt)
{
  int startpl,endpl,starts,ends,startl,endl,startv,endv;
  int whichpl,limits,scale,vartype;

  if(DOIMAGEDUMP==0) return(0);

  ////////////////////////////
  //
  // Image Loop
  //
  ////////////////////////////

  // STARTI=0 is normal
  // ENDI=NPRDUMP is normal, but can be NPRDUMP+2 to get linear RHO/UU
  if(PRODUCTION==0){
    if(GAMMIEIMAGE){
      startpl=GAMMIESTARTI;
      endpl=GAMMIEENDI;
      starts=0;
      ends=1;
      startl=0;
      endl=0;
      startv=0;
      endv=0;
    }
    else{
      startpl=JONSTARTI;
      endpl=JONENDI;
      starts=0;
      ends=NUMSCALE-1;
      startl=0;
      endl=NUMLIMITS-1;
      startv=0;
      //    endv=1;
      endv=2; // including failures now
    }
  }
  else{
    // no special images for production mode, just basic log density
    if(GAMMIEIMAGE){
      if(DOEVOLVERHO){
        startpl=0;
        endpl=0;
      }
      else{
        startpl=B1;
        endpl=B1;
      }
      starts=0;
      ends=0;
      startl=0;
      endl=0;
      startv=0;
      endv=0;
    }
    else{
      if(DOEVOLVERHO){
        startpl=0;
        endpl=0;
      }
      else{
        startpl=B1;
        endpl=B1;
      }
      starts=0;
      ends=0;
      startl=0;
      endl=0;
      startv=0;
      endv=0;
    }
  }



  for(vartype=startv;vartype<=endv;vartype++){
    for(limits=startl;limits<=endl;limits++){
      for(scale=starts;scale<=ends;scale++){
        for (whichpl = startpl; whichpl <= endpl; whichpl++) {
          if(
             ((vartype<=1)||((vartype==2)&&(DODEBUG)&&(whichpl<=NUMFAILFLOORFLAGS-1)))||
             ((limits==0)&&(scale==LINEAR)&&(vartype==2)&&(whichpl==4)&&DO_VORTICITY_IMAGE)
             ){
            if(image(dump_cnt,whichpl,scale,limits,vartype)>=1) return(1);
          }
        }
      }
    }
  }

  return(0);
}

/// define image parameters and data range
int imagedefs(int whichpl, int scale, int limits, int vartype)
{
  int i = 0, j = 0, k = 0, l = 0, col = 0, floor;
  FILE *fp;
  // whichpl : whichpl primitive variable
  FTYPE pr,iq, liq, lmax;
  unsigned char liqb;
  FTYPE min,max,sum;
  FTYPE minptr[NPR], maxptr[NPR], sumptr[NPR];
  char truemyidtxt[MAXFILENAME];
  FTYPE U[NPR];
  struct of_state q;
  FTYPE X[NDIM],V[NDIM],r,th;
  FTYPE lmin,aa;
  int compute_vorticity(FTYPE (*p)[NSTORE2][NSTORE3][NPR],FTYPE (*pvort)[NSTORE2][NSTORE3][NPR],int whichpl);
  int pl,pliter;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"imagedefs() called in parallel region\n");
    myexit(3487346);
  }
#endif


  ////////////////////////////
  //
  // Image output setup/definition
  //
  // Purpose is to set pimage to correct variable type (primitive or conservative), limits, scale, and which k.
  // Then image() outputs that one thing to file
  //
  ////////////////////////////

  GLOBALPOINT(pimage)=GLOBALPOINT(dUgeomarray); // assume dUgeomarray is only used within each substep and not across substeps
  if(limits==ZOOM){ // zoom in on dynamic range of values to see fine details
    DUMPGENLOOP{ // diagnostic loop // OPENMPOPTMARK: Could optimize this, but not frequently done
      if(vartype==0){
        if(whichpl==RHO || whichpl==UU || pl==PRAD0 || SCALARPL(pl)){
          bl_coord_ijk_2(i,j,k,CENT,X,V);
          r=V[1];
          th=V[2];
          if(whichpl==0) GLOBALMACP0A1(pimage,i,j,k,whichpl)=GLOBALMACP0A1(pdump,i,j,k,whichpl)/(RHOMIN*pow(r,-1.5));
          if(whichpl==1) GLOBALMACP0A1(pimage,i,j,k,whichpl)=GLOBALMACP0A1(pdump,i,j,k,whichpl)/(UUMIN*pow(r,-2.5));
        }
        else{
          if(scale==LINEAR) GLOBALMACP0A1(pimage,i,j,k,whichpl)=GLOBALMACP0A1(pdump,i,j,k,whichpl);
          else if(scale==LOG) GLOBALMACP0A1(pimage,i,j,k,whichpl)=fabs(GLOBALMACP0A1(pdump,i,j,k,whichpl))+MINVECTOR;
        }
      }
      else if(vartype==1){// conserved quantity
        // computes too much (all conserved quantites every time)
        if(DOENOFLUX == NOENOFLUX){
          get_geometry(i,j,k,CENT,ptrgeom) ;
          if(!failed){
            if(get_state(GLOBALMAC(pdump,i,j,k),ptrgeom,&q)>=1) return(1);
            if(primtoU(UDIAG,GLOBALMAC(pdump,i,j,k),&q,ptrgeom,U, NULL)>=1) return(1);
          }
        }
        else{
          PALLLOOP(pl) U[pl]=GLOBALMACP0A1(udump,i,j,k,pl);
        }
        if(scale==LINEAR) GLOBALMACP0A1(pimage,i,j,k,whichpl)=U[whichpl];
        else if(scale==LOG) GLOBALMACP0A1(pimage,i,j,k,whichpl)=fabs(U[whichpl]/ptrgeom->gdet)+MINVECTOR;
      }
      else if(vartype==2){ // failure quantity (no diff from below right now -- could zoom in on single failure regions)
        if(whichpl<NUMFAILFLOORFLAGS){
          floor=whichpl;
          if(scale==LINEAR) GLOBALMACP0A1(pimage,i,j,k,whichpl)=(FTYPE)GLOBALMACP0A3(failfloorcount,i,j,k,0,IMAGETS,floor); // always finalstep==0
          else if(scale==LOG) GLOBALMACP0A1(pimage,i,j,k,whichpl)=fabs((FTYPE)GLOBALMACP0A3(failfloorcount,i,j,k,0,IMAGETS,floor)+1); // always finalstep==0
        }
      }
    }
  }
  else{ // normal image
    DUMPGENLOOP{ // diagnostic loop // OPENMPOPTMARK: Could optimize this, but not frequently done
      if(vartype==0){
        if(whichpl==RHO || whichpl==UU || pl==PRAD0 || SCALARPL(pl)) GLOBALMACP0A1(pimage,i,j,k,whichpl)=GLOBALMACP0A1(pdump,i,j,k,whichpl);
        else{
          if(scale==LINEAR) GLOBALMACP0A1(pimage,i,j,k,whichpl)=GLOBALMACP0A1(pdump,i,j,k,whichpl);
          else if(scale==LOG) GLOBALMACP0A1(pimage,i,j,k,whichpl)=fabs(GLOBALMACP0A1(pdump,i,j,k,whichpl))+MINVECTOR;
        }
      }
      else if(vartype==1){// conserved quantity
        // computes too much (all conserved quantites every time)
        if(DOENOFLUX == NOENOFLUX){
          get_geometry(i,j,k,CENT,ptrgeom) ;
          if(!failed){
            if(get_state(GLOBALMAC(pdump,i,j,k),ptrgeom,&q)>=1) return(1);
            if(primtoU(UDIAG,GLOBALMAC(pdump,i,j,k),&q,ptrgeom,U, NULL)>=1) return(1);
          }
        }
        else{
          PALLLOOP(pl) U[pl]=GLOBALMACP0A1(udump,i,j,k,pl);
        }
        if(scale==LINEAR) GLOBALMACP0A1(pimage,i,j,k,whichpl)=U[whichpl];
        else if(scale==LOG) GLOBALMACP0A1(pimage,i,j,k,whichpl)=fabs(U[whichpl]/ptrgeom->gdet)+MINVECTOR;
      }
      else if(vartype==2){ // failure quantity
        if(whichpl<NUMFAILFLOORFLAGS){
          floor=whichpl;
          if(scale==LINEAR) GLOBALMACP0A1(pimage,i,j,k,whichpl)=(FTYPE)GLOBALMACP0A3(failfloorcount,i,j,k,0,IMAGETS,floor); // finalstep==0
          else if(scale==LOG) GLOBALMACP0A1(pimage,i,j,k,whichpl)=fabs((FTYPE)GLOBALMACP0A3(failfloorcount,i,j,k,0,IMAGETS,floor)+1); // finalstep==0
        }
      }
    }
#if(DO_VORTICITY_IMAGE)
    if(vartype==2){
      // overwrite vartype==2 with vorticity for whichpl==4
      // e.g. !dr82 images/im4f1s0l0100.r8 for final vorticity
      compute_vorticity(pglobal,pimage,4);
    }
#endif
  }



  ////////////////////////////
  //
  // Image FILE open/initialize
  //
  ////////////////////////////



  ////////////////////
  //
  // Image paramters setup (whole purpose currently is to find lmin and aa)
  // 
  /////////////////////

  /* density mapping is logarithmic, in 255 steps between e^lmax and
     e^lmin */

#define ZOOMFACTOR (10000)

  prminmaxsum(GLOBALPOINT(pimage),whichpl,1,maxptr,minptr,sumptr);
  if(limits==NORMAL){
    max=maxptr[whichpl];
    min=minptr[whichpl];
    sum=sumptr[whichpl];
#ifdef TESTNUMBER 
#if( DO_VORTICITY_IMAGE && (TESTNUMBER == 26 || TESTNUMBER == 27) )
    if( vartype==2 && whichpl==4 && limits == 0 && scale == LINEAR ) {
      //The exact vorticity should be changing between -5 & 10, so make the limits 
      //a bit larger than that to avoid "black" spots on the image
      if( debugfail >= 1)  trifprintf( "Vorticity: old: min = %lg, max = %lg", min, max );
      min = -5.5 / coordparams.timescalefactor;
      max = 10.5 / coordparams.timescalefactor;
      if( debugfail >= 1)  trifprintf( ", new: min = %lg, max = %lg\n", min, max );
    }
#endif
#endif
  }
  else{
    if(whichpl==RHO || whichpl==UU || pl==PRAD0 || SCALARPL(pl)){
      max=maxptr[whichpl]/ZOOMFACTOR;
      min=minptr[whichpl];
    }
    else{
      if(scale==LINEAR){
        max=maxptr[whichpl]/ZOOMFACTOR;
        min=minptr[whichpl]/ZOOMFACTOR;
      }
      else{
        max=maxptr[whichpl]/ZOOMFACTOR;
        min=minptr[whichpl];
      }
    }
  }
  sum=sumptr[whichpl];
  logsfprintf("whichpl: %d scale: %d limits: %d : min,max,avg: %21.15g %21.15g %21.15g\n",whichpl,scale,limits,min,max,sum/totalzones);

  if(scale==LOG){
    lmax = log(max);
    lmin = log(min);
  } else if(scale==LINEAR) {
    lmax = max;
    lmin = min;
  }
  else{
    dualfprintf(fail_file,"no such scale=%d\n",scale);
    myexit(1);
  }

  if (lmax != lmin)
    aa = 256. / (lmax - lmin);
  else
    aa = 0;



  // set the only paramters needed to dump image
  imageparms[ORIGIN]=aa;
  imageparms[LMIN]=lmin;
  // extra:
  imageparms[LMAX]=lmax;

  return(0);
}



/// create image dump
int image(long dump_cnt, int whichpl, int scale, int limits, int vartype)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


#if(USEOPENMP)
  if(omp_in_parallel()){
    dualfprintf(fail_file,"image() called in parallel region\n");
    myexit(3487347);
  }
#endif


  trifprintf("begin dumping image# %ld whichpl: %d scale: %d limits: %d vartype: %d\n ",dump_cnt,whichpl,scale,limits,vartype);

  // global vars so can avoid passing through general functions
  imagescale=scale;
  imagevartype=vartype;
  imagelimits=limits;
  imagewhichpl=whichpl;



  whichdump=IMAGEDUMPTYPE;
  datatype=MPI_UNSIGNED_CHAR;

  // actual prefix
  if(GAMMIEIMAGE&&(GAMMIESTARTI==GAMMIEENDI)&&(GAMMIESTARTI==0)){
    sprintf(fileprefix, "images/im");
  }
  else{
    if(vartype==0) sprintf(fileprefix, "images/im%1dp%1ds%1dl", whichpl, scale, limits);
    else if(vartype==1)  sprintf(fileprefix, "images/im%1dc%1ds%1dl", whichpl, scale, limits);
    else if(vartype==2)  sprintf(fileprefix, "images/im%1df%1ds%1dl", whichpl, scale, limits);
  }
  strcpy(fileformat,"%04ld"); // actual format
  strcpy(filesuffix,".r8"); // actual suffix

  // setup the image definitions (min,max, what data, etc.)
  if(imagedefs(whichpl,scale,limits,vartype)>=1) return(1);


  // MIXEDOUTPUT tells dump_gen() to treat as .r8 file, with text header and binary dump
  if(dump_gen(WRITEFILE,dump_cnt,MIXEDOUTPUT,whichdump,datatype,fileprefix,fileformat,filesuffix,image_header,image_content)>=1) return(1);

  trifprintf("end dumping image# %ld whichpl: %d scale: %d limits: %d vartype: %d\n ",dump_cnt,whichpl,scale,limits,vartype);

  return(0);

}

/// image header
int image_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{ 
  int realtotalsize[NDIM];


  realtotalsize[1]=totalsize[1]+2*EXTRADUMP1;
  realtotalsize[2]=totalsize[2]+2*EXTRADUMP2;
  realtotalsize[3]=totalsize[3]+2*EXTRADUMP3;
  
  ////////////////////////////
  //
  // HEADER file open/initialize
  //
  ////////////////////////////

  // write header
  if(bintxt==TEXTOUTPUT){

    fprintf(headerptr, "RAW\n# t=%21.15g nstep=%ld vartype=%2d  whichpl=%2d scale=%2d limits=%2d ",t, nstep, imagevartype, imagewhichpl, imagescale, imagelimits);
    fprintf(headerptr, "aa=%g lmin=%g lmax=%g ",imageparms[ORIGIN],imageparms[LMIN],imageparms[LMAX]);
    fprintf(headerptr,"\n");
    
    if( (realtotalsize[1]>1)&&(realtotalsize[2]>1)&&(realtotalsize[3]>1)){
      fprintf(headerptr, "%i %i\n255\n", realtotalsize[1],realtotalsize[2]*realtotalsize[3]);
    }
    else{
      if( (realtotalsize[2]>1)&&(realtotalsize[3]>1)){
        fprintf(headerptr, "%i %i\n255\n", realtotalsize[2],realtotalsize[3]);
      }
      else if( (realtotalsize[1]>1)&&(realtotalsize[3]>1)){
        fprintf(headerptr, "%i %i\n255\n", realtotalsize[1],realtotalsize[3]);
      }
      else if( (realtotalsize[1]>1)&&(realtotalsize[2]>1)){
        fprintf(headerptr, "%i %i\n255\n", realtotalsize[1],realtotalsize[2]);
      }
      else{
        if(realtotalsize[1]>1){
          fprintf(headerptr, "%i %i\n255\n", realtotalsize[1],1);
        }
        else if(realtotalsize[2]>1){
          fprintf(headerptr, "%i %i\n255\n", realtotalsize[2],1);
        }
        else if(realtotalsize[3]>1){
          fprintf(headerptr, "%i %i\n255\n", realtotalsize[3],1);
        }
        else{
          dualfprintf(fail_file,"Shouldn't reach here: %d %d %d\n",realtotalsize[1],realtotalsize[2],realtotalsize[3]);
          myexit(248742);
        }
      }
    }
  }
  else{
    dualfprintf(fail_file,"Shouldn't be trying to write binary header to image file\n");
    myexit(1);
  }
  fflush(headerptr);

  return(0);
}





/// image contents number
extern void set_image_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  // now setup the data output/input organization for chunking method for each number of columns
  *numcolumns=1;


  *numversion=0;
}


/// image content
/// uses global vars to get aa and lmin
int image_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  unsigned char liqb;
  FTYPE pr,iq,liq;
  FTYPE aa;
  FTYPE lmin;

  aa=imageparms[ORIGIN];
  lmin=imageparms[LMIN];

  pr=GLOBALMACP0A1(pimage,i,j,k,imagewhichpl);
  if (imagescale==LOG) iq = log(pr);
  else if(imagescale==LINEAR) iq = pr;
  
  liq = aa * (iq - lmin);
  if (liq > 255.)
    liq = 255.;
  if (liq < 0.)
    liq = 0.;
  
  liqb=(unsigned char)(liq);
  
  myset(datatype,&liqb,0,1,writebuf);

  return(0);
}
