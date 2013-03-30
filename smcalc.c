// compile with:
//   icc -Wall -O3 -axiMKW -unroll -ipo -rcd -tpp7 -march=pentium4 -mcpu=pentium4 -lm -o smcalc smcalc.c
//
// or just:
// gcc -O3 -Wall -o smcalc smcalc.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define sign(a) (copysign(1.0,a))

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#define RHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7
#define SMALL 1E-300
#define THIRD          0.3333333333333333333333333333333333L  /* 1/3 */

#include "global.nondepmnemonics.h"
#define WHICHPARA PARA4
#include "para_and_paraenohybrid.h"

#define NDIM 4
#define NPR 8
#define FTYPE double

#define PLOOP(pl) for(pl=0;pl<NPR;pl++)

// memory mapping
#define MAP(pl,max) (k*nx*ny*max+j*nx*max+i*max+(pl))
#define MAPGEN(i,j,k,pl,max) ((k)*nx*ny*max + (j)*nx*max + (i)*(max) + (pl))

// +-1 mappings
#define MAPIP1(pl,max) (k*nx*ny*max+j*nx*max+(i+1)*max+pl)
#define MAPIM1(pl,max) (k*nx*ny*max+j*nx*max+(i-1)*max+pl)
#define MAPJP1(pl,max) (k*nx*ny*max+(j+1)*nx*max+i*max+pl)
#define MAPJM1(pl,max) (k*nx*ny*max+(j-1)*nx*max+i*max+pl)
#define MAPKP1(pl,max) ((k+1)*nx*ny*max+j*nx*max+i*max+pl)
#define MAPKM1(pl,max) ((k-1)*nx*ny*max+j*nx*max+i*max+pl)

// variable mapping
#define MAP1 MAP(0,1)
#define MAPDER(pl) MAP(pl,3) // 3 dimensions
#define MAPJ(pl) MAP(pl,4) // 4 current components
#define MAPF(pl) MAP(pl,6) // 6 faraday components
#define MAPG(pl) MAP(pl,numgenvars)
#define MAPFICALC(pl) MAP(pl,3) // 3 dimensions

#define MAPP(pl) MAP(pl,NPR)
#define MAPPGEN(i,j,k,pl) MAPGEN(i,j,k,pl,NPR)
#define MAPPIP1(pl) MAPIP1(pl,NPR)
#define MAPPIM1(pl) MAPIM1(pl,NPR)
#define MAPPJP1(pl) MAPJP1(pl,NPR)
#define MAPPJM1(pl) MAPJM1(pl,NPR)
#define MAPPKP1(pl) MAPKP1(pl,NPR)
#define MAPPKM1(pl) MAPKM1(pl,NPR)



#define READHEADER 1 // should always be 1 generally
#define WRITEHEADER 1

// must keep up to date unless reading in numcolumns
#define NUMDUMPCOL1 36 // number of columns including i,j,k
#define NUMDUMPCOL2 56 // number of columns including i,j,k
// could be 56 too for jcon/jcov/fcon/fcov stuff if included in dump file

// DUMPTYPE 
#define REALDUMP 0
#define SINGLECOLUMN 1
#define ANYCOLUMN 2   // 2 or greater number of columns
// 0: real dump file
// 1: forfline dump file as in gammie.m (forfldump)
// 2: number of columns given for arbitrary size


//CALCTYPE
#define DOFARADAY 0
#define DOFLINE 1
#define DODER 2 // only use with DUMPTYPE==SINGLECOLUMN
#define DOAVG 3
#define DOFICALC 4


//#define MAXFILENAME (200)

static FTYPE Ficalc(int dimen, FTYPE *V, FTYPE *P);
static FTYPE ftilde( int dimen, int shift, FTYPE *Vabs, FTYPE *Pabs);


FTYPE SP0user;


int main(
         int argc,
         char *argv[],
         char *envp[]
         )
{
  int pl;
  FILE * dumpin;
  FILE *dumpout;
  FILE *colin;
  char filenamein[MAXFILENAME];
  char filenameout[MAXFILENAME];
  char dumpnamein[MAXFILENAME];
  char dumpnameout[MAXFILENAME];
  char systemcall[MAXFILENAME];
  int index;
  int i,j,k;
  int ii,jj,kk;
  int si,sj,sk;
  int iii,jjj,kkk;
  FTYPE X[NDIM],r,th,ph,divb,ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM],vmin1,vmax1,vmin2,vmax2,vmin3,vmax3;
  int totalsize[NDIM];
  long realnstep;
  FTYPE readnstep;
  FTYPE t,startx[NDIM],dx[NDIM],gam,a,R0,Rin,Rout,hslope,dt;
  int defcoord;
  FTYPE *p,*gdet,*aphi,*da,*db,*result1,*faradaydd,*faradayuu,*jcon,*jcov,*faradaydotj,*der;
  FTYPE *Ptot,*ficalc;
  int nx,ny,nz;
  int indexi,indexj,indexk;
  FTYPE result2;
  int DUMPTYPE,CALCTYPE,DERTYPE;
  FTYPE *genvar;
  FTYPE ftemp;
  int numgenvars;
  int dumpnum,startdump,enddump;
  int numcol;
  int inow,jnow;
  FTYPE (*Ypl)[20],*V,*P;
  FTYPE a_Ypl[NPR][20],a_V[20],a_P[20];
  int dimen,ismall;
  int numargs,argi;






  // shift sufficiently
  Ypl=(FTYPE (*) [20]) (&(a_Ypl[0][10]));
  V=&a_V[10];
  P=&a_P[10];


  numargs=1+8;
  if(argc<numargs){
    fprintf(stderr,"./smcalc DUMPTYPE CALCTYPE DERTYPE <Sp0> nx ny nz inname outname <startdump> <enddump>\n");
    fprintf(stderr,"DUMPTYPE: 0=realdump 1=single column >=2 any number of columns specified\n");
    fprintf(stderr,"CALCTYPE: 0=faraday 1=fline 2=derivatives 3=average 4=ficalc1,2,3\n");
    fprintf(stderr,"DERTYPE=0->centered when possible 1->numerical backward if possible 2=forward if possible 10+ of 0 11=+ of 1 12=+ of 2\n");
    fprintf(stderr,"Sp0=Shock strength. (e.g. 0.75 for flash)\n");
    exit(1);
  }


  argi=0;
  // always needed:
  argi++; DUMPTYPE=atoi(argv[argi]);
  argi++; CALCTYPE=atoi(argv[argi]);
  argi++; DERTYPE=atoi(argv[argi]);
  // optionals:
  if(CALCTYPE==4) argi++; SP0user=(FTYPE)atof(argv[argi]);

  // always needed:
  argi++; nx=atoi(argv[argi]);
  argi++; ny=atoi(argv[argi]);
  argi++; nz=atoi(argv[argi]);
  argi++; strcpy(filenamein,argv[argi]);
  argi++; strcpy(filenameout,argv[argi]);

  if(CALCTYPE==DOAVG){
    numargs=1+10;
    if(argc<numargs){
      fprintf(stderr,"CALCTYPE==DOAVG\n");
      fprintf(stderr,"./smcalc DUMPTYPE CALCTYPE nx ny nz inname outname <startdump> <enddump>\n");
      fprintf(stderr,"DUMPTYPE: 0=realdump 1=single column >=2 any number of columns specified\n");
      fprintf(stderr,"CALCTYPE: 0=faraday 1=fline 2=derivatives 3=average\n");
      fprintf(stderr,"DERTYPE=0->centered when possible 1->numerical backward if possible 2=forward if possible 10+ of 0 11=+ of 1 12=+ of 2\n");
      exit(1);
    }
    argi++; startdump=atoi(argv[argi]);
    argi++; enddump=atoi(argv[argi]);
  }
  else{
    startdump=0;
    enddump=0;
  }
  
  fprintf(stderr,"nx=%d ny=%d nz=%d dumptype=%d calctype=%d dertype=%d\n",nx,ny,nz,DUMPTYPE,CALCTYPE,DERTYPE); fflush(stderr);

  if(CALCTYPE==DOAVG){
    if(DUMPTYPE>=ANYCOLUMN) genvar=(FTYPE*)malloc(sizeof(FTYPE)*nx*ny*nz*DUMPTYPE);
    else if(DUMPTYPE==SINGLECOLUMN) genvar=(FTYPE*)malloc(sizeof(FTYPE)*nx*ny*nz*(DUMPTYPE+2));
    else if(DUMPTYPE==REALDUMP) genvar=(FTYPE*)malloc(sizeof(FTYPE)*nx*ny*nz*NUMDUMPCOL2); // bigger of 2

    if(genvar==NULL){
      fprintf(stderr,"problems allocating genvar memory\n"); fflush(stderr); exit(1);
    }

    if((DUMPTYPE>=ANYCOLUMN)||(DUMPTYPE==SINGLECOLUMN)){
      numgenvars=DUMPTYPE;
    }
    else if(DUMPTYPE==REALDUMP){
      numgenvars=NUMDUMPCOL2;
    }
    fprintf(stderr,"numgenvars=%d\n",numgenvars);

    for(index=0;index<nx*ny*nz*numgenvars;index++)  genvar[index]=0;

    sprintf(dumpnameout,"%s%s",filenameout,"avg");
    if((dumpout=fopen(dumpnameout,"wt"))==NULL){
      fprintf(stderr,"cannot open %s\n",dumpnameout);
      exit(1);
    }
  }
  else{
    if((dumpout=fopen(filenameout,"wt"))==NULL){
      fprintf(stderr,"cannot open %s\n",filenameout);
      exit(1);
    }
  }

  p=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz*NPR);
  gdet=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz);
  if(CALCTYPE==DOFLINE){
    aphi=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz);
    da=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz);
    db=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz);
    result1=(FTYPE *)malloc(sizeof(FTYPE)*MAX(nx,ny)); // used for either x-line or y-line processing of fieldline
  }
  else if(CALCTYPE==DOFARADAY){
    faradaydd=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz*6);
    faradayuu=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz*6);
  }
  else if(CALCTYPE==DODER){
    der=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz*3); // up to 3 derivatives since 3 dimensions (i.e. no mixed derivatives, only straight)
  }
  else if(CALCTYPE==DOFICALC){
    ficalc=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz*3); // for 3 dimensions
    Ptot=(FTYPE *)malloc(sizeof(FTYPE)*nx*ny*nz*1);
  }



  /////////////////
  //
  // check columns in dump files
  //
  /////////////////
  for(dumpnum=startdump;dumpnum<=enddump;dumpnum++){

    if((CALCTYPE!=DOAVG)||(startdump==enddump)){
      if((dumpin=fopen(filenamein,"rt"))==NULL){
        fprintf(stderr,"problem with file opening %s\n",filenamein);
        exit(1);
      }
      sprintf(systemcall,"cat %s | head -2 | tail -1 | wc -w > numcol.txt",filenamein);
      system(systemcall);
    }
    else{
      sprintf(dumpnamein,"%s%04d",filenamein,dumpnum);
      if((dumpin=fopen(dumpnamein,"rt"))==NULL){
        fprintf(stderr,"problem with file opening %s\n",dumpnamein);
        exit(1);
      }
      sprintf(systemcall,"cat %s | head -2 | tail -1 | wc -w > numcol.txt",dumpnamein);
      system(systemcall);
    }
    colin=fopen("numcol.txt","rt");
    fscanf(colin,"%d",&numcol);
    fclose(colin);
    if(CALCTYPE!=DOFLINE){
      if((DUMPTYPE==SINGLECOLUMN)&&(numcol!=4)){
        fprintf(stderr,"numcol=%d for singlecolumn, should be 4\n",numcol);
        //      exit(1);
      }
    }
    else{
      if((DUMPTYPE==SINGLECOLUMN)&&(numcol!=7)){
        fprintf(stderr,"numcol=%d for singlecolumn/dofline, should be 7\n",numcol);
        //      exit(1);
      }
    }
    if(CALCTYPE==DOAVG){
      if(numcol==NUMDUMPCOL1) numgenvars=NUMDUMPCOL1;
      else if(numcol==NUMDUMPCOL2) numgenvars=NUMDUMPCOL2;
      else{
        fprintf(stderr,"numcol=%d for realdump\n",numcol);
        exit(1);
      }
    }



    
       
    ////////////////////
    //
    // Process header
    //
    ////////////////////

#if(READHEADER)
    fprintf(stderr,"start reading header\n"); fflush(stderr);
    fscanf(dumpin, "%lf %d %d %d %lf %lf %lf %lf %lf %lf %ld %lf %lf %lf %lf %lf %lf %lf %d", &t,&totalsize[1],&totalsize[2],&totalsize[3],&startx[1],&startx[2],&startx[3],&dx[1],&dx[2],&dx[3],&realnstep,&gam,&a,&R0,&Rin,&Rout,&hslope,&dt,&defcoord);
    
    //realnstep=(long)readnstep;
    while(fgetc(dumpin)!='\n');
#endif
#if(WRITEHEADER)
    if(dumpnum==startdump){
      fprintf(stderr,"start writing header"); fflush(stderr);
      //      fprintf(stderr,"dt=%g defcoord=%d\n",dt,defcoord);
      fprintf(dumpout, "%10.5g %d %d %d %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %ld %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %d\n", t,totalsize[1], totalsize[2], totalsize[3], startx[1], startx[2], startx[3], dx[1],dx[2],dx[3],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord); fflush(dumpout);
    }
#endif

    if(CALCTYPE==DOAVG){
      fprintf(stderr,"start reading file: dumpnum=%d dumpname=%s\n",dumpnum,dumpnamein); fflush(stderr);
    }
    else{
      fprintf(stderr,"start reading file: dumpname=%s\n",filenamein); fflush(stderr);
    }



  


    ////////////////////
    //
    // LOOP over variables
    //
    ////////////////////

    for(iii=0;iii<nx*ny*nz;iii++){
      
      // all files assumed to have i and j and k
      fscanf(dumpin,"%d",&ii); //i
      fscanf(dumpin,"%d",&jj); //j
      fscanf(dumpin,"%d",&kk); //k

      // below accounts for if put in boundary zones.  Assumes now iii==0 will be first zone!  This wasn't assumed previously
      if(iii==0){
        si=ii;
        sj=jj;
        sk=kk;
      }
      i=ii-si;
      j=jj-sj;
      k=kk-sk;

      if((i<0)||(j<0)||(k<0)||(i>=nx)||(j>=ny)||(k>=nz)){
        fprintf(stderr,"data requested larger than said size: i=%d j=%d k=%d\n",i,j,k); fflush(stderr);
        exit(1);
      }


      /////////////////////////////////
      //      
      // read in different file types
      //
      /////////////////////////////////
      if(DUMPTYPE==REALDUMP){

        if(CALCTYPE!=DOAVG){
          // 3+3+NPR+1+4*4+4+1=36
          for(pl=1;pl<=3;pl++) fscanf(dumpin,"%lf",&X[pl]);
          fscanf(dumpin,"%lf",&r);
          fscanf(dumpin,"%lf",&th);
          fscanf(dumpin,"%lf",&ph);
          PLOOP(pl) fscanf(dumpin,"%lf", &p[MAPP(pl)]);
          fscanf(dumpin,"%lf",&divb);
   
          for (pl = 0; pl < NDIM; pl++)
            fscanf(dumpin,"%lf",&(ucon[pl]));
          for (pl = 0; pl < NDIM; pl++)
            fscanf(dumpin,"%lf",&(ucov[pl]));
          for (pl = 0; pl < NDIM; pl++)
            fscanf(dumpin,"%lf",&(bcon[pl]));
          for (pl = 0; pl < NDIM; pl++)
            fscanf(dumpin,"%lf",&(bcov[pl]));
   
          fscanf(dumpin,"%lf",&vmin1);
          fscanf(dumpin,"%lf",&vmax1);
          fscanf(dumpin,"%lf",&vmin2);
          fscanf(dumpin,"%lf",&vmax2);
          fscanf(dumpin,"%lf",&vmin3);
          fscanf(dumpin,"%lf",&vmax3);
   
          fscanf(dumpin,"%lf",&gdet[MAP1]);
        }
        else{
          genvar[MAPG(0)]=i+genvar[MAPG(0)];
          genvar[MAPG(1)]=j+genvar[MAPG(1)];
          genvar[MAPG(2)]=k+genvar[MAPG(2)];
          for(pl=3;pl<numgenvars;pl++){
            fscanf(dumpin,"%lf", &ftemp);
            genvar[MAPG(pl)]=ftemp+genvar[MAPG(pl)];
          }
   
        }
      }
      else if(DUMPTYPE==SINGLECOLUMN){
        if(CALCTYPE==DOFLINE){
          for(pl=5;pl<=7;pl++) fscanf(dumpin,"%lf", &p[MAPP(pl)]);
          //   for(pl=5;pl<=7;pl++) fprintf(stderr,"%lf ", p[MAPP(pl)]);
          fscanf(dumpin,"%lf",&gdet[MAP1]);
          //   fprintf(stderr,"%lf \n",gdet[MAP1]);
        }
        else if(CALCTYPE==DODER){ // dummy variable position
          for(pl=0;pl<=0;pl++) fscanf(dumpin,"%lf", &p[MAPP(pl)]);
        }
        else if(CALCTYPE==DOAVG){ // dummy variable position
          genvar[MAPG(0)]=i+genvar[MAPG(0)];
          genvar[MAPG(1)]=j+genvar[MAPG(1)];
          genvar[MAPG(2)]=k+genvar[MAPG(2)];
          for(pl=3;pl<numgenvars;pl++){
            fscanf(dumpin,"%lf", &ftemp);
            genvar[MAPG(pl)]=ftemp+genvar[MAPG(pl)];
          }
        }     
      }
      else if(DUMPTYPE>=ANYCOLUMN){
        if(CALCTYPE==DOAVG){ // dummy variable position
          //fprintf(stderr,"i=%d j=%d mapg(0)=%d\n",i,j,MAPG(0));
          genvar[MAPG(0)]=i+genvar[MAPG(0)];
          genvar[MAPG(1)]=j+genvar[MAPG(1)];
          genvar[MAPG(2)]=k+genvar[MAPG(2)];
          for(pl=3;pl<numgenvars;pl++){
            fscanf(dumpin,"%lf", &ftemp);
            //fprintf(stderr,"mapg(%d)=%d\n",pl,MAPG(pl));
            genvar[MAPG(pl)]=ftemp+genvar[MAPG(pl)];
          }
        }
        else if(CALCTYPE==DOFICALC){ // dummy variable position
          for(pl=0;pl<NPR;pl++){
            fscanf(dumpin,"%lf", &p[MAPP(pl)]);
            //     fprintf(stderr,"iii=%d ii=%d jj=%d i=%d j=%d :: p[%d]=%21.15g\n",iii,ii,jj,i,j,pl,p[MAPP(pl)]);
          }
          fscanf(dumpin,"%lf",&Ptot[MAP1]);
          fscanf(dumpin,"%lf",&gdet[MAP1]);
        }
      }
      

      // perform read-time calculations
      if(CALCTYPE==DOFLINE){
        // field line stuff
        da[MAP1]=-p[MAPP(6)]*gdet[MAP1]*dx[1];
        db[MAP1]=p[MAPP(5)]*gdet[MAP1]*dx[2];
        //fprintf(stdout,"%g %g :: %g %g :: %g %g\n",p[MAPP(6)],p[MAPP(5)],dx[1],dx[2],da[MAP1],db[MAP1]);
      }
      else if(CALCTYPE==DOFARADAY){
        // faraday (down down)
        ftemp=gdet[MAP1];
        faradaydd[MAPF(0)]=ftemp*(bcon[3]*ucon[2]-bcon[2]*ucon[3]);
        faradaydd[MAPF(1)]=ftemp*(bcon[1]*ucon[3]-bcon[3]*ucon[1]);
        faradaydd[MAPF(2)]=ftemp*(bcon[2]*ucon[1]-bcon[1]*ucon[2]);
        faradaydd[MAPF(3)]=ftemp*(bcon[3]*ucon[0]-bcon[0]*ucon[3]);
        faradaydd[MAPF(4)]=ftemp*(bcon[0]*ucon[2]-bcon[2]*ucon[0]);
        faradaydd[MAPF(5)]=ftemp*(bcon[1]*ucon[0]-bcon[0]*ucon[1]);
 
        // faraday (up up)
        ftemp=-1.0/ftemp;
        faradayuu[MAPF(0)]=ftemp*(bcov[3]*ucov[2]-bcov[2]*ucov[3]);
        faradayuu[MAPF(1)]=ftemp*(bcov[1]*ucov[3]-bcov[3]*ucov[1]);
        faradayuu[MAPF(2)]=ftemp*(bcov[2]*ucov[1]-bcov[1]*ucov[2]);
        faradayuu[MAPF(3)]=ftemp*(bcov[3]*ucov[0]-bcov[0]*ucov[3]);
        faradayuu[MAPF(4)]=ftemp*(bcov[0]*ucov[2]-bcov[2]*ucov[0]);
        faradayuu[MAPF(5)]=ftemp*(bcov[1]*ucov[0]-bcov[0]*ucov[1]);
      }

    }// end over variable loop
    
    fclose(dumpin);
    fprintf(stderr,"done reading file\n"); fflush(stderr);
  }// end over dump loop


  fprintf(stderr,"start post-read calculations\n"); fflush(stderr);






  ///////////////////////////////////
  //
  // perform post-read calculations
  //
  ///////////////////////////////////

  if(CALCTYPE==DOAVG){
    fprintf(stderr,"done reading all files\n"); fflush(stderr);

    for(pl=0;pl<nx*ny*nz*numgenvars;pl++){
      genvar[pl]=genvar[pl]/(enddump-startdump+1);
    }
  }
  else if(CALCTYPE==DOFLINE){

#if(0) // new Sasha-Jon algorithm
    // result1 is of size nx always

    // first set 0 point
    for(iii=0;iii<nx*ny*nz;iii++)    aphi[iii]=0.0;

    for(kkk=0;kkk<nz;kkk++){ // same procedure for each kkk out of nz

      // generate x-direction integral since 1-D
      for(iii=0;iii<nx;iii++){
        result1[iii]=-da[kkk*nx*ny+0]*0.5;
        // jjj here is index in i.  We are doing integral over dx1
        for(jjj=0;jjj<=iii;jjj++){
          if(jjj==iii) result1[iii]+=da[kkk*nx*ny+jjj]*0.5;
          else result1[iii]+=da[kkk*nx*ny+jjj];
        }
      }
      // now find i,j'th aphi
      for(iii=0;iii<nx*ny;iii++){
        i=indexi=(int)(iii%nx);
        j=indexj=(int)((iii%(nx*ny))/nx);
        k=indexk=(int)(iii/(nx*ny));
      
        // here jjj is over dx2
        result2=-db[kkk*nx*ny+0*nx+indexi]*0.5;
        for(jjj=0;jjj<=indexj;jjj++){
          if(jjj==indexj) result2+=db[kkk*nx*ny+jjj*nx+indexi]*0.5;
          else result2+=db[kkk*nx*ny+jjj*nx+indexi];
        }
      
        aphi[kkk*nx*ny+indexj*nx+indexi]=result1[indexi]+result2;
      }
    }// end over kkk out of nz

#elif(1) // new Sasha-Jon algorithm 2

    // first set 0 point
    for(iii=0;iii<nx*ny*nz;iii++)    aphi[iii]=0.0;
    
    for(kkk=0;kkk<nz;kkk++){ // same procedure for each kkk out of nz
     
      //      int fakekkk=0; // forced
      int fakekkk=kkk;

      indexi=0;
      // generate x-direction integral since 1-D
      for(jnow=0;jnow<ny;jnow++){
        result1[jnow]=-db[fakekkk*nx*ny+0]*0.5;
        // jjj here is index in i.  We are doing integral over dx1
        for(indexj=0;indexj<=jnow;indexj++){
          if(indexj==jnow) result1[jnow]+=db[fakekkk*nx*ny+nx*indexj+indexi]*0.5;
          else result1[jnow]+=db[fakekkk*nx*ny+nx*indexj+indexi];
        }
      }
      // now find i,j'th aphi
      for(iii=0;iii<nx*ny;iii++){
        i=indexi=(int)(iii%nx);
        j=indexj=(int)((iii%(nx*ny))/nx);
        k=indexk=(int)(iii/(nx*ny));
      
        // here jjj is over dx2
        result2=-da[fakekkk*nx*ny+indexj*nx+0]*0.5;
        for(inow=0;inow<=indexi;inow++){
          if(inow==indexi) result2+=da[fakekkk*nx*ny+indexj*nx+inow]*0.5;
          else result2+=da[fakekkk*nx*ny+indexj*nx+inow];
        }
      
        aphi[fakekkk*nx*ny+indexj*nx+indexi]=result1[indexj]+result2;
      }
    }// end over kkk out of nz

#endif

  }
  else if(CALCTYPE==DODER){

    int signit;
    if(DERTYPE<10) signit=-1;
    else signit=1;

    for(iii=0;iii<nx*ny*nz;iii++){
      i=indexi=(int)(iii%nx);
      j=indexj=(int)((iii%(nx*ny))/nx);
      k=indexk=(int)(iii/(nx*ny));


      // x1
      if(nx>1){
        if( ((DERTYPE==2 || DERTYPE==12) || (indexi==0))&&(indexi!=nx-1) ){// forward difference
          if(signit==-1) der[MAPDER(0)]=(p[MAPPIP1(0)]-p[MAPP(0)])/dx[1];
          else  der[MAPDER(0)]=(fabs(p[MAPPIP1(0)])+fabs(p[MAPP(0)]))/dx[1];
        }
        else if( ((DERTYPE==1 || DERTYPE==11)||(indexi==nx-1))&&(indexi!=0) ){// backward difference
          if(signit==-1) der[MAPDER(0)]=(p[MAPP(0)]-p[MAPPIM1(0)])/dx[1];
          else der[MAPDER(0)]=(fabs(p[MAPP(0)])+fabs(p[MAPPIM1(0)]))/dx[1];
        }
        else if(DERTYPE==0 || DERTYPE==10) {// centered difference
          if(signit==-1) der[MAPDER(0)]=0.5*(p[MAPPIP1(0)]-p[MAPPIM1(0)])/dx[1];
          else der[MAPDER(0)]=0.5*(fabs(p[MAPPIP1(0)])+fabs(p[MAPPIM1(0)]))/dx[1];
        }
        else{
          der[MAPDER(0)]=0.0;
        }
      }
      else{
        der[MAPDER(0)]=0.0;
      }


      // x2
      if(ny>1){
        if( ((DERTYPE==2 || DERTYPE==12) || (indexj==0))&&(indexj!=ny-1) ){// forward difference
          if(signit==-1) der[MAPDER(1)]=(p[MAPPJP1(0)]-p[MAPP(0)])/dx[2];
          else der[MAPDER(1)]=(fabs(p[MAPPJP1(0)])+fabs(p[MAPP(0)]))/dx[2];
        }
        else if( ((DERTYPE==1 || DERTYPE==11)||(indexj==ny-1))&&(indexj!=0) ){// backward difference
          if(signit==-1) der[MAPDER(1)]=(p[MAPP(0)]-p[MAPPJM1(0)])/dx[2];
          else der[MAPDER(1)]=(fabs(p[MAPP(0)])+fabs(p[MAPPJM1(0)]))/dx[2];
        }
        else if(DERTYPE==0 || DERTYPE==10){// centered difference
          if(signit==-1) der[MAPDER(1)]=0.5*(p[MAPPJP1(0)]+signit*p[MAPPJM1(0)])/dx[2];
          else der[MAPDER(1)]=0.5*(fabs(p[MAPPJP1(0)])+fabs(p[MAPPJM1(0)]))/dx[2];
        }
        else{
          der[MAPDER(1)]=0.0;
        }
      }
      else{
        der[MAPDER(1)]=0.0;
      }

      // x3
      if(nz>1){
        if( ((DERTYPE==2 || DERTYPE==12) || (indexk==0))&&(indexk!=nz-1) ){// forward difference
          if(signit==-1) der[MAPDER(2)]=(p[MAPPKP1(0)]-p[MAPP(0)])/dx[3];
          else der[MAPDER(2)]=(fabs(p[MAPPKP1(0)])+fabs(p[MAPP(0)]))/dx[3];
        }
        else if( ((DERTYPE==1 || DERTYPE==11)||(indexk==nz-1))&&(indexk!=0) ){// backward difference
          if(signit==-1) der[MAPDER(2)]=(p[MAPP(0)]-p[MAPPKM1(0)])/dx[3];
          else der[MAPDER(2)]=(fabs(p[MAPP(0)])+fabs(p[MAPPKM1(0)]))/dx[3];
        }
        else if(DERTYPE==0 || DERTYPE==10){// centered difference
          if(signit==-1) der[MAPDER(2)]=0.5*(p[MAPPKP1(0)]+signit*p[MAPPKM1(0)])/dx[3];
          else der[MAPDER(2)]=0.5*(fabs(p[MAPPKP1(0)])+fabs(p[MAPPKM1(0)]))/dx[3];
        }
        else{
          der[MAPDER(2)]=0.0;
        }
      }
      else{
        der[MAPDER(2)]=0.0;
      }

    }// end over iii
  }
  else if(CALCTYPE==DOFICALC){


    for(dimen=1;dimen<=3;dimen++){
      for(iii=0;iii<nx*ny*nz;iii++){
        i=indexi=(int)(iii%nx);
        j=indexj=(int)((iii%(nx*ny))/nx);
        k=indexk=(int)(iii/(nx*ny));
#define NUMBC 3
        if(i>=NUMBC && j>=NUMBC && k>=NUMBC && i<=nx-1-NUMBC && j<=ny-1-NUMBC && k<=nz-1-NUMBC){
          // create V,P,Y
          for(ismall=-NUMBC;ismall<=+NUMBC;ismall++){

            // get V = gdet*rho0*u^i
            V[ismall]  = p[MAPPGEN(i + ismall*(dimen==1) , j + ismall*(dimen==2) , k + ismall*(dimen==3) , U1+dimen-1)]; // u^i
            V[ismall] *= p[MAPPGEN(i + ismall*(dimen==1) , j + ismall*(dimen==2) , k + ismall*(dimen==3) , RHO)]; // rho0
            V[ismall] *= gdet[MAPGEN(i + ismall*(dimen==1) , j + ismall*(dimen==2)  , k + ismall*(dimen==3) , 0, 1)]; // gdet

            //     P[ismall] = p[MAPPGEN(i + ismall*(dimen==1),j + ismall*(dimen==2) , k + ismall*(dimen==3) , UU)]; // assumes ideal gas and only taking ratios of pressures
            P[ismall] = Ptot[MAPGEN(i + ismall*(dimen==1) , j + ismall*(dimen==2)  , k + ismall*(dimen==3) , 0, 1)]; // true Ptot


            //     PLOOP(pl) Ypl[pl][ismall] = p[MAPPGEN(i + ismall*(dimen==1),j + ismall*(dimen==2)  , k + ismall*(dimen==3) , pl)];

            //     fprintf(stderr,"i=%d j=%d k=%d ismall=%d V=%21.15g P=%21.15g\n",i,j,k,ismall,V[ismall],P[ismall]);
          }

          // now compute Ficalc for all pl
          //   ficalc[MAPFICALC(dimen-1)]=Ficalc(dimen, V, P, Ypl);

          ficalc[MAPFICALC(dimen-1)]=Ficalc(dimen, V, P);

          //   fprintf(stderr,"dimen=%d ficalc=%21.15g\n",dimen,ficalc[MAPFICALC(dimen-1)]);
        }
        else{
          // then set to 0
          ficalc[MAPFICALC(dimen-1)]=0.0;
        }
      }// end over iii
    }// end over dimen
  }
  fprintf(stderr,"done post-read calculations\n"); fflush(stderr);







  //////////////////////////////
  //
  // write calculations to file
  //
  //////////////////////////////

  if(CALCTYPE==DOAVG){
    fprintf(stderr,"start writing file: dumpname=%s\n",dumpnameout); fflush(stderr);
  }
  else{
    fprintf(stderr,"start writing file: dumpname=%s\n",filenameout); fflush(stderr);
  }

  for(iii=0;iii<nx*ny*nz;iii++){
    i=indexi=(int)(iii%nx);
    j=indexj=(int)((iii%(nx*ny))/nx);
    k=indexk=(int)(iii/(nx*ny));

    //    if(DUMPTYPE==REALDUMP){
    //  fprintf(dumpout,"%d %d %d ",i,j,k);
    // }

    if(CALCTYPE==DOFLINE){
      fprintf(dumpout,"%21.15g ",aphi[indexk*nx*ny+indexj*nx+indexi]);
    }
    else if(CALCTYPE==DOFARADAY){
      fprintf(dumpout,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "
              ,faradaydd[indexk*nx*ny*6+indexj*nx*6+indexi*6+0]
              ,faradaydd[indexk*nx*ny*6+indexj*nx*6+indexi*6+1]
              ,faradaydd[indexk*nx*ny*6+indexj*nx*6+indexi*6+2]
              ,faradaydd[indexk*nx*ny*6+indexj*nx*6+indexi*6+3]
              ,faradaydd[indexk*nx*ny*6+indexj*nx*6+indexi*6+4]
              ,faradaydd[indexk*nx*ny*6+indexj*nx*6+indexi*6+5]
              ,faradayuu[indexk*nx*ny*6+indexj*nx*6+indexi*6+0]
              ,faradayuu[indexk*nx*ny*6+indexj*nx*6+indexi*6+1]
              ,faradayuu[indexk*nx*ny*6+indexj*nx*6+indexi*6+2]
              ,faradayuu[indexk*nx*ny*6+indexj*nx*6+indexi*6+3]
              ,faradayuu[indexk*nx*ny*6+indexj*nx*6+indexi*6+4]
              ,faradayuu[indexk*nx*ny*6+indexj*nx*6+indexi*6+5]
              );
    }
    else if(CALCTYPE==DODER){
      fprintf(dumpout,"%21.15g %21.15g %21.15g "
              ,der[MAPDER(0)]
              ,der[MAPDER(1)]
              ,der[MAPDER(2)]
              );
    }
    else if(CALCTYPE==DOFICALC){
      fprintf(dumpout,"%21.15g %21.15g %21.15g ",ficalc[MAPFICALC(0)],ficalc[MAPFICALC(1)],ficalc[MAPFICALC(2)]);
    }
    else if(CALCTYPE==DOAVG){
      for(pl=0;pl<numgenvars;pl++){
        fprintf(dumpout,"%21.15g ",genvar[MAPG(pl)]);
      }
    }
    fprintf(dumpout,"\n");
  }



  //////////////////////////////
  //
  // free memory
  //
  //////////////////////////////
  if(CALCTYPE==DOAVG){
    free(genvar);
  }
  free(p);
  free(gdet);
  if(CALCTYPE==DOFLINE){
    free(aphi);
    free(da);
    free(db);
    free(result1);
  }
  else if(CALCTYPE==DOFARADAY){
    free(faradaydd);
    free(faradayuu);
  }
  else if(CALCTYPE==DODER){
    free(der);
  }
  else if(CALCTYPE==DOFICALC){
    free(ficalc);
  }
  fprintf(stderr,"done with all smcalc: dumpout=%lld\n",(long long int) dumpout); fflush(stderr);
  fclose(dumpout);

}



// SP0->SP0user so smcalc can choose shock strength to get

// PPM FLATTENER parameter
static FTYPE ftilde( int dimen, int shift, FTYPE *Vabs, FTYPE *Pabs)
{
  FTYPE Ftilde,Ftilde1,Ftilde2;
  FTYPE Sp;
  FTYPE *V, *P;
  FTYPE P2diff,Pbottom;


  // shift as needed
  P = Pabs + shift;
  V = Vabs + shift;

  // FLASH Equation 43
  P2diff=P[2]-P[-2];
  Pbottom=sign(P2diff)/(fabs(P2diff)+SMALL); // singularity avoidance but keeps signature
  Sp = (P[1] - P[-1]) * Pbottom ;



  // FLASH Equation 45
  Ftilde = max( 0, min( 1.0, 10.0 * (Sp - SP0user) ) );

  //  Ftilde*=Ftilde*Ftilde*Ftilde;

  // FLASH Equation 46
  Ftilde1 = fabs(P[1] - P[-1]) / (min(fabs(P[1]), fabs(P[-1]))+ SMALL );
  Ftilde *= ( (FTYPE)(Ftilde1>=THIRD) );
  //  if(Ftilde1<THIRD) Ftilde=0.0;

  // FLASH Equation 47
  Ftilde2 = V[1] - V[-1];
  Ftilde *= ( (FTYPE)(Ftilde2<=0.0) );
  //  if(Ftilde2>0.0) Ftilde=0.0;
  
  //  if(Sp>0.8){
  //    fprintf(stderr,"Pbottom=%21.15g Sp=%21.15g Ftilde=%21.15g\n",Pbottom,Sp,Ftilde);
  //  }
  
  return( Ftilde );
}


// PPM FLATTENERS (final formula)
static FTYPE  Ficalc(int dimen, FTYPE *V, FTYPE *P)
{
  FTYPE ftilde( int dimen, int shift, FTYPE *P, FTYPE *V);
  int signdP;
  FTYPE Fi;

  signdP = (P[1] - P[-1] > 0) * 2 - 1;
  // FLASH Equation 48
  Fi = max( ftilde(dimen, 0, V,P), ftilde(dimen, -signdP, V,P) );

  return(Fi);
}

