#include "decs.h"

// mpi.h has following datatypes corresponding to the C types
// pick one per dump file. (no per column types yet)
// same for image.c
// same for restart.c
/*
  #define MPI_CHAR           ((MPI_Datatype)1)
  #define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
  #define MPI_BYTE           ((MPI_Datatype)3)
  #define MPI_SHORT          ((MPI_Datatype)4)
  #define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
  #define MPI_INT            ((MPI_Datatype)6)
  #define MPI_UNSIGNED       ((MPI_Datatype)7)
  #define MPI_LONG           ((MPI_Datatype)8)
  #define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
  #define MPI_FLOAT          ((MPI_Datatype)10)
  #define MPI_DOUBLE         ((MPI_Datatype)11)
  #define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
  #define MPI_LONG_LONG_INT  ((MPI_Datatype)13)
*/



/* Follow these steps to create a new dump file

   1) defs.???.h: create the storage variable if needed

   2) set_array???.c : shift variable if necessary (follow examples)

   3) global.nondepmnemonics.h : change NUMDUMPTYPES and add label and add to MYDUMPNAMELIST

   4) global.???.h change NUMDUMPTYPES and add another entry if want separate timing of output

   5) dump.c : add dnumcolumns[LABEL]=NUMCOLUMNS where NUMCOLUMNS is number of entries in dump file

   6) dump.c : follow examples from here (dump() uses dump_header() and dump_content()).  One must define the header and content function and the wrapper (3 functions) or use an existing header function

   7) diag.c : follow example of "dump", dumpc, tlastdump, etc.

   8) init.c : DTdumpgen, dumpcntgen, and other things.

   9) global.dump.h : add global prototypes

*/







int init_dumps(void)
{

  trifprintf("begin: init_dumps\n");

  ///////////////////////////
  //
  // Output nprlist information to SM-readable file
  //
  ///////////////////////////
  output_nprlist_info();

  ///////////////////////////
  //
  // setup number of columns per dump file
  // (see dumpgen.c or dump.c for how used)
  //
  ///////////////////////////
  init_dnumcolumns_dnumversion();


  if(mpicombine==1 && mpicombinetype==MPICOMBINEMINMEM){
    ///////////////////////////
    //
    // setup link list (only used for MINMEM method)
    //
    ///////////////////////////
    init_linklists();
  }

  trifprintf("end: init_dumps\n");


  return(0);
}


void output_nprlist_info(void)
{
  int pliter,pl;
  int numversion;
  int numlines;
  FILE *out;
  
  // only CPU=0
  if(myid==0){

    out=fopen("nprlistinfo.dat","wt");
    if(out==NULL){
      dualfprintf(fail_file,"Couldn't open nprlistinfo.dat\n");
      myexit(12358235);
    }

    numlines=7;
    numversion=0; // version number of this file
    
    myfprintf(out,"%d %d\n",numlines,numversion);

    // NPR: (conserved quantities: U0-U? in SM)
    PLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");

    // NPR2INTERP:
    PINTERPLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");

    // NPR2NOTINTERP:
    PNOTINTERPLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");

    // NPRBOUND:
    PBOUNDLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");

    // NPRFLUXBOUND:
    PFLUXBOUNDLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");

    // NPRDUMP:
    PDUMPLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");

    // NPRINVERT:
    PINVERTLOOP(pliter,pl){
      myfprintf(out,"%d ",pl);
    }
    if(pliter==0) myfprintf(out,"-1"); // nothing in this list
    myfprintf(out,"\n");


    fclose(out);

  }// end if CPU==0


}



// setup number of columns per dump file (see dumpgen.c or dump.c for how used)
void init_dnumcolumns_dnumversion(void)
{
  char tempdumpnamelist[NUMDUMPTYPES][MAXFILENAME]=MYDUMPNAMELIST;
  int i;

  extern void set_image_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_dump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_gdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_avg_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_avg2_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_debug_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_enodebug_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_fieldline_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_dissdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_dumpother_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_fluxdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_eosdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_vpotdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_failfloordudump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);

  extern void set_rupperpoledump_read_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_rupperpoledump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);

  extern void set_rdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);

  extern void set_rmetricdump_read_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_rmetricdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);



  // assign local to global -- do this since global not easily assigned to value with defs and decs approach
  int dumpiter;
  for(dumpiter=0;dumpiter<NUMDUMPTYPES;dumpiter++){
    strcpy(dumpnamelist[dumpiter],tempdumpnamelist[dumpiter]);
  }



  // always numcolumns=0 for fake dump
  // version=0 shouldn't matter for fake dump
  dnumcolumns[FAKEDUMPTYPE]=0; dnumversion[FAKEDUMPTYPE]=0;
  

  // image
  set_image_content_dnumcolumns_dnumversion(&dnumcolumns[IMAGEDUMPTYPE],&dnumversion[IMAGEDUMPTYPE]);
  // dump
  set_dump_content_dnumcolumns_dnumversion(&dnumcolumns[MAINDUMPTYPE],&dnumversion[MAINDUMPTYPE]);
  // gdump
  set_gdump_content_dnumcolumns_dnumversion(&dnumcolumns[GRIDDUMPTYPE],&dnumversion[GRIDDUMPTYPE]);
  // avg
  set_avg_content_dnumcolumns_dnumversion(&dnumcolumns[AVG1DUMPTYPE],&dnumversion[AVG1DUMPTYPE]);
  // avg2
  set_avg2_content_dnumcolumns_dnumversion(&dnumcolumns[AVG2DUMPTYPE],&dnumversion[AVG2DUMPTYPE]);
  // debug
  set_debug_content_dnumcolumns_dnumversion(&dnumcolumns[DEBUGDUMPTYPE],&dnumversion[DEBUGDUMPTYPE]);
  // enodebug
  set_enodebug_content_dnumcolumns_dnumversion(&dnumcolumns[ENODEBUGDUMPTYPE],&dnumversion[ENODEBUGDUMPTYPE]);
  // fieldline
  set_fieldline_content_dnumcolumns_dnumversion(&dnumcolumns[FIELDLINEDUMPTYPE],&dnumversion[FIELDLINEDUMPTYPE]);
  // dissdump
  set_dissdump_content_dnumcolumns_dnumversion(&dnumcolumns[DISSDUMPTYPE],&dnumversion[DISSDUMPTYPE]);
  // dumpother
  set_dumpother_content_dnumcolumns_dnumversion(&dnumcolumns[OTHERDUMPTYPE],&dnumversion[OTHERDUMPTYPE]);
  // fluxdump
  set_fluxdump_content_dnumcolumns_dnumversion(&dnumcolumns[FLUXDUMPTYPE],&dnumversion[FLUXDUMPTYPE]);
  // eosdump
  set_eosdump_content_dnumcolumns_dnumversion(&dnumcolumns[EOSDUMPTYPE],&dnumversion[EOSDUMPTYPE]);
  // vpotdump
  set_vpotdump_content_dnumcolumns_dnumversion(&dnumcolumns[VPOTDUMPTYPE],&dnumversion[VPOTDUMPTYPE]);
  // failfloordudump
  set_failfloordudump_content_dnumcolumns_dnumversion(&dnumcolumns[FAILFLOORDUDUMPTYPE],&dnumversion[FAILFLOORDUDUMPTYPE]);

  // rdump (must come after all normal dumps since dnumcolumns used by restart to store other things needed up restart that are dealt with also above)
  // rupperpoledump
  set_rupperpoledump_content_dnumcolumns_dnumversion(&dnumcolumns[RESTARTUPPERPOLEDUMPTYPE],&dnumversion[RESTARTUPPERPOLEDUMPTYPE]);
  set_rdump_content_dnumcolumns_dnumversion(&dnumcolumns[RESTARTDUMPTYPE],&dnumversion[RESTARTDUMPTYPE]);
  // rmetricdump
  set_rmetricdump_content_dnumcolumns_dnumversion(&dnumcolumns[RESTARTMETRICDUMPTYPE],&dnumversion[RESTARTMETRICDUMPTYPE]);



  trifprintf("dump number of columns(see global.nondepmnemonics.h)\n");
  for(i=0;i<NUMDUMPTYPES;i++){
    trifprintf("%s dnumcolumns[%d]=%d\n",dumpnamelist[i],i,dnumcolumns[i]);
  }
  trifprintf("\n");


}







int dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dump# %ld ... ",dump_cnt);

  whichdump=MAINDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dump_content)>=1) return(1);

  /////// output the symmetry information to the fail file
  //writesyminfo();
  ///////

  trifprintf("end dumping dump# %ld ... ",dump_cnt);


  return(0);
  
}



/////// output the symmetry information to the fail file; symmetrizes w.r.t. i == j
void writesyminfo( void )
{
  int i, j;

  for( i = 0; i < N1; i++ ) {
  }

}


int dump_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{
  int dump_header_general(int whichdump, int whichdumpversion, int numcolumns, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr);
  int retval;

  retval=dump_header_general(whichdump, whichdumpversion, numcolumns, realnstep,dt, bintxt, headerptr);

  return(retval);

}

int fluxdump_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{
  int dump_header_general(int whichdump, int whichdumpversion, int numcolumns, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr);
  int retval;

  retval=dump_header_general(whichdump, whichdumpversion, numcolumns, fluxdumprealnstep,fluxdumpdt, bintxt, headerptr);

  return(retval);

}

int dump_header_general(int whichdump, int whichdumpversion, int numcolumns, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr)
{
  int realtotalsize[NDIM];
  FTYPE realstartx[NDIM];
  FTYPE X[NDIM];
  int is,ie,js,je,ks,ke;

  // set computational domain for output header
  is=Uconsevolveloop[FIS]+SHIFTX1DN;
  ie=Uconsevolveloop[FIE]+SHIFTX1UP;
  js=Uconsevolveloop[FJS]+SHIFTX2DN;
  je=Uconsevolveloop[FJE]+SHIFTX2UP;
  ks=Uconsevolveloop[FKS]+SHIFTX3DN;
  ke=Uconsevolveloop[FKE]+SHIFTX3UP;


  // get real startx's (assumes rectangular grid)
  realtotalsize[1]=totalsize[1]+2*EXTRADUMP1;
  if(EXTRADUMP1!=0){
    coord(0-EXTRADUMP1,0,0,FACE1,X);
    realstartx[1]=X[1];
  }
  else realstartx[1]=startx[1];

  realtotalsize[2]=totalsize[2]+2*EXTRADUMP2;
  if(EXTRADUMP2!=0){
    coord(0,0-EXTRADUMP2,0,FACE2,X);
    realstartx[2]=X[2];
  }
  else realstartx[2]=startx[2];

  realtotalsize[3]=totalsize[3]+2*EXTRADUMP3;
  if(EXTRADUMP3!=0){
    coord(0,0,0-EXTRADUMP3,FACE3,X);
    realstartx[3]=X[3];
  }
  else realstartx[3]=startx[3];
  
  // dx is the same (constant)

  // 15+3=18 elements total
  if(bintxt==BINARYOUTPUT){
    fwrite(&tsteppartf,sizeof(FTYPE),1,headerptr);
    fwrite(&realtotalsize[1],sizeof(int),1,headerptr);
    fwrite(&realtotalsize[2],sizeof(int),1,headerptr);
    fwrite(&realtotalsize[3],sizeof(int),1,headerptr);
    fwrite(&realstartx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&realstartx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&realstartx[3],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[3],sizeof(FTYPE),1,headerptr);
    fwrite(&localrealnstep,sizeof(long),1,headerptr);
    fwrite(&gam,sizeof(FTYPE),1,headerptr);
    fwrite(&a,sizeof(FTYPE),1,headerptr);
    fwrite(&R0,sizeof(FTYPE),1,headerptr);
    fwrite(&Rin,sizeof(FTYPE),1,headerptr);
    fwrite(&Rout,sizeof(FTYPE),1,headerptr);
    fwrite(&hslope,sizeof(FTYPE),1,headerptr);
    fwrite(&localdt,sizeof(FTYPE),1,headerptr);
    fwrite(&defcoord,sizeof(int),1,headerptr);
    fwrite(&MBH,sizeof(FTYPE),1,headerptr);
    fwrite(&QBH,sizeof(FTYPE),1,headerptr);
    fwrite(&EP3,sizeof(FTYPE),1,headerptr);
    fwrite(&THETAROT,sizeof(FTYPE),1,headerptr);
    fwrite(&is,sizeof(int),1,headerptr);
    fwrite(&ie,sizeof(int),1,headerptr);
    fwrite(&js,sizeof(int),1,headerptr);
    fwrite(&je,sizeof(int),1,headerptr);
    fwrite(&ks,sizeof(int),1,headerptr);
    fwrite(&ke,sizeof(int),1,headerptr);
    fwrite(&whichdump,sizeof(int),1,headerptr);
    fwrite(&whichdumpversion,sizeof(int),1,headerptr);
    fwrite(&numcolumns,sizeof(int),1,headerptr);
  }
  else{
#if(REALTYPE==DOUBLETYPE)
    fprintf(headerptr, "%21.15g %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %ld %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d %21.15g %21.15g %21.15g %21.15g %d %d %d %d %d %d %d %d %d\n", tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2], dx[3], localrealnstep,gam,a,R0,Rin,Rout,hslope,localdt,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "%31.25Lg %d %d %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %ld %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d %d %d %d %d %d %d %d %d\n", tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2],dx[3],localrealnstep,gam,a,R0,Rin,Rout,hslope,localdt,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns);
#endif
  }
  fflush(headerptr);
  return(0);
}       



// number of columns for dump file
void set_dump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{


  // always NPRDUMP
  if(GAMMIEDUMP)  *numcolumns=2*3 + NPRDUMP+NPR + 3 + 1 + NDIM * NDIM + 6 + 1
#if(CALCFARADAYANDCURRENTS)
                    + NDIM*2
                    + 2*6
#endif
                    ;
  else{
    *numcolumns=3*3 + NPRDUMP + 3 + (nprend+1) + 1 + NDIM * NDIM + 6 + 1  //replace NPR -> (nprend+1) since nprend, not NPR, controls dumping.  Fixes: DOEXTRAINTERP=1 case
#if(CALCFARADAYANDCURRENTS)
      + NDIM*2
      + 2*6
#endif
      ;    // 61 total if also doing currents and faraday, 41 otherwise

    if(FLUXB==FLUXCTSTAG && 0){ // DEBUG (change corresponding code in dump.c)
      *numcolumns+= NPR2INTERP*COMPDIM*2 + NPR + COMPDIM*3*2 + COMPDIM*3*2*2;
    }
  }


  // Version number:
  *numversion=0;

}



int dump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl,pliter;
  FTYPE r, th, vmin[NDIM], vmax[NDIM];
  int ignorecourant;
  struct of_state qdontuse;
  struct of_state *qptr=&qdontuse;
  FTYPE X[NDIM],V[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  FTYPE ftemp;
  FTYPE jcov[NDIM];
  FTYPE fcov[NUMFARADAY];
  FTYPE rho,u,pressure,cs2,Sden;
  int dir,l,m,n,o;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int loc=CENT;

  //////////////
  //
  // some calculations
  //

  bl_coord_ijk_2(i, j, k, loc, X,V);
  // if failed, then data output for below invalid, but columns still must exist    

  get_geometry(i, j, k, loc, ptrgeom);

  if (!failed) {
    if (get_state(GLOBALMAC(pdump,i,j,k), ptrgeom, qptr) >= 1)
      FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
    if (vchar(GLOBALMAC(pdump,i,j,k), qptr, 1, ptrgeom, &vmax[1], &vmin[1],&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 1);
    if (vchar(GLOBALMAC(pdump,i,j,k), qptr, 2, ptrgeom, &vmax[2], &vmin[2],&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 2);
    if (vchar(GLOBALMAC(pdump,i,j,k), qptr, 3, ptrgeom, &vmax[3], &vmin[3],&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 3);
  }
  else {// do a per zone check, otherwise set to 0
    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
    if (get_state(GLOBALMAC(pdump,i,j,k), ptrgeom, qptr) >= 1){
      for (pl = 0; pl < NDIM; pl++)
        qptr->ucon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
        qptr->ucov[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
        qptr->bcon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
        qptr->bcov[pl]=0;
    }
    if (vchar(GLOBALMAC(pdump,i,j,k), qptr, 1, ptrgeom, &vmax[1], &vmin[1],&ignorecourant) >= 1){
      vmax[1]=vmin[1]=0;
    }
    
    if (vchar(GLOBALMAC(pdump,i,j,k), qptr, 2, ptrgeom, &vmax[2], &vmin[2],&ignorecourant) >= 1){
      vmax[2]=vmin[2]=0;
    }

    if (vchar(GLOBALMAC(pdump,i,j,k), qptr, 3, ptrgeom, &vmax[3], &vmin[3],&ignorecourant) >= 1){
      vmax[3]=vmin[3]=0;
    }

    whocalleducon=0; // return to normal state
    
  }


  // GODMARK: use of globals: Ok for dumps since probably always globally dumped, not dumped per section computed or something similar
  setfdivb(&divb, GLOBALPOINT(pdump), GLOBALPOINT(pstagdump), GLOBALPOINT(udump), GLOBALPOINT(Bhatdump), i, j, k); // udump also set externally GODMARK

  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns


  //static
  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);  //ti
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);  //tj
    ftemp=(FTYPE)(k+startpos[3]);
    myset(datatype,&ftemp,0,1,writebuf);  //tk
  }
  myset(datatype,X,1,3,writebuf);  //x1, x2, x3
  myset(datatype,V,1,3,writebuf);  //r, h, ph
  // 9

  ////////////////////////
  //
  // rest dynamic

  // primitives
  // must use PDUMPLOOP() since may be any order unlike NPR loop
  PDUMPLOOP(pliter,pl) myset(datatype,&(GLOBALMACP0A1(pdump,i,j,k,pl)),0,1,writebuf); // NPRDUMP  //rho u v1 v2 v3 B1 B2 B3 ??

  ////////////
  //
  // output some EOS stuff since in general not simple function of rho0,u
  rho = GLOBALMACP0A1(pdump,i,j,k,RHO);
  u = GLOBALMACP0A1(pdump,i,j,k,UU);


  pressure = pressure_rho0_u_simple(i,j,k,loc,rho,u);
  cs2 = cs2_compute_simple(i,j,k,loc,rho,u);
  Sden = compute_entropy_simple(i,j,k,loc,rho,u);
  
  myset(datatype,&pressure,0,1,writebuf); // 1 //p
  myset(datatype,&cs2,0,1,writebuf); // 1 //cs2
  myset(datatype,&Sden,0,1,writebuf); // 1 //Sden

  //////////////////////
  //
  // output the conserved quantities since not easily inverted and at higher order aren't invertable from point primitives
  // PLOOP() used since conserved quantities always fill full PLOOP, while PDUMPLOOP is for primitives that may be duplicate among conserved quantities
  PLOOP(pliter,pl) myset(datatype,&(GLOBALMACP0A1(udump,i,j,k,pl)),0,1,writebuf); // NPR //U0 U1 U2 U3 U4 U5 U6 U7 ??
  myset(datatype,&divb,0,1,writebuf); // 1 //divb

  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(qptr->ucon[pl]),0,1,writebuf); //uu0 uu1 uu2 uu3 
  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(qptr->ucov[pl]),0,1,writebuf); //ud0 ud1 ud2 ud3
  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(qptr->bcon[pl]),0,1,writebuf); //bu0 bu1 bu2 bu3 
  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(qptr->bcov[pl]),0,1,writebuf); //bd0 bd1 bd2 bd3 
  // 4*4
    
  myset(datatype,&vmin[1],0,1,writebuf);  //v1m v1p v2m v2p v3m v3p
  myset(datatype,&vmax[1],0,1,writebuf);
  myset(datatype,&vmin[2],0,1,writebuf);
  myset(datatype,&vmax[2],0,1,writebuf);
  myset(datatype,&vmin[3],0,1,writebuf);
  myset(datatype,&vmax[3],0,1,writebuf);
  // 6

  // one static term
  myset(datatype,&(ptrgeom->gdet),0,1,writebuf); // 1 //gdet  //end of default read


  if(CALCFARADAYANDCURRENTS){ // NIM*2+6*2 = 8+12=20
    // updated 11/16/2003
    // new 10/23/2003
    // current density 
    lower_vec(GLOBALMAC(jcon,i,j,k),ptrgeom,jcov); 
    myset(datatype,GLOBALMAC(jcon,i,j,k),0,NDIM,writebuf); // (NDIM)
    myset(datatype,jcov,0,NDIM,writebuf);// (NDIM)
    // faraday (2*6)
    lowerf(GLOBALMAC(fcon,i,j,k),ptrgeom,fcov);
    myset(datatype,GLOBALMAC(fcon,i,j,k),0,NUMFARADAY,writebuf); //  (6)
    myset(datatype,fcov,0,NUMFARADAY,writebuf); // (6)
  }


  // DEBUG: Also add +3 to numcolumns for this to work
  if(0){
    if(FLUXB==FLUXCTSTAG) myset(datatype,GLOBALMAC(pstagdump,i,j,k),B1,3,writebuf);
    else{
      FTYPE plblob[NPR]={0};
      myset(datatype,plblob,B1,3,writebuf);
    }
  }


  if(FLUXB==FLUXCTSTAG && 0){ // DEBUG (change corresponding code in dump.c)
    // uses jrdp3dudebug in gtwod.m that assumes CALCFARADAYANDCURRENTS==0
    for(l=1;l<=COMPDIM;l++) myset(datatype,GLOBALMACP1A0(pl_ctdump,l,i,j,k),0,NPR2INTERP,writebuf); // 3*8 = 24
    for(l=1;l<=COMPDIM;l++) myset(datatype,GLOBALMACP1A0(pr_ctdump,l,i,j,k),0,NPR2INTERP,writebuf); // 3*8 = 24
    myset(datatype,GLOBALMAC(pstagdump,i,j,k),0,NPR,writebuf); // 8 // GODMARK: use of globals
    //    for(dir=1;dir<=COMPDIM;dir++) for(pl=B1;pl<=B3;pl++) for(n=0;n<=1;n++) myset(datatype,&GLOBALMACP3A0(pbcorninterp,dir,pl,n,i,j,k),0,1,writebuf); // 3*3*2 = 18
    for(dir=1;dir<=COMPDIM;dir++) for(pl=1;pl<=COMPDIM;pl++) for(n=0;n<NUMCS+1;n++) for(o=0;o<NUMCS;o++) myset(datatype,&GLOBALMACP1A3(pvbcorninterp,dir,i,j,k,pl,n,o),0,1,writebuf); // 3*3*(2+1)*2 = 54 (was 36)
  }

  return (0);
}






int debugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping debug dump# %ld ... ",dump_cnt);

  whichdump=DEBUGDUMPTYPE;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/debug");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,debug_content)>=1) return(1);

  trifprintf("end dumping debug# %ld ... ",dump_cnt);

  return(0);

}


extern void set_debug_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(DODEBUG){
    *numcolumns=2*NUMFAILFLOORFLAGS*NUMTSCALES;
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;

}

int debug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  // could also make everything FTYPE and convert like for normal i,j dump file


  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,GLOBALMAC(failfloorcount,i,j,k),0,2*NUMTSCALES*NUMFAILFLOORFLAGS,writebuf);
    
  return(0);
}





int enodebugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping enodebug dump# %ld ... ",dump_cnt);

  whichdump=ENODEBUGDUMPTYPE;
  //  datatype=MPI_FTYPE;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/enodebug");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,eno_dump_header,enodebug_content)>=1) return(1);

  trifprintf("end dumping enodebug# %ld ... ",dump_cnt);

  return(0);

}



void set_enodebug_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(DOENODEBUG){
    //*numcolumns=NUMENODEBUGS;
    *numcolumns=(3-1)* NUMENOINTERPTYPES * (NPR-4) * NUMENODEBUGS;  //SASMARK2
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;


}


int enodebug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  // could also make everything FTYPE and convert like for normal i,j dump file
  //myset(datatype,GLOBALMAC(enodebugarray,i,j,k),0,3*NUMENOINTERPTYPES*NPR*NUMENODEBUGS,writebuf);
  myset(datatype,GLOBALMAC(enodebugarray,i,j,k),0,(3-1)*NUMENOINTERPTYPES*(NPR-4)*NUMENODEBUGS,writebuf);  //atch corrected
    
  return(0);
}


int eno_dump_header(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE *headerptr)
{
  int dump_header_general(int whichdump, int whichdumpversion, int numcolumns, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr);
  int retval;

  retval=dump_header_general(whichdump, whichdumpversion, numcolumns, realnstep,dt, bintxt, headerptr);

  return(retval);
}       








int avgdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
        


  trifprintf("begin dumping avgdump# %ld ... ",dump_cnt);

  whichdump=AVG1DUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/avg");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,avg_content)>=1) return(1);

  trifprintf("end dumping avgdump# %ld ... ",dump_cnt);


  return(0);

}


void set_avg_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  // 36+29+8*2+4*2+2+12*2+96*2=339
  *numcolumns=3*3 + 1 + NUMNORMDUMP  // (6+1+29=36)
    + NUMNORMDUMP // |normal terms| (29)
#if(CALCFARADAYANDCURRENTS)
    + NDIM*2 // jcon/jcov (8)
    + NDIM*2 // |jcon|/|jcov| (8)
#endif
    + NDIM*2 // massflux/|massflux|
    + NUMOTHER*2 // other stuff and fabs of each
#if(CALCFARADAYANDCURRENTS)
    +6*2 // fcon/fcov (12)
    +6*2 // |fcon|,|fcov| (12)
#endif
    +7*16 // Tud all 7 parts, all 4x4 terms (112)
    +7*16 // |Tud| all 7 parts, all 4x4 terms (112)
    ;


  if(DOAVG2){
    *numcolumns-=224;
  }

  // Version number:
  *numversion=0;



}



int avg_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl = 0, l = 0, col = 0;
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int loc=CENT;


  bl_coord_ijk_2(i, j, k, loc, X,V);
  get_geometry(i, j, k, loc, ptrgeom);

  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(k+startpos[3]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);


  myset(datatype,&(ptrgeom->gdet),0,1,writebuf);

  // now do time average stuff
  myset(datatype,GLOBALMAC(normalvarstavg,i,j,k),0,NUMNORMDUMP,writebuf);
  myset(datatype,GLOBALMAC(anormalvarstavg,i,j,k),0,NUMNORMDUMP,writebuf);

#if(CALCFARADAYANDCURRENTS)
  myset(datatype,GLOBALMAC(jcontavg,i,j,k),0,NDIM,writebuf);
  myset(datatype,GLOBALMAC(jcovtavg,i,j,k),0,NDIM,writebuf);
  myset(datatype,GLOBALMAC(ajcontavg,i,j,k),0,NDIM,writebuf);
  myset(datatype,GLOBALMAC(ajcovtavg,i,j,k),0,NDIM,writebuf);
#endif
  myset(datatype,GLOBALMAC(massfluxtavg,i,j,k),0,NDIM,writebuf);
  myset(datatype,GLOBALMAC(amassfluxtavg,i,j,k),0,NDIM,writebuf);

  myset(datatype,GLOBALMAC(othertavg,i,j,k),0,NUMOTHER,writebuf);
  myset(datatype,GLOBALMAC(aothertavg,i,j,k),0,NUMOTHER,writebuf);

#if(CALCFARADAYANDCURRENTS)
  myset(datatype,GLOBALMAC(fcontavg,i,j,k),0,NUMFARADAY,writebuf);
  myset(datatype,GLOBALMAC(fcovtavg,i,j,k),0,NUMFARADAY,writebuf);
  myset(datatype,GLOBALMAC(afcontavg,i,j,k),0,NUMFARADAY,writebuf);
  myset(datatype,GLOBALMAC(afcovtavg,i,j,k),0,NUMFARADAY,writebuf);
#endif

#if(DOAVG2==0)
  myset(datatype,GLOBALMAC(tudtavg,i,j,k),0,NUMSTRESSTERMS,writebuf);
  myset(datatype,GLOBALMAC(atudtavg,i,j,k),0,NUMSTRESSTERMS,writebuf);
#endif

  return(0);

}


int avg2dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];



  trifprintf("begin dumping avg2dump# %ld ... ",dump_cnt);

  whichdump=AVG2DUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/avg2");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");

  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,avg2_content)>=1) return(1);

  trifprintf("end dumping avg2dump# %ld ... ",dump_cnt);


  return(0);

}

void set_avg2_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{


  if(DOAVG2){
    *numcolumns=10 + 224; // otherwise doesn't exist so don't need to set
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;



}

int avg2_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl = 0, l = 0, col = 0;
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int loc=CENT;

  bl_coord_ijk_2(i, j, k, loc, X,V);
  get_geometry(i, j, k, loc, ptrgeom);
  // if you change # of outputted vars, remember to change numcolumns above

  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(k+startpos[3]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);

  myset(datatype,&(ptrgeom->gdet),0,1,writebuf);
  // 10

  myset(datatype,GLOBALMAC(tudtavg,i,j,k),0,NUMSTRESSTERMS,writebuf);
  myset(datatype,GLOBALMAC(atudtavg,i,j,k),0,NUMSTRESSTERMS,writebuf);
  // 112*2

  // total=10+112*2=234

  return(0);
}


int gdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
        

  trifprintf("begin dumping gdump# %ld ... ",dump_cnt);

  whichdump=GRIDDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/gdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");

  // dump_cnt==-1 means no file number
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,gdump_content)>=1) return(1);

  trifprintf("end dumping gdump# %ld ... ",dump_cnt);

  return(0);
}



extern void set_gdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{


  // 205+4+4*4 currently
  //*numcolumns=3*3+NDIM*NDIM*NDIM+NPG*NDIM*NDIM*2+NPG+4+4*4;
  //NPG was replaced with unity in order to avoid excessive dumping of info (only center info now available)
  *numcolumns=3*3  +   NDIM*NDIM*NDIM  +   1*NDIM*NDIM*2   +   1  +  NDIM   +   NDIM*NDIM;
  //t^i x^i V^i,     \Gamma^\mu_{\nu\tau},     g^{\mu\nu} g_{\mu\nu}, \sqrt{-g}, \gamma_\mu, dx^\mu/dx^\nu

  // Version number:
  *numversion=0;

}



int gdump_content(int i, int j, int k, MPI_Datatype datatype, void *writebuf)
{
  int pl = 0, l = 0, m = 0, n = 0, col = 0;
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;
  FTYPE *ptrftemp;
  FTYPE dxdxp[NDIM][NDIM];
  int myii,myjj,mykk;
  LOCALMETRICTEMPVARS;
  int loc=CENT;
  FTYPE generalmatrixlower[NDIM][NDIM];
  FTYPE generalmatrixupper[NDIM][NDIM];
  int jj,kk;



  bl_coord_ijk_2(i, j, k, loc, X, V);
  dxdxprim_ijk(i, j, k, loc, dxdxp);



  ftemp=(FTYPE)(i+startpos[1]);
  myset(datatype,&ftemp,0,1,writebuf);
  ftemp=(FTYPE)(j+startpos[2]);
  myset(datatype,&ftemp,0,1,writebuf);
  ftemp=(FTYPE)(k+startpos[3]);
  myset(datatype,&ftemp,0,1,writebuf);
  // 3
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);
  // 6




#if(MCOORD!=CARTMINKMETRIC)
  myii=i;
  myjj=j;
  mykk=k;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif


  ptrftemp=(FTYPE*)(&GLOBALMETMACP0A3(conn,myii,myjj,mykk,0,0,0));
  myset(datatype,ptrftemp,0,NDIM*NDIM*NDIM,writebuf);

  // get local metric quantities for this loc,i,j,k
  GETLOCALMETRIC(loc,myii,myjj,mykk);

  DLOOP(jj,kk){
    generalmatrixlower[jj][kk]=localgcov[GIND(jj,kk)];
    generalmatrixupper[jj][kk]=localgcon[GIND(jj,kk)];
  }

    
  ptrftemp=&generalmatrixupper[0][0];
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);
  ptrftemp=&generalmatrixlower[0][0];
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);
  ptrftemp=&localgdet[0];
  myset(datatype,ptrftemp,0,1,writebuf);
  //ptrftemp=(FTYPE*)(&localgdetvol); // can take a peek if GDETVOLDIFF==1
  //  myset(datatype,ptrftemp,0,1,writebuf);

    
  ptrftemp=&GLOBALMETMACP0A1(conn2,myii,myjj,mykk,0);
  myset(datatype,ptrftemp,0,NDIM,writebuf);

  // 4*4
  ptrftemp=&dxdxp[0][0];
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);


  return(0);

}



int fieldlinedump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
        


  trifprintf("begin dumping fieldlinedump# %ld ... ",dump_cnt);

  whichdump=FIELDLINEDUMPTYPE;
  datatype=MPI_FLOAT; // don't need good precision
  strcpy(fileprefix,"dumps/fieldline");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // MIXEDOUTPUT means text header and forced binary data
  if(dump_gen(WRITEFILE,dump_cnt,MIXEDOUTPUT,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,fieldline_content)>=1) return(1);

  trifprintf("end dumping fieldlinedump# %ld ... ",dump_cnt);


  return(0);

}


extern void set_fieldline_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{
  if(DOFIELDLINE){
    *numcolumns=NUMFIELDLINEQUANTITIES;
  }
  else *numcolumns=0;


  // Version number:
  *numversion=0;

}


int fieldline_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl = 0, l = 0, col = 0;
  struct of_state q;
  //FTYPE U[NPR];
  FTYPE FL[NPR];
  // must be same precision as written content
  float ftemp;
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  int loc=CENT;


  //////////////
  //
  // some calculations
  //

  // if failed, then data output for below invalid, but columns still must exist    
  get_geometry(i, j, k, loc, ptrgeom);
  if (!failed) {
    if (get_state(GLOBALMAC(pdump,i,j,k), ptrgeom, &q) >= 1)
      FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
  }
  else {// do a per zone check, otherwise set to 0
    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
    if (get_state(GLOBALMAC(pdump,i,j,k), ptrgeom, &q) >= 1){
      for (pl = 0; pl < NDIM; pl++)
        q.ucon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
        q.ucov[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
        q.bcon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
        q.bcov[pl]=0;
    }
    whocalleducon=0; // return to normal state
    
  }

  MYFUN(primtoflux(UDIAG,GLOBALMAC(pdump,i,j,k), &q, RR, ptrgeom, FL),"step_ch.c:fluxcalc()", "primtoflux_calc() dir=1/2 l", RR);


  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns



  ////////////////////
  //
  // 2 various things

  // rho (for various things)
  ftemp=(float)GLOBALMACP0A1(pdump,i,j,k,RHO);
  myset(datatype,&ftemp,0,1,writebuf);

  // u (for various things)
  ftemp=(float)GLOBALMACP0A1(pdump,i,j,k,UU);
  myset(datatype,&ftemp,0,1,writebuf);


  //////////////////////
  //
  // 2 things for jet/energy per baryon at infinity

  // -u_t (-hu_t can be found from this and rho/u/p above)
  ftemp=(float)(-q.ucov[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // -T^t_t/(rho u^t)
  //  ftemp=(float)(-U[UU]/(ptrgeom->gdet * GLOBALMACP0A1(pdump,i,j,k,RHO)*q.ucon[TT]));
  //myset(datatype,&ftemp,0,1,writebuf);

  // -T^r_t/(rho u^r)
  if(q.ucon[RR]!=0.0){
    ftemp=(float)(-FL[UU]/(ptrgeom->gdet * GLOBALMACP0A1(pdump,i,j,k,RHO)*q.ucon[RR]));
  }
  else ftemp=0.0;
  myset(datatype,&ftemp,0,1,writebuf);


  // 1 extra thing

  // u^t
  ftemp=(float)(q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);


  ///////////////////////////
  //
  // 6 things for the field line stuff

  // v^r [ in grid frame]
  ftemp=(float)(q.ucon[1]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // v^\theta
  ftemp=(float)(q.ucon[2]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // v^\phi
  ftemp=(float)(q.ucon[3]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // B^r
  ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,B1));
  myset(datatype,&ftemp,0,1,writebuf);

  // B^\theta
  ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,B2));
  myset(datatype,&ftemp,0,1,writebuf);

  // B^\phi
  ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,B3));
  myset(datatype,&ftemp,0,1,writebuf);

#if( FIELDLINEGDETB == 1)
  //it is useful to have access to gdet*B^i at cell faces directly for plotting field lines
  ftemp=(float)(GLOBALMACP0A1(udump,i,j,k,B1));
  myset(datatype,&ftemp,0,1,writebuf);

  ftemp=(float)(GLOBALMACP0A1(udump,i,j,k,B2));
  myset(datatype,&ftemp,0,1,writebuf);

  ftemp=(float)(GLOBALMACP0A1(udump,i,j,k,B3));
  myset(datatype,&ftemp,0,1,writebuf);
#endif
  
  // see grmhd-dualfcon2omegaf.nb
  // below can be obtained from above set of v and B
  // \Omega_F_1
  //  ftemp=(float)(v3-B3*v2/(B2+SMALL));
  //myset(datatype,&ftemp,0,1,writebuf);

  // \Omega_F_2
  //ftemp=(float)(v3-B3*v1/(B1+SMALL));
  // myset(datatype,&ftemp,0,1,writebuf);

  return(0);

}




int dissdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dissdump# %ld ... ",dump_cnt);

  whichdump=DISSDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dissdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dissdump_content)>=1) return(1);

  trifprintf("end dumping dissdump# %ld ... ",dump_cnt);


  return(0);
  
}


void set_dissdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(DODISS){
    *numcolumns=NUMDISSFUNPOS;
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;

}



int dissdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,&GLOBALMAC(dissfunpos,i,j,k),0,NUMDISSFUNPOS,writebuf);

  return (0);
}



int dumpother(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dumpother# %ld ... ",dump_cnt);

  whichdump=OTHERDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dumpother");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dumpother_content)>=1) return(1);

  /////// output the symmetry information to the fail file
  //writesyminfo();
  ///////

  trifprintf("end dumping dumpother# %ld ... ",dump_cnt);


  return(0);
  
}


void set_dumpother_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(DODUMPOTHER){ // panalytic + numpother quantities
    *numcolumns=NPR+NUMPOTHER;
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;


}


int dumpother_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl,pliter;

#if(ANALYTICMEMORY==0)
  dualfprintf(fail_file,"Need to set ANALYTICMEMORY==1 for using analytical memory in dumpother_content()\n");
  myexit(9857652);
#endif


  myset(datatype,&GLOBALMAC(panalytic,i,j,k),0,NPR,writebuf);

  for(pl=0;pl<NUMPOTHER;pl++){
    myset(datatype,&GLOBALMACP1A0(pother,pl,i,j,k),0,1,writebuf);
  }

  return (0);
}




int fluxdumpdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping fluxdump# %ld ... ",dump_cnt);

  whichdump=FLUXDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/fluxdump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,fluxdump_header,fluxdump_content)>=1) return(1);

  trifprintf("end dumping fluxdump# %ld ... ",dump_cnt);


  return(0);
  
}

void set_fluxdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(FLUXDUMP){ // dU, flux, and ppprimitives for flux
    *numcolumns=NUMFLUXDUMP;
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;


}



int fluxdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl,pliter;

  myset(datatype,&GLOBALMAC(fluxdump,i,j,k),0,NUMFLUXDUMP,writebuf);

  return (0);
}



// dump stuff related specially to EOS
int eosdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping eosdump# %ld ... ",dump_cnt);

  whichdump=EOSDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/eosdump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,eosdump_content)>=1) return(1);

  /////// output the symmetry information to the fail file
  //writesyminfo();
  ///////

  trifprintf("end dumping eosdump# %ld ... ",dump_cnt);


  return(0);
  
}


void set_eosdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  //#if(WHICHEOS==KAZFULL)
  // all EOSs output same size data so uniform format
  // otherwise have to also put this condition in dump.c when outputting so don't overwrite memory!
  *numcolumns=MAXPARLIST+1+MAXNUMEXTRAS+MAXPROCESSEDEXTRAS; // 1 is temperature
  //#else
  //  *numcolumns=0;
  //#endif

  // Version number:
  *numversion=0;


}


int eosdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  FTYPE rho,u;
  FTYPE height,ye,ynu,temp;
  FTYPE extras[MAXNUMEXTRAS];
  FTYPE processed[MAXPROCESSEDEXTRAS];
  int numextras;
  int extraiter;
  int numparms;
  FTYPE parlist[MAXPARLIST];
  int loc=CENT;

  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns


  ////////////
  //
  // output some EOS stuff since in general not simple function of rho0,u
  rho = GLOBALMACP0A1(pdump,i,j,k,RHO);
  u = GLOBALMACP0A1(pdump,i,j,k,UU);

  // compute EOS stuff
  get_EOS_parms_simple(&numparms, i,j,k,loc, parlist);
  temp = compute_temp_simple(i,j,k,loc,rho,u);
  // get extra EOS stuff
  get_extrasprocessed_simple(1, i, j, k, loc, GLOBALMAC(pdump,i,j,k), extras, processed);
  //  compute_allextras(0,rho,u,&numextras,extras);
  

  // write EOS stuff
  // write out all stuff for all EOSs so uniform format to read for any run (otherwise need version information in file)
  myset(datatype,&parlist,0,MAXPARLIST,writebuf);
  myset(datatype,&temp,0,1,writebuf); // 1
  // write extras EOS stuff
  myset(datatype,extras,0,MAXNUMEXTRAS,writebuf); // numextras
  myset(datatype,processed,0,MAXPROCESSEDEXTRAS,writebuf); // numprocessedextras


  return (0);
}





int vpotdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping vpotdump# %ld ... ",dump_cnt);

  whichdump=VPOTDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/vpotdump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,vpotdump_content)>=1) return(1);

  trifprintf("end dumping vpotdump# %ld ... ",dump_cnt);


  return(0);
  
}



void set_vpotdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(DOVPOTDUMP){
    *numcolumns=NUMVPOTDUMP;
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;



}


int vpotdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int jj;

  // NOTEMARK: see also restart.c since this is added to restart 
  for(jj=0;jj<NUMVPOTDUMP;jj++){
    myset(datatype,&GLOBALMACP1A0(vpotarraydump,jj,i,j,k),0,1,writebuf); // 1 each
  }

  return (0);
}








int failfloordudump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping failfloordudump# %ld ... ",dump_cnt);

  whichdump=FAILFLOORDUDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/failfloordudump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,failfloordudump_content)>=1) return(1);

  trifprintf("end dumping failfloordudump# %ld ... ",dump_cnt);


  return(0);
  
}



void set_failfloordudump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion)
{

  if(DOFLOORDIAG){
    *numcolumns=NPR;
  }
  else *numcolumns=0;

  // Version number:
  *numversion=0;



}


int failfloordudump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl;

  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,&GLOBALMAC(failfloordu,i,j,k),0,NPR,writebuf); // NPR
  
  return (0);
}









// fake dump so can push out data in case still in MPI=2 delayed writing buffer
int fakedump(long dump_cnt)// arg not used
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping fakedump");

  whichdump=FAKEDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/fakedump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,-1,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,fakedump_header,fakedump_content)>=1) return(1);

  trifprintf("end dumping fakedump");


  return(0);
  
}



int fakedump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int jj;

  // blank

  return (0);
}


int fakedump_header(int whichdump, int whichdumpversion, int numcolumns,int bintxt, FILE *headerptr)
{
  int fake;

  fake=0;


  if(bintxt==BINARYOUTPUT){
    fwrite(&fake,sizeof(int),1,headerptr);
  }
  else{
#if(REALTYPE==DOUBLETYPE)
    fprintf(headerptr, "DONE: %d\n",fake);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "DONE: %d\n",fake);
#endif
  }
  fflush(headerptr);
  return(0);
}       

