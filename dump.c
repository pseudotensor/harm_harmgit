#include "decs.h"

/*! \file dump.c
     \brief file dumping

     // mpi.h has following datatypes corresponding to the C types
     // pick one per dump file. (no per column types yet)
     // same for image.c
     // same for restart.c


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



// Follow these steps to create a new dump file

   1) defs.???.h: create the storage variable if needed

   2) set_array???.c : shift variable if necessary (follow examples)

   3) global.nondepmnemonics.h : change NUMDUMPTYPES and add label and add to MYDUMPNAMELIST

   4) global.???.h change NUMDUMPTYPES and MYDUMPNAMELIST and add another entry if want separate timing of output

   5) dump.c : add dnumcolumns[LABEL]=NUMCOLUMNS where NUMCOLUMNS is number of entries in dump file

   6) dump.c : follow examples from here (dump() uses dump_header() and dump_content()).  One must define the header and content function and the wrapper (3 functions) or use an existing header function

   7) diag.c : follow example of "dump", dumpc, tlastdump, etc.  E.g. follow DISSDUMPTYPE and repeat for new entry.  Two parts of code.

   8) init.c : DTdumpgen, dumpcntgen, and other things.  Unless default over array works.

   9) global.dump.h : add global prototypes

*/






/// initialize/get-ready for dumping
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

/// output all npr-type info to files
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



/// setup number of columns per dump file (see dumpgen.c or dump.c for how used)
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
  extern void set_raddump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_vpotdump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_failfloordudump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_dissmeasuredump_content_dnumcolumns_dnumversion(int *numcolumns, int *numversion);
  extern void set_fluxsimpledump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion);

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
  // raddump
  set_raddump_content_dnumcolumns_dnumversion(&dnumcolumns[RADDUMPTYPE],&dnumversion[RADDUMPTYPE]);
  // vpotdump
  set_vpotdump_content_dnumcolumns_dnumversion(&dnumcolumns[VPOTDUMPTYPE],&dnumversion[VPOTDUMPTYPE]);
  // failfloordudump
  set_failfloordudump_content_dnumcolumns_dnumversion(&dnumcolumns[FAILFLOORDUDUMPTYPE],&dnumversion[FAILFLOORDUDUMPTYPE]);
  // dissmeasuredump
  set_dissmeasuredump_content_dnumcolumns_dnumversion(&dnumcolumns[DISSMEASUREDUMPTYPE],&dnumversion[DISSMEASUREDUMPTYPE]);
  // fluxsimpledump
  set_fluxsimpledump_content_dnumcolumns_dnumversion(&dnumcolumns[FLUXSIMPLEDUMPTYPE],&dnumversion[FLUXSIMPLEDUMPTYPE]);

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






/// Primary full dump file outputted usually not too frequently
int dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


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



/// output the symmetry information to the fail file; symmetrizes w.r.t. i == j
void writesyminfo( void )
{
  int i, j;

  for( i = 0; i < N1; i++ ) {
  }

}

/// Primary full dump file header
int dump_header(int whichdump, int whichdumpversion, int numcolumnsvar, int bintxt, FILE *headerptr)
{
  int dump_header_general(int whichdump, int whichdumpversion, int numcolumnsvar, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr);
  int retval;

  retval=dump_header_general(whichdump, whichdumpversion, numcolumnsvar, realnstep,dt, bintxt, headerptr);

  return(retval);

}

/// fluxdump header
int fluxdump_header(int whichdump, int whichdumpversion, int numcolumnsvar, int bintxt, FILE *headerptr)
{
  int dump_header_general(int whichdump, int whichdumpversion, int numcolumnsvar, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr);
  int retval;

  retval=dump_header_general(whichdump, whichdumpversion, numcolumnsvar, fluxdumprealnstep,fluxdumpdt, bintxt, headerptr);

  return(retval);

}

/// general dump header (used for full dumps and other dumps)
int dump_header_general(int whichdump, int whichdumpversion, int numcolumnsvar, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr)
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
    // SUPERGODMARK: For jon_interp to not need init.h specific stuff, need to at least add MCOORD so knows if spherical polar coordinates or not.  defcoord is insufficient in general.
    // Maybe also: COMPDIM, WHICHVEL, WHICHEOM, REMOVERESTMASSFROMUU, EOMTYPE, NPR, DOENTROPY, 
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
    fwrite(&numcolumnsvar,sizeof(int),1,headerptr);

    // init.h type stuff related to physical settings or what data is outputted
    int var;
    var=TRACKVPOT; fwrite(&var,sizeof(int),1,headerptr);
    var=MCOORD; fwrite(&var,sizeof(int),1,headerptr);
    var=DODISS; fwrite(&var,sizeof(int),1,headerptr);
    var=DOEVOLVEMETRIC; fwrite(&var,sizeof(int),1,headerptr);
    var=WHICHVEL; fwrite(&var,sizeof(int),1,headerptr);
    var=WHICHEOM; fwrite(&var,sizeof(int),1,headerptr);
    var=REMOVERESTMASSFROMUU; fwrite(&var,sizeof(int),1,headerptr);
    var=RELTYPE; fwrite(&var,sizeof(int),1,headerptr);
    var=EOMTYPE; fwrite(&var,sizeof(int),1,headerptr);
    var=WHICHEOS; fwrite(&var,sizeof(int),1,headerptr);
    var=DOENTROPY; fwrite(&var,sizeof(int),1,headerptr);
    var=WHICHENTROPYEVOLVE; fwrite(&var,sizeof(int),1,headerptr);
    var=CALCFARADAYANDCURRENTS; fwrite(&var,sizeof(int),1,headerptr);
    var=DOPOLEDEATH; fwrite(&var,sizeof(int),1,headerptr);
    var=DOPOLESMOOTH; fwrite(&var,sizeof(int),1,headerptr);
    var=DOPOLEGAMMADEATH; fwrite(&var,sizeof(int),1,headerptr);
    var=IF3DSPCTHENMPITRANSFERATPOLE; fwrite(&var,sizeof(int),1,headerptr);
    var=EOMRADTYPE; fwrite(&var,sizeof(int),1,headerptr);
    var=WHICHRADSOURCEMETHOD; fwrite(&var,sizeof(int),1,headerptr);
    var=OUTERDEATH; fwrite(&var,sizeof(int),1,headerptr);
    var=OUTERDEATHRADIUS; fwrite(&var,sizeof(int),1,headerptr);
    

       
  }
  else{

#define DUMPHEADERLIST tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2], dx[3], localrealnstep,gam,a,R0,Rin,Rout,hslope,localdt,defcoord,MBH,QBH,EP3,THETAROT,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumnsvar,  TRACKVPOT    ,MCOORD    ,DODISS    ,DOEVOLVEMETRIC    ,WHICHVEL    ,WHICHEOM    ,REMOVERESTMASSFROMUU    ,RELTYPE    ,EOMTYPE    ,WHICHEOS    ,DOENTROPY    ,WHICHENTROPYEVOLVE    ,CALCFARADAYANDCURRENTS    ,DOPOLEDEATH    ,DOPOLESMOOTH    ,DOPOLEGAMMADEATH    ,IF3DSPCTHENMPITRANSFERATPOLE    ,EOMRADTYPE    ,WHICHRADSOURCEMETHOD    ,OUTERDEATH    ,OUTERDEATHRADIUS

#if(REALTYPE==DOUBLETYPE)
    fprintf(headerptr, "%21.15g %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %ld %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d %21.15g %21.15g %21.15g %21.15g %d %d %d %d %d %d %d %d %d  %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %21.15g\n", DUMPHEADERLIST);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "%31.25Lg %d %d %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %ld %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d %d %d %d %d %d %d %d %d  %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %31.25Lg\n", DUMPHEADERLIST);
#endif
  }
  fflush(headerptr);
  return(0);
} 



/// number of columns for dump file
void set_dump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{



  if(DOMAINDUMPDIAG){
    // always NPRDUMP
    if(GAMMIEDUMP)  *numcolumnsvar=2*3 + NPRDUMP+NPR + 3 + 1 + 4 * NDIM + 6 + 1
#if(CALCFARADAYANDCURRENTS)
                      + NDIM*2
                      + 2*6
#endif
                      ;
    else{
      *numcolumnsvar=3*3 + NPRDUMP + 3 + (nprend+1) + 1 + 4 * NDIM + 6 + 1  //replace NPR -> (nprend+1) since nprend, not NPR, controls dumping.  Fixes: DOEXTRAINTERP=1 case
#if(CALCFARADAYANDCURRENTS)
        + NDIM*2
        + 2*6
#endif
        ;    // 61 total if also doing currents and faraday, 41 otherwise

      if(FLUXB==FLUXCTSTAG && 0){ // DEBUG (change corresponding code in dump.c)
        *numcolumnsvar+= NPR2INTERP*COMPDIM*2 + NPR + COMPDIM*3*2 + COMPDIM*3*2*2;
      }
    }
  }
  else{
    *numcolumnsvar=0;
  }

  // Version number:
  *numversion=0;

}



/// Primary full dump file contents
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
    if(N1NOT1){
      if (vchar_all(GLOBALMAC(pdump,i,j,k), qptr, 1, ptrgeom, &vmax[1], &vmin[1],&ignorecourant) >= 1) FAILSTATEMENT("dump.c:dump()", "vchar_all() dir=1or2", 1);
    }
    else{
      vmax[1]=vmin[1]=0;
    }
    if(N2NOT1){
      if (vchar_all(GLOBALMAC(pdump,i,j,k), qptr, 2, ptrgeom, &vmax[2], &vmin[2],&ignorecourant) >= 1) FAILSTATEMENT("dump.c:dump()", "vchar_all() dir=1or2", 2);
    }
    else{
      vmax[2]=vmin[2]=0;
    }
    if(N3NOT1){
      if (vchar_all(GLOBALMAC(pdump,i,j,k), qptr, 3, ptrgeom, &vmax[3], &vmin[3],&ignorecourant) >= 1) FAILSTATEMENT("dump.c:dump()", "vchar_all() dir=1or2", 3);
    }
    else{
      vmax[3]=vmin[3]=0;
    }

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
    if (vchar_all(GLOBALMAC(pdump,i,j,k), qptr, 1, ptrgeom, &vmax[1], &vmin[1],&ignorecourant) >= 1){
      vmax[1]=vmin[1]=0;
    }
    
    if (vchar_all(GLOBALMAC(pdump,i,j,k), qptr, 2, ptrgeom, &vmax[2], &vmin[2],&ignorecourant) >= 1){
      vmax[2]=vmin[2]=0;
    }

    if (vchar_all(GLOBALMAC(pdump,i,j,k), qptr, 3, ptrgeom, &vmax[3], &vmin[3],&ignorecourant) >= 1){
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
  // if you change # of outputted vars, remember to change numcolumnsvar


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


  // DEBUG: Also add +3 to numcolumnsvar for this to work
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





/// debug dump
int debugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


  trifprintf("begin dumping debug dump# %ld ... ",dump_cnt);

  whichdump=DEBUGDUMPTYPE;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/debugdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,debug_content)>=1) return(1);

  trifprintf("end dumping debug# %ld ... ",dump_cnt);

  return(0);

}

/// debug dump content number
extern void set_debug_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DODEBUG && DODEBUGDUMP){ //different than doing just DODEBUG
    *numcolumnsvar=2*NUMFAILFLOORFLAGS*NUMTSCALES;
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;

}

/// debug dump contents
int debug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  // could also make everything FTYPE and convert like for normal i,j dump file


  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,GLOBALMAC(failfloorcount,i,j,k),0,2*NUMTSCALES*NUMFAILFLOORFLAGS,writebuf);

  
  return(0);
}




/// eno debug dump
int enodebugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


  trifprintf("begin dumping enodebug dump# %ld ... ",dump_cnt);

  whichdump=ENODEBUGDUMPTYPE;
  //  datatype=MPI_FTYPE;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/enodebugdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,eno_dump_header,enodebug_content)>=1) return(1);

  trifprintf("end dumping enodebug# %ld ... ",dump_cnt);

  return(0);

}


/// eno debug dump content number
void set_enodebug_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOENODEBUG){
    //*numcolumnsvar=NUMENODEBUGS;
    *numcolumnsvar=(3-1)* NUMENOINTERPTYPES * (NPR-4) * NUMENODEBUGS;  //SASMARK2
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;


}

/// eno debug dump content
int enodebug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  // could also make everything FTYPE and convert like for normal i,j dump file
  //myset(datatype,GLOBALMAC(enodebugarray,i,j,k),0,3*NUMENOINTERPTYPES*NPR*NUMENODEBUGS,writebuf);
  myset(datatype,GLOBALMAC(enodebugarray,i,j,k),0,(3-1)*NUMENOINTERPTYPES*(NPR-4)*NUMENODEBUGS,writebuf);  //atch corrected
    
  return(0);
}

/// eno debug dump header
int eno_dump_header(int whichdump, int whichdumpversion, int numcolumnsvar, int bintxt, FILE *headerptr)
{
  int dump_header_general(int whichdump, int whichdumpversion, int numcolumnsvar, long localrealnstep, SFTYPE localdt, int bintxt, FILE *headerptr);
  int retval;

  retval=dump_header_general(whichdump, whichdumpversion, numcolumnsvar, realnstep,dt, bintxt, headerptr);

  return(retval);
} 







/// time-average data dump
int avgdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
 


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

/// time-average data dump content list
void set_avg_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOAVG){
    // 36+29+8*2+4*2+2+12*2+96*2=339
    *numcolumnsvar=3*3 + 1 + NUMNORMDUMP  // (6+1+29=36)
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
      *numcolumnsvar-=224;
    }
  }
  else{
    *numcolumnsvar=0;
  }

  // Version number:
  *numversion=0;



}


/// time-average data dump contents
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

/// 2nd time-average dump
int avg2dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};



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

/// 2nd time-average dump contents number
void set_avg2_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{


  if(DOAVG2){
    *numcolumnsvar=10 + 224; // otherwise doesn't exist so don't need to set
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;



}

/// 2nd time-average dump contents
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
  // if you change # of outputted vars, remember to change numcolumnsvar above

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



/// grid dump
int gdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
 

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


/// grid dump contents number
extern void set_gdump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOGDUMP){
    // 205+4+4*4 currently
    //*numcolumnsvar=3*3+NDIM*NDIM*NDIM+NPG*NDIM*NDIM*2+NPG+4+4*4;
    //NPG was replaced with unity in order to avoid excessive dumping of info (only center info now available)
    *numcolumnsvar=3*3  +   NDIM*NDIM*NDIM  +   1*NDIM*NDIM*2   +   1  +  NDIM   +   NDIM*NDIM;
    //t^i x^i V^i,     \Gamma^\mu_{\nu\tau},     g^{\mu\nu} g_{\mu\nu}, \sqrt{-g}, \gamma_\mu, dx^\mu/dx^\nu
  }
  else{
    *numcolumnsvar=0;
  }
  // Version number:
  *numversion=0;

}


/// grid dump contents
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


/// Primary quick(but sufficient) dump (originally used for showing field lines move)
int fieldlinedump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};
 


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

/// fieldline data dump content number
extern void set_fieldline_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{


/// dump.c's fieldlinedump()
/// CHANGES alot, make sure # is correct!
/// Add 4 radiation terms if doing radiation
#if( FIELDLINEGDETB == 1)
#define NUMFIELDLINEQUANTITIES (14-2 + NUMYFL*(DOYFL!=0) + (DOYL!=0) + (DOYNU!=0) + (1+NDIM+10)*(EOMRADTYPE!=EOMRADNONE))
/// rho, u, <yfl,yl,ynu>, u^t, v1,v2,v3,B1,B2,B3 <gdetB1,gdetB2,gdetB3>
/// radiation adds: vrad1,vrad2,vrad3
#else
#define NUMFIELDLINEQUANTITIES (11-2 + NUMYFL*(DOYFL!=0) + (DOYL!=0) + (DOYNU!=0) + (1+NDIM+10)*(EOMRADTYPE!=EOMRADNONE))
/// rho, u, <yfl,yl,ynu>, u^t, v1,v2,v3,B1,B2,B3
/// radiation adds: vrad1,vrad2,vrad3
#endif



  if(DOFIELDLINE){
    *numcolumnsvar=NUMFIELDLINEQUANTITIES;
  }
  else *numcolumnsvar=0;




  // Version number:
  *numversion=0;

}

/// fieldline data dump contents
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


  FTYPE *pr=GLOBALMAC(pdump,i,j,k);

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
      for (pl = 0; pl < NDIM; pl++) q.ucon[pl]=0;
      for (pl = 0; pl < NDIM; pl++) q.ucov[pl]=0;
      for (pl = 0; pl < NDIM; pl++) q.bcon[pl]=0;
      for (pl = 0; pl < NDIM; pl++) q.bcov[pl]=0;

      if(EOMRADTYPE!=EOMRADNONE){
        for (pl = 0; pl < NDIM; pl++) q.uradcon[pl]=0;
        for (pl = 0; pl < NDIM; pl++) q.uradcov[pl]=0;
      }
    }
    whocalleducon=0; // return to normal state
    
  }

  MYFUN(primtoflux(UDIAG,pr, &q, RR, ptrgeom, FL, NULL),"step_ch.c:fluxcalc()", "primtoflux_calc() dir=1/2 l", RR);


  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumnsvar



  ////////////////////
  //
  // 2 various things

  // rho (for various things)
  ftemp=(float)GLOBALMACP0A1(pdump,i,j,k,RHO);
  myset(datatype,&ftemp,0,1,writebuf); // 1

  // u (for various things)
  ftemp=(float)GLOBALMACP0A1(pdump,i,j,k,UU);
  myset(datatype,&ftemp,0,1,writebuf); // 1

  if(DOYFL){
    // Y_fl
    if(YFL1>=0){
      ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YFL1));
      myset(datatype,&ftemp,0,1,writebuf); // 1
    }
    if(YFL2>=0){
      ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YFL2));
      myset(datatype,&ftemp,0,1,writebuf); // 1
    }
    if(YFL3>=0){
      ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YFL3));
      myset(datatype,&ftemp,0,1,writebuf); // 1
    }
    if(YFL4>=0){
      ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YFL4));
      myset(datatype,&ftemp,0,1,writebuf); // 1
    }
    if(YFL5>=0){
      ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YFL5));
      myset(datatype,&ftemp,0,1,writebuf); // 1
    }
  }

  if(DOYL){
    // Y_l
    ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YL));
    myset(datatype,&ftemp,0,1,writebuf); // 1
  }

  if(DOYNU){
    // Y_nu
    ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,YNU));
    myset(datatype,&ftemp,0,1,writebuf); // 1
  }

  //////////////////////
  //
  // 2 things for jet/energy per baryon at infinity

  // -u_t (-hu_t can be found from this and rho/u/p above)
  //  ftemp=(float)(-q.ucov[0]);
  //  myset(datatype,&ftemp,0,1,writebuf); // 1

  // -T^t_t/(rho u^t)
  //  ftemp=(float)(-U[UU]/(ptrgeom->gdet * GLOBALMACP0A1(pdump,i,j,k,RHO)*q.ucon[TT]));
  //myset(datatype,&ftemp,0,1,writebuf);

  // -T^r_t/(rho u^r)
  //  if(q.ucon[RR]!=0.0){
  //    ftemp=(float)(-FL[UU]/(ptrgeom->gdet * GLOBALMACP0A1(pdump,i,j,k,RHO)*q.ucon[RR]));
  //  }
  //  else ftemp=0.0;
  //  myset(datatype,&ftemp,0,1,writebuf); // 1


  // 1 extra thing

  // u^t
  ftemp=(float)(q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf); // 1


  ///////////////////////////
  //
  // 6 things for the field line stuff

  // v^r [ in grid frame]
  ftemp=(float)(q.ucon[1]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf); // 1

  // v^\theta
  ftemp=(float)(q.ucon[2]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf); // 1

  // v^\phi
  ftemp=(float)(q.ucon[3]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf); // 1

  // B^r
  ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,B1));
  myset(datatype,&ftemp,0,1,writebuf); // 1

  // B^\theta
  ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,B2));
  myset(datatype,&ftemp,0,1,writebuf); // 1

  // B^\phi
  ftemp=(float)(GLOBALMACP0A1(pdump,i,j,k,B3));
  myset(datatype,&ftemp,0,1,writebuf); // 1

#if( FIELDLINEGDETB == 1)
  //it is useful to have access to gdet*B^i at cell faces directly for plotting field lines
  ftemp=(float)(GLOBALMACP0A1(udump,i,j,k,B1));
  myset(datatype,&ftemp,0,1,writebuf); // 1

  ftemp=(float)(GLOBALMACP0A1(udump,i,j,k,B2));
  myset(datatype,&ftemp,0,1,writebuf); // 1

  ftemp=(float)(GLOBALMACP0A1(udump,i,j,k,B3));
  myset(datatype,&ftemp,0,1,writebuf); // 1
#endif
  
  // see grmhd-dualfcon2omegaf.nb
  // below can be obtained from above set of v and B
  // \Omega_F_1
  //  ftemp=(float)(v3-B3*v2/(B2+SMALL));
  //myset(datatype,&ftemp,0,1,writebuf); // 1

  // \Omega_F_2
  //ftemp=(float)(v3-B3*v1/(B1+SMALL));
  // myset(datatype,&ftemp,0,1,writebuf); // 1


  if(EOMRADTYPE!=EOMRADNONE){

    // Erf (for various things)
    ftemp=(float)GLOBALMACP0A1(pdump,i,j,k,PRAD0);
    myset(datatype,&ftemp,0,1,writebuf); // 1

    // urad^t [ in grid frame]
    ftemp=(float)(q.uradcon[0]);
    myset(datatype,&ftemp,0,1,writebuf); // 1

    // vrad^r [ in grid frame]
    ftemp=(float)(q.uradcon[1]/q.uradcon[0]);
    myset(datatype,&ftemp,0,1,writebuf); // 1
    
    // vrad^\theta
    ftemp=(float)(q.uradcon[2]/q.uradcon[0]);
    myset(datatype,&ftemp,0,1,writebuf); // 1
    
    // vrad^\phi
    ftemp=(float)(q.uradcon[3]/q.uradcon[0]);
    myset(datatype,&ftemp,0,1,writebuf); // 1






    // add kappaes, kappa, kappan, kappaemit, kappanemit so don't have to keep python scripts up to date with form or how opacity is dealth with.

    //radiative stress tensor in the lab frame
    FTYPE Rij[NDIM][NDIM];

    //this call returns R^i_j, i.e., the first index is contra-variant and the last index is co-variant
    mhdfull_calc_rad(pr, ptrgeom, &q, Rij);

    //the four-velocity of fluid in lab frame
    FTYPE *ucon,*ucov;
    ucon = q.ucon;
    ucov = q.ucov;

    //Eradff = R^a_b u_a u^b
    FTYPE Ruu=0.; DLOOP(i,j) Ruu+=Rij[i][j]*ucov[i]*ucon[j];
    // get relative Lorentz factor between gas and radiation
    FTYPE gammaradgas = 0.0;
    int jj;
    DLOOPA(jj) gammaradgas += - (q.ucov[jj] * q.uradcon[jj]);

    // get B
    FTYPE bsq = dot(q.bcon, q.bcov);
    FTYPE B=sqrt(bsq);

    FTYPE nradff=0;
    FTYPE kappa,kappan=0;
    FTYPE kappaemit,kappanemit=0;
    FTYPE kappaes;
    FTYPE lambda,nlambda=0;
    FTYPE Tgas=0,Tradff=0;
    calc_Tandopacityandemission(pr,ptrgeom,&q,Ruu,gammaradgas,B,&Tgas,&Tradff,&nradff,&kappa,&kappan,&kappaemit,&kappanemit,&kappaes, &lambda, &nlambda);

    myset(datatype,&Tgas,0,1,writebuf); // 1
    myset(datatype,&Tradff,0,1,writebuf); // 1
    myset(datatype,&nradff,0,1,writebuf); // 1
    myset(datatype,&kappa,0,1,writebuf); // 1
    myset(datatype,&kappan,0,1,writebuf); // 1
    myset(datatype,&kappaemit,0,1,writebuf); // 1
    myset(datatype,&kappanemit,0,1,writebuf); // 1
    myset(datatype,&kappaes,0,1,writebuf); // 1
    myset(datatype,&lambda,0,1,writebuf); // 1
    myset(datatype,&nlambda,0,1,writebuf); // 1

  }



  return(0);

}



/// dissipation dump
int dissdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


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

/// dissipation dump contents number
void set_dissdump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DODISS){
    *numcolumnsvar=NUMDISSFUNPOS;
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;

}


/// dissipation dump contents
int dissdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,&GLOBALMAC(dissfunpos,i,j,k),0,NUMDISSFUNPOS,writebuf);

  return (0);
}


/// other dump stuff
int dumpother(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


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

/// other dump stuff contents number
void set_dumpother_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DODUMPOTHER){ // panalytic + numpother quantities
    *numcolumnsvar=NPR+NUMPOTHER;
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;


}

/// other dump stuff contents
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



/// flux (usually for debug) dump
int fluxdumpdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


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

/// flux dump contents number
void set_fluxdump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(FLUXDUMP==1){ // dU, flux, and ppprimitives for flux
    *numcolumnsvar=NUMFLUXDUMP;
  }
  else *numcolumnsvar=0;



  // Version number:
  *numversion=0;


}


/// flux dump contents
int fluxdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl,pliter;

  myset(datatype,&GLOBALMAC(fluxdump,i,j,k),0,NUMFLUXDUMP,writebuf);

  return (0);
}



/// dump stuff related specially to EOS
int eosdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};

  
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


/// dump stuff related specially to EOS -- contents number
void set_eosdump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOEOSDIAG==1){
    // all EOSs output same size data so uniform format
    // otherwise have to also put this condition in dump.c when outputting so don't overwrite memory!
    *numcolumnsvar=MAXPARLIST+1+MAXNUMEXTRAS+MAXPROCESSEDEXTRAS; // 1 is temperature
  }
  else{
    *numcolumnsvar=0;
  }
    
  // Version number:
  *numversion=0;


}

/// dump stuff related specially to EOS -- contents
int eosdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  FTYPE rho,u;
  FTYPE height,ye,ynu,temp;
  // NOTEMARK: have to fill with zero or else might not be number and if ldouble is used then seg faults
  FTYPE extras[MAXNUMEXTRAS]={0};
  FTYPE processed[MAXPROCESSEDEXTRAS]={0};
  int numextras;
  int extraiter;
  int numparms;
  FTYPE parlist[MAXPARLIST]={0};
  int loc=CENT;


  //  dualfprintf(fail_file,"NUMEOSGLOBALS=%g LASTEOSGLOBAL=%g FIRSTEOSGLOBAL=%g\n",NUMEOSGLOBALS,LASTEOSGLOBAL,FIRSTEOSGLOBAL);

  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumnsvar


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














/// dump stuff related specially to RAD
int raddump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};

  
  trifprintf("begin dumping raddump# %ld ... ",dump_cnt);

  whichdump=RADDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/raddump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,raddump_content)>=1) return(1);

  /////// output the symmetry information to the fail file
  //writesyminfo();
  ///////

  trifprintf("end dumping raddump# %ld ... ",dump_cnt);


  return(0);
  
}

/// rad dump contents number
void set_raddump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(EOMRADTYPE!=EOMRADNONE && DORADDIAG){
    *numcolumnsvar=NDIM*2 + NDIM+1 + NDIM + NDIM*2 + 1*14 + 4*3;
  }
  else{
    *numcolumnsvar=0;
  }

  // Version number:
  *numversion=0;


}

/// rad dump contents
int raddump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int loc=CENT;
  FTYPE *pr=&GLOBALMACP0A1(pdump,i,j,k,0);
  

  // leave if shouldn't write anything
  if(EOMRADTYPE==EOMRADNONE) return(0);

  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumnsvar


  // get X,V
  FTYPE X[NDIM],V[NDIM];
  bl_coord_ijk_2(i, j, k, loc, X,V);
  // get geometry
  struct of_geom geomdontuse;
  struct of_geom *ptrgeom=&geomdontuse;
  get_geometry(i, j, k, loc, ptrgeom);

  // get full state
  struct of_state q;
  get_state(pr,ptrgeom,&q);

  myset(datatype,q.uradcon,0,NDIM,writebuf); // NDIM
  myset(datatype,q.uradcov,0,NDIM,writebuf); // NDIM
 
  // Transform these lab frame coordinate basis primitives to fluid frame E,F^i
  int whichvel=VEL4; // desired final ff version
  int whichcoord=MCOORD; // desired final ff version
  FTYPE pradffortho[NPR]={-1};
  FTYPE prff[NPR]; // not new compared to pr
  prad_fforlab(&whichvel, &whichcoord, LAB2FF, i,j,k,CENT,ptrgeom,pradffortho, pr, prff);
  myset(datatype,pradffortho,PRAD0,NDIM,writebuf); // NDIM

  // get 4-force in lab and fluid frame
  FTYPE Gdpl[NPR]={0},Gdabspl[NPR]={0},chi,Tgas,Trad;
  int computestate=0;// already computed above
  int computeentropy=1;
  koral_source_rad_calc(computestate,computeentropy,pr, ptrgeom, Gdpl, Gdabspl, &chi, &Tgas, &Trad, &q);
  myset(datatype,Gdpl,PRAD0,NDIM,writebuf); // NDIM
  myset(datatype,Gdabspl,PRAD0,NDIM,writebuf); // NDIM

  // get gas T
  myset(datatype,&Tgas,0,1,writebuf); // 1

  // get radiation T from actual Erf
  FTYPE Tradlte=calc_LTE_TfromE(pr[PRAD0]);
  myset(datatype,&Tradlte,0,1,writebuf); // 1

  // get radiation's fluid frame T from actual Erf in fluid frame
  FTYPE Tradff=calc_LTE_TfromE(pradffortho[PRAD0]);
  myset(datatype,&Tradff,0,1,writebuf); // 1

  // get Trad -- full fluid frame T
  myset(datatype,&Trad,0,1,writebuf); // 1


  //radiative stress tensor in the lab frame
  FTYPE Rij[NDIM][NDIM];

  //this call returns R^i_j, i.e., the first index is contra-variant and the last index is co-variant
  mhdfull_calc_rad(pr, ptrgeom, &q, Rij);

  //the four-velocity of fluid in lab frame
  FTYPE *ucon,*ucov;
  ucon = q.ucon;
  ucov = q.ucov;

  //Eradff = R^a_b u_a u^b
  FTYPE Ruu=0.; DLOOP(i,j) Ruu+=Rij[i][j]*ucov[i]*ucon[j];
  // get relative Lorentz factor between gas and radiation
  FTYPE gammaradgas = 0.0;
  int jj;
  DLOOPA(jj) gammaradgas += - (q.ucov[jj] * q.uradcon[jj]);

  // get Erf in fluid frame
  myset(datatype,&Ruu,0,1,writebuf); // 1

  // get Erf [assuming LTE]
  FTYPE Erf;
  Erf=calc_LTE_Efromurho(pr[UU],pr[RHO]);
  myset(datatype,&Erf,0,1,writebuf); // 1

  // get B
  FTYPE bsq = dot(q.bcon, q.bcov);
  FTYPE B=sqrt(bsq);

  // get lambda
  //  FTYPE Tgas=calc_PEQ_Tfromurho(pr[UU],pr[RHO]);
  //  calc_rad_lambda(pr, ptrgeom, Tgas, &lambda, &nlambda, &kappaemit, &kappanemit);
  // get absorption opacities
  //  FTYPE Tgas;
  FTYPE nradff=0;
  FTYPE kappa,kappan=0;
  FTYPE kappaemit,kappanemit=0;
  FTYPE kappaes;
  FTYPE lambda,nlambda=0;
  calc_Tandopacityandemission(pr,ptrgeom,&q,Ruu,gammaradgas,B,&Tgas,&Tradff,&nradff,&kappa,&kappan,&kappaemit,&kappanemit,&kappaes, &lambda, &nlambda);

  myset(datatype,&nradff,0,1,writebuf); // 1
  myset(datatype,&kappa,0,1,writebuf); // 1
  myset(datatype,&kappan,0,1,writebuf); // 1
  myset(datatype,&kappaemit,0,1,writebuf); // 1
  myset(datatype,&kappanemit,0,1,writebuf); // 1
  myset(datatype,&kappaes,0,1,writebuf); // 1
  myset(datatype,&lambda,0,1,writebuf); // 1
  myset(datatype,&nlambda,0,1,writebuf); // 1

  // get tau
  FTYPE tautot[NDIM]={0},tautotmax;
  calc_tautot(pr, ptrgeom, &q, tautot, &tautotmax);
  myset(datatype,tautot,0,NDIM,writebuf); // NDIM
  myset(datatype,&tautotmax,0,1,writebuf); // 1



  // get radiative vchar
  int ignorecourant;
  FTYPE vmin1=0,vmax1=0,vmax21=0,vmin21=0;
  vchar_rad(pr, &q, 1, ptrgeom, &vmax1, &vmin1, &vmax21, &vmin21, &ignorecourant);
  myset(datatype,&vmin1,0,1,writebuf); // 1
  myset(datatype,&vmax1,0,1,writebuf); // 1
  myset(datatype,&vmin21,0,1,writebuf); // 1
  myset(datatype,&vmax21,0,1,writebuf); // 1
  FTYPE vmin2=0,vmax2=0,vmax22=0,vmin22=0;
  vchar_rad(pr, &q, 2, ptrgeom, &vmax2, &vmin2, &vmax22, &vmin22, &ignorecourant);
  myset(datatype,&vmin2,0,1,writebuf); // 1
  myset(datatype,&vmax2,0,1,writebuf); // 1
  myset(datatype,&vmin22,0,1,writebuf); // 1
  myset(datatype,&vmax22,0,1,writebuf); // 1
  FTYPE vmin3=0,vmax3=0,vmax23=0,vmin23=0;
  vchar_rad(pr, &q, 3, ptrgeom, &vmax3, &vmin3, &vmax23, &vmin23, &ignorecourant);
  myset(datatype,&vmin3,0,1,writebuf); // 1
  myset(datatype,&vmax3,0,1,writebuf); // 1
  myset(datatype,&vmin23,0,1,writebuf); // 1
  myset(datatype,&vmax23,0,1,writebuf); // 1

  return (0);
}











/// Vector potential dump
int vpotdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


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


/// Vector potential dump contents number
void set_vpotdump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOVPOTDUMP){
    *numcolumnsvar=NUMVPOTDUMP;
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;



}

/// Vector potential dump contents
int vpotdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int jj;

  // NOTEMARK: see also restart.c since this is added to restart 
  for(jj=0;jj<NUMVPOTDUMP;jj++){
    myset(datatype,&GLOBALMACP1A0(vpotarraydump,jj,i,j,k),0,1,writebuf); // 1 each
  }

  return (0);
}







/// failure and floor dU dump
int failfloordudump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


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


/// failfloor dump content number
void set_failfloordudump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOFLOORDIAG){
    *numcolumnsvar=NPR; // normal dU from floor
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;



}

/// failfloor dump contents
int failfloordudump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl;

  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,&GLOBALMAC(failfloordu,i,j,k),0,NPR,writebuf); // NPR
  
  return (0);
}




/// failure and floor dU dump
int fluxsimpledump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


  trifprintf("begin dumping fluxsimpledump# %ld ... ",dump_cnt);

  whichdump=FLUXSIMPLEDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/fluxsimpledump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,fluxsimpledump_content)>=1) return(1);

  trifprintf("end dumping fluxsimpledump# %ld ... ",dump_cnt);


  return(0);
  
}


/// simple fluxes dump content number
void set_fluxsimpledump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DOFLOORDIAG && N1NOT1==1){
    *numcolumnsvar=NPR; // radial fluxes only
  }
  else *numcolumnsvar=0;


  if(FLUXDUMP==2) *numcolumnsvar+= (NUMFLUXESTOSAVE*NUMPHYSICALFLUXTERMS);


  // Version number:
  *numversion=0;



}

/// fluxsimple dump contents
int fluxsimpledump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl;

  myset(datatype,&GLOBALMACP0A1(F1,i,j,k,0),0,NPR,writebuf); // NPR


  FTYPE ftemp;
  if(FLUXDUMP==2){ // can reduce this alot (rho not in en or yflx's, etc.)
    int fluxterm;
    int pliter;
    for(fluxterm=0;fluxterm<NUMPHYSICALFLUXTERMS;fluxterm++){
      PLOOP(pliter,pl){
        if(FLUXESTOSAVEPL(pl)){
          ftemp=(GLOBALMACP0A1(fluxdump,i,j,k,fluxterm*NPR + pl)); // (float)
          myset(datatype,&ftemp,0,1,writebuf);
        }
      }
    }
  }
  
  return (0);
}




/// dissipation measure dump
int dissmeasuredump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


  trifprintf("begin dumping dissmeasuredump# %ld ... ",dump_cnt);

  whichdump=DISSMEASUREDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dissmeasuredump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dissmeasuredump_content)>=1) return(1);

  trifprintf("end dumping dissmeasuredump# %ld ... ",dump_cnt);


  return(0);
  
}

/// dissipation measure dump contents number
void set_dissmeasuredump_content_dnumcolumns_dnumversion(int *numcolumnsvar, int *numversion)
{

  if(DODISSMEASURE && DODISSMEASUREDIAG){
    *numcolumnsvar=NSPECIAL+1; // dissmeasurepl and dissmeasure
    *numcolumnsvar+=3*2; // Fi for each direction and gas/rad
  }
  else *numcolumnsvar=0;

  // Version number:
  *numversion=0;

}


/// dissipation measure dump contents
int dissmeasuredump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // NOTEMARK: see also restart.c since this is added to restart 
  myset(datatype,&GLOBALMAC(dissmeasurearray,i,j,k),0,NSPECIAL+1+3*2,writebuf);

  return (0);
}



/// fake dump so can push out data in case still in MPI=2 delayed writing buffer
/// Kraken had problems with this when computing and dumping currents.
int fakedump(long dump_cnt)// arg not used
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME]={'\0'};
  char filesuffix[MAXFILENAME]={'\0'};
  char fileformat[MAXFILENAME]={'\0'};


  trifprintf("begin dumping fakedump");

  //  return(0); // Kraken had problems with this when computing and dumping currents.

  whichdump=FAKEDUMPTYPE;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/fakedump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,-1,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,fakedump_header,fakedump_content)>=1) return(1);

  trifprintf("end dumping fakedump");


  return(0);
  
}


// fake dump content number
int fakedump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int jj;

  // blank

  return (0);
}


// fake dump header
int fakedump_header(int whichdump, int whichdumpversion, int numcolumnsvar,int bintxt, FILE *headerptr)
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

