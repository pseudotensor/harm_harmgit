#if(SENSITIVE==LONGDOUBLETYPE)
// assume sensitive>=realtype in precision
#if(REALTYPE==LONGDOUBLETYPE) // was FLOATTYPE==REALTYPE and SENS=DOUBLETYPE
#define HEADEROLDONEIN "%Lf"
#define HEADEROLD2IN "%Lf %Lf"
#define HEADEROLD3IN "%Lf %Lf %Lf"
#define HEADEROLD4IN "%Lf %Lf %Lf %Lf"
#define HEADEROLD5IN "%Lf %Lf %Lf %Lf %Lf"
#define HEADEROLD6IN "%Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADEROLD7IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADEROLD8IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADEROLD9IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define RESTARTHEADEROLD "%d %d %d "                                    \
  "%Lf %Lf %ld %Lf %Lf %Lf "                                            \
  "%Lf %Lf %Lf %Lf %ld %Lf %ld %ld %ld %ld %ld "                        \
  "%Lf %d %d %d %d %d %d %d %d "                                        \
  "%Lf %Lf %Lf %Lf %d "                                                 \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#elif(REALTYPE==DOUBLETYPE)
#define HEADEROLDONEIN "%lf"
#define HEADEROLD2IN "%lf %lf"
#define HEADEROLD3IN "%lf %lf %lf"
#define HEADEROLD4IN "%lf %lf %lf %lf"
#define HEADEROLD5IN "%lf %lf %lf %lf %lf"
#define HEADEROLD6IN "%lf %lf %lf %lf %lf %lf"
#define HEADEROLD7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADEROLD8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADEROLD9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADEROLD "%d %d %d "                                    \
  "%Lf %Lf %ld %lf %lf %lf "                                            \
  "%Lf %Lf %Lf %Lf %ld %Lf %ld %ld %ld %ld %ld "                        \
  "%Lf %d %d %d %d %d %d %d %d "                                        \
  "%lf %lf %lf %lf %d "                                                 \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADEROLDONEIN "%f"
#define HEADEROLD2IN "%f %f"
#define HEADEROLD3IN "%f %f %f"
#define HEADEROLD4IN "%f %f %f %f"
#define HEADEROLD5IN "%f %f %f %f %f"
#define HEADEROLD6IN "%f %f %f %f %f %f"
#define HEADEROLD7IN "%f %f %f %f %f %f %f"
#define HEADEROLD8IN "%f %f %f %f %f %f %f %f"
#define HEADEROLD9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADEROLD "%d %d %d "                                    \
  "%lf %lf %ld %f %f %f "                                               \
  "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "                        \
  "%lf %d %d %d %d %d %d %d %d "                                        \
  "%f %f %f %f %d "                                                     \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==DOUBLETYPE)
// assume sensitive>=realtype in precision
#if(REALTYPE==DOUBLETYPE)
#define HEADEROLDONEIN "%lf"
#define HEADEROLD2IN "%lf %lf"
#define HEADEROLD3IN "%lf %lf %lf"
#define HEADEROLD4IN "%lf %lf %lf %lf"
#define HEADEROLD5IN "%lf %lf %lf %lf %lf"
#define HEADEROLD6IN "%lf %lf %lf %lf %lf %lf"
#define HEADEROLD7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADEROLD8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADEROLD9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADEROLD "%d %d %d "                                    \
  "%lf %lf %ld %lf %lf %lf "                                            \
  "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "                        \
  "%lf %d %d %d %d %d %d %d %d "                                        \
  "%lf %lf %lf %lf %d "                                                 \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADEROLDONEIN "%f"
#define HEADEROLD2IN "%f %f"
#define HEADEROLD3IN "%f %f %f"
#define HEADEROLD4IN "%f %f %f %f"
#define HEADEROLD5IN "%f %f %f %f %f"
#define HEADEROLD6IN "%f %f %f %f %f %f"
#define HEADEROLD7IN "%f %f %f %f %f %f %f"
#define HEADEROLD8IN "%f %f %f %f %f %f %f %f"
#define HEADEROLD9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADEROLD "%d %d %d "                                    \
  "%lf %lf %ld %f %f %f "                                               \
  "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "                        \
  "%lf %d %d %d %d %d %d %d %d "                                        \
  "%f %f %f %f %d "                                                     \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==FLOATTYPE)

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADEROLD "" // dumb, so crash on compile
#elif(REALTYPE==FLOATTYPE)
#define HEADEROLDONEIN "%f"
#define HEADEROLD2IN "%f %f"
#define HEADEROLD3IN "%f %f %f"
#define HEADEROLD4IN "%f %f %f %f"
#define HEADEROLD5IN "%f %f %f %f %f"
#define HEADEROLD6IN "%f %f %f %f %f %f"
#define HEADEROLD7IN "%f %f %f %f %f %f %f"
#define HEADEROLD8IN "%f %f %f %f %f %f %f %f"
#define HEADEROLD9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADEROLD "%d %d %d "                                    \
  "%f %f %ld %f %f %f "                                                 \
  "%f %f %f %f %ld %f %ld %ld %ld %ld %ld "                             \
  "%f %d %d %d %d %d %d %d %d "                                         \
  "%f %f %f %f %d "                                                     \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#endif


// 2
// 6
// 10
// 3
// 5
// 4
// SUM=30
// 23+10=33
// total=63
#define WRITERESTARTHEADEROLD "%lld %lld %lld "                         \
  "%21.15g %21.15g %ld %21.15g %21.15g %21.15g "                        \
  "%21.15g %21.15g %21.15g %21.15g %ld %21.15g %ld %ld %ld %ld %ld"     \
  "%21.15g %d %d %d %d %d %d %d %d "                                    \
  "%21.15g %21.15g %21.15g %21.15g %d "                                 \
  "%d %d %d %d %d %d "                                                  \
  "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "

#define HEADEROLDONEOUT "%21.15g "




// extra calls only for reading header.  If want output header written as old form, then don't erase data just not outputted since have DISS&&0 and DISSVSR&&0 in restart.old.c
int extrarestartfunction_old(void)
{
  int enerregion;
  int dissloop;
  int i;
  int pl,pliter;


  gamideal=gam;
  MBH=1.0;
  QBH=0.0;
  EP3=0.0;
  THETAROT=0.0;

  
  if(DODISS) ENERREGIONLOOP(enerregion) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) dissreg_tot[enerregion][dissloop]=0.0;
   
  if(DODISSVSR) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) for(i=0;i<ncpux1*N1;i++) dissvsr_tot[dissloop][i]=0.0;


  // must set some things manually then
  DTdumpgen[DISSDUMPTYPE] = DTdumpgen[MAINDUMPTYPE]/5.0;
  dumpcntgen[DISSDUMPTYPE]=dumpcntgen[MAINDUMPTYPE]*5; // not exactly right, so should manually set


  DTdumpgen[FLUXDUMPTYPE]=DTdumpgen[OTHERDUMPTYPE]=DTdumpgen[EOSDUMPTYPE]=DTdumpgen[VPOTDUMPTYPE]=DTdumpgen[MAINDUMPTYPE];
  DTdumpgen[FIELDLINEDUMPTYPE]=DTdumpgen[ENERDUMPTYPE];

  dumpcntgen[FLUXDUMPTYPE]=dumpcntgen[OTHERDUMPTYPE]=dumpcntgen[EOSDUMPTYPE]=dumpcntgen[VPOTDUMPTYPE]=dumpcntgen[MAINDUMPTYPE];
  dumpcntgen[FIELDLINEDUMPTYPE]=dumpcntgen[ENERDUMPTYPE];


  // these things were not stored in old restart file
  // below floor model is only used if rescaletype!=4
  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prfloorcoef[RHO]=RHOMIN/100.0;
    prfloorcoef[UU]=UUMIN/100.0;
  }
  else{
    prfloorcoef[RHO]=RHOMIN;
    prfloorcoef[UU]=UUMIN;
  }


  PALLLOOP(pl) do_transverse_flux_integration[pl]=do_transverse_flux_integration[RHO];
  PALLLOOP(pl) do_source_integration[pl]=do_source_integration[RHO];
  PALLLOOP(pl) do_conserved_integration[pl]=do_conserved_integration[RHO];



  return(0);
}


// headerptr created and only used here OR passed a given pointer
int read_restart_header_old(int whichdump, int whichdumpversion, int numcolumns, int bintxt, FILE*headerptr)
{
  int ii;
  int dir,pl,pliter;
  int enerregion, floor, tscale;
  int idum1,idum2,idum3;

  trifprintf("begin reading header of restart file\n");

  if(bintxt==BINARYOUTPUT){

    /* read in global variables, in binary */
    fread(&idum1, sizeof(int), 1, headerptr);
    fread(&idum2, sizeof(int), 1, headerptr);
    fread(&idum3, sizeof(int), 1, headerptr); // extra 3D thing

    // all cpus read the rest of header the same
    fread(&t, sizeof(SFTYPE), 1, headerptr);
    fread(&tf, sizeof(SFTYPE), 1, headerptr);
    fread(&nstep, sizeof(long), 1, headerptr);
    fread(&a, sizeof(FTYPE), 1, headerptr);
    fread(&gam, sizeof(FTYPE), 1, headerptr);
    fread(&cour, sizeof(FTYPE), 1, headerptr);
      
    fread(&DTdumpgen[MAINDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fread(&DTdumpgen[AVG1DUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fread(&DTdumpgen[ENERDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fread(&DTdumpgen[IMAGEDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    // fread(&DTr, sizeof(SFTYPE), 1, headerptr) ;
    fread(&DTr, sizeof(long), 1, headerptr);
    fread(&DTdumpgen[DEBUGDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fread(&dumpcntgen[MAINDUMPTYPE], sizeof(long), 1, headerptr);
    fread(&dumpcntgen[IMAGEDUMPTYPE], sizeof(long), 1, headerptr);
    fread(&rdump_cnt, sizeof(long), 1, headerptr);
    fread(&dumpcntgen[AVG1DUMPTYPE], sizeof(long), 1, headerptr);
    fread(&dumpcntgen[DEBUGDUMPTYPE], sizeof(long), 1, headerptr);
      
      
    fread(&dt, sizeof(SFTYPE), 1, headerptr);
    fread(&lim[1], sizeof(int), 1, headerptr);
    fread(&lim[2], sizeof(int), 1, headerptr);
    fread(&lim[3], sizeof(int), 1, headerptr);
    fread(&TIMEORDER, sizeof(int), 1, headerptr);
    fread(&fluxmethod, sizeof(int), 1, headerptr);
    fread(&FLUXB, sizeof(int), 1, headerptr);
    fread(&UTOPRIMVERSION, sizeof(int), 1, headerptr);
    fread(&failed, sizeof(int), 1, headerptr);
      
    fread(&R0, sizeof(FTYPE), 1, headerptr);
    fread(&Rin, sizeof(FTYPE), 1, headerptr);
    fread(&Rout, sizeof(FTYPE), 1, headerptr);
    fread(&hslope, sizeof(FTYPE), 1, headerptr);
    fread(&defcoord, sizeof(FTYPE), 1, headerptr);

    fread(&BCtype[X1UP],sizeof(int), 1, headerptr);
    fread(&BCtype[X1DN],sizeof(int), 1, headerptr);
    fread(&BCtype[X2UP],sizeof(int), 1, headerptr);
    fread(&BCtype[X2DN],sizeof(int), 1, headerptr);
    fread(&BCtype[X3UP],sizeof(int), 1, headerptr); // 3D thing
    fread(&BCtype[X3DN],sizeof(int), 1, headerptr);

    // new May 6, 2003
    fread(&realnstep, sizeof(long), 1, headerptr);
    fread(&debugfail,sizeof(int), 1, headerptr);
    fread(&whichrestart,sizeof(int), 1, headerptr);
    fread(&cooling,sizeof(int), 1, headerptr);
    fread(&restartsteps[0],sizeof(long), 1, headerptr);
    fread(&restartsteps[1],sizeof(long), 1, headerptr);
    fread(&GAMMIEDUMP,sizeof(int), 1, headerptr);
    fread(&GAMMIEIMAGE,sizeof(int), 1, headerptr);
    fread(&GAMMIEENER,sizeof(int), 1, headerptr);
    fread(&DODIAGS,sizeof(int), 1, headerptr);
    fread(&DOENERDIAG,sizeof(int), 1, headerptr);
    fread(&DOGDUMPDIAG,sizeof(int), 1, headerptr);
    fread(&DORDUMPDIAG,sizeof(int), 1, headerptr);
    fread(&DODUMPDIAG,sizeof(int), 1, headerptr);
    fread(&DOAVGDIAG,sizeof(int), 1, headerptr);
    fread(&DOIMAGEDIAG,sizeof(int), 1, headerptr);
    fread(&DOAREAMAPDIAG,sizeof(int), 1, headerptr);
    fread(&POSDEFMETRIC,sizeof(int), 1, headerptr);
    fread(&periodicx1,sizeof(int), 1, headerptr);
    fread(&periodicx2,sizeof(int), 1, headerptr);
    fread(&periodicx3,sizeof(int), 1, headerptr); // 3D thing
    fread(&binaryoutput,sizeof(int), 1, headerptr);
    fread(&sortedoutput,sizeof(int), 1, headerptr);
    fread(&defcon,sizeof(FTYPE), 1, headerptr);
    fread(&SAFE,sizeof(FTYPE), 1, headerptr);
    fread(&RHOMIN,sizeof(FTYPE), 1, headerptr);
    fread(&UUMIN,sizeof(FTYPE), 1, headerptr);
    fread(&RHOMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&UUMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&BSQORHOLIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&BSQOULIMIT,sizeof(FTYPE), 1, headerptr);
    fread(&GAMMAMAX,sizeof(FTYPE), 1, headerptr);
    fread(&GAMMADAMP,sizeof(FTYPE), 1, headerptr);
    fread(&GAMMAFAIL,sizeof(FTYPE), 1, headerptr);      
    // end new May 6, 2003

    PDUMPLOOP(pliter,pl) fread(&prMAX[pl],sizeof(FTYPE),1,headerptr);

    fread(&rescaletype,sizeof(int),1,headerptr);

    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fread(&avgscheme[1],sizeof(int),1,headerptr);
    fread(&avgscheme[2],sizeof(int),1,headerptr);
    fread(&avgscheme[3],sizeof(int),1,headerptr);
    fread(&do_transverse_flux_integration[RHO],sizeof(int),1,headerptr);
    fread(&do_source_integration[RHO],sizeof(int),1,headerptr);
    fread(&do_conserved_integration[RHO],sizeof(int),1,headerptr);
    fread(&INVERTFROMAVERAGEIFFAILED,sizeof(int),1,headerptr);
    fread(&LIMIT_AC_PRIM_FRAC_CHANGE,sizeof(int),1,headerptr);
    fread(&MAX_AC_PRIM_FRAC_CHANGE,sizeof(FTYPE),1,headerptr);
    fread(&DOENOFLUX,sizeof(int),1,headerptr);


    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(pliter,pl) fread(&pcumreg_tot[enerregion][dir][pl],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fread(&fladdreg_tot[enerregion][pl],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fread(&sourceaddreg_tot[enerregion][pl],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fread(&Ureg_init_tot[enerregion][pl],sizeof(FTYPE),1,headerptr);
    TSCALELOOP(tscale) FLOORLOOP(floor) fread(&failfloorcountlocal_tot[0][tscale][floor],sizeof(CTYPE),1,headerptr); // keep old style with finalstep==0
    // end new June 6,2003

    //    if(DODISS&&0) ENERREGIONLOOP(enerregion) fread(&dissreg_tot[enerregion],sizeof(FTYPE),1,headerptr);
    
    if(DOLUMVSR) fread(&lumvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);
    //    if(DODISSVSR&&0) fread(&dissvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);


    fread(&UORHOLIMIT,sizeof(FTYPE), 1, headerptr);

 
  }
  // assumes headerptr is normal file pointer already defined
  else{
    fscanf(headerptr,RESTARTHEADEROLD,
           &idum1,&idum2,&idum3,
           &t,&tf,&nstep,&a,&gam,&cour,
           &(DTdumpgen[MAINDUMPTYPE]),&(DTdumpgen[AVG1DUMPTYPE]),&(DTdumpgen[ENERDUMPTYPE]),&(DTdumpgen[IMAGEDUMPTYPE]),&DTr,&(DTdumpgen[DEBUGDUMPTYPE]),&(dumpcntgen[MAINDUMPTYPE]),&(dumpcntgen[IMAGEDUMPTYPE]),&rdump_cnt,&(dumpcntgen[AVG1DUMPTYPE]),&(dumpcntgen[DEBUGDUMPTYPE]),
           &dt,&lim[1],&lim[2],&lim[3],&TIMEORDER,&fluxmethod,&FLUXB,&UTOPRIMVERSION,&failed,
           &R0,&Rin,&Rout,&hslope,&defcoord,
           &BCtype[X1UP],&BCtype[X1DN],&BCtype[X2UP],&BCtype[X2DN],&BCtype[X3UP],&BCtype[X3DN],
           &realnstep,&debugfail,&whichrestart,&cooling,&restartsteps[0],&restartsteps[1],&GAMMIEDUMP,&GAMMIEIMAGE,&GAMMIEENER,&DODIAGS,&DOENERDIAG,&DOGDUMPDIAG,&DORDUMPDIAG,&DODUMPDIAG,&DOAVGDIAG,&DOIMAGEDIAG,&DOAREAMAPDIAG,&POSDEFMETRIC,&periodicx1,&periodicx2,&periodicx3,&binaryoutput,&sortedoutput,&defcon,&SAFE,&RHOMIN,&UUMIN,&RHOMINLIMIT,&UUMINLIMIT,&BSQORHOLIMIT,&BSQOULIMIT,&GAMMAMAX,&GAMMADAMP,&GAMMAFAIL
           );
    PDUMPLOOP(pliter,pl) fscanf(headerptr,HEADEROLDONEIN,&prMAX[pl]);

    fscanf(headerptr,"%d",&rescaletype);

    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fscanf(headerptr,"%d",&avgscheme[1]);
    fscanf(headerptr,"%d",&avgscheme[2]);
    fscanf(headerptr,"%d",&avgscheme[3]);
    fscanf(headerptr,"%d",&do_transverse_flux_integration[RHO]);
    fscanf(headerptr,"%d",&do_source_integration[RHO]);
    fscanf(headerptr,"%d",&do_conserved_integration[RHO]);
    fscanf(headerptr,"%d",&INVERTFROMAVERAGEIFFAILED);
    fscanf(headerptr,"%d",&LIMIT_AC_PRIM_FRAC_CHANGE);
    fscanf(headerptr,HEADEROLDONEIN,&MAX_AC_PRIM_FRAC_CHANGE);
    fscanf(headerptr,"%d",&DOENOFLUX);


    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(pliter,pl) fscanf(headerptr,HEADEROLDONEIN,&pcumreg_tot[enerregion][dir][pl]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fscanf(headerptr,HEADEROLDONEIN,&fladdreg_tot[enerregion][pl]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fscanf(headerptr,HEADEROLDONEIN,&sourceaddreg_tot[enerregion][pl]);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fscanf(headerptr,HEADEROLDONEIN,&Ureg_init_tot[enerregion][pl]);
#if(COUNTTYPE==LONGLONGINTTYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fscanf(headerptr,"%lld",&failfloorcountlocal_tot[0][tscale][floor]);
#elif(COUNTTYPE==DOUBLETYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fscanf(headerptr,"%lf",&failfloorcountlocal_tot[0][tscale][floor]);
#endif
    // end new June 6,2003
    
    //    if(DODISS&&0) ENERREGIONLOOP(enerregion) fscanf(headerptr,HEADEROLDONEIN,&dissreg_tot[enerregion][0]);
    
    // assumes same CPU geometry during restart
    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) fscanf(headerptr,HEADEROLDONEIN,&lumvsr_tot[ii]);

    // assumes same CPU geometry during restart
    //    if(DODISSVSR&&0) for(ii=0;ii<ncpux1*N1;ii++) fscanf(headerptr,HEADEROLDONEIN,&dissvsr_tot[ii]);


    fscanf(headerptr,HEADEROLDONEIN,&UORHOLIMIT);


    // below handled in dump_gen()
    // flush to just after the header line in case binary read of data
    //    while(fgetc(headerptr)!='\n');
  }

  /////////////////
  //
  // some checks
  if (idum1 != totalsize[1]) {
    dualfprintf(fail_file, "error reading restart file; N1 differs\n");
    dualfprintf(fail_file, "got totalsize[1]=%d needed totalsize[1]=%d\n",idum1,totalsize[1]);
    myexit(3);
  }
  if (idum2 != totalsize[2]) {
    dualfprintf(fail_file, "error reading restart file; N2 differs\n");
    dualfprintf(fail_file, "got totalsize[2]=%d needed totalsize[2]=%d\n",idum2,totalsize[2]);
    myexit(4);
  }
  if (idum3 != totalsize[3]) {
    dualfprintf(fail_file, "error reading restart file; N3 differs\n");
    dualfprintf(fail_file, "got totalsize[3]=%d needed totalsize[3]=%d\n",idum3,totalsize[3]);
    myexit(4);
  }


  trifprintf("end reading header of restart file\n");

  return(0);
}


// all cpus should do this
int restart_read_defs_old(void)
{
  int enerregion;
  int floor,tscale;
  int dir,pl,pliter;
  int ii;

  if(myid==0){
    ////////////
    //
    // Define for cpu=0 only, which will continue to keep track of the total after restart
    //
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(pliter,pl) pcumreg[enerregion][dir][pl]=pcumreg_tot[enerregion][dir][pl];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fladdreg[enerregion][pl]=fladdreg_tot[enerregion][pl];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) sourceaddreg[enerregion][pl]=sourceaddreg_tot[enerregion][pl];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) Ureg_init[enerregion][pl]=Ureg_init_tot[enerregion][pl];

    TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcountlocal[0][tscale][floor]=failfloorcountlocal_tot[0][tscale][floor];


    //    if(DODISS&&0) ENERREGIONLOOP(enerregion) dissreg[enerregion][0]=dissreg_tot[enerregion][0];
    // assume (*dissfunpos)[][] not restored since zeroed out each dump

    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=lumvsr_tot[ii];
    //    if(DODISSVSR&&0) for(ii=0;ii<ncpux1*N1;ii++) dissvsr[ii]=dissvsr_tot[ii];

  }



  // now that CPU=0 has the restart header, pass to other CPUs
#if(USEMPI)
  MPI_Bcast(&t, 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&tf, 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&nstep, 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&a, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&gam, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&cour, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  
  MPI_Bcast(&DTdumpgen[MAINDUMPTYPE], 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DTdumpgen[AVG1DUMPTYPE], 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DTdumpgen[ENERDUMPTYPE], 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DTdumpgen[IMAGEDUMPTYPE], 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DTr, 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DTdumpgen[DEBUGDUMPTYPE], 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&dumpcntgen[MAINDUMPTYPE], 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&dumpcntgen[IMAGEDUMPTYPE], 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&rdump_cnt, 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&dumpcntgen[AVG1DUMPTYPE], 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&dumpcntgen[DEBUGDUMPTYPE], 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  
  MPI_Bcast(&dt, 1, MPI_SFTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&lim[1], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&lim[2], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&lim[3], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&TIMEORDER, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&fluxmethod, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&FLUXB, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&UTOPRIMVERSION, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&failed, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  
  MPI_Bcast(&R0, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&Rin, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&Rout, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&hslope, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&defcoord, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&BCtype[X1UP], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&BCtype[X1DN], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&BCtype[X2UP], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&BCtype[X2DN], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&BCtype[X3UP], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD); // 3D thing
  MPI_Bcast(&BCtype[X3DN], 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&realnstep, 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&debugfail, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&whichrestart, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&cooling, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&restartsteps[0], 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&restartsteps[1], 1, MPI_LONG, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&GAMMIEDUMP, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&GAMMIEIMAGE, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&GAMMIEENER, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DODIAGS, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DOENERDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DOGDUMPDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DORDUMPDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DODUMPDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DOAVGDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DOIMAGEDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DOAREAMAPDIAG, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&POSDEFMETRIC, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&periodicx1, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&periodicx2, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&periodicx3, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD); // 3D thing
  MPI_Bcast(&binaryoutput, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&sortedoutput, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&defcon, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&SAFE, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&RHOMIN, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&UUMIN, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&RHOMINLIMIT, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&UUMINLIMIT, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&BSQORHOLIMIT, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&BSQOULIMIT, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&GAMMAMAX, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&GAMMADAMP, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&GAMMAFAIL, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);

  // new June 6, 2003 (cumulatives)
  // don't broadcast cumulatives or initial conditions, only needed by cpu=0  
  // end new June 6,2003

  PDUMPLOOP(pliter,pl) MPI_Bcast(&prMAX[pl], 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);

  MPI_Bcast(&rescaletype, 1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);

  // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
  MPI_Bcast(&avgscheme[1],1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&avgscheme[2],1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&avgscheme[3],1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&do_transverse_flux_integration,1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&do_source_integration,1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&do_conserved_integration,1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&INVERTFROMAVERAGEIFFAILED,1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&LIMIT_AC_PRIM_FRAC_CHANGE,1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&MAX_AC_PRIM_FRAC_CHANGE,1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);
  MPI_Bcast(&DOENOFLUX,1, MPI_INT, MPIid[0], MPI_COMM_GRMHD);


  MPI_Bcast(&UORHOLIMIT, 1, MPI_FTYPE, MPIid[0], MPI_COMM_GRMHD);


#endif



  return(0);
}



int write_restart_header_old(int whichdump, int whichdumpversion, int numcolumns, int bintxt,FILE*headerptr)
{
  int dir,pl,pliter,floor,tscale;
  int enerregion;
  int ii;

  trifprintf("begin writing header of restart file\n");

  if(bintxt==BINARYOUTPUT){

    /* write out key global variables, in binary */
    fwrite(&totalsize[1], sizeof(int), 1, headerptr);
    fwrite(&totalsize[2], sizeof(int), 1, headerptr);
    fwrite(&totalsize[3], sizeof(int), 1, headerptr); // extra 3D thing
      
    fwrite(&t, sizeof(SFTYPE), 1, headerptr);
    fwrite(&tf, sizeof(SFTYPE), 1, headerptr);
    fwrite(&nstep, sizeof(long), 1, headerptr);
    fwrite(&a, sizeof(FTYPE), 1, headerptr);
    fwrite(&gam, sizeof(FTYPE), 1, headerptr);
    fwrite(&cour, sizeof(FTYPE), 1, headerptr);
      
    fwrite(&DTdumpgen[MAINDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fwrite(&DTdumpgen[AVG1DUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fwrite(&DTdumpgen[ENERDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fwrite(&DTdumpgen[IMAGEDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    // fwrite(&DTr, sizeof(SFTYPE),1, headerptr) ;
    fwrite(&DTr, sizeof(long), 1, headerptr);
    fwrite(&DTdumpgen[DEBUGDUMPTYPE], sizeof(SFTYPE), 1, headerptr);
    fwrite(&dumpcntgen[MAINDUMPTYPE], sizeof(long), 1, headerptr);
    fwrite(&dumpcntgen[IMAGEDUMPTYPE], sizeof(long), 1, headerptr);
    fwrite(&rdump_cnt, sizeof(long), 1, headerptr);
    fwrite(&dumpcntgen[AVG1DUMPTYPE], sizeof(long), 1, headerptr);
    fwrite(&dumpcntgen[DEBUGDUMPTYPE], sizeof(long), 1, headerptr);
      
    fwrite(&dt, sizeof(SFTYPE), 1, headerptr);
    fwrite(&lim[1], sizeof(int), 1, headerptr);
    fwrite(&lim[2], sizeof(int), 1, headerptr);
    fwrite(&lim[3], sizeof(int), 1, headerptr);
    fwrite(&TIMEORDER, sizeof(int), 1, headerptr);
    fwrite(&fluxmethod, sizeof(int), 1, headerptr);
    fwrite(&FLUXB, sizeof(int), 1, headerptr);
    fwrite(&UTOPRIMVERSION, sizeof(int), 1, headerptr);
    fwrite(&failed, sizeof(int), 1, headerptr);
      
    fwrite(&R0, sizeof(FTYPE), 1, headerptr);
    fwrite(&Rin, sizeof(FTYPE), 1, headerptr);
    fwrite(&Rout, sizeof(FTYPE), 1, headerptr);
    fwrite(&hslope, sizeof(FTYPE), 1, headerptr);
    fwrite(&defcoord, sizeof(FTYPE), 1, headerptr);

    fwrite(&BCtype[X1UP],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X1DN],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X2UP],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X2DN],sizeof(int), 1, headerptr);
    fwrite(&BCtype[X3UP],sizeof(int), 1, headerptr); // 3d thing
    fwrite(&BCtype[X3DN],sizeof(int), 1, headerptr);

    fwrite(&realnstep, sizeof(long), 1, headerptr);
    fwrite(&debugfail,sizeof(int), 1, headerptr);
    fwrite(&whichrestart,sizeof(int), 1, headerptr);
    fwrite(&cooling,sizeof(int), 1, headerptr);
    fwrite(&restartsteps[0],sizeof(long), 1, headerptr);
    fwrite(&restartsteps[1],sizeof(long), 1, headerptr);
    fwrite(&GAMMIEDUMP,sizeof(int), 1, headerptr);
    fwrite(&GAMMIEIMAGE,sizeof(int), 1, headerptr);
    fwrite(&GAMMIEENER,sizeof(int), 1, headerptr);
    fwrite(&DODIAGS,sizeof(int), 1, headerptr);
    fwrite(&DOENERDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOGDUMPDIAG,sizeof(int), 1, headerptr);
    fwrite(&DORDUMPDIAG,sizeof(int), 1, headerptr);
    fwrite(&DODUMPDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOAVGDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOIMAGEDIAG,sizeof(int), 1, headerptr);
    fwrite(&DOAREAMAPDIAG,sizeof(int), 1, headerptr);
    fwrite(&POSDEFMETRIC,sizeof(int), 1, headerptr);
    fwrite(&periodicx1,sizeof(int), 1, headerptr);
    fwrite(&periodicx2,sizeof(int), 1, headerptr);
    fwrite(&periodicx3,sizeof(int), 1, headerptr); // 3d thing
    fwrite(&binaryoutput,sizeof(int), 1, headerptr);
    fwrite(&sortedoutput,sizeof(int), 1, headerptr);
    fwrite(&defcon,sizeof(FTYPE), 1, headerptr);
    fwrite(&SAFE,sizeof(FTYPE), 1, headerptr);
    fwrite(&RHOMIN,sizeof(FTYPE), 1, headerptr);
    fwrite(&UUMIN,sizeof(FTYPE), 1, headerptr);
    fwrite(&RHOMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&UUMINLIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&BSQORHOLIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&BSQOULIMIT,sizeof(FTYPE), 1, headerptr);
    fwrite(&GAMMAMAX,sizeof(FTYPE), 1, headerptr);
    fwrite(&GAMMADAMP,sizeof(FTYPE), 1, headerptr);
    fwrite(&GAMMAFAIL,sizeof(FTYPE), 1, headerptr);

    PDUMPLOOP(pliter,pl) fwrite(&prMAX[pl],sizeof(FTYPE),1,headerptr);

    fwrite(&rescaletype,sizeof(int),1,headerptr);

    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fwrite(&avgscheme[1],sizeof(int),1,headerptr);
    fwrite(&avgscheme[2],sizeof(int),1,headerptr);
    fwrite(&avgscheme[3],sizeof(int),1,headerptr);
    fwrite(&do_transverse_flux_integration,sizeof(int),1,headerptr);
    fwrite(&do_source_integration,sizeof(int),1,headerptr);
    fwrite(&do_conserved_integration,sizeof(int),1,headerptr);
    fwrite(&INVERTFROMAVERAGEIFFAILED,sizeof(int),1,headerptr);
    fwrite(&LIMIT_AC_PRIM_FRAC_CHANGE,sizeof(int),1,headerptr);
    fwrite(&MAX_AC_PRIM_FRAC_CHANGE,sizeof(FTYPE),1,headerptr);
    fwrite(&DOENOFLUX,sizeof(int),1,headerptr);


    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(pliter,pl) fwrite(&pcumreg_tot[enerregion][dir][pl],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fwrite(&fladdreg_tot[enerregion][pl],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fwrite(&sourceaddreg_tot[enerregion][pl],sizeof(FTYPE),1,headerptr);
    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fwrite(&Ureg_init_tot[enerregion][pl],sizeof(FTYPE),1,headerptr);

    TSCALELOOP(tscale) FLOORLOOP(floor) fwrite(&failfloorcountlocal_tot[0][tscale][floor],sizeof(CTYPE),1,headerptr);
    // end new June 6,2003

    //    if(DODISS&&0) ENERREGIONLOOP(enerregion) fwrite(&dissreg_tot[enerregion][0],sizeof(FTYPE),1,headerptr);

    if(DOLUMVSR) fwrite(lumvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);
    //    if(DODISSVSR&&0) fwrite(dissvsr_tot,sizeof(SFTYPE),ncpux1*N1,headerptr);


    fwrite(&UORHOLIMIT,sizeof(FTYPE), 1, headerptr);


  }
  else if(bintxt==TEXTOUTPUT){
    fprintf(headerptr,WRITERESTARTHEADEROLD,
            totalsize[1],totalsize[2], totalsize[3], //3
            t,tf,nstep,a,gam,cour, // 6
            DTdumpgen[MAINDUMPTYPE],DTdumpgen[AVG1DUMPTYPE],DTdumpgen[ENERDUMPTYPE],DTdumpgen[IMAGEDUMPTYPE],DTr,DTdumpgen[DEBUGDUMPTYPE],dumpcntgen[MAINDUMPTYPE],dumpcntgen[IMAGEDUMPTYPE],rdump_cnt,dumpcntgen[AVG1DUMPTYPE],dumpcntgen[DEBUGDUMPTYPE], // 11
            dt,lim[1],lim[2],lim[3],TIMEORDER,fluxmethod,FLUXB,UTOPRIMVERSION,failed,  // 3
            R0,Rin,Rout,hslope,defcoord, // 5
            BCtype[X1UP],BCtype[X1DN],BCtype[X2UP],BCtype[X2DN],BCtype[X3UP],BCtype[X3DN], // 6
            realnstep,debugfail,whichrestart,cooling,restartsteps[0],restartsteps[1],GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG,DOIMAGEDIAG,DOAREAMAPDIAG,POSDEFMETRIC,periodicx1,periodicx2,periodicx3,binaryoutput,sortedoutput,defcon,SAFE,RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT,BSQORHOLIMIT,BSQOULIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL // 23+10+1
            // total = 65 (not updated)
            );

    // 74 total so far

    // end new June 6,2003
    PDUMPLOOP(pliter,pl) fprintf(headerptr,HEADEROLDONEOUT,prMAX[pl]);
    fprintf(headerptr,"%d ",rescaletype);

    // 83 total so far


    // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
    fprintf(headerptr,"%d ",avgscheme[1]);
    fprintf(headerptr,"%d ",avgscheme[2]);
    fprintf(headerptr,"%d ",avgscheme[3]);
    fprintf(headerptr,"%d ",do_transverse_flux_integration[RHO]);
    fprintf(headerptr,"%d ",do_source_integration[RHO]);
    fprintf(headerptr,"%d ",do_conserved_integration[RHO]);
    fprintf(headerptr,"%d ",INVERTFROMAVERAGEIFFAILED);
    fprintf(headerptr,"%d ",LIMIT_AC_PRIM_FRAC_CHANGE);
    fprintf(headerptr,HEADEROLDONEOUT,MAX_AC_PRIM_FRAC_CHANGE);
    fprintf(headerptr,"%d ",DOENOFLUX);

    // 93 total so far

    // new June 6, 2003 (cumulatives)
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(pliter,pl) fprintf(headerptr,HEADEROLDONEOUT,pcumreg_tot[enerregion][dir][pl]);

    // 93 + 3*6*8 = 237 so far

    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fprintf(headerptr,HEADEROLDONEOUT,fladdreg_tot[enerregion][pl]);

    // 237 + 3*8 = 261 so far

    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fprintf(headerptr,HEADEROLDONEOUT,sourceaddreg_tot[enerregion][pl]);

    // 261 + 3*8 = 285 so far

    ENERREGIONLOOP(enerregion) PDUMPLOOP(pliter,pl) fprintf(headerptr,HEADEROLDONEOUT,Ureg_init_tot[enerregion][pl]);

    // 285 + 3*8 = 309 so far


#if(COUNTTYPE==LONGLONGINTTYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fprintf(headerptr,"%lld ",failfloorcountlocal_tot[0][tscale][floor]);
#elif(COUNTTYPE==DOUBLETYPE)
    TSCALELOOP(tscale) FLOORLOOP(floor) fprintf(headerptr,"%21.15g ",failfloorcountlocal_tot[0][tscale][floor]);
#endif

    // 309 + 4*9 = 345
    
    //    if(DODISS&&0) ENERREGIONLOOP(enerregion) fprintf(headerptr,HEADEROLDONEOUT,dissreg_tot[enerregion][0]);

    // assumes same CPU geometry during restart
    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) fprintf(headerptr,HEADEROLDONEOUT,lumvsr_tot[ii]);

    // 345 + 64 = 409

    // assumes same CPU geometry during restart
    //    if(DODISSVSR&&0) for(ii=0;ii<ncpux1*N1;ii++) fprintf(headerptr,HEADEROLDONEOUT,dissvsr_tot[ii]);


    fprintf(headerptr,"%21.15g ",UORHOLIMIT);

    // 409 + 1 = 410


    fprintf(headerptr,"\n");
  }
  fflush(headerptr);


  trifprintf("end writing header of restart file\n");
  return(0);
}


