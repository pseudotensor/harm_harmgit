
///////////////////////////////////////////////////////////////////
//
// debug stuff
//
///////////////////////////////////////////////////////////////////



int faildebug1(int numnormterms, int whichcons, FTYPE *U_target, FTYPE *EOSextra, FTYPE *pr0, struct of_geom *ptrgeom)
{
  FTYPE tolx, tolf;
  int ntrial, mnewtfailed;
  int pl,pliter;
  FTYPE Ustart[NPR];
  FTYPE dUriemann[NPR];
  FTYPE dU[NPR],dUcomp[NUMSOURCES][NPR];
  FTYPE bsq;
  FTYPE mhd[NDIM][NDIM];
  FTYPE ucon[NDIM];
  FTYPE flux[NDIM][NPR];
  FTYPE **alpha;
  struct of_state q;
  int i,j,k;
  int myii,myjj,mykk;
  FTYPE Uwithgeom[NPR];
  FTYPE norm;


  PLOOP(pliter,pl) dUriemann[pl]=0.0; // fudge


  if (get_state(pr0, ptrgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 2);
  if (primtoU(UNOTHING,pr0, &q, ptrgeom, Ustart, NULL) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 2);
  
  mhd_calc(pr0,0, ptrgeom, &q, mhd[0], NULL);
  mhd_calc(pr0,1, ptrgeom, &q, mhd[1], NULL);
  mhd_calc(pr0,2, ptrgeom, &q, mhd[2], NULL);
  mhd_calc(pr0,3, ptrgeom, &q, mhd[3], NULL);
  primtoflux(UNOTHING,pr0, &q, 1,ptrgeom, flux[0], NULL);
  primtoflux(UNOTHING,pr0, &q, 2,ptrgeom, flux[1], NULL);
  primtoflux(UNOTHING,pr0, &q, 3,ptrgeom, flux[2], NULL);
  PLOOP(pliter,pl) Uwithgeom[pl]=Ustart[pl]*ptrgeom->EOMFUNCMAC(pl);
  FTYPE CUf[NUMDTCUFS]={1.0};
  int timeorder=0;
  FTYPE *CUimp=&CUf[NUMPREDTCUFS+timeorder];
  int didreturnpf=0;
  int eomtype=EOMDEFAULT;
  source(pr0, pr0, pr0, &didreturnpf, &eomtype, ptrgeom, &q, Uwithgeom, Uwithgeom, CUf, CUimp, 0.0, dUriemann, dUcomp,dU);
  bsq_calc(pr0, ptrgeom, &bsq);
  
  
  alpha=dmatrix(1, 5, 1, 5);

  if (dudp_calc(WHICHEOS,whichcons,EOSextra, pr0, &q, ptrgeom, alpha) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "dudp_calc()", 1);
  // more general normalization
  // copied from usrfun
#if(NORMMETHOD==-1)
  norm=1.0;
#elif(NORMMETHOD==0)
  norm=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
  norm=0.0;
  numnormterms=0;
  for (j = 0; j < INVERTNPR ; j++){
    for (k = 0; k < INVERTNPR ; k++){
      if(alpha[j + 1][k + 1]>NUMEPSILON){
        norm+=fabs(alpha[j + 1][k + 1]);
        numnormterms++;
      }
    }
  }
  norm=(FTYPE)(numnormterms)/(norm); // (i.e. inverse of average)
#endif    
  for (j = 0; j < INVERTNPR; j++)
    for (k = 0; k < INVERTNPR; k++)
      alpha[j + 1][k + 1] *= (norm);
  
  dualfprintf(fail_file,"myTud={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file, "%21.15g``20 ",mhd[i][j]);
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  
  dualfprintf(fail_file,"myF={");
  for(i=0;i<=2;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NPR;j++){
      dualfprintf(fail_file, "%21.15g``20 ",flux[i][j]);
      if(j<NPR-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<2)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  
  dualfprintf(fail_file,"myucon={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3]);
  dualfprintf(fail_file,"myucov={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3]);
  dualfprintf(fail_file,"mybcon={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3]);
  dualfprintf(fail_file,"mybcov={%21.15g, %21.15g, %21.15g, %21.15g}\n",q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3]);
  dualfprintf(fail_file,"myg=%21.15g\n",ptrgeom->gdet);
  
  //    PLOOP(pliter,pl) dualfprintf(fail_file,"pr0[%d]=%21.15g Ustart[%d]=%21.15g U_target[%d]=%21.15g\n",pl,pr0[pl],pl,Ustart[pl],pl,U_target[pl]);
  dualfprintf(fail_file,"pr0={");
  PLOOP(pliter,pl){
    dualfprintf(fail_file,"%21.15g``20 ",pr0[pl]);
    if(pl<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"myU0={");
  PLOOP(pliter,pl){
    dualfprintf(fail_file,"%21.15g``20 ",Ustart[pl]);
    if(pl<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"mydU={");
  PLOOP(pliter,pl){
    dualfprintf(fail_file,"%21.15g``20 ",dU[pl]);
    if(pl<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"Utarget={");
  PLOOP(pliter,pl){
    dualfprintf(fail_file,"%21.15g``20 ",U_target[pl]);
    if(pl<NPR-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  dualfprintf(fail_file,"norm=%21.15g``20\n",norm);
  dualfprintf(fail_file,"mybsq=%21.15g``20\n",bsq);
  
  dualfprintf(fail_file,"alpha={");
  for(i=1;i<=5;i++){
    dualfprintf(fail_file,"{");
    for(j=1;j<=5;j++){
      dualfprintf(fail_file, "%21.15g``20 ",alpha[i][j]);
      if(j<5)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<5)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  

#if(MCOORD!=CARTMINKMETRIC)
  myii=ptrgeom->i;
  myjj=ptrgeom->j;
  mykk=ptrgeom->k;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif

  
  dualfprintf(fail_file,"myconn={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file,"{");
      for(k=0;k<NDIM;k++){
        dualfprintf(fail_file, "%21.15g``20 ",GLOBALMETMACP0A3(conn,myii,myjj,mykk,i,j,k));
        if(k<NDIM-1)       dualfprintf(fail_file,", ");
      }
      dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  dualfprintf(fail_file,"mygcon={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file, "%21.15g``20 ",ptrgeom->gcon[GIND(i,j)]);
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");
  
  dualfprintf(fail_file,"mygcov={");
  for(i=0;i<NDIM;i++){
    dualfprintf(fail_file,"{");
    for(j=0;j<NDIM;j++){
      dualfprintf(fail_file, "%21.15g``20 ",ptrgeom->gcov[GIND(i,j)]);
      if(j<NDIM-1)       dualfprintf(fail_file,", ");
    }
    dualfprintf(fail_file,"}\n"); // \n breaks line so can copy easily to mathematica
    if(i<NDIM-1)       dualfprintf(fail_file,", ");
  }
  dualfprintf(fail_file,"}\n");

  return(0);
}




// same as faildebug1 but with column formatted output, instead of copy/paste into mathematica format
int faildebug2(int numnormterms, int whichcons, FTYPE *U_target, FTYPE *EOSextra, FTYPE *pr0, struct of_geom *ptrgeom)
{
  FTYPE tolx, tolf;
  int ntrial, mnewtfailed;
  int pl,pliter;
  FTYPE Ustart[NPR];
  FTYPE dUriemann[NPR];
  FTYPE dU[NPR],dUcomp[NUMSOURCES][NPR];
  FTYPE bsq;
  FTYPE mhd[NDIM][NDIM];
  FTYPE ucon[NDIM];
  FTYPE flux[NDIM][NPR];
  FTYPE **alpha;
  struct of_state q;
  int i,j,k;
  FILE *out;
  FTYPE X[NDIM],V[NDIM];
  int myii,myjj,mykk;
  FTYPE Uwithgeom[NPR];
  FTYPE norm;


  coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,X);
  bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,V);

  PLOOP(pliter,pl) dUriemann[pl]=0.0; // fudge


  if (get_state(pr0, ptrgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 2);
  if (primtoU(UNOTHING,pr0, &q, ptrgeom, Ustart, NULL) >= 1)
    FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 2);
  
  mhd_calc(pr0,0, ptrgeom, &q, mhd[0], NULL);
  mhd_calc(pr0,1, ptrgeom, &q, mhd[1], NULL);
  mhd_calc(pr0,2, ptrgeom, &q, mhd[2], NULL);
  mhd_calc(pr0,3, ptrgeom, &q, mhd[3], NULL);
  primtoflux(UNOTHING,pr0, &q, 1,ptrgeom, flux[0], NULL);
  primtoflux(UNOTHING,pr0, &q, 2,ptrgeom, flux[1], NULL);
  primtoflux(UNOTHING,pr0, &q, 3,ptrgeom, flux[2], NULL);
  PLOOP(pliter,pl) Uwithgeom[pl]=Ustart[pl]*ptrgeom->EOMFUNCMAC(pl);
  FTYPE CUf[NUMDTCUFS]={1.0};
  int timeorder=0;
  FTYPE *CUimp=&CUf[NUMPREDTCUFS+timeorder];
  int didreturnpf=0;
  int eomtype=EOMDEFAULT;
  source(pr0, pr0, pr0, &didreturnpf, &eomtype, ptrgeom, &q, Uwithgeom, Uwithgeom, CUf, CUimp, 0.0, dUriemann, dUcomp,dU);
  bsq_calc(pr0, ptrgeom, &bsq);
  
  
  alpha=dmatrix(1, 5, 1, 5);

  if (dudp_calc(WHICHEOS,whichcons,EOSextra, pr0, &q, ptrgeom, alpha) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "dudp_calc()", 1);
  // more general normalization
  // copied from usrfun
#if(NORMMETHOD==-1)
  norm=1.0;
#elif(NORMMETHOD==0)
  norm=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
  norm=0.0;
  numnormterms=0;
  for (j = 0; j < INVERTNPR; j++){
    for (k = 0; k < INVERTNPR; k++){
      if(alpha[j + 1][k + 1]>NUMEPSILON){
        norm+=fabs(alpha[j + 1][k + 1]);
        numnormterms++;
      }
    }
  }
  norm=(FTYPE)(numnormterms)/(norm); // (i.e. inverse of average)
#endif    
  for (j = 0; j < INVERTNPR; j++)
    for (k = 0; k < INVERTNPR; k++)
      alpha[j + 1][k + 1] *= (norm);



  out=fopen("myconsts.txt","wt");  
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",X[1],X[2],X[3],V[1],V[2],V[3],gam,a,hslope,R0);
  fclose(out);

  
  out=fopen("myTud.txt","wt");
  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      fprintf(out,"%21.15g ",mhd[i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  out=fopen("myF.txt","wt");
  for(i=0;i<=2;i++){ // mathematica rows
    for(j=0;j<NPR;j++){ // mathematica columns
      fprintf(out, "%21.15g ",flux[i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  out=fopen("myub.txt","wt");  
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3]);
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3]);
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3]);
  fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3]);
  fclose(out);

  out=fopen("mygnbsq.txt","wt");  
  fprintf(out,"%21.15g %21.15g %21.15g\n",ptrgeom->gdet,norm,bsq);
  fclose(out);
  
  out=fopen("mypU.txt","wt");  
  PLOOP(pliter,pl) fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",pr0[pl],Ustart[pl],U_target[pl],dU[pl]);
  fclose(out);

  
  out=fopen("alpha.txt","wt");  
  for(i=1;i<=5;i++){
    for(j=1;j<=5;j++){
      fprintf(out, "%21.15g ",alpha[i][j]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

#if(MCOORD!=CARTMINKMETRIC)
  myii=ptrgeom->i;
  myjj=ptrgeom->j;
  mykk=ptrgeom->k;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif
  
   
  out=fopen("myconn.txt","wt");  
  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      for(k=0;k<NDIM;k++){
        fprintf(out, "%21.15g ",GLOBALMETMACP0A3(conn,myii,myjj,mykk,i,j,k));
      }
      fprintf(out,"\n");
    }
  }
  fclose(out);
  
  out=fopen("mygcon.txt","wt");  

  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      fprintf(out, "%21.15g ",ptrgeom->gcon[GIND(i,j)]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  out=fopen("mygcov.txt","wt");  

  for(i=0;i<NDIM;i++){
    for(j=0;j<NDIM;j++){
      fprintf(out, "%21.15g ",ptrgeom->gcov[GIND(i,j)]);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  return(0);
}
