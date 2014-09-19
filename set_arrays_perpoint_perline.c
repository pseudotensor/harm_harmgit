/*! \file set_arrays_perpoint_perline.c
     \brief Sets allocation, pointer shift, and dummy assignment for per-point quantities or 1D arrays
*/


#include "decs.h"



/// all arrays that are not multi-dimensional grid arrays (that have no particular physical direction associated with that direction)
/// set_arrays() includes 1D arrays that are always (say) related to physical dimension like radius.  For example, N1 is always related to radius for SPC.
/// The distinction is made in order to LOOP faster over multi-dimensional arrays
void set_arrays_perpoint_perline()
{
  int dissloop;
  extern int init_selfgrav(void);



  interporder = (int (*)) (&(a_interporder[NUMNEGINTERPS]));
  usedqarray = (int (*)) (&(a_usedqarray[NUMNEGINTERPS]));
  
  
  
#if(DOLUMVSR)
  // yes, for each cpu
  lumvsr=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
  if(lumvsr==NULL){
    dualfprintf(fail_file,"Couldn't open lumvsr memory\n");
    myexit(1);
  }
  
  lumvsr_tot=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
  if(lumvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open lumvsr_tot memory\n");
    myexit(1);
  }
#endif
  
#if(DODISSVSR)
  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){
    // yes, for each cpu
    dissvsr[dissloop]=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(dissvsr[dissloop]==NULL){
      dualfprintf(fail_file,"Couldn't open dissvsr memory: %d\n",dissloop);
      myexit(1);
    }
    
    dissvsr_tot[dissloop]=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(dissvsr_tot[dissloop]==NULL){
      dualfprintf(fail_file,"Couldn't open dissvsr_tot memory: %d\n",dissloop);
      myexit(1);
    }
  }
  //for(ii=0;ii<ncpux1*N1;ii++) dissvsr[ii]=0;
  //for(ii=0;ii<ncpux1*N1;ii++) dissvsr_tot[ii]=0;
#endif
  
#if(DOSELFGRAVVSR)
  // yes, for each cpu
  dMvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dMvsr+=N1BND;
  if(dMvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dMvsr memory\n");
    myexit(1);
  }
  
  dMvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dMvsr_tot+=N1BND;
  if(dMvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dMvsr_tot memory\n");
    myexit(1);
  }

  dVvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dVvsr+=N1BND;
  if(dVvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dVvsr memory\n");
    myexit(1);
  }
  
  dVvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dVvsr_tot+=N1BND;
  if(dVvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dVvsr_tot memory\n");
    myexit(1);
  }

  vrsqvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  vrsqvsr+=N1BND;
  if(vrsqvsr==NULL){
    dualfprintf(fail_file,"Couldn't open vrsqvsr memory\n");
    myexit(1);
  }
  
  vrsqvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  vrsqvsr_tot+=N1BND;
  if(vrsqvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open vrsqvsr_tot memory\n");
    myexit(1);
  }

  // yes, for each cpu
  dTrrvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dTrrvsr+=N1BND;
  if(dTrrvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dTrrvsr memory\n");
    myexit(1);
  }
  
  dTrrvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dTrrvsr_tot+=N1BND;
  if(dTrrvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dTrrvsr_tot memory\n");
    myexit(1);
  }

  Mvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Mvsr_tot+=N1BND;
  if(Mvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Mvsr_tot memory\n");
    myexit(1);
  }

  Mvsrface1_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Mvsrface1_tot+=N1BND;
  if(Mvsrface1_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Mvsrface1_tot memory\n");
    myexit(1);
  }

  MOrvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  MOrvsr_tot+=N1BND;
  if(MOrvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open MOrvsr_tot memory\n");
    myexit(1);
  }

  // the potentials are located in ghost cells so outer boundary CPUs have sufficient data
  // +1 corresponds to last outer face, so not used for center phi but allocated to keep arrays same size
  phivsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  phivsr_tot+=N1BND;
  if(phivsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open phivsr_tot memory\n");
    myexit(1);
  }

  rcent=(FTYPE*)malloc(NUMGRAVPOS*sizeof(FTYPE));
  rcent+=N1BND;
  if(rcent==NULL){
    dualfprintf(fail_file,"Couldn't open rcent memory\n");
    myexit(1);
  }

  rcent_tot=(FTYPE*)malloc(NUMGRAVPOS*sizeof(FTYPE));
  rcent_tot+=N1BND;
  if(rcent_tot==NULL){
    dualfprintf(fail_file,"Couldn't open rcent_tot memory\n");
    myexit(1);
  }


  // yes, for each cpu
  dJvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dJvsr+=N1BND;
  if(dJvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dJvsr memory\n");
    myexit(1);
  }
  
  dJvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dJvsr_tot+=N1BND;
  if(dJvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dJvsr_tot memory\n");
    myexit(1);
  }

  Jvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Jvsr_tot+=N1BND;
  if(Jvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Jvsr_tot memory\n");
    myexit(1);
  }

  Jvsrface1_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Jvsrface1_tot+=N1BND;
  if(Jvsrface1_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Jvsrface1_tot memory\n");
    myexit(1);
  }


#endif
  
  if(DOSELFGRAVVSR){
    // initialize self-gravity functions since not set yet and will otherwise be used to set first metric
    init_selfgrav();
  }

}
