
/*! \file mpi_grmhd_grray_liaison.c
     \brief Genreal MPI hooks to allow GRMHD or GRRAY+GRMHD to run together
*/


#include "decs.h"





#if(USEMPIANY)
/////////////////////////////
//
// Include these declarations and functions only if any code that uses this function is doing MPI
//
/////////////////////////////


///////////////
/// declarations
static void get_processtypelist(int processtype, MPI_Comm localcomm, int **local_processtypelistlocal, int *local_sizeproclist);
static void liaison_init_mpi_processtypes(MPI_Comm localcomm, int **processtypelistlocal_world, int *local_sizeproclist);
static void grmhd_init_mpi_processtypes(MPI_Comm localcomm, int **processtypelistlocal_world, int *local_sizeproclist);
static void grray_init_mpi_processtypes(MPI_Comm localcomm, int **processtypelistlocal_world, int *local_sizeproclist);


/// wrapper for liaison code version of get_processtypelist()
void liaison_init_mpi_processtypes(MPI_Comm localcomm, int **processtypelistlocal_world, int *local_sizeproclist)
{
  get_processtypelist(BINARYLIAISONTYPE,localcomm,processtypelistlocal_world,local_sizeproclist);

}

/// wrapper for grmhd code version of get_processtypelist()
void grmhd_init_mpi_processtypes(MPI_Comm localcomm, int **processtypelistlocal_world, int *local_sizeproclist)
{
  get_processtypelist(BINARYGRMHDTYPE,localcomm,processtypelistlocal_world,local_sizeproclist);

}

/// wrapper for grray code version of get_processtypelist()
void grray_init_mpi_processtypes(MPI_Comm localcomm, int **processtypelistlocal_world, int *local_sizeproclist)
{
  get_processtypelist(BINARYGRRAYTYPE,localcomm,processtypelistlocal_world,local_sizeproclist);

}



/// get processor list for a single input processtype per process
/// assumes all processes call this function together
void get_processtypelist(int processtype, MPI_Comm localcomm, int **processtypelistlocal, int *sizeproclistlocal)
{
  int i;


  // may already have this, but get again to keep code modular
  MPI_Comm_size(localcomm, sizeproclistlocal);
  stderrfprintf("processtype=%d sizeproclistlocal=%d\n",processtype,*sizeproclistlocal); fflush(stderr);

  //  return;

  // allocate processor list 
  *processtypelistlocal=(int*) malloc((*sizeproclistlocal)*sizeof(int));
  if(*processtypelistlocal==NULL){
    stderrfprintf("Could not allocate processtypelistlocal\n");
    exit(1);
  }
  stderrfprintf("After malloc of processtypelistlocal: %d %d\n",processtype,myid_world); fflush(stderr);

  // Gather process type list
  MPI_Allgather(&processtype,1,MPI_INT,*processtypelistlocal,1,MPI_INT,localcomm);
  stderrfprintf("After MPI_Allgather: %d\n",myid_world); fflush(stderr);
  //  MPI_Barrier(localcomm);

  
  //  stderrfprintf("processtype=%d sizeproclistlocal0=%d\n",processtype,*sizeproclistlocal);
  //  for(i=0;i<*sizeproclistlocal;i++){
  //    stderrfprintf("it0[%d]=%d\n",i,(*processtypelistlocal)[i]); fflush(stderr);
  //  }


}


#endif // end if any code is doing MPI










#if(USEMPILIAISON)
/////////////////////////////
//
// Below functions are used if really doing liaison method with MPI
//
/////////////////////////////


/// declaration
static void init_MPI_group_grmhd_grray_liaison(
                                               int *processtypelistlocal_world,int sizeproclistlocal_world,
                                               MPI_Group *MPI_GROUP_WORLD,
                                               MPI_Group *MPI_GROUP_GRMHD_LIAISON, MPI_Comm *MPI_COMM_GRMHD_LIAISON,
                                               MPI_Group *MPI_GROUP_GRRAY_LIAISON, MPI_Comm *MPI_COMM_GRRAY_LIAISON);




/// After MPI_init() this must be called by liaison code
/// assumes global variables used
/// just a wrapper
void liaison_init_mpi_liaisonmode_globalset(void)
{
  void liaison_init_mpi_liaisonmode(
                                    int **processtypelistlocal_world,int *sizeproclistlocal_world,
                                    int **processtypelistlocal_grmhd_liaison,int *sizeproclistlocal_grmhd_liaison,
                                    int **processtypelistlocal_grray_liaison,int *sizeproclistlocal_grray_liaison,
                                    int **processtypelistlocal_grmhd,int *sizeproclistlocal_grmhd,
                                    int **processtypelistlocal_grray,int *sizeproclistlocal_grray,
                                    int **processtypelistlocal_liaison_from_grmhd,int *sizeproclistlocal_liaison_from_grmhd,
                                    int **processtypelistlocal_liaison_from_grray,int *sizeproclistlocal_liaison_from_grray,
                                    MPI_Group *MPI_GROUP_WORLD,
                                    MPI_Group *MPI_GROUP_GRMHD_LIAISON, MPI_Comm *MPI_COMM_GRMHD_LIAISON,
                                    MPI_Group *MPI_GROUP_GRRAY_LIAISON, MPI_Comm *MPI_COMM_GRRAY_LIAISON,
                                    MPI_Group *MPI_GROUP_GRMHD, MPI_Comm *MPI_COMM_GRMHD,
                                    MPI_Group *MPI_GROUP_GRRAY, MPI_Comm *MPI_COMM_GRRAY,
                                    MPI_Group *MPI_GROUP_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LIAISON_FROM_GRMHD, 
                                    MPI_Group *MPI_GROUP_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LIAISON_FROM_GRRAY);


  liaison_init_mpi_liaisonmode(
                               &processtypelist_world,&sizeproclist_world,
                               &processtypelist_grmhd_liaison,&sizeproclist_grmhd_liaison,
                               &processtypelist_grray_liaison,&sizeproclist_grray_liaison,
                               &processtypelist_grmhd,&sizeproclist_grmhd,
                               &processtypelist_grray,&sizeproclist_grray,
                               &processtypelist_liaison_from_grmhd,&sizeproclist_grmhd,
                               &processtypelist_liaison_from_grray,&sizeproclist_grray,
                               &MPI_GROUP_WORLD,
                               &MPI_GROUP_GRMHD_LIAISON, &MPI_COMM_GRMHD_LIAISON,
                               &MPI_GROUP_GRRAY_LIAISON, &MPI_COMM_GRRAY_LIAISON,
                               &MPI_GROUP_GRMHD, &MPI_COMM_GRMHD,
                               &MPI_GROUP_GRRAY, &MPI_COMM_GRRAY,
                               &MPI_GROUP_LIAISON_FROM_GRMHD, &MPI_COMM_LIAISON_FROM_GRMHD, 
                               &MPI_GROUP_LIAISON_FROM_GRRAY, &MPI_COMM_LIAISON_FROM_GRRAY);
}





/// After MPI_init() this must be called by grmhd code
/// assumes global variables used
/// just a wrapper
void grmhd_init_mpi_liaisonmode_globalset(void)
{
  void grmhd_init_mpi_liaisonmode(
                                  int **processtypelistlocal_world,int *sizeproclistlocal_world,
                                  int **processtypelistlocal_grmhd_liaison,int *sizeproclistlocal_grmhd_liaison,
                                  int **processtypelistlocal_grray_liaison,int *sizeproclistlocal_grray_liaison,
                                  int **processtypelistlocal_grmhd,int *sizeproclistlocal_grmhd,
                                  int **processtypelistlocal_liaison_from_grmhd,int *sizeproclistlocal_liaison_from_grmhd,
                                  MPI_Group *MPI_GROUP_WORLD,
                                  MPI_Group *MPI_GROUP_GRMHD_LIAISON, MPI_Comm *MPI_COMM_GRMHD_LIAISON,
                                  MPI_Group *MPI_GROUP_GRRAY_LIAISON, MPI_Comm *MPI_COMM_GRRAY_LIAISON,
                                  MPI_Group *MPI_GROUP_GRMHD, MPI_Comm *MPI_COMM_GRMHD,
                                  MPI_Group *MPI_GROUP_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LIAISON_FROM_GRMHD);

  grmhd_init_mpi_liaisonmode(
                             &processtypelist_world,&sizeproclist_world,
                             &processtypelist_grmhd_liaison,&sizeproclist_grmhd_liaison,
                             &processtypelist_grray_liaison,&sizeproclist_grray_liaison,
                             &processtypelist_grmhd,&sizeproclist_grmhd,
                             &processtypelist_liaison_from_grmhd,&sizeproclist_grmhd,
                             &MPI_GROUP_WORLD,
                             &MPI_GROUP_GRMHD_LIAISON, &MPI_COMM_GRMHD_LIAISON,
                             &MPI_GROUP_GRRAY_LIAISON, &MPI_COMM_GRRAY_LIAISON,
                             &MPI_GROUP_GRMHD, &MPI_COMM_GRMHD,
                             &MPI_GROUP_LIAISON_FROM_GRMHD, &MPI_COMM_LIAISON_FROM_GRMHD);
}







/// After MPI_init() this must be called by grray code
/// assumes global variables used
/// just a wrapper
void grray_init_mpi_liaisonmode_globalset(void)
{
  void grray_init_mpi_liaisonmode(
                                  int **processtypelistlocal_world,int *sizeproclistlocal_world,
                                  int **processtypelistlocal_grmhd_liaison,int *sizeproclistlocal_grmhd_liaison,
                                  int **processtypelistlocal_grray_liaison,int *sizeproclistlocal_grray_liaison,
                                  int **processtypelistlocal_grray,int *sizeproclistlocal_grray,
                                  int **processtypelistlocal_liaison_from_grray,int *sizeproclistlocal_liaison_from_grray,
                                  MPI_Group *MPI_GROUP_WORLD,
                                  MPI_Group *MPI_GROUP_GRMHD_LIAISON, MPI_Comm *MPI_COMM_GRMHD_LIAISON,
                                  MPI_Group *MPI_GROUP_GRRAY_LIAISON, MPI_Comm *MPI_COMM_GRRAY_LIAISON,
                                  MPI_Group *MPI_GROUP_GRRAY, MPI_Comm *MPI_COMM_GRRAY,
                                  MPI_Group *MPI_GROUP_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LIAISON_FROM_GRRAY);

  grray_init_mpi_liaisonmode(
                             &processtypelist_world,&sizeproclist_world,
                             &processtypelist_grmhd_liaison,&sizeproclist_grmhd_liaison,
                             &processtypelist_grray_liaison,&sizeproclist_grray_liaison,
                             &processtypelist_grray,&sizeproclist_grray,
                             &processtypelist_liaison_from_grray,&sizeproclist_liaison_from_grray,
                             &MPI_GROUP_WORLD,
                             &MPI_GROUP_GRMHD_LIAISON, &MPI_COMM_GRMHD_LIAISON,
                             &MPI_GROUP_GRRAY_LIAISON, &MPI_COMM_GRRAY_LIAISON,
                             &MPI_GROUP_GRRAY, &MPI_COMM_GRRAY,
                             &MPI_GROUP_LIAISON_FROM_GRRAY, &MPI_COMM_LIAISON_FROM_GRRAY);

}






/// After MPI_init() this must be called by liaison code if want to avoid global variables
void liaison_init_mpi_liaisonmode(
                                  int **processtypelistlocal_world,int *sizeproclistlocal_world,
                                  int **processtypelistlocal_grmhd_liaison,int *sizeproclistlocal_grmhd_liaison,
                                  int **processtypelistlocal_grray_liaison,int *sizeproclistlocal_grray_liaison,
                                  int **processtypelistlocal_grmhd,int *sizeproclistlocal_grmhd,
                                  int **processtypelistlocal_grray,int *sizeproclistlocal_grray,
                                  int **processtypelistlocal_liaison_from_grmhd,int *sizeproclistlocal_liaison_from_grmhd,
                                  int **processtypelistlocal_liaison_from_grray,int *sizeproclistlocal_liaison_from_grray,
                                  MPI_Group *MPI_GROUP_LOCAL_WORLD,
                                  MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                  MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                  MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                  MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                  MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD, 
                                  MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY)
{
  void init_MPI_group_grmhd_grray_liaison_split(
                                                int *processtypelistlocal_grmhd_liaison,int sizeproclistlocal_grmhd_liaison,
                                                int *processtypelistlocal_grray_liaison,int sizeproclistlocal_grray_liaison,
                                                MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                                MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                                MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                                MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                                MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD, 
                                                MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY);

#if(USEMPILIAISON)


  // get processtypes
  liaison_init_mpi_processtypes(MPI_COMM_WORLD,processtypelistlocal_world,sizeproclistlocal_world);

  // global group sets:
  init_MPI_group_grmhd_grray_liaison(
                                     *processtypelistlocal_world,*sizeproclistlocal_world,
                                     MPI_GROUP_LOCAL_WORLD,
                                     MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_COMM_LOCAL_GRMHD_LIAISON,
                                     MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_COMM_LOCAL_GRRAY_LIAISON);

  // get processtypes for newly created subgroups
  // 1) Get GRMHD+LIAISON
  liaison_init_mpi_processtypes(*MPI_COMM_LOCAL_GRMHD_LIAISON,processtypelistlocal_grmhd_liaison,sizeproclistlocal_grmhd_liaison);
  // 2) Get GRRAY+LIAISON
  liaison_init_mpi_processtypes(*MPI_COMM_LOCAL_GRRAY_LIAISON,processtypelistlocal_grray_liaison,sizeproclistlocal_grray_liaison);

  // liaison split from liaison+grmhd and liaison+grray groups
  init_MPI_group_grmhd_grray_liaison_split(
                                           *processtypelistlocal_grmhd_liaison,*sizeproclistlocal_grmhd_liaison,
                                           *processtypelistlocal_grray_liaison,*sizeproclistlocal_grray_liaison,
                                           MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_COMM_LOCAL_GRMHD_LIAISON,
                                           MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_COMM_LOCAL_GRRAY_LIAISON,
                                           MPI_GROUP_LOCAL_GRMHD, MPI_COMM_LOCAL_GRMHD,
                                           MPI_GROUP_LOCAL_GRRAY, MPI_COMM_LOCAL_GRRAY,
                                           MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_COMM_LOCAL_LIAISON_FROM_GRMHD,
                                           MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_COMM_LOCAL_LIAISON_FROM_GRRAY);

  // ensure NULL
  if(*MPI_COMM_LOCAL_GRMHD!=MPI_COMM_NULL){
    stderrfprintf("MPI_COMM_LOCAL_GRMHD still has members of LIAISON!\n"); fflush(stderr);
    exit(1);
  }
  // ensure NULL
  if(*MPI_COMM_LOCAL_GRRAY!=MPI_COMM_NULL){
    stderrfprintf("MPI_COMM_LOCAL_GRRAY still has members of LIAISON!\n"); fflush(stderr);
    exit(1);
  }

  // get processtypes for newly created subgroups
  // 1) Do LIAISON from GRMHD+LIAISON
  liaison_init_mpi_processtypes(*MPI_COMM_LOCAL_LIAISON_FROM_GRMHD,processtypelistlocal_liaison_from_grmhd,sizeproclistlocal_liaison_from_grmhd);

  // 2) Do LIAISON from GRRAY+LIAISON
  liaison_init_mpi_processtypes(*MPI_COMM_LOCAL_LIAISON_FROM_GRRAY,processtypelistlocal_liaison_from_grray,sizeproclistlocal_liaison_from_grray);


#endif


}






/// After MPI_init() this must be called by grmhd code if want to avoid global variables
void grmhd_init_mpi_liaisonmode(
                                int **processtypelistlocal_world,int *sizeproclistlocal_world,
                                int **processtypelistlocal_grmhd_liaison,int *sizeproclistlocal_grmhd_liaison,
                                int **processtypelistlocal_grray_liaison,int *sizeproclistlocal_grray_liaison,
                                int **processtypelistlocal_grmhd,int *sizeproclistlocal_grmhd,
                                int **processtypelistlocal_liaison_from_grmhd,int *sizeproclistlocal_liaison_from_grmhd,
                                MPI_Group *MPI_GROUP_LOCAL_WORLD,
                                MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD)
{
  void init_MPI_group_grmhd_liaison_split(
                                          int *processtypelistlocal_grmhd_liaison,int sizeproclistlocal_grmhd_liaison,
                                          MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                          MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                          MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD);


  stderrfprintf("USEMPILIAISON=%d\n",USEMPILIAISON);

#if(USEMPILIAISON)



  // get processtypes
  grmhd_init_mpi_processtypes(MPI_COMM_WORLD,processtypelistlocal_world,sizeproclistlocal_world);


  // global group sets:
  init_MPI_group_grmhd_grray_liaison(
                                     *processtypelistlocal_world,*sizeproclistlocal_world,
                                     MPI_GROUP_LOCAL_WORLD,
                                     MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_COMM_LOCAL_GRMHD_LIAISON,
                                     MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_COMM_LOCAL_GRRAY_LIAISON);

  // ensure really NULL
  if(*MPI_COMM_LOCAL_GRRAY_LIAISON!=MPI_COMM_NULL){
    stderrfprintf("GRRAY_LIAISON still has members of GRMHD!\n"); fflush(stderr);
    exit(1);
  }

  // get processtypes for newly created subgroups
  // 1) Get GRMHD+LIAISON proc type list
  grmhd_init_mpi_processtypes(*MPI_COMM_LOCAL_GRMHD_LIAISON,processtypelistlocal_grmhd_liaison,sizeproclistlocal_grmhd_liaison);

  // liaison split from liaison+grmhd group
  init_MPI_group_grmhd_liaison_split(
                                     *processtypelistlocal_grmhd_liaison,*sizeproclistlocal_grmhd_liaison,
                                     MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_COMM_LOCAL_GRMHD_LIAISON,
                                     MPI_GROUP_LOCAL_GRMHD, MPI_COMM_LOCAL_GRMHD,
                                     MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_COMM_LOCAL_LIAISON_FROM_GRMHD);

  // ensure really NULL
  if(*MPI_COMM_LOCAL_LIAISON_FROM_GRMHD!=MPI_COMM_NULL){
    stderrfprintf("MPI_COMM_LOCAL_LIAISON_FROM_GRMHD still has members of GRMHD!\n"); fflush(stderr);
    exit(1);
  }

  // get processtypes for newly created subgroups
  // 1) Get GRMHD proc type list
  grmhd_init_mpi_processtypes(*MPI_COMM_LOCAL_GRMHD,processtypelistlocal_grmhd,sizeproclistlocal_grmhd);



#endif


}



/// After MPI_init() this must be called by grray code if want to avoid global variables
void grray_init_mpi_liaisonmode(
                                int **processtypelistlocal_world,int *sizeproclistlocal_world,
                                int **processtypelistlocal_grmhd_liaison,int *sizeproclistlocal_grmhd_liaison,
                                int **processtypelistlocal_grray_liaison,int *sizeproclistlocal_grray_liaison,
                                int **processtypelistlocal_grray,int *sizeproclistlocal_grray,
                                int **processtypelistlocal_liaison_from_grray,int *sizeproclistlocal_liaison_from_grray,
                                MPI_Group *MPI_GROUP_LOCAL_WORLD,
                                MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY)
{
  void init_MPI_group_grray_liaison_split(
                                          int *processtypelistlocal_grray_liaison,int sizeproclistlocal_grray_liaison,
                                          MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                          MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                          MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY);

#if(USEMPILIAISON)


  // get processtypes
  grray_init_mpi_processtypes(MPI_COMM_WORLD,processtypelistlocal_world,sizeproclistlocal_world);


  // global group sets:
  init_MPI_group_grmhd_grray_liaison(
                                     *processtypelistlocal_world,*sizeproclistlocal_world,
                                     MPI_GROUP_LOCAL_WORLD,
                                     MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_COMM_LOCAL_GRMHD_LIAISON,
                                     MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_COMM_LOCAL_GRRAY_LIAISON);

  // ensure really NULL
  if(*MPI_COMM_LOCAL_GRMHD_LIAISON!=MPI_COMM_NULL){
    stderrfprintf("GRMHD_LIAISON still has members of GRRAY!\n"); fflush(stderr);
    exit(1);
  }

  // get processtypes for newly created subgroups
  // 1) Get GRRAY+LIAISON proc type list
  grray_init_mpi_processtypes(*MPI_COMM_LOCAL_GRRAY_LIAISON,processtypelistlocal_grray_liaison,sizeproclistlocal_grray_liaison);

  // liaison split from liaison+grray group
  init_MPI_group_grray_liaison_split(
                                     *processtypelistlocal_grray_liaison,*sizeproclistlocal_grray_liaison,
                                     MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_COMM_LOCAL_GRRAY_LIAISON,
                                     MPI_GROUP_LOCAL_GRRAY, MPI_COMM_LOCAL_GRRAY,
                                     MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_COMM_LOCAL_LIAISON_FROM_GRRAY);

  // ensure really NULL
  if(*MPI_COMM_LOCAL_LIAISON_FROM_GRRAY!=MPI_COMM_NULL){
    stderrfprintf("MPI_COMM_LOCAL_LIAISON_FROM_GRRAY still has members of GRRAY!\n"); fflush(stderr);
    exit(1);
  }

  // get processtypes for newly created subgroups
  // 1) Get GRRAY proc type list
  grray_init_mpi_processtypes(*MPI_COMM_LOCAL_GRRAY,processtypelistlocal_grray,sizeproclistlocal_grray);



#endif


}










/// sub function (that all processes end up calling) in order to determine groups and communicators that split MPI_COMM_WORLD into sub groups
/// create sub groups:
/// 1) GRMHD+LIAISON from MPI_COMM_WORLD
/// 2) GRRAY+LIAISON from MPI_COMM_WORLD
void init_MPI_group_grmhd_grray_liaison(
                                        int *processtypelistlocal,int sizeproclistlocal,
                                        MPI_Group *MPI_GROUP_LOCAL_WORLD,
                                        MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                        MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON)
{
  int *ranks; 
  int i,j,k,numranks;
  int procsdone;
  int sizetemp;

  // allocate things that are truenumprocs in size
  ranks=(int*)malloc(sizeof(int)*truenumprocs);
  if(ranks==NULL){
    stderrfprintf("Problem allocating memory for ranks with truenumprocs=%d\n",truenumprocs); fflush(stderr);
    myexit(3876252356);
  }
  for(i=0;i<truenumprocs;i++) ranks[i]=0;


#if(USEMPILIAISON)  

  ////////////////////////////////
  // get group for MPI_COMM_WORLD
  MPI_Comm_group(MPI_COMM_WORLD, MPI_GROUP_LOCAL_WORLD);


  //  stderrfprintf("sizeproclistlocal=%d\n",sizeproclistlocal);
  //  for(i=0;i<sizeproclistlocal;i++){
  //    stderrfprintf("it[%d]=%d\n",i,processtypelistlocal[i]); fflush(stderr);
  //  }


  ////////////////////////////////
  // create GRMHD+LIAISON GROUP ensuring GRMHD's are ordered as originally (in case ordered at startup) and LIAISON's at end
  j=0;
  for(i=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYGRMHDTYPE){
      ranks[j]=i;
      j++;
    }
  }
  for(i=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYLIAISONTYPE){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  

  // now create group and communicator
  // Below "ranks" is more primitive than GRMHD code ranks, so below shouldn't be mapped
  MPI_Group_incl(*MPI_GROUP_LOCAL_WORLD, numranks, ranks, MPI_GROUP_LOCAL_GRMHD_LIAISON);
  MPI_Comm_create(MPI_COMM_WORLD, *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_COMM_LOCAL_GRMHD_LIAISON); 



  ////////////////////////////////
  // create GRRAY+LIAISON GROUP ensuring GRRAY's are ordered as originally (in case ordered at startup) and LIAISON's at end
  j=0;
  for(i=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYGRRAYTYPE){
      ranks[j]=i;
      j++;
    }
  }
  for(i=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYLIAISONTYPE){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;

  // now create group and communicator
  MPI_Group_incl(*MPI_GROUP_LOCAL_WORLD, numranks, ranks, MPI_GROUP_LOCAL_GRRAY_LIAISON);
  MPI_Comm_create(MPI_COMM_WORLD, *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_COMM_LOCAL_GRRAY_LIAISON);

  //  stderrfprintf("numranks=%d\n",numranks);
  //  for(i=0;i<numranks;i++){
  //    stderrfprintf("ranks[%d]=%d\n",i,ranks[i]); fflush(stderr);
  //  }


#endif // end if USEMPILIAISON

  free(ranks);

}





/// create sub groups (only liaison code calls this full function)
/// 1) GRMHD from GRMHD+LIAISON
/// 2) LIAISON_FROM_GRMHD from GRMHD+LIAISON
/// 3) GRRAY from GRRAY+LIAISON
/// 4) LIAISON_FROM_GRRAY from GRRAY+LIAISON
void init_MPI_group_grmhd_grray_liaison_split(
                                              int *processtypelistlocal_grmhd_liaison,int sizeproclistlocal_grmhd_liaison,
                                              int *processtypelistlocal_grray_liaison,int sizeproclistlocal_grray_liaison,
                                              MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                              MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                              MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                              MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                              MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD, 
                                              MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY)
{

  void init_MPI_group_grmhd_liaison_split(
                                          int *processtypelistlocal_grmhd_liaison,int sizeproclistlocal_grmhd_liaison,
                                          MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                          MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                          MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD);
  void init_MPI_group_grray_liaison_split(
                                          int *processtypelistlocal_grray_liaison,int sizeproclistlocal_grray_liaison,
                                          MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                          MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                          MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY);


#if(USEMPILIAISON)

  // split GRMHD+LIAISON (uses same function used by GRMHD code)
  init_MPI_group_grmhd_liaison_split(processtypelistlocal_grmhd_liaison,sizeproclistlocal_grmhd_liaison,
                                     MPI_GROUP_LOCAL_GRMHD_LIAISON,  MPI_COMM_LOCAL_GRMHD_LIAISON,
                                     MPI_GROUP_LOCAL_GRMHD,  MPI_COMM_LOCAL_GRMHD,
                                     MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD,  MPI_COMM_LOCAL_LIAISON_FROM_GRMHD);


  // split GRRAY+LIAISON (uses same function used by GRRAY code)
  init_MPI_group_grray_liaison_split(processtypelistlocal_grray_liaison,sizeproclistlocal_grray_liaison,
                                     MPI_GROUP_LOCAL_GRRAY_LIAISON,  MPI_COMM_LOCAL_GRRAY_LIAISON,
                                     MPI_GROUP_LOCAL_GRRAY,  MPI_COMM_LOCAL_GRRAY,
                                     MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY,  MPI_COMM_LOCAL_LIAISON_FROM_GRRAY);
  


#endif


}



/// create sub groups (grmhd code and liaison code will call this)
/// 1,2) GRMHD and LIAISON_FROM_GRMHD from GRMHD+LIAISON
void init_MPI_group_grmhd_liaison_split(
                                        int *processtypelistlocal,int sizeproclistlocal,
                                        MPI_Group *MPI_GROUP_LOCAL_GRMHD_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRMHD_LIAISON,
                                        MPI_Group *MPI_GROUP_LOCAL_GRMHD, MPI_Comm *MPI_COMM_LOCAL_GRMHD,
                                        MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRMHD)
{
  int *ranks;
  int i,j,k,numranks;
  int procsdone;


  // allocate things that are truenumprocs in size
  ranks=(int*)malloc(sizeof(int)*truenumprocs);
  if(ranks==NULL){
    stderrfprintf("Problem allocating memory for ranks with truenumprocs=%d\n",truenumprocs); fflush(stderr);
    myexit(3876252356);
  }
  for(i=0;i<truenumprocs;i++) ranks[i]=0;


#if(USEMPILIAISON)


  ////////////////////////////////
  // decompose GRMHD+LIAISON GROUP
  for(i=0,j=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYGRMHDTYPE){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;

  // now create group and communicator
  MPI_Group_incl(*MPI_GROUP_LOCAL_GRMHD_LIAISON, numranks, ranks, MPI_GROUP_LOCAL_GRMHD);
  MPI_Comm_create(*MPI_COMM_LOCAL_GRMHD_LIAISON, *MPI_GROUP_LOCAL_GRMHD, MPI_COMM_LOCAL_GRMHD);



  ////////////////////////////////
  // decompose GRRAY+LIAISON GROUP
  for(i=0,j=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYLIAISONTYPE){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;

  // now create group and communicator
  MPI_Group_incl(*MPI_GROUP_LOCAL_GRMHD_LIAISON, numranks, ranks, MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD);
  MPI_Comm_create(*MPI_COMM_LOCAL_GRMHD_LIAISON, *MPI_GROUP_LOCAL_LIAISON_FROM_GRMHD, MPI_COMM_LOCAL_LIAISON_FROM_GRMHD); 


#endif

  free(ranks);

}



/// create sub groups (grray code and liaison code will call this)
/// 1,2) GRRAY and LIAISON_FROM_GRRAY from GRRAY+LIAISON
void init_MPI_group_grray_liaison_split(
                                        int *processtypelistlocal,int sizeproclistlocal,
                                        MPI_Group *MPI_GROUP_LOCAL_GRRAY_LIAISON, MPI_Comm *MPI_COMM_LOCAL_GRRAY_LIAISON,
                                        MPI_Group *MPI_GROUP_LOCAL_GRRAY, MPI_Comm *MPI_COMM_LOCAL_GRRAY,
                                        MPI_Group *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_Comm *MPI_COMM_LOCAL_LIAISON_FROM_GRRAY)
{
  int *ranks;
  int i,j,k,numranks;
  int procsdone;


  // allocate things that are truenumprocs in size
  ranks=(int*)malloc(sizeof(int)*truenumprocs);
  if(ranks==NULL){
    stderrfprintf("Problem allocating memory for ranks with truenumprocs=%d\n",truenumprocs); fflush(stderr);
    myexit(3876252356);
  }
  for(i=0;i<truenumprocs;i++) ranks[i]=0;



#if(USEMPILIAISON)


  ////////////////////////////////
  // decompose GRRAY+LIAISON GROUP
  for(i=0,j=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYGRRAYTYPE){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;

  // now create group and communicator
  MPI_Group_incl(*MPI_GROUP_LOCAL_GRRAY_LIAISON, numranks, ranks, MPI_GROUP_LOCAL_GRRAY);
  MPI_Comm_create(*MPI_COMM_LOCAL_GRRAY_LIAISON, *MPI_GROUP_LOCAL_GRRAY, MPI_COMM_LOCAL_GRRAY);



  ////////////////////////////////
  // decompose GRRAY+LIAISON GROUP
  for(i=0,j=0;i<sizeproclistlocal;i++){
    if(processtypelistlocal[i]==BINARYLIAISONTYPE){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;

  // now create group and communicator
  MPI_Group_incl(*MPI_GROUP_LOCAL_GRRAY_LIAISON, numranks, ranks, MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY);
  MPI_Comm_create(*MPI_COMM_LOCAL_GRRAY_LIAISON, *MPI_GROUP_LOCAL_LIAISON_FROM_GRRAY, MPI_COMM_LOCAL_LIAISON_FROM_GRRAY); 

#endif


  free(ranks);

}


// when using these communicators, must make sure the call to a communication using it isn't done by a non-member cpu!
// (above: stupid, I know, should just skip if non-member cpu tries a function)

//void init_MPI_group_liaison_free(void)
//{
//      MPI_Comm_free(&combound[i]); // messy since makes nonmember have NULL comm.  Should make as is, and if not member, then skip
//    MPI_Group_free(&grprem[i]);
//}






/// method to exit MPI run
int final_myexit(void)
{
#if(USEMPI)
  // must abort since no clear to communicate to other cpus
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif

  stderrfprintf( "END\n");
  fflush(stderr);
  exit(0);
    
  return(0);
}















#else
/////////////////////////////
//
// NOT normal liaison mode (for normal GRMHD or GRRAY run or debugging)
//
/////////////////////////////





/// method to exit MPI run
/// can exit normally for normal uncoupled run
int final_myexit(void)
{
#if(USEMPI)
  // finish up MPI
  // Barrier required
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  
  stderrfprintf( "END\n");
  fflush(stderr);
  exit(0);
  
  
  return(0);
}







#if(DOINGLIAISONTYPECODE)


/// what to do when no MPI or liaison mode but want to set-up liaison
void liaison_init_mpi_liaisonmode_globalset(void)
{


#if(USEMPILIAISON)

  ////////////////////////////////
  // get group for MPI_COMM_WORLD
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  MPI_COMM_LIAISON_FROM_GRMHD=MPI_COMM_WORLD;
  MPI_GROUP_LIAISON_FROM_GRMHD=MPI_GROUP_WORLD;

  liaison_init_mpi_processtypes(MPI_COMM_WORLD, &processtypelist_world, &sizeproclist_world);
  liaison_init_mpi_processtypes(MPI_COMM_LIAISON_FROM_GRMHD, &processtypelist_liaison_from_grmhd, &sizeproclist_liaison_from_grmhd);

#else
  sizeproclist_world=sizeproclist_liaison_from_grmhd=1;
  processtypelist_world=(int*)malloc(sizeof(int)*1);
  processtypelist_liaison_from_grmhd=(int*)malloc(sizeof(int)*1);
  if(processtypelist_world==NULL || processtypelist_liaison_from_grmhd==NULL){
    stderrfprintf("Couldn't allocate memory for processtypelist_world,liaison_from_grmhd\n");
    exit(1);
  }
  processtypelist_world[0]=BINARYLIAISONTYPE;
  processtypelist_liaison_from_grmhd[0]=BINARYLIAISONTYPE;
#endif

  // needed for parts of liaison code that uses set-up of grid in grmhd code
  sizeproclist_grmhd=sizeproclist_world;

#if(USEMPILIAISON||USEMPIGRMHD)
  MPI_GROUP_GRMHD=MPI_GROUP_LIAISON_FROM_GRMHD;
  MPI_COMM_GRMHD=MPI_COMM_LIAISON_FROM_GRMHD;
#endif

  
}




#endif // end if liaison type code



#if(DOINGGRMHDTYPECODE)

/// NON-LIAISON GRMHD VERSION
/// Only those things needed by GRMHD code without liaison
void grmhd_init_mpi_liaisonmode_globalset(void)
{


#if(USEMPIGRMHD)

  ////////////////////////////////
  // get group for MPI_COMM_WORLD
  stderrfprintf("MPICOMM1\n");
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  stderrfprintf("MPICOMM2\n");
  MPI_COMM_GRMHD=MPI_COMM_WORLD;
  stderrfprintf("MPICOMM3\n");
  MPI_GROUP_GRMHD=MPI_GROUP_WORLD;

  stderrfprintf("MPICOMM4\n");
  grmhd_init_mpi_processtypes(MPI_COMM_WORLD, &processtypelist_world, &sizeproclist_world);
  stderrfprintf("MPICOMM5\n");
  grmhd_init_mpi_processtypes(MPI_COMM_GRMHD, &processtypelist_grmhd, &sizeproclist_grmhd);
  stderrfprintf("MPICOMM6\n");

#else
  sizeproclist_world=sizeproclist_grmhd=1;
  processtypelist_world=(int*)malloc(sizeof(int)*1);
  processtypelist_grmhd=(int*)malloc(sizeof(int)*1);
  if(processtypelist_world==NULL || processtypelist_grmhd==NULL){
    stderrfprintf("Couldn't allocate memory for processtypelist_world,grmhd\n");
    exit(1);
  }
  processtypelist_world[0]=BINARYGRMHDTYPE;
  processtypelist_grmhd[0]=BINARYGRMHDTYPE;
#endif
  
}


#endif // end if grmhd type code



#if(DOINGGRRAYTYPECODE)

/// NON-LIAISON GRRAY VERSION
/// Only those things needed by GRRAY code without liaison
void grray_init_mpi_liaisonmode_globalset(void)
{


#if(USEMPIGRRAY)

  ////////////////////////////////
  // get group for MPI_COMM_WORLD
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  MPI_COMM_GRRAY=MPI_COMM_WORLD;
  MPI_GROUP_GRRAY=MPI_GROUP_WORLD;

  grray_init_mpi_processtypes(MPI_COMM_WORLD, &processtypelist_world, &sizeproclist_world);
  grray_init_mpi_processtypes(MPI_COMM_GRRAY, &processtypelist_grray, &sizeproclist_grray);

#else
  sizeproclist_world=sizeproclist_grray=1;
  processtypelist_world=(int*)malloc(sizeof(int)*1);
  processtypelist_grray=(int*)malloc(sizeof(int)*1);
  if(processtypelist_world==NULL || processtypelist_grray==NULL){
    stderrfprintf("Couldn't allocate memory for processtypelist_world,grray\n");
    exit(1);
  }
  processtypelist_world[0]=BINARYGRRAYTYPE;
  processtypelist_grray[0]=BINARYGRRAYTYPE;
#endif
  
}



#endif // end if grray type code


#endif // end if not doing any liaison mode



