
/*! \file mpidefs.mpi_grmhd_grray_liaison.h
     \brief Variable definitions for mpi_grmhd_grray_liaison.c
*/

#if(DOINGLIAISON)

#if(USEMPILIAISON || USEMPIGRMHD)
//////////
//
// group and communicators needed by all codes
//
//////////
MPI_Group MPI_GROUP_WORLD;
//MPI_Comm MPI_COMM_WORLD;
MPI_Group MPI_GROUP_GRMHD_LIAISON;
MPI_Comm MPI_COMM_GRMHD_LIAISON;

//////////
//
// liaison and GRMHD specific groups
//
/////////
MPI_Group MPI_GROUP_GRMHD;
MPI_Comm MPI_COMM_GRMHD;
MPI_Group MPI_GROUP_LIAISON_FROM_GRMHD;
MPI_Comm MPI_COMM_LIAISON_FROM_GRMHD;
#endif

#if(USEMPILIAISON || USEMPIGRRAY)
MPI_Group MPI_GROUP_GRRAY_LIAISON;
MPI_Comm MPI_COMM_GRRAY_LIAISON;

/////////
//
// liaison and GRRAY specific groups
//
/////////
MPI_Group MPI_GROUP_GRRAY;
MPI_Comm MPI_COMM_GRRAY;

MPI_Group MPI_GROUP_LIAISON_FROM_GRRAY;
MPI_Comm MPI_COMM_LIAISON_FROM_GRRAY;
#endif



// process types in each group
int *processtypelist_world;
int *processtypelist_grmhd_liaison;
int *processtypelist_grray_liaison;
int *processtypelist_grmhd;
int *processtypelist_grray;

int *processtypelist_liaison_from_grmhd;
int *processtypelist_liaison_from_grray;

// number of CPUs in group
int sizeproclist_world;
int sizeproclist_grmhd_liaison;
int sizeproclist_grray_liaison;
int sizeproclist_grmhd;
int sizeproclist_grray;

int sizeproclist_liaison_from_grmhd;
int sizeproclist_liaison_from_grray;


#elif(DOINGLIAISONTYPECODE==1)


#if(USEMPILIAISON)
MPI_Group MPI_GROUP_WORLD;
MPI_Group MPI_GROUP_LIAISON_FROM_GRMHD;
MPI_Comm MPI_COMM_LIAISON_FROM_GRMHD;
#endif
int *processtypelist_world;
int *processtypelist_liaison_from_grmhd;
int sizeproclist_world;
int sizeproclist_liaison_from_grmhd;

// needed for parts of liaison code that uses set-up of grid in grmhd code
int sizeproclist_grmhd;
#if(USEMPILIAISON || USEMPIGRMHD)
MPI_Group MPI_GROUP_GRMHD;
MPI_Comm MPI_COMM_GRMHD;
#endif


#elif(DOINGGRMHDTYPECODE==1)


#if(USEMPIGRMHD)
MPI_Group MPI_GROUP_WORLD;
MPI_Group MPI_GROUP_GRMHD;
MPI_Comm MPI_COMM_GRMHD;
#endif
int *processtypelist_world;
int *processtypelist_grmhd;
int sizeproclist_world;
int sizeproclist_grmhd;

#elif(DOINGGRRAYTYPECODE==1)

#if(USEMPIGRRAY)
MPI_Group MPI_GROUP_WORLD;
MPI_Group MPI_GROUP_GRRAY;
MPI_Comm MPI_COMM_GRRAY;
#endif
int *processtypelist_world;
int *processtypelist_grray;
int sizeproclist_world;
int sizeproclist_grray;

#endif
