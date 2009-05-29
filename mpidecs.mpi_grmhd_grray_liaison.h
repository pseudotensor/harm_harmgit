
#if(DOINGLIAISON)

#if(USEMPILIAISON || USEMPIGRMHD)
//////////
//
// group and communicators needed by all codes
//
//////////
extern MPI_Group MPI_GROUP_WORLD;
//MPI_Comm MPI_COMM_WORLD;
extern MPI_Group MPI_GROUP_GRMHD_LIAISON;
extern MPI_Comm MPI_COMM_GRMHD_LIAISON;

//////////
//
// liaison and GRMHD specific groups
//
/////////
extern MPI_Group MPI_GROUP_GRMHD;
extern MPI_Comm MPI_COMM_GRMHD;
extern MPI_Group MPI_GROUP_LIAISON_FROM_GRMHD;
extern MPI_Comm MPI_COMM_LIAISON_FROM_GRMHD;
#endif

#if(USEMPILIAISON || USEMPIGRRAY)
extern MPI_Group MPI_GROUP_GRRAY_LIAISON;
extern MPI_Comm MPI_COMM_GRRAY_LIAISON;

/////////
//
// liaison and GRRAY specific groups
//
/////////
extern MPI_Group MPI_GROUP_GRRAY;
extern MPI_Comm MPI_COMM_GRRAY;

extern MPI_Group MPI_GROUP_LIAISON_FROM_GRRAY;
extern MPI_Comm MPI_COMM_LIAISON_FROM_GRRAY;
#endif



// process types in each group
extern int *processtypelist_world;
extern int *processtypelist_grmhd_liaison;
extern int *processtypelist_grray_liaison;
extern int *processtypelist_grmhd;
extern int *processtypelist_grray;

extern int *processtypelist_liaison_from_grmhd;
extern int *processtypelist_liaison_from_grray;

// number of CPUs in group
extern int sizeproclist_world;
extern int sizeproclist_grmhd_liaison;
extern int sizeproclist_grray_liaison;
extern int sizeproclist_grmhd;
extern int sizeproclist_grray;

extern int sizeproclist_liaison_from_grmhd;
extern int sizeproclist_liaison_from_grray;


#elif(DOINGLIAISONTYPECODE==1)


#if(USEMPILIAISON)
extern MPI_Group MPI_GROUP_WORLD;
extern MPI_Group MPI_GROUP_LIAISON_FROM_GRMHD;
extern MPI_Comm MPI_COMM_LIAISON_FROM_GRMHD;
#endif
extern int *processtypelist_world;
extern int *processtypelist_liaison_from_grmhd;
extern int sizeproclist_world;
extern int sizeproclist_liaison_from_grmhd;

// needed for parts of liaison code that uses set-up of grid in grmhd code
extern int sizeproclist_grmhd;
#if(USEMPILIAISON || USEMPIGRMHD)
extern MPI_Group MPI_GROUP_GRMHD;
extern MPI_Comm MPI_COMM_GRMHD;
#endif


#elif(DOINGGRMHDTYPECODE==1)


#if(USEMPIGRMHD)
extern MPI_Group MPI_GROUP_WORLD;
extern MPI_Group MPI_GROUP_GRMHD;
extern MPI_Comm MPI_COMM_GRMHD;
#endif
extern int *processtypelist_world;
extern int *processtypelist_grmhd;
extern int sizeproclist_world;
extern int sizeproclist_grmhd;

#elif(DOINGGRRAYTYPECODE==1)

#if(USEMPIGRRAY)
extern MPI_Group MPI_GROUP_WORLD;
extern MPI_Group MPI_GROUP_GRRAY;
extern MPI_Comm MPI_COMM_GRRAY;
#endif
extern int *processtypelist_world;
extern int *processtypelist_grray;
extern int sizeproclist_world;
extern int sizeproclist_grray;

#endif
