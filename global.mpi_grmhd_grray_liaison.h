
/*! \file global.mpi_grmhd_grray_liaison.h
    \brief global definitions/macros that does LIAISON to connect Broderick ray tracing code and HARM



//////////////
//
// To use mpi_grmhd_grray_liaison.c code:
//
//////////////
//
// 0) Define DOINGLIAISON 0 or 1 in your code that gets included in decs.h
// 1) Need to define a DOINGGRMHDTYPECODE or DOINGGRRAYTYPECODE to 1 in your GRMHD+LIAISON or GRRAY+LIAISON codes
// 2) Use: mpidefs.mpi_grmhd_grray_liaison.h also contains definitions
// 3) Use: global.mpi_grmhd_grray_liaison.h contains global settings for macros and function declarations
// 4) Use one of 3 *_globalset() functions below or use local versions as inside grmhd_grray_liaison.c in order to initialize use of code when doing or not doing liaisonmode
// 5) 

*/


/// assumes local code sets USEMPIGRMHD and USEMPIGRRAY
#define USEMPINONLIAISON (USEMPIGRMHD&&DOINGGRMHDTYPECODE || USEMPIGRRAY&&DOINGGRRAYTYPECODE) // no choice
#define USEMPILIAISON (DOINGLIAISON) // can turn off USEMPILIAISON if DOINGLIAISON==1 but want no MPI code for liaison.  However, if DOINGLIAISON==0 then USEMPILIAISON must be 0

#define USEMPIANY (USEMPINONLIAISON || USEMPILIAISON) // no choice



//////// each cpu identifies its binary type in an array
#define BINARYNOTYPE 0
#define BINARYGRMHDTYPE 1
#define BINARYGRRAYTYPE 2
#define BINARYLIAISONTYPE 3

#if(USEMPIANY==0)
#define MPI_Group int
#define MPI_Comm int
#endif

/////////////
///
/// functions to call if using global variables
/// If want to use local functions, then see inside mpi_grmhd_grray_liaison.c for function declarations
///
/////////////
extern void liaison_init_mpi_liaisonmode_globalset(void);
extern void grmhd_init_mpi_liaisonmode_globalset(void);
extern void grray_init_mpi_liaisonmode_globalset(void);

