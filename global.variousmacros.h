

/*! \file global.variousmacros.h
    \brief macros and definitions related to various things
    // Various macros

*/

/*! \file gslincludes.h
    \brief macros and definitions for GSL library

*/


#define MYDMIN(a,b) (mydminarg1=(a),mydminarg2=(b),(mydminarg1) < (mydminarg2) ? \
                     (mydminarg1) : (mydminarg2))

#define delta(i,j) ((i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define mink(I,J) (I != J ? (0.) : (I == 0 ? (-1.) : (1.)))

#define pfixupeach(pr,i,j,k,which,min) {if(pr[which]<min){ fladd[which]+=dV*MYGDET(i,j,k,CENT)*(min-pr[which]); pr[which]=min;}}

#define pfixup(pr,i,j,k) {pfixupeach(pr,i,j,k,RHO,RHOMIN); pfixupeach(pr,i,j,k,UU,UUMIN); }

// #define FAILSTATEMENT(file,function,number) {fprintf(fail_file,"%s
// %d-%s(): failure\n",file,number,function); fflush(fail_file);
// dualfprintf(fail_file,"MAC(rho,i,j,k): %21.15g MAC(uu,i,j,k): %21.15g MAC(rho2,i,j,k):
// %21.15g MAC(uu2,i,j,k): %21.15g i: %d j: %d pl:
// %d\n",MACP0A1(p,i,j,k,RHO),MACP0A1(p,i,j,k,UU),MACP0A1(ph,i,j,k,RHO),MACP0A1(ph,i,j,k,UU),i,j,pl);
// return(1);}

#define FAILSTATEMENTVOID(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function);} }

#if(USEOPENMP==0)
#define FAILSTATEMENT(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function);} return(1);}
#else
// can't have return in OpenMP parallel section, so use this:
#define FAILSTATEMENT(file,function,number) FAILSTATEMENTVOID(file,function,number)
#endif



#if(JONCHECKS2 && PRODUCTION==0 && (USEOPENMP==0))
#define MYFUN(fun,one,two,three) if(fun>=1){ FAILSTATEMENT(one,two,three);}
#else
// if PRODUCTION>0 then avoid if statement
#define MYFUN(fun,one,two,three) {fun;}
#endif


#if(PRODUCTION>1)
/// blank-out the function call since error_check() slows things down in MPI
#define error_check(wherefrom) (0)
#endif


//https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
#if(PRODUCTION>0)
#define prod0dualfprintf(cond,...) if(cond){ dualfprintf(__VA_ARGS__); }
#else
#define prod0dualfprintf(cond,...)
#endif



