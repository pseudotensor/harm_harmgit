//KORAL - problem.rad.h
//choice of the problem plus some definitions



//available problems:

//1 RADBEAM2D - beam of light 
//2 RADINFALL - radial inflow
//3 DONUT - 2d Polish donut
//4 GEODESICINFALL - geodesic infall with blobs or not
//5 EDDINFALL - infall with flux from inside
//6 RADTUBE- radiative shock tubes as in Farris et al 09
//7 BONDI - like in Fragile's paper
//8 HDTUBE - relativistic shock tube
//9 HDTUBE2D - in 2d
//10 RADPULSE - radiative blob spreading around
//11 RADSHADOW - radiative shadow
//12 RADATM - atmosphere enlighted
//13 DONUTOSC - 2d Polish donut oscillating
//14 RADWAVEBC - 1d linear rad wave imposed on boundary
//15 RADWAVE - 1d linear rad wave with periodic BC
//16 RADPULSE3D - radiative blob spreading around
//17 RADDBLSHADOW - radiative shadow with two beams inclined

#define PROBLEM 15

#if(PROBLEM==1)

#define PR_DEFINE "RADPROBLEMS/RADBEAM2D/define.h"
#define PR_BC "RADPROBLEMS/RADBEAM2D/bc.c"
#define PR_INIT "RADPROBLEMS/RADBEAM2D/init.c"
#define PR_KAPPA "RADPROBLEMS/RADBEAM2D/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADBEAM2D/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADBEAM2D/anasol.c"

#endif

#if(PROBLEM==2)

#define PR_DEFINE "RADPROBLEMS/RADINFALL/define.h"
#define PR_BC "RADPROBLEMS/RADINFALL/bc.c"
#define PR_INIT "RADPROBLEMS/RADINFALL/init.c"
#define PR_KAPPA "RADPROBLEMS/RADINFALL/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADINFALL/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADINFALL/anasol.c"

#endif

#if(PROBLEM==3)

#define PR_DEFINE "RADPROBLEMS/DONUT/define.h"
#define PR_BC "RADPROBLEMS/DONUT/bc.c"
#define PR_INIT "RADPROBLEMS/DONUT/init.c"
#define PR_KAPPA "RADPROBLEMS/DONUT/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/DONUT/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/DONUT/anasol.c"

#endif

#if(PROBLEM==4)

#define PR_DEFINE "RADPROBLEMS/GEODESICINFALL/define.h"
#define PR_BC "RADPROBLEMS/GEODESICINFALL/bc.c"
#define PR_INIT "RADPROBLEMS/GEODESICINFALL/init.c"
#define PR_KAPPA "RADPROBLEMS/GEODESICINFALL/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/GEODESICINFALL/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/GEODESICINFALL/anasol.c"

#endif

#if(PROBLEM==5)

#define PR_DEFINE "RADPROBLEMS/EDDINFALL/define.h"
#define PR_BC "RADPROBLEMS/EDDINFALL/bc.c"
#define PR_INIT "RADPROBLEMS/EDDINFALL/init.c"
#define PR_KAPPA "RADPROBLEMS/EDDINFALL/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/EDDINFALL/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/EDDINFALL/anasol.c"

#endif

#if(PROBLEM==6)

#define PR_DEFINE "RADPROBLEMS/RADTUBE/define.h"
#define PR_BC "RADPROBLEMS/RADTUBE/bc.c"
#define PR_INIT "RADPROBLEMS/RADTUBE/init.c"
#define PR_KAPPA "RADPROBLEMS/RADTUBE/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADTUBE/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADTUBE/anasol.c"

#endif

#if(PROBLEM==7)

#define PR_DEFINE "RADPROBLEMS/BONDI/define.h"
#define PR_BC "RADPROBLEMS/BONDI/bc.c"
#define PR_INIT "RADPROBLEMS/BONDI/init.c"
#define PR_KAPPA "RADPROBLEMS/BONDI/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/BONDI/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/BONDI/anasol.c"

#endif

#if(PROBLEM==8)

#define PR_DEFINE "RADPROBLEMS/HDTUBE/define.h"
#define PR_BC "RADPROBLEMS/HDTUBE/bc.c"
#define PR_INIT "RADPROBLEMS/HDTUBE/init.c"
#define PR_KAPPA "RADPROBLEMS/HDTUBE/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/HDTUBE/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/HDTUBE/anasol.c"

#endif

#if(PROBLEM==9)

#define PR_DEFINE "RADPROBLEMS/HDTUBE2D/define.h"
#define PR_BC "RADPROBLEMS/HDTUBE2D/bc.c"
#define PR_INIT "RADPROBLEMS/HDTUBE2D/init.c"
#define PR_KAPPA "RADPROBLEMS/HDTUBE2D/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/HDTUBE2D/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/HDTUBE2D/anasol.c"

#endif

#if(PROBLEM==10)

#define PR_DEFINE "RADPROBLEMS/RADPULSE/define.h"
#define PR_BC "RADPROBLEMS/RADPULSE/bc.c"
#define PR_INIT "RADPROBLEMS/RADPULSE/init.c"
#define PR_KAPPA "RADPROBLEMS/RADPULSE/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADPULSE/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADPULSE/anasol.c"

#endif

#if(PROBLEM==11)

#define PR_DEFINE "RADPROBLEMS/RADSHADOW/define.h"
#define PR_BC "RADPROBLEMS/RADSHADOW/bc.c"
#define PR_INIT "RADPROBLEMS/RADSHADOW/init.c"
#define PR_KAPPA "RADPROBLEMS/RADSHADOW/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADSHADOW/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADSHADOW/anasol.c"

#endif

#if(PROBLEM==12)

#define PR_DEFINE "RADPROBLEMS/RADATM/define.h"
#define PR_BC "RADPROBLEMS/RADATM/bc.c"
#define PR_INIT "RADPROBLEMS/RADATM/init.c"
#define PR_KAPPA "RADPROBLEMS/RADATM/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADATM/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADATM/anasol.c"

#endif

#if(PROBLEM==13)

#define PR_DEFINE "RADPROBLEMS/DONUTOSC/define.h"
#define PR_BC "RADPROBLEMS/DONUTOSC/bc.c"
#define PR_INIT "RADPROBLEMS/DONUTOSC/init.c"
#define PR_KAPPA "RADPROBLEMS/DONUTOSC/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/DONUTOSC/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/DONUTOSC/anasol.c"

#endif

#if(PROBLEM==14)

#define PR_DEFINE "RADPROBLEMS/RADWAVEBC/define.h"
#define PR_BC "RADPROBLEMS/RADWAVEBC/bc.c"
#define PR_INIT "RADPROBLEMS/RADWAVEBC/init.c"
#define PR_KAPPA "RADPROBLEMS/RADWAVEBC/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADWAVEBC/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADWAVEBC/anasol.c"

#endif

#if(PROBLEM==15)

#define PR_DEFINE "RADPROBLEMS/RADWAVE/define.h"
#define PR_BC "RADPROBLEMS/RADWAVE/bc.c"
#define PR_INIT "RADPROBLEMS/RADWAVE/init.c"
#define PR_KAPPA "RADPROBLEMS/RADWAVE/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADWAVE/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADWAVE/anasol.c"

#endif

#if(PROBLEM==16)

#define PR_DEFINE "RADPROBLEMS/RADPULSE3D/define.h"
#define PR_BC "RADPROBLEMS/RADPULSE3D/bc.c"
#define PR_INIT "RADPROBLEMS/RADPULSE3D/init.c"
#define PR_KAPPA "RADPROBLEMS/RADPULSE3D/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADPULSE3D/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADPULSE3D/anasol.c"

#endif

#if(PROBLEM==17)

#define PR_DEFINE "RADPROBLEMS/RADDBLSHADOW/define.h"
#define PR_BC "RADPROBLEMS/RADDBLSHADOW/bc.c"
#define PR_INIT "RADPROBLEMS/RADDBLSHADOW/init.c"
#define PR_KAPPA "RADPROBLEMS/RADDBLSHADOW/kappa.c"
#define PR_KAPPAES "RADPROBLEMS/RADDBLSHADOW/kappaes.c"
#define PR_ANASOL "RADPROBLEMS/RADDBLSHADOW/anasol.c"

#endif

/*********************/
//including problem specific definitions from RADPROBLEMS/XXX/define.h
/*********************/
#include PR_DEFINE

