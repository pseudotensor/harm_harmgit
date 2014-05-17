
/*! \file global.stepch.h
    \brief macros and definitions related to step_ch (RK stepping)

// how much of Ui, dU, and Uf to keep for final solution
// ultimately ucum is actual solution used to find final pf
// Cunew[1]*dt*dU + Cunew[2]*CUf[2]*dt*dU + \Sum_allpriorsubsteps CUf[1]*CUF[2 prior]*dt*dU[prior]

*/


#if(MAXTIMEORDER==4)
#define UFSET(theCUf,thedt,theUi,theUf,thedUriemann,thedUgeom,thedUnongeom) ((theCUf[0])*(theUi) + (theCUf[1])*(theUf) + (theCUf[2])*(thedt)*((thedUriemann)+(thedUgeom))  + (theCUf[NUMPREDTCUFS+0])*(thedt)*(thedUnongeom[0]) + (theCUf[NUMPREDTCUFS+1])*(thedt)*(thedUnongeom[1]) + (theCUf[NUMPREDTCUFS+2])*(thedt)*(thedUnongeom[2]) + (theCUf[NUMPREDTCUFS+3])*(thedt)*(thedUnongeom[3]) )
#define UCUMUPDATE(theCunew,thedt,theUi,theUf,thedUriemann,thedUgeom,thedUnongeom) ((theCunew[0])*(theUi) + (theCunew[1])*(thedt)*((thedUriemann)+(thedUgeom)) + (theCunew[2])*(theUf) + (theCunew[NUMPREDTCUFS+0])*(thedt)*(thedUnongeom[0]) + (theCunew[NUMPREDTCUFS+1])*(thedt)*(thedUnongeom[1]) + (theCunew[NUMPREDTCUFS+2])*(thedt)*(thedUnongeom[2]) + (theCunew[NUMPREDTCUFS+3])*(thedt)*(thedUnongeom[3]) )
#elif(MAXTIMEORDER==5)
#define UFSET(theCUf,thedt,theUi,theUf,thedUriemann,thedUgeom,thedUnongeom) ((theCUf[0])*(theUi) + (theCUf[1])*(theUf) + (theCUf[2])*(thedt)*((thedUriemann)+(thedUgeom))  + (theCUf[NUMPREDTCUFS+0])*(thedt)*(thedUnongeom[0]) + (theCUf[NUMPREDTCUFS+1])*(thedt)*(thedUnongeom[1]) + (theCUf[NUMPREDTCUFS+2])*(thedt)*(thedUnongeom[2]) + (theCUf[NUMPREDTCUFS+3])*(thedt)*(thedUnongeom[3]) + (theCUf[NUMPREDTCUFS+4])*(thedt)*(thedUnongeom[4]) )
#define UCUMUPDATE(theCunew,thedt,theUi,theUf,thedUriemann,thedUgeom,thedUnongeom) ((theCunew[0])*(theUi) + (theCunew[1])*(thedt)*((thedUriemann)+(thedUgeom)) + (theCunew[2])*(theUf) + (theCunew[NUMPREDTCUFS+0])*(thedt)*(thedUnongeom[0]) + (theCunew[NUMPREDTCUFS+1])*(thedt)*(thedUnongeom[1]) + (theCunew[NUMPREDTCUFS+2])*(thedt)*(thedUnongeom[2]) + (theCunew[NUMPREDTCUFS+3])*(thedt)*(thedUnongeom[3]) + (theCunew[NUMPREDTCUFS+4])*(thedt)*(thedUnongeom[4]) )
#else
#error "Need to setup UFSET for other MAXTIMEORDER"
#endif


/// inverse of above UFSET() (assuming only 1 dU)
#define dUfromUFSET(theCUf,thedt,theUi,theUf,theUfnew) ( ((theUfnew) - ((theCUf[0])*(theUi) + (theCUf[1])*(theUf)))/((theCUf[2])*(thedt)) )


