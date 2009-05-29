



#define UFSET(theCUf,thedt,theUi,theUf,thedUriemann,thedUgeom) (theCUf[0]*theUi + theCUf[1]*theUf + theCUf[2]*thedt*(thedUriemann+thedUgeom))

  // how much of Ui, dU, and Uf to keep for final solution
  // ultimately ucum is actual solution used to find final pf
// Cunew[1]*dt*dU + Cunew[2]*CUf[2]*dt*dU + \Sum_allpriorsubsteps CUf[1]*CUF[2 prior]*dt*dU[prior]
#define UCUMUPDATE(theCunew,thedt,theUi,theUf,thedUriemann,thedUgeom) (theCunew[0]*theUi + theCunew[1]*thedt*(thedUriemann+thedUgeom) + theCunew[2]*theUf)
