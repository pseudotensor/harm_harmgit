//KORAL - bc.c
//returns absorption opacity coefficient
//**********************
//included in:
//ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
//from rad.c
//**********************

if(NTUBE==1) return 0.4*rho;
if(NTUBE==2) return 0.2*rho;
if(NTUBE==3) return 0.3*rho;
if(NTUBE==31) return 25.*rho;
if(NTUBE==4) return 0.08*rho;
if(NTUBE==41) return 0.7*rho;
if(NTUBE==5) return 1000.*rho;
