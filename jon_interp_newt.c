#include "decs.h"
#include "NRUTIL.H"

#define NRANSI

// modifiable
#define MAXITS 200
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0


/*
// for debugging
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);fprintf(stderr,"its=%d test=%g\n",its,test);return;}
*/
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}

void newt(int useanalyticjac
	  ,FTYPE parms[]
	  ,FTYPE x[], int n, int *check
	  ,void (*vecfunc)(int n, FTYPE *parms, FTYPE v[], FTYPE f[])
		   ,int (*usrfun)(int n, FTYPE *parms, FTYPE *Xguess, FTYPE *spc_diff, FTYPE **alpha)
	  )
{
  //  void fdjac(int n, FTYPE parms[], FTYPE x[], FTYPE fvec[], FTYPE **df,
  //	     void (*vecfunc)(int n, FTYPE parms[], FTYPE v[], FTYPE f[]));
//	void lubksb(FTYPE **a, int n, int *indx, FTYPE b[]);
	//int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
	int i,its,j,*indx;
	FTYPE d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
	// DEBUG VARS
	int qq,pp;



	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	fvec=vector(1,n);
	nn=n;
	nrfuncv=vecfunc;
	f=nrfmin(parms,x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test < 0.01*TOLF) {
		*check=0;
		//fprintf(stderr,"end tolf0.01\n");
		FREERETURN
	}
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(FTYPE)n);
	for (its=1;its<=MAXITS;its++) {
	  if(useanalyticjac){
	    (*usrfun)(n,parms,x,fvec,fjac);
	  }
	  else{
	    fdjac(n,parms,x,fvec,fjac,vecfunc);
	  }
		///////////////////////////DEBUG
	  /*
		for(qq=1;qq<n+1;qq++) for(pp=1;pp<n+1;pp++){
		  fprintf(stderr,"fjac[%d][%d]=%g\n",qq,pp,fjac[qq][pp]);
		}
		fflush(stderr);
		//exit(0);
		*/		
		///////////////////////////DEBUG
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		ludcmp(fjac,n,indx,&d);
		//		if(its==5) exit(0);
		lubksb(fjac,n,indx,p);
		lnsrch(n,parms,xold,fold,g,p,x,&f,stpmax,check,nrfmin);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			//fprintf(stderr,"end tolf\n");
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 1 : 0);
			if(SIMPLEDEBUGINTERP) fprintf(stderr,"end check\n");
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX){
		  //		  fprintf(stderr,"end tolx\n");
		  FREERETURN
		    }
	}
	nrerror("MAXITS exceeded in newt");
}
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. */
