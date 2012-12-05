//int
//analytical_solution(ldouble t,int ix,int iy,int iz,ldouble *uu,ldouble *pp,ldouble *vv)
//{
/*
	      ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
	      ldouble xx,yy,zz;
	   
	      xx=get_x(ix,0);
	      yy=get_x(iy,1);
	      zz=get_x(iz,2);
	      ldouble gg[4][5],eup[4][4],elo[4][4];
	      pick_g(ix,iy,iz,gg);
	      calc_LNRFes(gg,eup,elo);


	      struct tempst1
	      {
		double x,t,k;
	      } tst1;



	      double f_int(double y, void * params) 
	      {
		struct tempst1 *ppp = (struct tempst1 *) params;
		double x=ppp->x;
		double k=ppp->k;
		double t=ppp->t;

		ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;
		
		Tgas=T_AMB*(1.+BLOBP*exp(-((xx)*(xx)+(yy)*(yy)+(zz)*(zz))/BLOBW/BLOBW));
		E=calc_LTE_EfromT(Tgas);

		double gy=E;
		return exp(-(x-y)*(x-y)/4./k/t)*gy;
	      }



ldouble chi= calc_kappaes(RHO_AMB,T_AMB,-1.,-1.,-1.);

	      
double k=1./3./(double)chi;


	      tst1.k=k;
	      tst1.x=xx;
	      tst1.t=t;

	      gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);
       
double result, error;
	      
	      gsl_function F;
	      F.function = &f_int;
	      F.params = &tst1;
	      
	      if(t>0.)
	      gsl_integration_qagi (&F,  1.e-5, 1e-10, 1000, w, &result, &error); 
	      
	      gsl_integration_workspace_free (w);
     
vv[3]=(ldouble)(result/sqrt(4.*Pi*k*(double)t));
if(ix==NX/2) printf("%Le %Le %e %Le\n",t,1./3./chi,k,vv[3]);
*/
 
