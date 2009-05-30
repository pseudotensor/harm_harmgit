/*
 * parabolic interpolation subroutin  
 * ref. Collella && Woodward's PPM paper
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "decs.h"




void para(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

         y[0]=x1;
         y[1]=x2;
         y[2]=x3;
         y[3]=x4;
         y[4]=x5;

         /*CW1.7 */
         for(i=1 ; i<4 ; i++) {
               Dqm = 2. *(y[i]-y[i-1]);
               Dqp = 2. *(y[i+1]-y[i]);
               Dqc = 0.5 *(y[i+1]-y[i-1]);
               aDqm = fabs(Dqm) ;
               aDqp = fabs(Dqp) ;
               aDqc = fabs(Dqc) ;
               s = Dqm*Dqp;
	       Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
	       aDqvanl=fabs(Dqvanl);
	       
	       if(lim == VANL) 
		 { 
		   if (s <=0.) dq[i]=0.;
		   else dq[i] = -aDqvanl*SIGN(Dqc);
		   //else dq[i]=MIN(MIN(aDqc,aDqvanl),MIN(aDqm,aDqp))*SIGN(Dqc);
		 }
	       else if(lim == MC)
		 {
		   if (s <=0.) dq[i]=0.;       //CW1.8
		   else dq[i]=-MIN(aDqc,MIN(aDqm,aDqp))*SIGN(Dqc);
		 }
	     
	       else if(lim == MINM)
		 {
		   if (s<=0.) dq[i] = 0.;
		   else if (aDqm<aDqp) dq[i] = -aDqm*SIGN(Dqc);
		   else dq[i]=-aDqp*SIGN(Dqc);
		 }
	       else if(lim == NLIM) //w/o slope limiter
		 {
		   //if(s<=0.) dq[i] = 0.;
		   dq[i] = Dqc;
		 }
         }

         /* CW1.6 */

         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;
	 
	 
	 l=MAX(MIN(y[2],y[1]),l);
	 l=MIN(MAX(y[2],y[1]),l);
	 r=MAX(MIN(y[2],y[3]),r);
	 r=MIN(MAX(y[2],y[3]),r);
	 
	 
         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));
	 
	 /*
         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

         else if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
         else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;
	 

         *lout=l;   //a_L,j
	 *rout=r;
	 */

	 if (qa <=0. ) {
	       *lout=y[2];
	       *rout=y[2];
         }  
	 else {
	        *lout = l;
		*rout = r;
	 }
		  
	 //2. at top/bottom of a steep gradient 
         if (qd*(qd-qe)<0.0) *lout=3.0*y[2]-2.0*r;
	 else *lout = l;
 
         if (qd*(qd+qe)<0.0) *rout=3.0*y[2]-2.0*l;
	 else *rout = r;
         //*dw=r-l;                      //CW1.5
         //*w6=6.0*(y[2]-0.5*(l+r));
}
