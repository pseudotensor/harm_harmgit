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


#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define SIGN(a) ( ((a) <0.) ? -1. : 1. )


void para(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

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

               if (s <=0.) dq[i]=0.;       //CW1.8
               else dq[i]=MIN(aDqc,MIN(aDqm,aDqp))*SIGN(Dqc);
         }

         /* CW1.6 */

         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;

         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));


         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

         if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
         else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;
	 

         *lout=l;   //a_L,j
	 *rout=r;
         //*dw=r-l;                      //CW1.5
         //*w6=6.0*(y[2]-0.5*(l+r));
}
