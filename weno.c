// ##############################################################
//
// weno.c
//
//
// ###############################################################
#include <stdio.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"


//double weno5(double a, double b, double c, double d, double e, double epsw);

void weno(double **f,  //fluxes at each edge
         double **ql, // left states
         double **qr, // right states
         int    is,   // start index of chain
         int    ie,   // end index of chain
         double epsw,  // epsilon in limiter
         int    imax, // chainSize
         int    nq)   // number of variables
{
 int i,n,im,k;
 double f0[imax-1],f1[imax-1],f2[imax-1];

  epsw = 1e-12;
  im = imax-2;

  // 1st order upwind
 
//  for (n=0;n<nq;n++)
//  {
//   for (i=is;i<im;i++) //loop through chain 
//   {
//     ql[i][n] = f[i][n]; //left and right states are the same as the face
//     qr[i][n] = f[i][n];
//   }
//  }
//  return;
 
  // 5th order weno scheme
  for (n=0;n<nq;n++) // number of flow variable
  {
   for (i=0;i<=ie;i++) // loop through the chain
   {   
     f0[i] = f[i][n];
   }
     
   for (i=is;i<im;i++)
   {
     ql[i][n]=weno5(f0[i-2],f0[i-1],f0[i],f0[i+1],f0[i+2],epsw);
     
     //printf("ql:%f\n",ql[i][n]); 
     qr[i][n]=weno5(f0[i+2],f0[i+1],f0[i],f0[i-1],f0[i-2],epsw);
     //printf("qr:%f\n",qr[i][n]); 
   }



 } //loop for number of flow variable

} //end function

//##########################################################
//
// weno5
//
//##########################################################
#include <stdio.h>
#include <stdlib.h>
//#include "ham2dtypes.h"
//#include "ham2dFunctionDefs.h"


double weno5(double a,  //fluxes at each edge
          double  b, // left states
          double  c, // right states
          double  d,   // start index of chain
          double  e,   // end index of chain
          double  epsw)
{
 double b1,b2,djm1,ejm1,dj,ej,djp1,ejp1;
 double dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0;
 double w0,w1,w2,wopt0,wopt1,wopt2;
 double sol,wtot,w;
 int    pexp,iweight;
    
      pexp = 2;
      
      b1 = 13./12.;
      b2 = 1./6.;
      djm1 = a-2.*b+c;
      ejm1 = a-4.*b+3.*c;
      dj   = b-2.*c+d;
      ej   = b-d;
      djp1 = c-2.*d+e;
      ejp1 = 3.*c-4.*d+e;
      dis0 = (b1*djm1*djm1+0.25*ejm1*ejm1)+epsw;
      dis1 = (b1*dj*dj+0.25*ej*ej)+epsw;
      dis2 = (b1*djp1*djp1+0.25*ejp1*ejp1)+epsw;
      
      //dis0 = pow(dis0,p);
      //dis1 = pow(dis1,p);
      //dis2 = pow(dis2,p);

      
      d01 = dis0/dis1;
      d02 = dis0/dis2;
      
      wopt0 = 0.1;
      wopt1 = 0.6;
      wopt2 = 0.3;

      a1ba0 = wopt1/wopt0*pow(abs(d01),pexp);
      a2ba0 = wopt2/wopt0*pow(abs(d02),pexp);


//      a1ba0 = 6.*d01;
//      a2ba0 = 3.*d02;


      w0 = 1.0/(1.0+a1ba0+a2ba0);
      w1 = a1ba0*w0;
      w2 = a2ba0*w0;
      wtot = w0+w1+w2;
      w0 = w0/wtot; w1 = w1/wtot; w2 = w2/wtot;

      iweight=0;
      if(iweight==0)
      {
        w=w0;
        w=w*(wopt0+wopt0*wopt0-w*3.0*wopt0+w*w)/(wopt0*wopt0+w*(1.0-2*wopt0));
        w0=w;
        w=w1;
        w=w*(wopt1+wopt1*wopt1-w*3.0*wopt1+w*w)/(wopt1*wopt1+w*(1.0-2*wopt1));
        w1=w;
        w=w2;
        w=w*(wopt2+wopt2*wopt2-w*3.0*wopt2+w*w)/(wopt2*wopt2+w*(1.0-2*wopt2));
        w2=w;

        wtot = w0+w1+w2;
        w0=w0/wtot; w1=w1/wtot; w2=w2/wtot;
      }

      q30 = 2.*a-7.*b+11.*c;
      q31 = -b+5.*c+2.*d;
      q32 = 2.*c+5.*d-e;

      sol = b2*(w0*q30+w1*q31+w2*q32);  //output
      
      //printf("sol:%f\n",sol);
      return sol;
        
} // end function
