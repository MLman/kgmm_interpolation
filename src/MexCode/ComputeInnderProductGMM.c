/* Om Shanti. Om Sai Ram. November 12+1, 2014. 7:30 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12+1, 2014. 5:40 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12+1, 2014. 4:35 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12+1, 2014. ~1:14 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12+1, 2014. 9:29 a.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12, 2014. 10:00 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12, 2014. 8:08 p.m. Madison, WI.
   Nagesh Adluru. */

#include <mex.h>
#include "math.h"
#include "pdflib.h"
#include "rnglib.h"
/*#include "vtmv.h"*/
/* mex -I"/home/adluru/GMMProject/JohnBurkardt/include" -largeArrayDims ComputeInnderProductGaussian.c /home/adluru/GMMProject/JohnBurkardt/libc/Linux/pdflib.o /home/adluru/GMMProject/JohnBurkardt/libc/Linux/rnglib.o */
/* mex -I"/home/adluru/GMMProject/JohnBurkardt/include" -largeArrayDims ComputeInnderProductGMM.c /home/adluru/GMMProject/JohnBurkardt/libc/Linux/pdflib.o /home/adluru/GMMProject/JohnBurkardt/libc/Linux/rnglib.o */

double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] );

/* Function to compute inner product of two Gaussian PDFs. */
double ComputeInnerProductGaussianDiagCov(double *mu1, double *mu2, double *Sigma1, double *Sigma2, int d)
{
  /* mexPrintf("Hello!"); */
  int x, y, idx, idx2;
  double twoPi = 6.283185307179586;
  double *r2, *xmu, *cInv, *Sigma, xcx, cDet;
  double innerProduct;
  
  /* x - mu. x = mu1, mu = mu2 */
  xmu = (double *)malloc(d * sizeof(double));
  for(x = 0; x < d; x++)
  {
    idx = x; 
    xmu[x] = *(mu1 + idx) - *(mu2 + idx);
  }
  
  /* Sigma = Sigma1 + Sigma2 .*/
  Sigma = (double *)malloc(1 * d * sizeof(double));
  for(y = 0; y < 1; y++)
    for(x = 0; x < d; x++)
    {
      /* Indexing [x][y]
      MATLAB stores arrays in column major ordering
      while C stores in row major ordering */
      idx = x * 1 + y; /* Essentially x but added stuff for readability. */
      *(Sigma + idx) = *(Sigma1 + idx) + *(Sigma2 + idx);
    }
  /* Computing the inverse of Sigma via Cholesky factor. */
  /* r2 = r8mat_pofac(d, Sigma); */
  /* cDet = r8mat_podet(d, r2); */
  /* cDet = 1; */
  cDet = 0;
  for(x = 0; x < d; x++)
  {
    idx = x;
    /* cDet *= *(Sigma + idx); */
    cDet += log(sqrt(*(Sigma + idx)));
  }
  /* cInv = r8mat_poinv(d, r2); */
  cInv = (double *)malloc(d * sizeof(double));
  for(y = 0; y < 1; y++)
    for(x = 0; x < d; x++)
    {
      idx = x * 1 + y;
      idx2 = x * 1;
      /* *(cInv + idx) = (x == y) ? 1 / *(Sigma + idx) : 0; */
      /* *(cInv + idx) = (x == y) ? 1 / *(Sigma + idx2) : 0; /* 1:25 p.m. Madison, WI. */ 
      /* *(cInv + idx) = (x == y) ? 1 / (sqrt(*(Sigma + idx2)) * sqrt(*(Sigma + idx2))) : 0;*/
      *(cInv + idx) = sqrt(*(Sigma + idx2));
      /* *(cInv + idx) = (x == y) ? sqrt(*(Sigma + idx2)) : 0; */
      
	/* mexPrintf("\n cInv[%d] = %e, Sigma = %e, logSqrtDet = %e", x, *(cInv + idx), *(Sigma + idx2), cDet);*/
    }
    
  /* xcx = r8mat_vtmv(d, d, xmu, cInv, xmu); */
  xcx = 0.0;
  for ( y = 0; y < 1; y++ )
  {
    for ( x = 0; x < d; x++ )
    {
      /* xcx += pow(xmu[x] * cInv[x+y*d], 2); */
       /*xcx += (xmu[x] * cInv[x+y*d]) * (xmu[x] * cInv[x+y*d]); */
      xcx += (xmu[x] / cInv[x+y*d]) * (xmu[x] / cInv[x+y*d]); 
    }
  }
  
  /* innerProduct = 1.0 / sqrt(pow(2.0 * pi, d))
      * 1.0 / sqrt(cDet)
      * exp(-0.5 * xcx); */
      
  /* innerProduct = exp(-log(sqrt(cDet)) - (0.5 * d * log(twoPi)) - (0.5 * xcx)); */
  innerProduct = exp(-cDet - (0.5 * d * log(twoPi)) - (0.5 * xcx)); /* cDet is logsqrtdet. */
  free(cInv);
  free(Sigma);
  free(xmu);
  return innerProduct;
}

/* Function to compute inner product of two Gaussian PDFs. */
double ComputeInnerProductGaussian(double *mu1, double *mu2, double *Sigma1, double *Sigma2, int d)
{
  /* mexPrintf("Hello!"); */
  int x, y, idx;
  /*double pi = 3.141592653589793;*/
  double twoPi = 6.283185307179586;
  double *r2, *xmu, *cInv, *Sigma, xcx, cDet;
  double innerProduct;
  
  /* x - mu. x = mu1, mu = mu2 */
  xmu = (double *)malloc(d * sizeof(double));
  for(x = 0; x < d; x++)
  {
    idx = x; 
    xmu[x] = *(mu1 + idx) - *(mu2 + idx);
  }
  
  /* Sigma = Sigma1 + Sigma2 .*/
  Sigma = (double *)malloc(d * d * sizeof(double));
  for(y = 0; y < d; y++)
    for(x = 0; x < d; x++)
    {
      /* Indexing [x][y]
      MATLAB stores arrays in column major ordering
      while C stores in row major ordering */
      idx = x + d * y;
      *(Sigma + idx) = *(Sigma1 + idx) + *(Sigma2 + idx);
    }
  /* Computing the inverse of Sigma via Cholesky factor. */
  r2 = r8mat_pofac(d, Sigma);
  cDet = r8mat_podet(d, r2);
 /* cDet = 0;
  for(x = 0; x < d; x++)
  {
    idx = x;
     cDet *= *(Sigma + idx); 
    cDet += log(*(r2 + idx));
  } */
  /* cInv = r8mat_poinv(d, r2); */
  cInv = r8mat_utsol (d, r2, xmu); /* From the function r8vec_multinormal_pdf in pdflib.c */
  
  /* xcx = r8mat_vtmv(d, d, xmu, cInv, xmu); */
  xcx = r8vec_dot_product (d, cInv, cInv); /*cInv now is actually xmu*r2 */
  
  /* mexPrintf("\n cDet = %e, xcx = %e", cDet, xcx); */
  
  /* innerProduct = 1.0 / sqrt(pow(2.0 * pi, d))
      * 1.0 / sqrt(cDet)
      * exp(-0.5 * xcx); */
      
  innerProduct = exp(-log(sqrt(cDet)) - (0.5 * d * log(twoPi)) - (0.5 * xcx));
  /*innerProduct = exp(-cDet - (0.5 * d * log(twoPi)) - (0.5 * xcx)); /* cDet is logsqrtdet. */
  free(cInv);
  free(Sigma);
  free(xmu);
  return innerProduct;
}

/* The MATLAB-C gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Declaring variables. */
  double *innerProduct, *mu1, *mu2, *Sigma1, *Sigma2;
  int i, numOfFields, numOfElements;
  int k1, k2, d1, d2, d10, d20;
  int j, jj;

  int idx1, idx2;
  int si, sj, mi;
  double d1Sqr, d2Sqr, currentInnerProduct;
  double *psi1, *mu1s, *Sigma1s, *psi2, *mu2s, *Sigma2s;

  /* Check for proper input and output. */
  if(nrhs != 2)
    mexErrMsgTxt("Two input structures are expected: gmm1 and gmm2.");
  else if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments. Function just returns one scalar.");
  for(i = 0; i < 1; i++)
  {
    if(!mxIsStruct(prhs[i]))
    {
      mexPrintf("\nInput %d is not a structure!", i + 1);
      mexErrMsgTxt("This input has to be a structure.");
    }
  }
  
  /* Doing some consistency checks. */
 /* numOfFields = mxGetNumberOfFields(prhs[0]);
  if(numOfFields != 4)
    mexErrMsgTxt("Number of fields in each structure must be 4: NComponents, mu, Sigma, PComponents");
  if(numOfFields != mxGetNumberOfFields(prhs[1]))
    mexErrMsgTxt("Both the input structures should have same number of fields!");*/
  
  /* Extracting values from the input parameters. */
  /* First GMM. */
  k1 = mxGetScalar(mxGetField(prhs[0], 0, "NComponents"));
  mu1s = (double *) mxGetPr(mxGetField(prhs[0], 0, "mu")); /* has to be k1 x d1 matrix */
  Sigma1s = (double *)mxGetPr(mxGetField(prhs[0], 0, "Sigma")); /* can be d1 x d1 x k1 matrix or 1 x d1 x k1 */
  psi1 = (double *)mxGetPr(mxGetField(prhs[0], 0, "PComponents")); /* has to be 1 x k1 */
  d1 = (mxGetDimensions(mxGetField(prhs[0], 0, "Sigma")))[1]; /* Just to be safe taking the second dimension */
  d10 = (mxGetDimensions(mxGetField(prhs[0], 0, "Sigma")))[0];
  
  /*mexPrintf("\n d1 = %d, d10 = %d, k1 = %d", d1, d20, k2);*/ /* Did not need to use. */
  /* mexPrintf("\n d1 = %d, d10 = %d, k1 = %d", d1, d10, k1); /* Used it :). */
  
  /* Second GMM. */
  k2 = mxGetScalar(mxGetField(prhs[1], 0, "NComponents"));
  mu2s = (double *)mxGetPr(mxGetField(prhs[1], 0, "mu")); /* has to be k2 x d2 matrix and d1 = d2*/
  Sigma2s = (double *)mxGetPr(mxGetField(prhs[1], 0, "Sigma")); /* can be d2 x d2 x k2 matrix or 1 x d2 x k2 */
  psi2 = (double *)mxGetPr(mxGetField(prhs[1], 0, "PComponents")); /* has to be 1 x k1 */
  d2 = (mxGetDimensions(mxGetField(prhs[1], 0, "Sigma")))[1]; /* Just to be safe taking the second dimension */
  /*d20 = (mxGetDimensions(mxGetField(prhs[0], 0, "Sigma")))[0];*/
  d20 = (mxGetDimensions(mxGetField(prhs[1], 0, "Sigma")))[0];
  
  /*mexPrintf("\n d2 = %d, d20 = %d, k2 = %d", d2, d20, k2);*/ /* Did not need to use. */
   /* mexPrintf("\n d2 = %d, d20 = %d, k2 = %d", d2, d20, k2); /* Used it :) */
  
  /* Doing some additional consistency checks. (Add more later!)*/
  if(d1 != d2 || d10 != d20)
    mexErrMsgTxt("The two GMMs must have same covariance dimensions!");
  
  /* Creating  pointers for output parameters. */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  innerProduct = (double *)mxGetPr(plhs[0]);
  
  /* Allocating memory to mu1, mu2, Sigma1, Sigma2. */
  d1Sqr = d1 * d10;
  d2Sqr = d2 * d20;
  mu1 = (double *)malloc(d1 * sizeof(double));
  mu2 = (double *)malloc(d1 * sizeof(double));
  Sigma1 = (double *)malloc(d1Sqr * sizeof(double));
  Sigma2 = (double *)malloc(d2Sqr * sizeof(double));
  
  /* Computing the inner product. */
  *(innerProduct) = 0;
  for(j = 0; j < k1; j++)
  {
    for(mi = 0; mi < d1; mi++)
    {
      idx1 = mi; /* Redudant but will help increase readability. */
      idx2 = mi * k1 + j;
      *(mu1 + idx1) = *(mu1s + idx2);
    }
    for(si = 0; si < d10; si++)
      for(sj = 0; sj < d1; sj++)
      {
	idx1 = sj * d10 + si; /* ~1:14 p.m. */
	idx2 = j * d1Sqr + sj * d10 + si;
	*(Sigma1 + idx1) = *(Sigma1s + idx2);
      }
    for(jj = 0; jj < k2; jj++)
    {
      for(mi = 0; mi < d2; mi++)
      {
	idx1 = mi; /* Redudant but will help increase readability. */
	idx2 = mi * k2 + jj;
	*(mu2 + idx1) = *(mu2s + idx2);
      }
    for(si = 0; si < d20; si++) /* Remember d2 = d1. */
	for(sj = 0; sj < d2; sj++)
	{
	  idx1 = sj * d20 + si; /* ~1:14 p.m. */
	  idx2 = jj * d2Sqr + sj * d20 + si;
	  *(Sigma2 + idx1) = *(Sigma2s + idx2);
	}
      currentInnerProduct = (d10 == 1) ? ComputeInnerProductGaussianDiagCov(mu1, mu2, Sigma1, Sigma2, d1) : ComputeInnerProductGaussian(mu1, mu2, Sigma1, Sigma2, d1);
      /* mexPrintf("\nInner product between component %d, %d = %f", j, jj, currentInnerProduct); */
      *(innerProduct) = *(innerProduct) + ((*(psi1 + j)) * (*(psi2 + jj))) * currentInnerProduct;
    }
  }
 
  /*j = 0;
  mexPrintf("\nSigma1 for component %d\n", j + 1);
  for(si = 0; si < d10; si++)
  {
      for(sj = 0; sj < d1; sj++)
      {
	idx1 = sj * d10 + si;
	idx2 = j * d1Sqr + sj * d10 + si;
	 mexPrintf("%e ", *(Sigma1s + idx2));
      }
	  mexPrintf("\n");
  }*/
  
  
  free(mu1);
  free(mu2);
  free(Sigma1);
  free(Sigma2);
}

/******************************************************************************/

double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_VTMV multiplies computes the scalar x' * A * y.

  Discussion:

    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of
    the matrix.

    Input, double X[N], the first vector factor.

    Input, double A[M*N], the M by N matrix.

    Input, double Y[M], the second vector factor.

    Output, double R8MAT_VTMV, the value of X' * A * Y.
*/
{
  int i;
  int j;
  double vtmv;

  vtmv = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      vtmv = vtmv + x[i] * a[i+j*m] * y[j];
    }
  }
  return vtmv;
}