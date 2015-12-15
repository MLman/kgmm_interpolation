/* Om Shanti. Om Sai Ram. November 12, 2014. 8:08 p.m. Madison, WI.
   Nagesh Adluru. */

#include <mex.h>
#include "math.h"
#include "pdflib.h"
#include "rnglib.h"
/*#include "vtmv.h"*/
/* mex -I"/home/adluru/GMMProject/JohnBurkardt/include" -largeArrayDims ComputeInnderProductGaussian.c /home/adluru/GMMProject/JohnBurkardt/libc/Linux/pdflib.o /home/adluru/GMMProject/JohnBurkardt/libc/Linux/rnglib.o */

double r8mat_vtmv ( int m, int n, double x[], double a[], double y[] );
   
/* Function to compute inner product of two Gaussian PDFs. */
void ComputeInnerProductGaussian(double *innerProduct, double *mu1, double *mu2, double *Sigma1, double *Sigma2, int d)
{
  /* mexPrintf("Hello!"); */
  int x, y, idx;
  double pi = 3.141592653589793;
  double *r2, *xmu, *cInv, *Sigma, xcx, cDet;
  
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
  cInv = r8mat_poinv(d, r2);
  xcx = r8mat_vtmv(d, d, xmu, cInv, xmu);
  *(innerProduct) = 1.0 / sqrt(pow(2.0 * pi, d))
      * 1.0 / sqrt(cDet)
      * exp(-0.5 * xcx);
  free(cInv);
  free(Sigma);
  free(xmu);
}

/* The MATLAB-C gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Declaring variables. */
  double *innerProduct, *mu1, *mu2, *Sigma1, *Sigma2;
  double d; /* Dimensionality of the Gaussian. */
  int i;

  /* Check for proper input and output. */
  if(nrhs != 4)
    mexErrMsgTxt("Four inputs are expected: mu1, mu2, Sigma1, Sigma2.");
  else if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments. Function just returns one scalar.");
  for(i = 0; i < 4; i++)
  {
    if(mxIsComplex(prhs[i]))
    {
      mexPrintf("\nInput %d is complex!", i + 1);
      mexErrMsgTxt("This input cannot be complex");
    }
  }
  
  /* Extracting values from the input parameters. */
  mu1 = mxGetPr(prhs[0]);
  mu2 = mxGetPr(prhs[1]);
  Sigma1 = mxGetPr(prhs[2]);
  Sigma2 = mxGetPr(prhs[3]);
  d = mxGetM(prhs[0]);
  
  /* Doing some additional consistency checks. */
  if(d != mxGetM(prhs[1]))
    mexErrMsgTxt("Both mu1 and mu2 should be same length column vectors.");
  for(i = 2; i <= 3; i++)
  {
    if(mxGetM(prhs[i]) != mxGetN(prhs[i]))
      mexErrMsgTxt("Sigma1 and Sigma2 need to be square matrices with number of rows matching to that of mu1 and mu2.");
    if(mxGetM(prhs[i]) != d)
      mexErrMsgTxt("Number of rows in Sigma1 and Sigma2 must match to that of mu1 and mu2.");
  }
  
  /* Creating  pointers for output parameters. */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  innerProduct = mxGetPr(plhs[0]);
  
  /* Computing the inner product. */
  ComputeInnerProductGaussian(innerProduct, mu1, mu2, Sigma1, Sigma2, d);
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