/* Om Shanti. Om Sai Ram. November 12, 2014. 11:04 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12, 2014. 10:00 p.m. Madison, WI.
   Om Shanti. Om Sai Ram. November 12, 2014. 8:08 p.m. Madison, WI.
   Nagesh Adluru. */

#include <mex.h>
/* mex -I"/home/adluru/GMMProject/JohnBurkardt/include" -largeArrayDims ComputeInnderProductGaussian.c /home/adluru/GMMProject/JohnBurkardt/libc/Linux/pdflib.o /home/adluru/GMMProject/JohnBurkardt/libc/Linux/rnglib.o */
/* mex Check3DArrayAccess.c */
 
/* The MATLAB-C gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Declaring variables. */
  double *array, *idxArray, *element;
  int i, j, k, d1, d2, d3, idx;
  
  /* Extracting the values from input parameters. */
  array = mxGetPr(mxGetField(prhs[0], 0, "array"));
  idxArray = mxGetPr(mxGetField(prhs[0], 0, "idxArray"));
  
  /* Creating  pointers for output parameters. */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  element = mxGetPr(plhs[0]);
  
  d1 = (mxGetDimensions(mxGetField(prhs[0], 0, "array")))[0];
  d2 = (mxGetDimensions(mxGetField(prhs[0], 0, "array")))[1];
  d3 = (mxGetDimensions(mxGetField(prhs[0], 0, "array")))[2];
  
  i = *(idxArray + 0);
  j = *(idxArray + 1);
  k = *(idxArray + 2);
  
  mexPrintf("\n d1 = %d, d2 = %d, d3 = %d", d1, d2, d3);
  mexPrintf("\n i = %d, j = %d, k = %d", i, j, k);
  
  idx = (k - 1) * (d2 * d1) + (j - 1) * (d1) + (i - 1);
  *(element) = *(array + idx);
}