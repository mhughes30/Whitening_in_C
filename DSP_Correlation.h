

#ifndef DSP_CORRELATION_H_
#define DSP_CORRELATION_H_

#include "DSP_Matrix.h"

// Calculates the correlation for a series of lags between two input arrays.
// Inputs: (1) x = input array 1
//		   (2) y = input array 2
//         (3) Rxy = preallocated correlation array between x and y
//         (4) lags = preallocated array containing lag values
//		   (5) length = the length of x or y
// Note:	x and y must be the same length
//          Also, rxy must be preallocated to zero. Otherwise, random values can occur!
//          This function can be used for autocorrelation as well.
void DSP_Correlation(float *x, float *y, float *Rxy, int *lags, unsigned int length);


// Calculates the Auto-Correlation for a series of lags. 
// Inputs: (1) x = input array 1
//         (2) Rxx = preallocated auto-correlation array
//		   (3) length = the length of x
// Note:	x and y must be the same length
//          Also, Rxx must be preallocated to zero. Otherwise, random values can occur!
//          This function is slightly more efficient than DSP_Correlation, since it requires less 
//          inputs and isn't concerned with the lag vector.
void DSP_AutoCorrelation(float *x, float *Rxx, unsigned int length);


// Calculates the sample offset between two data arrays
// Inputs: (1) rxy = the correlation array
//		   (2) lags = the lag array
//		   (3) length = the lenght of rxy of lags
int  DSP_GetSampleOffset(float *rxy, int *lags, unsigned int length);


// Produces the Auto-Correlation matrix of input x
// Inputs: (1) x = preallocated input array
//		   (2) lenX = length of x
//		   (3) dimA = the desired dimension of the square auto-correlation matrix
// Output: (1) the auto-correlation matrix of x
MATRIX *DSP_ProduceAutoCorrMat(float *x, unsigned int lenX, unsigned int dimA);


#endif

