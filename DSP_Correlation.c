
#include <stdlib.h>
#include "DSP_Correlation.h"


//----------------------------- DSP_Correlation -------------------------------------//
void DSP_Correlation(float *x, float *y, float *Rxy, int *lag, unsigned int length)
{
	int L;										// the lag index
	unsigned int n;							// the data index
	unsigned int curIdx;						// current Rxy index
	unsigned int halfLen = length - 1;			// half length of the data

	// Determine correlation for positive lags
	for (L = 0; L < length; L++)
	{
		curIdx = L + halfLen;
		lag[curIdx] = L;
		for (n = 0; n < length-L; n++)
			Rxy[curIdx] += x[n+L]*y[n];
	}
	// Determine correlation for negative lags
	for (L = -halfLen; L < 0; L++)
	{
		curIdx = L + halfLen;
		lag[curIdx] = L;
		for (n = -L; n < length; n++)
			Rxy[curIdx] += x[n+L]*y[n];
	}

	return;
}

//----------------------------- DSP_AutoCorrelation -------------------------------------//
void DSP_AutoCorrelation(float *x, float *Rxx, unsigned int length)
{
	int L;										// the lag index
	unsigned int n;							// the data index
	unsigned int curIdx;						// the current Rxx index
	unsigned int halfLen = length - 1;			// half length of the data

	// Determine correlation for positive lags
	for (L = 0; L < length; L++)
	{
		curIdx = L + halfLen;
		for (n = 0; n < length-L; n++)
			Rxx[curIdx] += x[n+L]*x[n];
	}
	// Determine correlation for negative lags
	for (L = -halfLen; L < 0; L++)
	{
		curIdx = L + halfLen;
		for (n = -L; n < length; n++)
			Rxx[curIdx] += x[n+L]*x[n];
	}

	return;
}

//----------------------------- DSP_GetSampleOffset -------------------------------------//
int DSP_GetSampleOffset(float *rxy, int *lags, unsigned int length)
{
	int sampOffset = 0;									// lag value associated with maxCorr
	unsigned int n;									// array index
	float maxCorr = DSP_MaxValueArray(rxy, length);		// max value o frxy

	for (n=0; n<length; n++)
	{
		if (rxy[n] == maxCorr)
		{
			sampOffset = lags[n];
			break;
		}
	}

	return sampOffset;
}


//----------------------------- DSP_ProduceAutoCorrMat -------------------------------------//
MATRIX *DSP_ProduceAutoCorrMat(float *x, unsigned int lenX, unsigned int dimA)
{
	MATRIX *A;				// autocorrelation matrix
	float *Rxx;				// autocorrelation vector
	unsigned int lenRxx;	// autocorrelation length
	unsigned int m, n;		// indexes
	unsigned int midIdx;	// center location of Rxx
	unsigned int curIdx;	// current index of Rxx for filling A
	float **a;				// double pointer to A

	//-----------Compute the AutoCorrelation of Input Vector x--------//
	lenRxx = 2*lenX-1;
	Rxx = (float*)calloc(lenRxx, sizeof(float));
	DSP_AutoCorrelation(x, Rxx, lenX);

	//---------------Compute the AutoCorrelation Matrix---------------//
	// build it up to the dimension dimA
	A = MatrixAllocate(dimA, dimA, sizeof(float));
	a = (float**)A->dataPtr;	
	curIdx = 0;		
	// Form the Toeplitz Matrix
	midIdx = lenX-1;
	for (m=0; m<dimA; m++)
	{
		for (n=0; n<dimA; n++)
		{
			curIdx = m + midIdx - n;
			a[m][n] = Rxx[curIdx];		
		}
	}

	// free dynamically-allocated memory
	free((float*)Rxx);

	return (A);
}
