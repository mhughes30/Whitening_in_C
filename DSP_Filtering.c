
#include <stdlib.h>
#include <math.h>
#include "DSP_Matrix.h"
#include "DSP_Filtering.h"


// A filter that behaves exactly the same as Matlab's FIR filter.
int DSP_Filter1D(float *x, float *h, float *y, int lenX, int lenH) 
{
	int n;					// the output index, y[n]
	int j, k;				// the input indexes, x[j] and h[k]
	int minK, maxK, maxJ;	// allowed index ranges for j and k
	double temp;

	// error checking code
	if(!x || !h || !y)
		return 0;
	if (lenX <=0 || lenH <= 0)
		return 0;

	for (n=0; n<lenX; ++n)
	{
		y[n] = 0;
		minK  = 0 > (n-lenX+1) ? 0:n-lenX+1;
		maxK = (lenH-1) < n ? lenH-1:n;
		maxJ = n > (lenX-1) ? lenX-1:n;
		for (j=maxJ, k=minK; k<=maxK; --j, ++k)
		{
			temp = x[j] * h[k];
			y[n] += (float)temp;
			//y[n] += x[j] * h[k];
		}
	}

	return 1;	
}

// A filter that ignores the pre- and post-lag values during computation
int DSP_FilterNoDelay(float *x, float *h, float *y, int lenX, int lenH) 
{
	int i, j;						// indexes
	int coefLen2;					// half length of filter h
	int accLen;						// length of accumulation (for convolution)
    float acc;						// accumulated convolution value
    float *xPtr, *dataPtr, *xEnd;	// pointers to input 
	float *coefStart, *coefPtr;		// pointers to filter coefficients

	// error checking code
	if(!x || !h || !y)
		return 0;
	if (lenX <=0 || lenH <= 0)
		return 0;

	// Set up the filter start and length
    coefStart = h;
    coefLen2 = (lenH + 1)/2;

	// Set up input data pointers
    xEnd = x + lenX - 1;
	xPtr = x + coefLen2 - 1;

	// initial value of accumulation length
    accLen = coefLen2;

	for (i=0; i<lenX; i++)
	{
		// set up input and filter pointers for accumulation 
	    dataPtr = xPtr;
        coefPtr = coefStart;
		// do accumulation and set result to filtered output
        acc = (*coefPtr++) * (*dataPtr--);
        for(j=1; j<accLen; j++)
		{
            acc += (*coefPtr++) * (*dataPtr--);
		}
		y[i] = acc;
		// check for end case
        if (xPtr == xEnd) 
		{
            accLen--;			// reduce accumulation length
            coefStart++;       // increment to next filter coefficient
        }
		// if not at end, then check for startup, add to input pointer 
        else 
		{
            if (accLen < lenH)
			{
				accLen++;
			}
            xPtr++;
        }
    }

	return 1;	
}

// Subtracts the DC component
void DSP_RemoveDC2(float *input, float *output, int inputLen)
{
	float meanDC;								// the mean of input
	unsigned int n;								// array index
	meanDC = DSP_MeanValueArray(input, inputLen);

	for (n=0; n<inputLen; n++)
	{
		output[n] = input[n] - meanDC;
	}

	return;
}

// Subtracts the DC component
void DSP_RemoveDC(float *input, int inputLen)
{
	float meanDC;									// the mean of input
	unsigned int n;								// array index
	meanDC = DSP_MeanValueArray(input, inputLen);

	for (n=0; n<inputLen; n++)
	{
		input[n] -= meanDC;
	}

	return;
}


// A decimation filter
void DSP_Decimation(float *input, float *output, unsigned int len, unsigned int decFactor) 
{
	unsigned int n, m;			// array indexes
	n = 0;
	m = 0;

	while(n<len)
	{
		output[m] = input[n];
		m++;
		n += decFactor;
	}

	return;
}

// compute the decimation length
unsigned int DSP_DecimationLength(unsigned int len, unsigned int decFactor)
{
	return (unsigned int)ceil((double)len/(double)decFactor);
}
