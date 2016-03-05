
#include "DSP_Whitening.h"
#include "DSP_Correlation.h"


void DSP_ComputeWhiteningFilter(float *x, float *hW, unsigned int lenX, unsigned int dimHW)
{
	MATRIX *RXX;				// autocorrelation matrix of x
	MATRIX *W;					// white-noise autocorrelation vector
	MATRIX *InvRXX;				// Inverse of RXX
	MATRIX *HW;					// Matrix form of whitening filter coefficients
	unsigned int dimRXX;		// size of RXX
	unsigned int m;			// index
	float **wPtr;				// data pointer to W
	float **rPtr;				// data pointer to RXX
	float **hPtr;				// data pointer to HW

	//------------------ Produce the Auto-Correlation Matrix, RXX --------------------//
	RXX = DSP_ProduceAutoCorrMat(x, lenX, dimHW);
	dimRXX = RXX->numRows;

	//---------------- Produce the White Noise Auto-Correlation Matrix ----------------//
	// (a column vector)
	// needs to be a MATRIX in order to perform matrix multiplication.
	W = MatrixAllocate(dimRXX, 1, sizeof(float));
	wPtr = (float**)W->dataPtr;
	rPtr = (float**)RXX->dataPtr;
	wPtr[0][0] = rPtr[0][0];				// set average power to that of the input x

	//------------------------- Compute the Inverse of RXX ---------------------------//
	InvRXX = MatrixInvert(RXX);
	MatrixFree(RXX);

	//------------------- Compute the Whitening Filter Coefficients ------------------//
	HW = MatrixMult(InvRXX, W);
	MatrixFree(InvRXX);
	MatrixFree(W);
	hPtr = (float**)HW->dataPtr;
	// convert MATRIX to a vector and scale the magnitude, so that the first coefficient is 1.
	for (m=0; m<dimRXX; m++)
	{
		hW[m] = -1*hPtr[m][0]/hPtr[0][0];
	}
	// now, set the first coefficient to zero.
	hW[0] = 0;

	return;
}

