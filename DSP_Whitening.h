

#ifndef DSP_WHITENING_H_
#define DSP_WHITENING_H_

// Produces the Whitening Filter Coefficients.
// Inputs: (1) x = preallocated input array
//		   (2) hW = the preallocated whitening filter array (the output)
//		   (3) lenX = length of x
//		   (4) dimHW = the length of the whitening filter hW
void DSP_ComputeWhiteningFilter(float *x, float *hW, unsigned int lenX, unsigned int dimHW);


#endif
