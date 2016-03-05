
#ifndef DSP_FILTERING_H_
#define DSP_FILTERING_H_

// Performs convolution based filtering, and has the same output as Matlab's filter() function
// Inputs: (1) x = the data array
//		   (2) h = the filter array
//         (3) y = the prealloated output array of the filter
//         (4) lenX = the length of input x
//         (5) lenH = the lenght of input h
// Outputs: Returns 1 if the inputs are as specified, and returns 0 otherwise
int DSP_Filter1D(float *x, float *h, float *y, int lenX, int lenH); 


// Performs convolution based filtering, as above but eliminates filter delay
// The filter delay is the first (lenH-1)/2 samples for odd length filters. So, for an
// odd length filter, it excludes the (lenH-1)/2 samples at both the beginning and end
// of the filtered sequence. Note that for even-length filters, the sample delay is not an integer value, so
// the filter output will have a half-sample delay (In this case, lenH/2 samples are excluded at the beginning 
// and end).
// Inputs: (1) x = the data array
//		   (2) h = the filter array
//         (3) y = the prealloated output array of the filter
//         (4) lenX = the length of input x
//         (5) lenH = the lenght of input h
// Outputs: Returns 1 if the inputs are as specified, and returns 0 otherwise
int DSP_FilterNoDelay(float *x, float *h, float *y, int lenX, int lenH); 


// Removes the DC term from the input array and outputs the result
// Inputs: (1) input = the data array
//		   (2) output = the preallocated output array
//         (3) inputLen = the length of the input
void DSP_RemoveDC2(float *input, float *output, int inputLen);

// Removes the DC term from the input array
// Inputs: (1) input = the data array
//         (2) inputLen = the length of the input
void DSP_RemoveDC(float *input, int inputLen);

// Performs decimation in time
// Inputs: (1) input = the data array
//		   (2) output = the preallocated output array
//         (4) len = the length of input
//         (5) decFactor = the integer value by which to perform decimation
// Note: One must insure that decimating doesn't lower the sampling frequency below Nyquist. 
void DSP_Decimation(float *input, float *output, unsigned int len, unsigned int decFactor) ;


// Calculates the length of the resulting decimated array.
// Inputs: (1) len = the length of the array to be decimated
//		   (2) decFactor = the decimation factor integer
// Note: This should be run prior to DSP_Decimation in order to easily determine the length 
//       of the output array for preallocation.
unsigned int DSP_DecimationLength(unsigned int len, unsigned int decFactor);


#endif

