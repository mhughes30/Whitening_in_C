
#ifndef DSP_MATRIX_H_
#define DSP_MATRIX_H_

/*------------------------- 1D Array Operations ---------------------------*/
// Multiplies two 1D arrays and returns the result.
// Inputs: (1) output = preallocated output array
//		   (2) input1
//         (3) input2
//         (4) N = the length of the arrays
// Note:	input1 and input2 must be the same length
void  DSP_MultiplyArrays(float *output, float *input1, float *input2, unsigned int N);

// Copies the contents of one 1D array and places it in another
// Inputs: (1) input = prealloated input array to be copied from
//		   (2) subInput = preallocated output array to be copied to
//         (3) inputLen = length of input
//         (4) subLen = length of output
// Note:	subInput can be shorter than or equal to the length of input.
void  DSP_ArrayCopy(float *input, float *subInput, unsigned int inputLen, unsigned int subLen);

// Subtract two 1D arrays and returns the result.
// Inputs: (1) output = preallocated output array
//		   (2) input1
//         (3) input2
//         (4) N = the length of the arrays
// Note:	input1 and input2 must be the same length
//          Performs (input1 - input2)
void DSP_SubtractArrays2(float *output, float *input1, float *input2, unsigned int N);

// Subtract two 1D arrays and returns the result.
// Inputs: (1) input1 = the input and output array
//         (2) input2 = the array to subtract from input1
//         (3) N = the length of the arrays
// Note:	input1 and input2 must be the same length
//          Performs input1 = (input1 - input2)
void DSP_SubtractArrays(float *input1, float *input2, unsigned int N);

// Copies the contents of one 1D array and places it in another
// Inputs: (1) dataIn  = prealloated input array to be copied from
//		   (2) dataOut = preallocated output array to be copied to
//         (3) numCopy = number of data points to copy
//         (4) offset  = starting element in dataIn from which to begin copying to dataOut
// Note:	dataOut length equals offset + dataIn length
void DSP_ArrayCopy2(float *dataIn, float *dataOut, unsigned int numCopy, unsigned int offset);

// Determines minimum value of a 1D array.
// Inputs: (1) input
//		   (2) N = length of input
float DSP_MinValueArray(float *input, unsigned int N);

// Determines maximum value of a 1D array.
// Inputs: (1) input
//		   (2) N = length of input
float DSP_MaxValueArray(float *input, unsigned int N);

// Determines mean value of a 1D array.
// Inputs: (1) input
//		   (2) N = length of input
float DSP_MeanValueArray(float *input, unsigned int N);
/*----------------------------------------------------------------------------*/

/*--------------------------- 2D Matrix Operations ---------------------------*/
// a 2D MATRIX Structure
typedef struct{
	int elemSize;				// data type size: sizeof(float), sizeof(int), etc
	unsigned int numRows;		// number of rows in matrix
	unsigned int numCols;		// number of columns in matrix
	char **dataPtr;				// 2D pointer to the matrix data
}	MATRIX;

// Macro for coying and scaling a matirx
// The backslashes are required.
// Inputs:	(1)	a = pointer to input MATRIX structure.
//			(2)	b = pointer to resultant MATRIX structure.
//			(3) s = scale factor variable to be multiplied.
//			(4) rows = number of rows in matrix b
//			(5) cols = number of columns in matrix b
//			(6) rowoff = number of rows to offset matrix b
//			(7) coloff = number of columns to offset matrix b
//			(8) typea  = legal C type describing the type of a
//			(9) typeb  = legal C type describing the type of b
#define SCALE_MAT(a,b,s,rows,cols,rowoff,coloff,typea,typeb)\
{ \
	typea **_AMX = (typea **)a->dataPtr;  \
    typeb **_BMX = (typeb **)b->dataPtr;  \
	typea *_PTA;  \
    typeb *_PTB;  \
    int _IX,_JX;  \
    for(_IX = 0 ; _IX < rows ; _IX++) \
	{  \
		_PTB = _BMX[_IX];  \
        _PTA = _AMX[_IX + rowoff] + coloff;  \
        for(_JX = 0 ; _JX < cols ; _JX++)  \
			*_PTB++ = (typeb)(s * (*_PTA++));  \
	}  \
}    

// Macro for coying and scaling a matirx
// The backslashes are required.
// Inputs:	(1)	a = pointer to first MATRIX structure.
//			(2)	b = pointer to second MATRIX structure.
//			(3) c = pointer to result MATRIX structure.
//			(4) rowsa = number of rows in matrix a
//			(5) colsa = number of columns in matrix a
//			(6) colsb = number of columns in matrix b
//			(7) typea = legal C type describing the type of a
//			(8) typeb = legal C type describing the type of b
//			(9) typec = legal C type describing the type of c
#define MULT_MAT(a,b,c,rowsa,colsa,colsb,typea,typeb,typec) \
{  \
	typea **_AMX = (typea **)a->dataPtr;  \
    typeb **_BMX = (typeb **)b->dataPtr;  \
    typec **_CMX = (typec **)c->dataPtr;  \
    typea *_PTA;  \
    typeb *_PTB;  \
    typec *_PTC;  \
    int _IX,_JX,_KX;  \
    for(_IX = 0 ; _IX < rowsa ; _IX++) \
	{  \
		_PTC = _CMX[_IX];  \
		_PTB = _BMX[0];  \
        for(_JX = 0 ; _JX < colsb ; _JX++) \
		{  \
			_PTA = _AMX[_IX];  \
            *_PTC = (*_PTA++) * (*_PTB++);  \
            for(_KX = 1 ; _KX < colsa ; _KX++)  \
				*_PTC += (*_PTA++)* _BMX[_KX][_JX];  \
            _PTC++;  \
        }  \
    }  \
}    

// Allocates Memory for a 2D Matrix
// Inputs: (1) numRows = # of rows
//		   (2) numCols = # of columns
//		   (3) elementSize = data type size: sizeof(float), etc
// Output: (1) the allocated matrix
MATRIX *MatrixAllocate(unsigned int numRows, unsigned int numCols, int elementSize);

// Uses the Gauss-Jordan Matrix Inversion Algorithm
// Inputs: (1) A = preallocated square matrix
// Output: (1) The inverted matrix
MATRIX *MatrixInvert(MATRIX *A);

// Performs matrix multiplication (m by n)*(n by x)
// Inputs: (1) A = preallocated matrix
//		   (2) B = preallocated matrix
// Output: (1) the product of A*B
MATRIX *MatrixMult(MATRIX *A, MATRIX *B);

// Frees memory associated with matrix A
// Inputs: (1) A = preallocated matrix
void MatrixFree(MATRIX *A);
/*----------------------------------------------------------------------------*/

#endif

