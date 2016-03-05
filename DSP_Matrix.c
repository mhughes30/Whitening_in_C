
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DSP_Matrix.h"


/*--------------------------- 1D Array Operations ---------------------------------*/
void DSP_MultiplyArrays(float *output, float *input1, float *input2, unsigned int N)
{
	unsigned int i;

	for (i=0; i<N; i++)
		output[i] = input1[i]*input2[i];

	return;
}


void DSP_SubtractArrays2(float *output, float *input1, float *input2, unsigned int N)
{
	unsigned int i;

	for (i=0; i<N; i++)
	{
		output[i] = input1[i] - input2[i];
	}

	return;
}

void DSP_SubtractArrays(float *input1, float *input2, unsigned int N)
{
	unsigned int i;

	for (i=0; i<N; i++)
	{
		input1[i] -= input2[i];
	}

	return;
}


float DSP_MinValueArray(float *input, unsigned int N)
{
	unsigned int i;
	float minValue = input[0];

	for (i=1; i<N; i++) 
	{
		if (input[i] < minValue) 
			minValue = input[i];
	}
	
	return minValue;
}


float DSP_MaxValueArray(float *input, unsigned int N)
{	
	float maxValue = input[0];
	unsigned int i;

	for (i=1; i<N; i++) 
	{
		if (input[i] > maxValue) 
			maxValue = input[i];
	}
	
	return maxValue;
}


float DSP_MeanValueArray(float *input, unsigned int N)
{	
	double sum = input[0];
	unsigned int i;

	for (i=1; i<N; i++)
	{
		sum += input[i];
	}
	
	return (sum/(double)N);
}

void  DSP_ArrayCopy(float *input, float *subInput, unsigned int inputLen, unsigned int subLen)
{
	unsigned int m;

	for (m=0; m<subLen; m++)
	{
		subInput[m] = input[m];
	}
	return;
}

// Copy an array given a certain offset and numCopy
void DSP_ArrayCopy2(float *dataIn, float *dataOut, unsigned int numCopy, unsigned int offset)
{
	unsigned int m;
	unsigned int mStart = offset;
	unsigned int mEnd = numCopy + offset;
	unsigned int outIdx = 0;

	for (m=mStart; m<mEnd; m++)
	{
		outIdx = m - offset;
		dataOut[outIdx] = dataIn[m];
	}

	return;
}
/*----------------------------------------------------------------------------------*/


/*------------------------------ 2D Matrix Operations ------------------------------*/
MATRIX *MatrixAllocate(unsigned int numRows, unsigned int numCols, int elementSize)
{
	int i;
	MATRIX *A;
	float **floatMatrix;
	
	A = (MATRIX *)calloc(1,sizeof(MATRIX));
	if (!A)
	{
		printf("\nError allocating matrix!");
		exit(1);
	}
	A->numRows = numRows;
	A->numCols = numCols;
	A->elemSize = elementSize;
	
	//switch(elementSize)
	//{
	//	case sizeof(float):
	//	{
	//		float **floatMatrix;
			floatMatrix = (float**)calloc(numRows,sizeof(float *));
			if (!floatMatrix)
			{
				printf("\nError allocating matrix!");
				exit(1);
			}
			for (i=0; i<numRows; i++)
			{
				floatMatrix[i] = (float*)calloc(numCols,sizeof(float));
				if (!floatMatrix[i])
				{
					printf("\nError allocating matrix!");
					exit(1);
				}
			}
			A->dataPtr = (char**)floatMatrix;
		//	break;
		//}
		//case sizeof(double):
		//{
		//	double **doubleMatrix;
		//	doubleMatrix = (double**)calloc(numRows,sizeof(double *));
		//	if (!doubleMatrix)
		//	{
		//		printf("\nError allocating matrix!");
		//		exit(1);
		//	}
		//	for (i=0; i<numRows; i++)
		//	{
		//		doubleMatrix[i] = (double*)calloc(numCols,sizeof(double));
		//		if (!doubleMatrix[i])
		//		{
		//			printf("\nError allocating matrix!");
		//			exit(1);
		//		}
		//	}
		//	A->dataPtr = (char**)doubleMatrix;
		//	break;
		//}
		//default:
		//{
		//	printf("\nError allocating matrix, unssupported type!");
		//	exit(1);
		//}
	//}
	
	return (A);
}

MATRIX *MatrixInvert(MATRIX *A)
{
	MATRIX *Ai;
	float **a;
	float big, pivotInverse, temp, absElement;
	int *pivotFlag, *swapCol, *swapRow;
	int i, n, row, col, swap, irow, icol;
//	unsigned int m;
	
	// make sure the matrix is square
	if (A->numRows != A->numCols)
	{
		printf("Error: Non-Square Matrix!\n");
		exit(1);
	}
	if (!A->dataPtr)
	{
		printf("Error: Non-Square Matrix!\n");
		exit(1);
	}
	n = A->numRows;

	Ai = MatrixAllocate(n, n, sizeof(float));
	SCALE_MAT(A,Ai,1,n,n,0,0,float,float)
	a = (float**)Ai->dataPtr;
	
	//printf("Test\n");
	//for (i=0; i<n; i++)
	//{
	//	for (m=0; m<n; m++)
	//	{
	//		printf(" %f", a[i][m]);
	//	}
	//	printf("\n");
	//}
	//
	pivotFlag = (int*)calloc(n, sizeof(int));
	swapRow = (int*)calloc(n,sizeof(int));
	swapCol = (int*)calloc(n,sizeof(int));
	if (!pivotFlag || !swapRow || !swapCol)
	{
		printf("Allocation Error!\n");
		exit(1);
	}

	for (i=0; i<n; i++)
	{	
		big = 0.0;
		for (row=0; row<n; row++)
		{
			if (!pivotFlag[row])
			{
				for (col=0; col<n; col++)
				{
					if (!pivotFlag[col])
					{
						absElement = fabs(a[row][col]);
						if (absElement >= big)
						{
							big = absElement;
							irow = row;
							icol = col;
						}
					}
				}
			}
		}
		pivotFlag[icol]++;
		if (irow != icol)
		{
			for (col=0; col<n; col++)
			{
				temp = a[irow][col];
				a[irow][col] = a[icol][col];
				a[icol][col] = temp;
			}
		}
		swapRow[i] = irow;
		swapCol[i] = icol;
		if (a[icol][icol] == 0.0)
		{
			printf("Error: Singular Matrix!\n");
			exit(1);
		}
		pivotInverse = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (col=0; col<n; col++)
		{
			a[icol][col] = a[icol][col]*pivotInverse;
		}
		for (row=0; row<n; row++)
		{
			if (row != icol)
			{
				temp = a[row][icol];
				a[row][icol] = 0.0;
				for (col=0; col<n; col++)
				{
					a[row][col] = a[row][col] - a[icol][col]*temp;
				}
			}
		}
	}

	for (swap=n-1; swap>=0; swap--)
	{
		if (swapRow[swap] != swapCol[swap])
		{
			for (row=0; row <n; row++)
			{
				temp = a[row][swapRow[swap]];
				a[row][swapRow[swap]] = a[row][swapCol[swap]];
				a[row][swapCol[swap]] = temp;
			}
		}
	}
	
	free((char*)pivotFlag);
	free((char*)swapRow);
	free((char*)swapCol);

	return (Ai);
}

MATRIX *MatrixMult(MATRIX *A, MATRIX *B)
{
    MATRIX *C;
    int aR, aC, bC;

    if(B->numRows != A->numCols) 
	{
        printf("\nERROR in matrix_mult: B row and A col sizes must agree");
        printf("\nA is %dx%d and B is %dx%d\n", A->numRows, A->numCols, B->numRows, B->numCols);
        exit(1);
    }

    if(!A->dataPtr) 
	{
        printf("\nERROR: No A matrix allocated for matrix_mult");
        exit(1);
    }
    if(!B->dataPtr) 
	{
        printf("\nERROR: No B matrix allocated for matrix_mult");
        exit(1);
    }

/* simplify the row and col variables */
    aR = A->numRows;
    aC = A->numCols;
    bC = B->numCols;

/* allocate C to be of the highest ranking type of A and B
   (largest element size gives the highest rank ) */

    if(A->elemSize > B->elemSize)
         C = MatrixAllocate(aR, bC, A->elemSize);
    else
         C = MatrixAllocate(aR, bC, B->elemSize);

    MULT_MAT(A, B, C, aR, aC, bC, float, float, float)
  
    return(C);
}

void MatrixFree(MATRIX *A)
{
	int i;
	char **a;
	
	if (!A || !A->dataPtr || !A->numRows || !A->numCols)
	{
		printf("\nError: No Matrix to free!");
		exit(1);
	}
	a = A->dataPtr;
	for (i=0; i< A->numRows; i++)
	{
		free(a[i]);
	}
	free((char*)a);
	a = 0;				//NULL
	free((char*)A);

	return;
}

