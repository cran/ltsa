// Dynamic vector and matrix memory allocation

#include "trenchR.h"

VECTOR Vector( long n )
{
	// allocate a double vector and set default values to 0
        
	VECTOR vector;
	vector=(VECTOR)Calloc(n, typeof(double));
	memset( vector, 0, n * sizeof(double) );

	return vector;
}

MATRIX Matrix( long n, long m )
{
	int i;
	// allocate a double matrix with subscript range m[n x m]
	MATRIX matrix;
   
	/* allocate pointers to rows */
	matrix = (MATRIX) Calloc( n, typeof(double*) );
   
	/* allocate rows and set pointers to them */
	matrix[0] = (double*) Calloc( n * m, typeof(double) );

   	memset( matrix[0], 0, n * m * sizeof(double) );
   
	for( i = 0; i < n; i++ ) matrix[i] = matrix[0] + m*i;

	/* return pointer to array of pointers to rows */
	return matrix;
}

void free_matrix( MATRIX matrix )
{
	Free( matrix[0] );
	Free( matrix );
} 

void free_vector( VECTOR vector )
{
	Free( vector );
}
 


