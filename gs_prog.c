#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gs_prog.h"

double * gs_prog(int n, int *ia, int *ja, double *a, double *b, double *u)
/*
   Goal
   ===
   Computation of the residual by the forward Gauss-Seidel method
   
   Arguments
   ===========
   (input)  n  - pointer of the number of unknowns
   (input)  ia - pointer of the vector 'ia' of the matrix A (fine grid)
   (input)  ja - pointer of the vector 'ja' of the matrix A (fine grid)
   (input)  a  - pointer of the vector 'a' of the matrix A (fine grid)
   (input)  b  - pointer of the vector 'b' of the matrix b (fine grid)
   (input)  u  - pointer of the solution vector (fine grid)
                
   (output) r_prog - pointer of the residual vector
*/
{
	int i,j;
	double *r = malloc(n * sizeof(double)); /* r = b - Au */
	double *y = malloc(n * sizeof(double)); /* y = Au */
	double *d = malloc(n * sizeof(double)); /* d= tril(A)^(-1) * r */
	double *sum = malloc(n * sizeof(double)); /* sum = A * r */
	double *r_prog = malloc(n * sizeof(double)); /* residual vector with the new approximation of the solution vector */
	
	/* Initialization of each element of the vectors to 0 */
	for (i=0; i<n;i++){
		sum[i] = 0.0;
		d[i] = 0.0; 
		y[i] = 0.0;
		r[i] = 0.0;
		r_prog[i] = 0.0;
	}
	
	/* Calculation of the vector y = Au and of the residual */
	for (i=0; i<n; i++){
		for (j=ia[i]; j< (ia[i+1]); j++){
			y[i] += a[j] * u[ja[j]];   
		}
		
		r[i] = b[i] - y[i];
	}
	
	/* Calculation of vector d  */ 
	d[0] = r[0] / a[0];
	
	for (i=1; i<n; i++){
		for(j=ia[i]; j<ia[i+1];j++){
			if (ja[j] < i){sum[i] += a[j] * d[ja[j]] ;}
		}
		
		for(j=ia[i]; j<ia[i+1];j++){
			if (ja[j] == i){d[i] = (r[i]-sum[i]) / (a[j]);}
		}
	}
	
	/* Calculation of the new approximation of u */
	for (i=0;i<n;i++){
		u[i] = u[i] + d[i];
	}
	
	/* Calculation of the new residual */
	for (i=0; i<n; i++){
		y[i] = 0.0;
		
		for (j=ia[i]; j< (ia[i+1]); j++){
			y[i] += a[j] * u[ja[j]];   
		}
		
		r_prog[i] = b[i] - y[i];
	}
	
	free(r); free(d); free(sum); free(y);
	
	return r_prog;
}
