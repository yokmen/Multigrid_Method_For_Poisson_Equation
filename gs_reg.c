#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gs_reg.h"

double * gs_reg(int n, int *ia, int *ja, double *a, double *b, double *u_corr)
/*
   Goal
   ===
   Computation of the residual by the backward Gauss-Seidel method
   
   Arguments
   ===========
   (input)  n  - pointer of the number of unknowns
   (input)  ia - pointer of the vector 'ia' of the matrix A (fine grid)
   (input)  ja - pointer of the vector 'ja' of the matrix A (fine grid)
   (input)  a  - pointer of the vector 'a' of the matrix A (fine grid)
   (input)  b  - pointer of the vector 'b' of the matrix b (fine grid)
   (input)  u_corr  - pointer of the solution vector (fine grid)
                
   (output) u - pointer of the new solution vector
*/
{
	int i,j;
	double *r = malloc(n * sizeof(double)); /* r = b - y */
	double *y = malloc(n * sizeof(double)); /* y = A * x */
	double *d = malloc(n * sizeof(double)); /* d= triu(A)^(-1) * r */
	double *sum = malloc(n * sizeof(double)); /* sum = A * r */	
	double *u = malloc(n * sizeof(double));
	
	/* Initialization of each element of the vectors to 0 */
	for (i=0; i<n;i++){
		y[i] = 0.0;
		sum[i] = 0.0;
		d[i] = 0.0; 
		r[i] = 0.0;
	}
	
	for (i=0; i<n;i++){
		for (j=ia[i]; j< (ia[i+1]); j++){
			y[i] += a[j] * u_corr[ja[j]];   
		}
	}
	
	/* Calculation of vector r = b - y */
	for (i=0; i<n; i++){
		r[i] = b[i] - y[i];
	}
	
	/* Calcul du vecteur d */ 
	d[n-1] = r[n-1] / a[ia[n]-1];
	
	for (i=n-2; i > -1; i--){
		for(j=ia[i]; j<ia[i+1];j++){
			if (ja[j] > i){sum[i] += a[j] * d[ja[j]] ;}
		}
		
		for(j=ia[i]; j<ia[i+1];j++){
			if (ja[j] == i){d[i] = (r[i]-sum[i]) / (a[j]);}
		}
	}
	
	for (i=0; i<n ; i++){
		u[i] = u_corr[i] + d[i];
	}
	
	free(r); free(d); free(y); free(sum);
	
	return u;
}
