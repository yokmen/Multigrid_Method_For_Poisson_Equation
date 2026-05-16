#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prob_grossier.h"
#include "indice.h"

int prob_grossier(double h_c, int **ia_c, int **ja_c, double **a_c,
         double Lx, double Ly, double x_i, double x_f, double y_i, double y_f, double u_i, double u_e)
/*           
   Goal
   ===
   Generate the linear system n_c x n_c for the coarse problem
                          
                             Au = b                                   

   which corresponds to the discretization on a regular m x m
   Cartesian grid of the two-dimensional Poisson equation
              
            d    d        d    d         
         - == ( == u ) - == ( == u )  =  0    on [0,Lx] x [0,Ly] \ [x_i,x_f] x [y_i,y_f]
           dx   dx       dy   dy           

  with Dirichlet boundary conditions
         
         u = u_i sur Γi et u = u_e sur Γe 
  
  The numbering of the unknowns is lexicographic, the direction x being traversed before that of y.
  Matrix A is returned in CSR format which is defined via three arrays: 'ia_c', 'ja_c' and 'a_c'.

  Arguments
  =========
  (input)  h_c  - discretization step [m] (coarse grid)
  (input)  Ly  - height of the outside capacitor [m]
  (input)  h  - discretization step [m]
  (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
  (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
  (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
  (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m] 
  (input)  u_i - potentiel électrique sur Γi [V]
  (input)  u_e - potentiel électrique sur Γe [V]
                
  (input)  ia_c - pointer of the vector 'ia_c' of the matrix A (coarse grid)
  (input)  ja_c - pointer of the vector 'ja_c' of the matrix A (coarse grid)
  (input)  a_c  - pointer of the vector 'a_c' of the matrix A (coarse grid)
*/
{
	int  nnz_c, ix, iy, ind, mx_c, my_c, nx_c, ny_c;
	double invh2_c;
	int X_c = (x_f - x_i)/ h_c +1;
	int Y_c = (y_f - y_i)/h_c + 1;
	
	mx_c = Lx / h_c + 1;
	my_c = Ly / h_c + 1;
	nx_c = mx_c - 2;
	ny_c = my_c - 2;
	
	
	if(nx_c < 0) {
		printf("\n ERREUR : mx_c = %d n'est pas une valeur valide\n\n",mx_c);
		return 1;
	}
	
	else if (ny_c < 0) {
		printf("\n ERREUR : my_c = %d n'est pas une valeur valide\n\n",my_c);
		return 1;
	}
	
	invh2_c = 1 / (h_c*h_c); /* invh2_c is the inverse square of h_c */
	
	int n_c  = nx_c*ny_c - X_c * Y_c ; /* number of unknowns */
	
	nnz_c = 5*nx_c*ny_c - 2*nx_c - 2*ny_c - 5*X_c*Y_c - 2*(X_c + Y_c); /* number of nonzero elements */
	
	/* array allocation */
	
	*ia_c  = malloc((n_c + 1) * sizeof(int));
	*ja_c  = malloc(nnz_c * sizeof(int));
	*a_c   = malloc(nnz_c * sizeof(double));
	
	if (*ia_c == NULL || *ja_c == NULL || *a_c == NULL) {
		printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
		return 1;
	}
	
	/* main part: matrix filling */
	
	ind = 0; /* in case nx = 0 or ny =0 */
	nnz_c = 0;
	for (iy = 0; iy < ny_c; iy++) {
		for (ix = 0; ix < nx_c; ix++) {
			
			/* If we are on ∂Γi or in Γi: we go to the next loop */
			if ( (ix >= x_i/h_c -1) && (ix <= x_f/h_c -1) && (iy >= y_i/h_c -1) && (iy <= y_f/h_c -1)){continue;}
			
			/* equation number */
			ind = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix, iy);
			
			/* mark start of next row in array 'ia' */
			(*ia_c)[ind] = nnz_c;
			
			/* Filling of matrices a, ia and ja */
		    if ( !( (iy ==0) || ( (iy==y_f/h_c) && (ix > x_i/h_c -2) && (ix < x_f/h_c ))  ) ){
				(*a_c)[nnz_c] = -invh2_c; 
				(*ja_c)[nnz_c] = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix, iy-1);
				nnz_c++;
			} 
		
			if ( !( (ix ==0) || ( (ix==x_f/h_c) && (iy > y_i/h_c -2) && (iy < y_f/h_c ) ) ) ){
				(*a_c)[nnz_c] = -invh2_c; 
				(*ja_c)[nnz_c] = ind -1;
				nnz_c++;
			} 
			
			
			(*a_c)[nnz_c] = 4*invh2_c;
			(*ja_c)[nnz_c] = ind;
			nnz_c++;
			
			if ( !((ix == nx_c - 1) || ( (ix==x_i/h_c - 2) && (iy > y_i/h_c -2) && (iy < y_f/h_c) ) ) ){
				(*a_c)[nnz_c] = -invh2_c; 
				(*ja_c)[nnz_c] = ind +1;
				nnz_c++;
			} 
			
			if ( !((iy == ny_c -1) || ( (iy==y_i/h_c - 2) && (ix > x_i/h_c -2) && (ix < x_f/h_c ) ))){
				(*a_c)[nnz_c] = -invh2_c; 
				(*ja_c)[nnz_c] = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix, iy+1);
				nnz_c++;
			}
		}
	}
	
	/* last element of array 'ia' */
	(*ia_c)[ind + 1] = nnz_c;
	
	return 0;
}
