#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prob.h"
#include "indice.h"

int prob(double h, int *n, int **ia, int **ja, double **a, double **b, 
         double Lx, double Ly, double x_i, double x_f, double y_i, double y_f, double u_i, double u_e)
/*           
   Goal
   ===
   Generate the n x n linear system
                          
                             Au = b                                   

   which corresponds to the discretization on a regular m x m Cartesian grid of the two-dimensional Poisson equation
              
            d    d        d    d         
         - == ( == u ) - == ( == u )  =  0    on [0,Lx] x [0,Ly] \ [x_i,x_f] x [y_i,y_f]
           dx   dx       dy   dy           

  with Dirichlet boundary conditions
         
         u = u_i on Γi and u = u_e on Γe 
  
  The numbering of the unknowns is lexicographic, the direction x being traversed before that of y.
  Matrix A is returned in CSR format which is defined via three arrays: 'ia', 'ja' and 'a'.

  Arguments
  =========
  (input)  Lx  - length of the outside capacitor [m]
  (input)  Ly  - height of the outside capacitor [m]
  (input)  h  - discretization step [m]
  (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
  (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
  (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
  (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m] 
  (input)  u_i - electric potential on Γi [V]
  (input)  u_e - electric potential on Γe [V]

  (input)  n  - pointer of the number of unknowns
  (input)  ia - pointer of the vector 'ia' of the matrix A (fine grid)
  (input)  ja - pointer of the vector 'ja' of the matrix A (fine grid)
  (input)  a  - pointer of the vector 'a' of the matrix A (fine grid)
  (input)  b  - pointer of the vector 'b' of the matrix b (fine grid)
*/
{
    int  nnz, ix, iy, ind, mx, my, nx, ny;
    double invh2;
    int X = (x_f - x_i)/ h +1;
	int Y = (y_f - y_i)/h + 1;
    
    mx = Lx / h + 1;
    my = Ly / h + 1;
    nx = mx - 2;
    ny = my - 2;
    
    if(nx < 0) {
        printf("\n ERREUR : mx = %d n'est pas une valeur valide\n\n",mx);
        return 1;
    }
    
    else if (ny < 0) {
		printf("\n ERREUR : my = %d n'est pas une valeur valide\n\n",my);
        return 1;
    }
    
    invh2 = 1 / (h*h); /* invh2 is the inverse square of h */
    
    *n  = nx*ny - X * Y ; /* number of unknowns  */
   
    nnz = 5*nx*ny - 2*nx - 2*ny - 5*X*Y - 2*(X + Y); /* number of nonzero elements */
     
    /* array allocation */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* main part: matrix filling */

    ind = 0; /* in case nx = 0 or ny =0  */
    nnz = 0;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
			
			/* If we are on ∂Γi or in Γi: we go to the next loop */
			if ( (ix >= x_i/h -1) && (ix <= x_f/h -1) && (iy >= y_i/h -1) && (iy <= y_f/h -1)){continue;}
			
			/* equation number */
			ind = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy);

			/* mark start of next row in array 'ia' */
			(*ia)[ind] = nnz;

			/* calculate the right-hand side */
			(*b)[ind] = 0; /* Right hand side initialization */
				
			/* Dirichlet conditions at the southern edge */
		    if ((iy ==0) || ( (iy==y_f/h) && (ix > x_i/h -2) && (ix < x_f/h ))){
				if (iy ==0){
					(*b)[ind] += u_e * invh2; /* Dirichlet conditions at the outer boundary */
				}
				else{
					(*b)[ind] += u_i * invh2; /* Dirichlet conditions at the inner boundary */
				} 
			} 
			    
		    else{
				(*a)[nnz] = -invh2; 
				(*ja)[nnz] = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy-1);;
				nnz++;
			}
		
			/* Dirichlet conditions at the western edge */
			if ((ix ==0) || ( (ix==x_f/h) && (iy > y_i/h -2) && (iy < y_f/h ) )){
		     	if (ix ==0){
						(*b)[ind] += u_e * invh2; /* Dirichlet conditions at the outer boundary */
				}
				else{
					(*b)[ind] += u_i * invh2;/* Dirichlet conditions at the inner boundary */
				} 
			} 
			    
			else{
				(*a)[nnz] = -invh2; 
				(*ja)[nnz] = ind -1;
				nnz++;
			}
				
			/* row fill: diagonal element */
			(*a)[nnz] = 4*invh2; 
			(*ja)[nnz] = ind;
			nnz++;
				
			/* Dirichlet conditions at the eastern boundary */
			if ((ix == nx - 1) || ( (ix==x_i/h - 2) && (iy > y_i/h -2) && (iy < y_f/h) )){
				if (ix == nx - 1){
					(*b)[ind] += u_e * invh2; /* Dirichlet conditions at the outer boundary */
				}
				else{
					(*b)[ind] += u_i * invh2; /* Dirichlet conditions at the inner boundary */
				} 
			} 
			    
		    else{
				(*a)[nnz] = -invh2; 
				(*ja)[nnz] = ind +1;
				nnz++;
			}
				
			/* Dirichlet conditions at the northern boundary */
			if ((iy == ny -1) || ( (iy==y_i/h - 2) && (ix > x_i/h -2) && (ix < x_f/h ) )){
				if (iy == ny-1){
					(*b)[ind] += u_e * invh2; /* Dirichlet conditions at the outer boundary */
				}
				else{
					(*b)[ind] += u_i * invh2; /* Dirichlet conditions at the inner boundary */
				}
			}
			    
		   else{
		    	(*a)[nnz] = -invh2; 
				(*ja)[nnz] = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy+1);
				nnz++;
			}
		}
	}
	/* 	last element of array 'ia'  */
	(*ia)[ind + 1] = nnz;
	
	return 0;
}
