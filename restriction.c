#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "restriction.h"
#include "indice.h"

double * restriction(double h, double Lx,double Ly, double x_i, double x_f, double y_i, double y_f, double *r_prog)
/*
  Goal
  ===
  Restrict the solution from the fine grid to the coarse grid
  
  Arguments
  =========
  (input)  Lx  - length of the outside capacitor [m]
  (input)  Ly  - height of the outside capacitor [m]
  (input)  h  - discretization step [m]
  (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
  (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
  (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
  (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m]
  (input)  r   - pointer to the residual vector calculated by forward Gauss-Seidel
                
  (output) r_c - pointer to the residual vector computed after restriction on the coarse grid
*/
{
	/* Parameters for the fine grid */
	int ix, iy, ind, mx, my, nx, ny, ind_c, mx_c, my_c, nx_c, ny_c, X_c, Y_c, n_c;
	int ind_sup_y, ind_inf_y, ind_inf_x, ind_sup_x, ind_supx_infy, ind_supx_supy, ind_infx_infy, ind_infx_supy;
	
	mx = Lx / h + 1;
	my = Ly / h + 1;
	nx = mx - 2;
	ny = my - 2;
	
	/* Parameters for coarse grid */
	double h_c = 2*h;
	X_c = (x_f - x_i)/ h_c +1;
	Y_c = (y_f - y_i)/h_c + 1;
	mx_c = Lx / h_c + 1; /* number of points in the grid along the x direction */
	my_c = Ly / h_c + 1; /* number of points in the grid along the y direction */
	nx_c = mx_c - 2;
	ny_c = my_c - 2;
	n_c  = nx_c*ny_c - X_c * Y_c ;
	
	double *r_c = malloc(n_c * sizeof(double));
	
	for (iy=0; iy<ny; iy++){
		for (ix=0; ix<nx ;ix++){
			
			/* If we are on ∂Γi or in Γi: we go to the next loop */
			if ( (ix >= x_i/h -1) && (ix <= x_f/h -1) && (iy >= y_i/h -1) && (iy <= y_f/h -1)){continue;}
			
			/* We only take into account the points on the coarse grid, i.e. ix and iy odd */
			if ( (ix %2 ==0) || (iy %2 == 0) ){continue;}
			
			/* equation number  */
			
			ind = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy);
			ind_sup_y = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy+1);
			ind_inf_y = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy-1);
			ind_inf_x = indice(h, Lx, x_i, x_f, y_i, y_f, ix-1, iy);
			ind_sup_x = indice(h, Lx, x_i, x_f, y_i, y_f, ix+1, iy);
			ind_supx_infy = indice(h, Lx, x_i, x_f, y_i, y_f, ix+1, iy-1);
			ind_supx_supy = indice(h, Lx, x_i, x_f, y_i, y_f, ix+1, iy+1);
			ind_infx_infy = indice(h, Lx, x_i, x_f, y_i, y_f, ix-1, iy-1);
			ind_infx_supy = indice(h, Lx, x_i, x_f, y_i, y_f, ix-1, iy+1);
			ind_c = indice(h_c, Lx, x_i, x_f, y_i, y_f, (ix-1)/2, (iy-1)/2);
			
			/* Dirichlet conditions at the southern edge  */
			if((iy==y_f/h) && (ix >= x_i/h -2) && (ix <= x_f/h)){
				if (ix == x_i/h - 2){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_supy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16 ;
					}
					
				else if (ix == x_f/h){
					r_c[ind_c] =  r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (ix == x_f/h-1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_sup_x] + r_prog[ind_sup_y])/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if(ix == x_i/h - 1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_supy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else {
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_supy] + r_prog[ind_infx_supy] )/16;
				}
			}
			
			/* Dirichlet conditions at the western edge */
			else if ((ix==x_f/h) && (iy >= y_i/h -2) && (iy <= y_f/h )){
				if (iy == y_i/h - 2){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_infy])/16;
				}
				
				else if (iy == y_f/h){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (iy == y_i/h -1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_infy])/16;
				}
				
				else if (iy == y_f/h - 1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_supy] )/16;
				}
				
				else{
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy])/16;
				}
			}
			
			/* Dirichlet conditions at the eastern edge */
			else if ((ix==x_i/h - 2) && (iy >= y_i/h -2) && (iy <= y_f/h)){
				if (iy == y_i/h-2){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (iy == y_f/h){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + (r_prog[ind_supx_supy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (iy == y_i/h-1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (iy == y_f/h -1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_y] )/8 + (r_prog[ind_supx_supy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else{
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_y] )/8 + (r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
			}
			
			/* Dirichlet conditions at the northern edge */
			else if ((iy==y_i/h - 2) && (ix >= x_i/h -2) && (ix <= x_f/h)){
				if (ix == x_i/h-2){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (ix == x_f/h){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_infy])/16;
				}
				
				else if (ix == x_i/h-1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x])/8 + ( r_prog[ind_supx_infy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
				
				else if (ix == x_f/h-1){
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x])/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_infy])/16;
				}
				
				else{
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x])/8 + ( r_prog[ind_supx_infy] + r_prog[ind_infx_infy])/16;
				}
			}
			
			else{
					r_c[ind_c] = r_prog[ind]/4 + (r_prog[ind_inf_x] + r_prog[ind_inf_y] + r_prog[ind_sup_x] + r_prog[ind_sup_y] )/8 + ( r_prog[ind_supx_infy] + r_prog[ind_supx_supy] + r_prog[ind_infx_infy] + r_prog[ind_infx_supy] )/16;
				}
		}
	}
	
	return r_c;
}
