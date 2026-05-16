#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prolongation.h"
#include "indice.h"

double *prolongation(double h, int n, double Lx, double Ly, double x_i, double x_f, double y_i, double y_f, double u_e, double u_i, double *u_c)
/*           
  But
  ===
   Prolonger la solution de la grille grossiere vers la grille fine
   
  Arguments
  =========
  (input) h  - pas de discrétisation (en m) (grille fine)
  (input) n  - nombre d'inconnus dans le systeme (grille fine)
  (input) Lx  - longueur de la grille (en m)
  (input) Ly  - largeur de la grille  (en m)
  (input) x_i - abscisse du début de Γi (=xl dans l'énoncé) (en m)
  (input) x_f - abscisse de fin de Γi (=xr dans l'énoncé) (en m)
  (input) y_i - ordonnée de début de Γi (=yb dans l'énoncé) (en m)
  (input) y_f - ordonnée de fin de Γi (=yt dans l'énoncé) (en m)
  (input) u_i - potentiel électrique sur Γi (en V)
  (input) u_e - potentiel électrique sur Γe (en V)
  (input) u_c - pointeur vers le vecteur solution de la grille grossiere
                
  (output) P_c - pointeur vers le tableau P_c  (grille fine)
*/
{
	int ix, iy, ind, ind_c, mx, my, nx, ny, ind_sup_x, ind_sup_y, ind_inf_x, ind_inf_y;
	int ind_supx_supy, ind_supx_infy, ind_infx_supy, ind_infx_infy;
	
	mx = Lx / h + 1; /* nombre de pts dans la grille selon la direction x */
	my = Ly / h + 1; /* nombre de pts dans la grille selon la direction y */
	nx = mx - 2; /* noeuds de Dirichlet extérieur ne sont pas pris en compte */
	ny = my - 2; /* noeuds de Dirichlet extérieur ne sont pas pris en compte */
	
	double h_c = 2*h;
	
	double *P_c = malloc(n * sizeof(double));
	
	for (iy=0; iy<ny; iy++){
		for (ix=0; ix<nx ;ix++){
			
			/* If we are on ∂Γi or in Γi: we go to the next loop */
			if ( (ix >= x_i/h -1) && (ix <= x_f/h -1) && (iy >= y_i/h -1) && (iy <= y_f/h -1)){continue;}
			
			/* fine grid equation number */
			ind = indice(h, Lx, x_i, x_f, y_i, y_f, ix, iy);
			ind_inf_x = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix/2-1, (iy-1)/2);
			ind_sup_x = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix/2, (iy-1)/2);
			ind_inf_y = indice(h_c, Lx, x_i, x_f, y_i, y_f, (ix-1)/2, iy/2-1);
			ind_sup_y = indice(h_c, Lx, x_i, x_f, y_i, y_f, (ix-1)/2, iy/2);
			ind_infx_infy = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix/2-1, iy/2-1);
			ind_infx_supy = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix/2-1, iy/2);
			ind_supx_infy = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix/2, iy/2-1);
			ind_supx_supy = indice(h_c, Lx, x_i, x_f, y_i, y_f, ix/2, iy/2);
			
			/* Filling of the prolongation matrix P_c */
			
			/* If ix and iy odd => direct transfer */
			if ((ix %2 != 0) && (iy %2 != 0)){
				ind_c = indice(h_c, Lx, x_i, x_f, y_i, y_f, (ix-1)/2, (iy-1)/2);
				P_c[ind] = u_c[ind_c];
			}
			
			/* If ix even and iy odd */
			else if ( (ix %2 == 0) && (iy %2 != 0)){
				if ((ix == 0) || ( (ix==x_f/h) && (iy > y_i/h -2) && (iy < y_f/h ) )){
					P_c[ind] = (u_c[ind_sup_x])/2; // Linear interpolation along x with outer west Dirichlet condition
				}
				
				else if ( (ix == nx-1) || ((ix==x_i/h - 2) && (iy > y_i/h -2) && (iy < y_f/h)) ){
					P_c[ind] = (u_c[ind_inf_x])/2; // Linear interpolation along x with Dirichlet condition is exterior
				}
				
				else{
					P_c[ind] = (u_c[ind_inf_x] + u_c[ind_sup_x])/2; // Linear interpolation along x
				}
			}
			
			/* If ix odd and iy even */
			else if ( (ix %2 != 0) && (iy %2 == 0) ){
				if ( (iy == 0) || ( (iy==y_f/h) && (ix > x_i/h -2) && (ix < x_f/h) ) ){
					P_c[ind] = (u_c[ind_sup_y])/2; // Linear interpolation along y with outer south Dirichlet condition
				}
				
				else if ( (iy == ny-1) || ( (iy==y_i/h - 2) && (ix > x_i/h - 2) && (ix < x_f/h) ) ){
					P_c[ind] = (u_c[ind_inf_y])/2; // Linear interpolation along y with Dirichlet condition is exterior
				}
				
				else{
					P_c[ind] = (u_c[ind_inf_y] + u_c[ind_sup_y])/2;
				}
			}
			
			/* If ix and iy even */
			else{
				/* Exterior Dirichlet conditions */
				if (ix==0){
					if (iy == 0){
						P_c[ind] = (u_c[ind_supx_supy])/4;
					}
					
					else if (iy==ny-1){
						P_c[ind] = (u_c[ind_supx_infy])/4;
					}
					
					else{
						P_c[ind] = (u_c[ind_supx_infy] + u_c[ind_supx_supy])/4;
					}
				}
				
				else if (ix == nx-1){
					if (iy == 0){
						P_c[ind] = (u_c[ind_infx_supy])/4;
					}
					
					else if (iy==ny-1){
						P_c[ind] = (u_c[ind_infx_infy])/4;
					}
					
					else{
						P_c[ind] = (u_c[ind_infx_infy] + u_c[ind_infx_supy])/4;
					}
				}
				
				else if (iy == 0){
					P_c[ind] = (u_c[ind_supx_supy] + u_c[ind_infx_supy])/4;
				}
				
				else if (iy == ny-1){
					P_c[ind] = (u_c[ind_infx_infy] + u_c[ind_supx_infy])/4;
				}
				
				/* Interior Dirichlet conditions */
				else if ( (ix==x_i/h - 2) && (iy >= y_i/h-2)  && (iy <= y_f/h)){
					if (iy == y_i/h - 2){
						P_c[ind] = (u_c[ind_supx_infy] + u_c[ind_infx_supy] + u_c[ind_infx_infy])/4;
					}
					
					else if (iy == y_f/h){
						P_c[ind] = (u_c[ind_infx_infy] + u_c[ind_infx_supy] + u_c[ind_supx_supy])/4;
					}
					
					else{
						P_c[ind] = (u_c[ind_infx_supy] + u_c[ind_infx_infy])/4;
					}
				}
				
				else if ( (ix == x_f/h) && (iy >= y_i/h-2)  && (iy <= y_f/h)){
					if (iy == y_i/h - 2){
						P_c[ind] = (u_c[ind_infx_infy] + u_c[ind_supx_supy] + u_c[ind_supx_infy])/4;
					}
					
					else if (iy == y_f/h){
						P_c[ind] = (u_c[ind_infx_supy] + u_c[ind_supx_supy] + u_c[ind_supx_infy])/4;
					}
					
					else{
						P_c[ind] = (u_c[ind_supx_supy] + u_c[ind_supx_infy])/4;
					}
				}
				
				else if ( (iy == y_i/h - 2) && (ix > x_i/h-2)  && (ix < x_f/h) ){
					if(ix == x_i/h -1){
						P_c[ind] = (u_c[ind_infx_supy] + u_c[ind_supx_infy] + u_c[ind_infx_infy])/4;
					}
					
					else if(ix == x_f/h -1){
						P_c[ind] = (u_c[ind_infx_infy] + u_c[ind_supx_infy] + u_c[ind_supx_supy])/4;
					}
					
					else{
						P_c[ind] = (u_c[ind_supx_infy] + u_c[ind_infx_infy])/4;
					}
				}
				
				else if ( (iy == y_f/h) && (ix > x_i/h-2)  && (ix < x_f/h) ){
					if (ix == x_i/h - 1){
						P_c[ind] = (u_c[ind_supx_supy] + u_c[ind_infx_supy] + u_c[ind_infx_infy])/4;
					}
					
					else if (ix == x_f/h - 1){
						P_c[ind] = (u_c[ind_infx_supy] + u_c[ind_supx_supy] + u_c[ind_supx_infy])/4; 
					}
					
					else {
						P_c[ind] = ( u_c[ind_supx_supy] + u_c[ind_infx_supy] )/4;
					}
				}
				
				/* General case */
				else {
					P_c[ind] = ( u_c[ind_supx_infy] + u_c[ind_supx_supy] + u_c[ind_infx_supy] + u_c[ind_infx_infy] )/4;
				}
			}
		}
	}
	
	return P_c;
}
