#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multigrid.h"
#include "gs_prog.h"
#include "gs_reg.h"
#include "restriction.h"
#include "prolongation.h"
#include "prob_grossier.h"
#include "prob.h"
#include "umfpk.h"
#include "norme.h"

double *multigrid(int cycle, int level, double h, double *b, double Lx, double Ly, double x_i, double x_f, double y_i, double y_f, double u_e, double u_i)
/*
   Goal
   ===
   Computation of the solution by the multigrid method
   
   Arguments
   ===========
   (input)  cycle  - number of cycle of the conjugate gradient method
   (input)  level  - number of level for the multigrid
   (input)  h  - discretization step [m]
   (input)  b  - pointer of the vector 'b' (fine grid)
   (input)  Lx  - length of the outside capacitor [m]
   (input)  Ly  - height of the outside capacitor [m]
   (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
   (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
   (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
   (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m] 
   (input) u_i - potentiel électrique sur Γi [V]
   (input) u_e - potentiel électrique sur Γe [V] 
                
   (output) u - Vector solution 
*/
{
	
	/* create data file */
	FILE *data3 = fopen("data3.txt", "w+");
	fprintf (data3, "%i %.16f\n",0, 0.0);
	
	int i,k,m;
	
	int  n, X, Y, nx, ny, n_c, nx_c, ny_c, X_c, Y_c, n_prov, nx_prov, ny_prov, X_prov, Y_prov,nnz_prov;
	double h_prov, h_c, norme_residu;
	
	/* Creation of the various data for the creation of the graphs */
	char pas[10];
	char nombre_cycle[10];
	char niveau[10];
	sprintf(nombre_cycle,"%i",cycle+1);
	sprintf(pas,"%f",h);
	sprintf(niveau,"%i",level);
	
	/* Initialization of the different tables */
	double ** u = malloc((level-1) * sizeof(double *));
	double ** r_prog = malloc((level-1) * sizeof(double *));
	double ** r_c = malloc((level-1) * sizeof(double *));
	int ** ia = malloc((level-1) * sizeof(int *));
	int ** ja = malloc((level-1) * sizeof(int *));
	double ** a  = malloc((level-1) * sizeof(double *));
	double ** P_c = malloc(level * sizeof(double *));
	double ** u_corr = malloc(level * sizeof(double *));
	
	/* Parameters for the fine grid */
	X = (x_f - x_i)/ h +1;
	Y = (y_f - y_i)/h + 1;
	nx = Lx / h - 1; /* number of points in the grid in the x-direction */
	ny = Ly / h - 1; /* number of points in the grid in the y-direction */
	n  = nx*ny - X * Y ; /* number of unknowns */
	
	/* Parameters for the coarsest grid */
	double h_coarsest = pow(2,level-1)*h;
	int X_coarsest = (x_f-x_i)/h_coarsest+1;
	int Y_coarsest = (y_f-y_i)/h_coarsest+1;
	int nx_coarsest = Lx / h_coarsest - 1; /* number of points in the coarsest grid in the x-direction */
	int ny_coarsest = Ly / h_coarsest - 1; /* number of points in the coarsest grid in the x-direction */
	int n_coarsest = nx_coarsest*ny_coarsest - X_coarsest * Y_coarsest ; /* number of unknowns */
	int nnz_coarsest = 5*nx_coarsest*ny_coarsest - 2*nx_coarsest - 2*ny_coarsest - 5*X_coarsest*Y_coarsest - 2*(X_coarsest + Y_coarsest); /* number of nonzero elements */
	
	double *a_coarsest = malloc(nnz_coarsest * sizeof(double));
	int *ja_coarsest = malloc(nnz_coarsest * sizeof(int));
	int *ia_coarsest = malloc((n_coarsest+1) * sizeof(int));
	double *u_coarsest = malloc(n_coarsest * sizeof(double));
	
	if (prob_grossier(h_coarsest, &ia_coarsest, &ja_coarsest, &a_coarsest, Lx, Ly, x_i, x_f, y_i, y_f, u_i, u_e))
	  return NULL;
	
	/* Initialisation */
	u[0] = malloc(n*sizeof(double));
	for (i=0 ; i<n ; i++){
		u[0][i] = 0.0;
	}
	
	for (m =1 ; m < cycle+1 ; m++){
		
		/* RESTRICTION PART */
		for (k = 0 ; k < level-1 ; k++){
			
			/* Parameters for the fine grid for the considered cycle */
			h_prov = pow(2,k) * h;
			X_prov = (x_f - x_i)/ h_prov +1;
			Y_prov = (y_f - y_i)/h_prov + 1;
			nx_prov = Lx / h_prov - 1; /* number of points in the fine grid along the x direction */
			ny_prov = Ly / h_prov - 1; /* number of points in the fine grid along the y direction */
			n_prov  = nx_prov*ny_prov - X_prov * Y_prov ; /* nombre d'inconnues */
			nnz_prov = 5*nx_prov*ny_prov - 2*nx_prov - 2*ny_prov - 5*X_prov*Y_prov - 2*(X_prov + Y_prov);
			
			/* Parameters for the coarse grid for the considered cycle */
			h_c = 2*h_prov;
			X_c = (x_f - x_i)/ h_c +1;
			Y_c = (y_f - y_i)/h_c + 1;
			nx_c = Lx / h_c - 1; /* number of points in the coarse grid along the x direction */
			ny_c = Ly / h_c - 1; /* number of points in the coarse grid along the y direction */
			n_c  = nx_c*ny_c - X_c * Y_c ; /* number of unknowns */
			
			if (m==1){
				r_prog[k] = malloc(n_prov * sizeof(double));
				r_c[k] = malloc(n_c * sizeof(double));
				a[k] = malloc(nnz_prov * sizeof(double));
				ja[k] = malloc(nnz_prov * sizeof(int));
				ia[k] = malloc((n_prov+1) * sizeof(int));
				
				if(k>0){u[k] = malloc(n_prov * sizeof(double));}
			}
			
			if (prob_grossier(h_prov, &ia[k], &ja[k], &a[k], Lx, Ly, x_i, x_f, y_i, y_f, u_i, u_e))
			  return NULL;
			
			if (k==0){
				r_prog[k] = gs_prog(n_prov, ia[k], ja[k], a[k], b, u[k]);
			}
			
			else{
				for (i=0; i<n_prov; i++){
					u[k][i] = 0.0;
				}
				
				r_prog[k] = gs_prog(n_prov, ia[k], ja[k], a[k], r_c[k-1], u[k]);
			}
			
			r_c[k] = restriction(h_prov, Lx, Ly, x_i, x_f, y_i, y_f, r_prog[k]);
			
			if (k == level - 2){
				if( solve_umfpack(n_c, ia_coarsest, ja_coarsest, a_coarsest, r_c[k], u_coarsest) ) // Solving with UMFPACK
					return NULL;
			}
		}
		
		/* PROLONGATION PART */
		for (k = level - 2 ; k > -1 ; k--){
			
			/* Parameters for the fine grid for the considered cycle */
			h_prov = pow(2,k) * h;
			X_prov = (x_f - x_i)/ h_prov +1;
			Y_prov = (y_f - y_i)/h_prov + 1;
			nx_prov = Lx / h_prov - 1; /* number of points in the fine grid along the x direction */
			ny_prov = Ly / h_prov - 1; /* number of points in the fine grid along the y direction */
			n_prov  = nx_prov*ny_prov - X_prov * Y_prov ; /* number of unknowns */
			nnz_prov = 5*nx_prov*ny_prov - 2*nx_prov - 2*ny_prov - 5*X_prov*Y_prov - 2*(X_prov + Y_prov);
			
			if(m==1){
				P_c[k] = malloc(n_prov * sizeof(double));
				u_corr[k] = malloc(n_prov * sizeof(double));
			}
			
			/* Prolongation */
			
			if (k==level-2){P_c[k] = prolongation(h_prov, n_prov, Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i, u_coarsest); }
			
			else{P_c[k] = prolongation(h_prov, n_prov, Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i, u[k+1]);}
			
			/* Correction */
			for (i=0; i<n_prov; i++){
				u_corr[k][i] = u[k][i] + P_c[k][i];
			}
			
			/* Post smoothing */
			if(k==0){u[k] = gs_reg(n_prov, ia[k], ja[k], a[k], b, u_corr[k]);}
			
			else{u[k] = gs_reg(n_prov, ia[k], ja[k], a[k], r_c[k-1], u_corr[k]);}
		}
		
		norme_residu = norme(n, ia[0], ja[0], a[0], b, u[0]);
		double log_norme_residu = log10(norme_residu);
		
		/* Writing data in the file data3 */
		fprintf (data3, "%i %.16e\n",m, log_norme_residu);
	}
	fclose(data3);
	
	/* Creation of the command file */
	FILE *cmd = fopen ("cmd3.txt","w+");
	fprintf(cmd,
            "set title 'Méthode %s-grilles pour un pas de discrétisation h=%s' \n" \
            "set xrange [0:%s]\n"
            "set tics out\n"
            "set xlabel 'Nombre de cycles'\n"
            "set ylabel 'norme du résidu en échelle logorithmique'\n"
            "unset key\n"
            "set grid\n"
            "plot 'data3.txt' linetype 7 linecolor 2\n",
            niveau, pas, nombre_cycle);
	fclose (cmd);
	
	/* execute commands in cmd file */
	
	//system("gnuplot -persistent cmd3.txt") ;
	
	free(r_prog); free(r_c); free(ia); free(ja); free(a); free(P_c); free(u_corr); free(a_coarsest); free(ja_coarsest); free(ia_coarsest); free(u_coarsest);
	
	return u[0];
}

