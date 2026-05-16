#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "umfpk.h"
#include "twogrid.h"
#include "gs_prog.h"
#include "twogrid.h"
#include "gs_reg.h"
#include "restriction.h"
#include "prolongation.h"
#include "prob_grossier.h"
#include "norme.h"
#include "prob.h"

double * twogrid(int cycle, double h, int *ia, int *ja, double *a, double *b, double Lx, double Ly, double x_i, double x_f, double y_i, double y_f, double u_e, double u_i)
/*           
  Goal
  ===
  Calcul du vecteur solution calculee par l'algorithme deux-grilles et afficher la norme du
  residu en fonction du nombre de cycle

  Arguments
  =========
  (input) cycle - Number of cycles of the two-grid algorithm
  (input)  Lx  - length of the outside capacitor [m]
  (input)  Ly  - height of the outside capacitor [m]
  (input)  h  - discretization step [m]
  (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
  (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
  (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
  (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m] 
  (input)  u_i - electric potential on Γi [V]
  (input)  u_e - electric potential on Γe [V]
                
  (output) u - Vector solution calculated by the two-grid algorithm
*/
{
	/* Creation of the various data for the creation of the graphs  */
	
	/* create data file  */
	FILE *data = fopen("data.txt", "w+");
	fprintf (data, "%i %.16e\n",0, 0.0);
	
	int mx, my, nx, ny, n, n_c, nnz_c, mx_c, my_c, nx_c, ny_c;
	
	/* Parameters for the fine grid */
	int X = (x_f - x_i)/ h +1;
	int Y = (y_f - y_i)/h + 1;
	mx = Lx / h + 1;
	my = Ly / h + 1;
	nx = mx - 2;
	ny = my - 2;
	n  = nx*ny - X * Y ; /* number of unknowns */
	
	char pas[10];
	char nombre_cycle[10];
	sprintf(nombre_cycle,"%i",cycle+1);
	sprintf(pas,"%f",h);
	
	/* Parameters for coarse grid */
	double h_c = 2*h;
	int X_c = (x_f-x_i)/h_c+1;
	int Y_c = (y_f-y_i)/h_c+1;
	mx_c = Lx / h_c + 1;
	my_c = Ly / h_c + 1;
	nx_c = mx_c - 2;
	ny_c = my_c - 2;
	n_c = nx_c*ny_c - X_c * Y_c ; /* number of unknowns */
	nnz_c = 5*nx_c*ny_c - 2*nx_c - 2*ny_c - 5*X_c*Y_c - 2*(X_c + Y_c); /* number of nonzero elements */
	
	int i,m;
	
	/* Array initialization */
	double norme_residu;
	double *u = malloc(n * sizeof(double));
	double *r_prog = malloc(n * sizeof(double));
	double *r_c = malloc(n_c * sizeof(double));
	double *P_c = malloc(n * sizeof(double));
	double *u_grossier = malloc(n_c * sizeof(double));
	
	int *ia_grossier = malloc((n_c + 1) * sizeof(int));
	int *ja_grossier = malloc(nnz_c * sizeof(int));
	double *a_grossier  = malloc(nnz_c * sizeof(double));
	double *u_corr = malloc(n * sizeof(double));
	
	/* Coarse problem initialization */
	if (prob_grossier(h_c, &ia_grossier, &ja_grossier, &a_grossier, Lx, Ly, x_i, x_f, y_i, y_f, u_i, u_e))
	  return NULL;
	
	/* Initialization of the solution vector to 0 */
	for (i=0; i < n; i++){
		u[i] = 0.0;
	}
	
	/* Two-grid algorithm */
	for (m=1; m<cycle+1; m++){
		/* Pre-smoothing of the problem */
		r_prog = gs_prog(n, ia, ja, a, b, u);
		
		/* Restriction of the problem */
		r_c = restriction(h, Lx, Ly, x_i, x_f, y_i, y_f, r_prog);
		
		if( solve_umfpack(n_c, ia_grossier, ja_grossier, a_grossier, r_c, u_grossier) ) // Solving with UMFPACK
		  return NULL;
		
		/* Prolongation */
		P_c = prolongation(h, n, Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i, u_grossier);
		
		/* Correction of the solution*/
		for (i=0; i<n; i++){
			u_corr[i] = u[i] + P_c[i];
		}
		
		/* Post-smoothing of the problem */
		u = gs_reg(n, ia, ja, a, b, u_corr);
		
		/* Calculation of the residual norm */
		norme_residu = norme(n, ia, ja, a, b, u);
		double log_norme_residu = log10(norme_residu);
		
		/* Writing data to the file data */
		fprintf(data, "%i %.16e\n",m, log_norme_residu);
	}
	
	fclose(data);
	
	/* Creation of the command file */	
	FILE *cmd = fopen ("cmd.txt","w+");
	fprintf(cmd,
            "set title 'Méthode deux-grilles pour un pas de discrétisation h = %s' \n" \
            "set xrange [0:%s]\n"
            "set tics out\n"
            "set xlabel 'Nombre de cycles'\n"
            "set ylabel 'norme du résidu en échelle logorithmique'\n"
            "unset key\n"
            "set grid\n"
            "plot 'data.txt' linetype 7 linecolor 2\n",
            pas, nombre_cycle);
	fclose (cmd);
	
	/* execute commands in cmd file */
	system("gnuplot -persistent cmd.txt") ;
	
	free(r_prog); free(r_c); free(P_c); free(u_grossier); free(ia_grossier); free(ja_grossier);free(a_grossier); free(u_corr);
	
	return u;
}
