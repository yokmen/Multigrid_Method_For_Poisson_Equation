#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "restriction.h"
#include "multigrid.h"
#include "norme.h"
#include "conjugate.h"

double *conjugate(int cycle, int level, double h, int *ia, int *ja, double *a, double *b, double Lx, double Ly, double x_i, double x_f, double y_i, double y_f, double u_e, double u_i)
/*
   Goal
   ===
   Computation of the solution by the conjugate gradient method with two-grid preconditioner
   
   Arguments
   ===========
   (input)  cycle  - number of cycle of the conjugate gradient method
   (input)  h  - discretization step [m]
   (input)  ia - pointer of the vector 'ia' of the matrix A (fine grid)
   (input)  ja - pointer of the vector 'ja' of the matrix A (fine grid)
   (input)  a  - pointer of the vector 'a' of the matrix A (fine grid)
   (input)  b  - pointer of the vector 'b' (fine grid)
   (input)  Lx  - length of the outside capacitor [m]
   (input)  Ly  - height of the outside capacitor [m]
   (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
   (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
   (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
   (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m] 
   (input) u_i - electrostatic potential on Γi [V]
   (input) u_e - electrostatic potential on Γe [V] 
                
   (output) u - Vector solution 
*/
{
	/* create data file */
	FILE *data4 = fopen("data4.txt", "w+");
	fprintf (data4, "%i %.16f\n",0, 0.0);
	
	int  m,i,j, mx, my, nx, ny, n;
	double norme_residu;
	
	/* Parameters for the fine grid */
	int X = (x_f - x_i)/ h +1;
	int Y = (y_f - y_i)/h + 1;
	mx = Lx / h + 1; /* number of points in the grid in the x-direction */
	my = Ly / h + 1; /* number of points in the grid in the y-direction */
	nx = mx - 2; 
	ny = my - 2; 
	n  = nx*ny - X * Y ; /* number of unknowns */
	
	/* Creation of the various data for the creation of the graphs */
	char pas[10];
	char nombre_cycle[10];
	sprintf(nombre_cycle,"%i",cycle+1);
	sprintf(pas,"%f",h);
	
	double *num_beta = malloc(cycle * sizeof(double));
	double *den_beta = malloc(cycle * sizeof(double));
	double *beta = malloc(cycle * sizeof(double));
	double *alpha = malloc(cycle * sizeof(double));
	double *num_alpha = malloc(cycle * sizeof(double));
	double *den_alpha = malloc(cycle * sizeof(double));
	
	double ** u = malloc((cycle+1) * sizeof(double *));
	double ** r = malloc((cycle+1) * sizeof(double *));
	double ** v = malloc(cycle * sizeof(double *));
	double ** d = malloc(cycle * sizeof(double *));
	double ** f = malloc(cycle * sizeof(double *));
	
	/* Initialisation */
	u[0] = malloc(n * sizeof(double));
	r[0] = malloc(n * sizeof(double));
	for (i=0; i<n ; i++){
		u[0][i] = 0.0;
		r[0][i] = b[i];
	}
	
	for (m = 0 ; m < cycle ; m++){
		
		v[m] = malloc(n * sizeof(double));
		d[m] = malloc(n * sizeof(double));
		f[m] = malloc(n * sizeof(double));
		u[m+1] = malloc(n * sizeof(double));
		r[m+1] = malloc(n * sizeof(double));
		
		/* Computation of the vector v = MG(0,r[m]) */
		v[m] = multigrid(1, level, h, r[m], Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i);
		
		/* Computation of beta*/
		num_beta[m] = 0.0;
		den_beta[m] = 0.0;
		
		if (m == 0){
			beta[m] = 0.0;
		}
		
		else{
			for (i=0; i<n ; i++){
				num_beta[m] += r[m][i] * v[m][i];
				den_beta[m] += r[m-1][i] * v[m-1][i];
			}
			beta[m] = num_beta[m]/den_beta[m];
		}
		
		/* Computation of d */
		if (m == 0){
			for (i=0 ; i<n ; i++){
				d[m][i] = v[m][i];
			}
		}
		
		else{
			for (i=0 ; i<n ; i++){
				d[m][i] = v[m][i] + beta[m] * d[m-1][i];
			}
		}
		
		/* Computation of alpha */
		for (i=0 ; i < n ; i++){
			for (j=ia[i] ; j < ia[i+1] ; j++){
				f[m][i] += a[j] * d[m][ja[j]];
			}
		}
		num_alpha[m] = 0.0;
		den_alpha[m] = 0.0;
		for (i=0; i<n ; i++){
			num_alpha[m] += r[m][i] * v[m][i];
			den_alpha[m] += d[m][i]*f[m][i];
		}
		
		alpha[m] = num_alpha[m] / den_alpha[m];
		
		/* Computation of the solution vector */
		for (i = 0; i < n; i ++){
			u[m+1][i] = u[m][i] + alpha[m]*d[m][i];
		}
		
		/* Computation of the new residual vector */
		for (i = 0; i < n; i ++){
			r[m+1][i] = r[m][i] - alpha[m]*f[m][i];
		}
		
		/* Computation of the residual norm */
		norme_residu = norme(n, ia, ja, a, b, u[m+1]);
		double log_norme_residu = log10(norme_residu);
		
		/* Writing data in the file data4 */
		fprintf (data4, "%i %.16e\n",m+1, log_norme_residu);
	}
	fclose(data4);
	
	/* Creation of the command file */	
	FILE *cmd = fopen ("cmd4.txt","w+");
	fprintf(cmd,
            "set title 'Méthode du gradient conjugué pour un pas de discrétisation h = %s' \n" \
            "set xrange [0:%s]\n"
            "set tics out\n"
            "set xlabel 'Nombre de cycles'\n"
            "set ylabel 'norme du résidu en échelle logorithmique'\n"
            "unset key\n"
            "set grid\n"
            "plot 'data4.txt' linetype 7 linecolor 2\n",
            pas, nombre_cycle);
	fclose (cmd);
	
	/* execute commands in cmd file */
	system("gnuplot -persistent cmd4.txt") ;
	
	free(num_beta); free(den_beta); free(beta); free(alpha); free(num_alpha); free(den_alpha); free(r); free(v); free(d); free(f); 
	
	return u[cycle];
}
