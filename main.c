#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "norme.h"
#include "twogrid.h"
#include "plot_potentiel.h"
#include "multigrid.h"
#include "conjugate.h"

/* Fonction main */

int main(int argc, char *argv[])
{
  /* Variable declaration */
  int cycle = 50;
  int level = 2;
  double h = pow(2,-level)/2000;
  double Lx = 0.005;
  double Ly = 0.004;
  double x_i = 0.003;
  double x_f = 0.004;
  double y_i = 0.001;
  double y_f = 0.003;
  double u_i = -5.0;
  double u_e = 5.0;
  int n, *ia, *ja; 
  double *a, *b;
  double tc1, tc2, tw1, tw2;
  
  /* Question 0*/ 
  /* Generating the probel in CSR format */
  
  printf("Question 0\n");
  if (prob(h, &n, &ia, &ja, &a, &b, Lx, Ly, x_i, x_f, y_i, y_f, u_i, u_e ))
     return 1;
  printf("\nProblème : ");
  printf("h=%f, n = %8d  nnz = %9d\n\n\n",h,n,ia[n]);
  
  /* Question 1 : 2-grid algorithm  */
  double *x_twogrid = malloc(n * sizeof(double)); 
  x_twogrid = twogrid(cycle, h, ia, ja, a, b, Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i);
  plot_potentiel(Lx, Ly, h, x_i, x_f, y_i, y_f, u_e, u_i, x_twogrid); 
  free(x_twogrid);
  
  /* Question 2 : multogrid algorithm */
  //double *x_multigrid = malloc(n * sizeof(double)); 
  //x_multigrid = multigrid(cycle, level, h, b, Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i);
  //plot_potentiel(Lx, Ly, h, x_i, x_f, y_i, y_f, u_e, u_i, x_multigrid); 
  //free(x_multigrid);
  
  /* Question 3 : Conjugate gradient algorithm with multigrid preconditioner */
  
  /* Note: Comment out line of code 209 of the multigrid function : system("gnuplot -persistent cmd3.txt") */
  
  //double *x_conjugate = malloc(n * sizeof(double));
  //x_conjugate = conjugate(cycle, level, h, ia, ja, a, b, Lx, Ly, x_i, x_f, y_i, y_f, u_e, u_i);
  //plot_potentiel(Lx, Ly, h, x_i, x_f, y_i, y_f, u_e, u_i, x_conjugate); 
  //free(x_conjugate);
  
  /* Solution time measurement */
  tc1 = mytimer_cpu(); tw1 = mytimer_wall();
  tc2 = mytimer_cpu(); tw2 = mytimer_wall();
  printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1);
  printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1);

  /* Freeing memory */
  free(ia); free(ja); free(a); free(b);
  
  
  return 0;
}
