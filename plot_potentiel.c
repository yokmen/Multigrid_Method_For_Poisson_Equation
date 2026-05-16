#include <stdio.h>
#include <stdlib.h>
#include "plot_potentiel.h"
#include "indice.h"

void plot_potentiel(double Lx, double Ly, double h, double x_i, double x_f, double y_i, double y_f, double u_e, double u_i, double *V)
/*
  Goal
  ===
  This function is used to visualize the electrical potential.
            
  Arguments
  =========                
   (input)  Lx  - length of the outside capacitor [m]
   (input)  Ly  - height of the outside capacitor [m]
   (input)  h  - discretization step [m]
   (input)  x_i - abscissa of the beginning of Γi (= xl in the statement) [m]
   (input)  x_f - abscissa of the end of Γi (=xr in the statement) [m]
   (input)  y_i - ordinate of the beginning Γi (=yb in the statement) [m]
   (input)  y_f - ordinate of the end of Γi (=yt in the statement) [m] 
   (input)  u_i - potentiel électrique sur Γi [V]
   (input)  u_e - potentiel électrique sur Γe [V]     
   (input)  V  - Electric potential [V] 
*/
{
	/* Creation of the various data for the creation of the graphs */
	char longueur[10];
	char largeur[10];
	char potentiel_exterieur[10];
	char potentiel_interieur[10];
	char nombre_points[10];
	
	int m=Lx/h+1;
	
	sprintf(longueur,"%f",Lx);
	sprintf(largeur,"%f",Ly);
	sprintf(potentiel_exterieur,"%f",u_e);
	sprintf(potentiel_interieur,"%f",u_i);
	sprintf(nombre_points,"%i",m);
	
	/* create data file */
	FILE *data2 = fopen("data2.txt", "w+");
    int ind;
	double x,y,u;
	for (y=0.0; y <= Ly; y += h){
		for (x=0.0; x < Lx+h; x += h){
			/* Conditions at the outer boundary */
			if (x==0.0 || x==Lx || y==0.0 || y==Ly ){u = u_e;}
			
			/* Conditions at the inner boundary */
			else if ((x >=x_i) && (x<=x_f+h) && (y>=y_i) && (y<=y_f+h)){u=u_i;}
		   
			/* Otherwise, the potential is equal to that calculated by the main function */
			else {
				ind = indice(h, Lx, x_i, x_f, y_i, y_f, (x/h-1), (y/h-1));
				u = V[ind];
			}
			
			/* Writing data to the file data2 */
			fprintf (data2, "%f %f %.16e\n",x,y,u); 
		}
		fprintf(data2, "\n");
	}
	
	fclose(data2);
	
	/* Creation of the command file */	
	FILE *cmd = fopen ("cmd2.txt","w+");
	fprintf(cmd,
            "set title 'Potentiel électrique [V] avec m=%s'\n" \
            "set xrange [0:%s]\n"
            "set yrange [0:%s]\n"
            "set cbrange [%s:%s]\n"
            "set tics out\n"
            "set xlabel 'Coordonées x [m]'\n"
            "set ylabel 'Coordonées y [m]'\n"
            "set tics nomirror \n"
            "plot 'data2.txt' with image notitle\n",
            nombre_points, longueur, largeur, potentiel_interieur, potentiel_exterieur);
	fclose (cmd);
	
	/* execute commands in cmd file */
	system("gnuplot -persistent cmd2.txt") ;
}
