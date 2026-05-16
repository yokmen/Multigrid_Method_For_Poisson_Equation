#include <stdio.h>
#include "indice.h"

int indice(double h,double Lx, double x_i, double x_f, double y_i, double y_f, int ix, int iy)
/*
   Goal
   ===
   This function calculates, in lexicographical order, the index of the equation according to the position in the grid.
   
   Arguments
   ===========
   (input)   h  - Discretization step (in m) 
   (input)   nx - Number of points in the grid by removing the outer Dirichlet nodes
   (input)  x_i - abscissa of the beginning of Γi (=xl in the statement) (in m)  
   (input)  x_f - abscissa of the end of Γi (=xl in the statement) (in m) 
   (input)  y_i - start ordinate of Γi (=yb in the statement) (in m)
   (input)  y_f - end ordinate of Γi (=yb in the statement) (in m) 
                
   (output) ind - Associated index
*/
{
	int ind;
	int nx = Lx / h -1 ;
	int X = (x_f - x_i)/ h +1;
	int Y = (y_f - y_i)/h + 1;
	
	/* Here the index does not change from the basic formula */
	if (iy < y_i/h -1){
		ind = ix + nx * iy;
	}
	
	/* Here the positioning in the grid is such that there is a hole and we are in the left side */ 
	else if ( (iy >= y_i/h -1) && (iy <= y_f/h -1) && (ix < x_i/h -1) ){
		ind = ix + nx * iy - X * (iy - y_i / h +1 );
    }
    
    /* Here the positioning in the grid is such that there is a hole and we are in the right side */ 
    else if ( (iy >= y_i/h -1) && (iy <= y_f/h -1) && (ix > x_f/h -1) ){
        ind = ix + nx * iy - X * (iy - y_i / h +2);
    }
    
    else{
	    ind = ix + nx * iy - X*Y;
	}
	
	return ind;
}
