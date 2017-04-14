/*
 * This program simulates a Sedov explosion in a cube with side 256m at 300K and 1atm
 * The blast liberates 10**10J from a cube of 2m (4?) in the center of the cube
 * Spatial resolution is 2m (1?) in each direction
 * Method used is finite volumes to solve Euler equations
 * This code STOPS when the shock wave is of 120m of radius
 * Prints radial density at r= 10,60,120m
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <time.h>

//------------------------------------------------------------------------------
//   DEFINITIONS & METHODS
//------------------------------------------------------------------------------
#define man "./sedov.x "
#define FLOAT float
#define atm 101325
void allocate2d(FLOAT* arr);
void index(int i, int j);

//------------------------------------------------------------------------------
//   GLOBAL VARIABLES
//------------------------------------------------------------------------------

// Box size in m 
int NN = 256;
// Spatial resolution in m
int resol = 2;
// Size of the tridimentional arrays
int N;
// Blast energy liberated in joules (from little square at the center)
FLOAT blast = pow(10.0,10.0);
// Gas temperature in Kelvin
FLOAT T = 300.0;
// Gas Pressure in Pascals
FLOAT P = 1.0*atm;
// Velocity v
// Pressure p, internal energy e, density rho
FLOAT** v;
FLOAT* p,e,rho;


//------------------------------------------------------------------------------
//   MAIN
//------------------------------------------------------------------------------


int main(int argc, char** argv){
	
	// N is defined 
	N = ((FLOAT)NN)/((FLOAT)resol);
	// Initialises variables
	v = malloc(3*sizeof(float));
	allocate3d(v[0]);
	allocate3d(v[1]);
	allocate3d(v[2]);
	allocate3d(p);
	allocate3d(rho);
	allocate3d(e);
	return 0;
}






//------------------------------------------------------------------------------
//  AUXILIAR METHODS
//------------------------------------------------------------------------------


/*
 * Allocates cubic nxnxn matrix into unidimentional FLOAT array arr
 *
 */
void allocate3d(FLOAT* arr){
	arr = malloc(N*N*N*sizeof(FLOAT));
}

/*
 * Gets unidimentional index given position i,j,k
 * i,j,k < N
 *
 */
void index(int i, int j, int k){
	return N*N*i + N*j + k;
}


/*
 * Initializes variables according to the corresponding initial conditions
 *
 */
void __init__(){
	int i;
	int j;
	int k;
	// Fills variables as if it was a normal ideal gas
	for( i = 0; i < N; i++ ){
		for( j = 0; j < N; j++ ){
			for( k = 0; k < N; k++ ){
				p[index(i,j,k)] = P;
				v[0][index(i,j,k)] = 0.0;
				v[1][index(i,j,k)] = 0.0;
`				v[2][index(i,j,k)] = 0.0;
				e[index(i,j,k)] = stateEq(v[0][index(i,j,k)],v[1][index(i,j,k)],v[2][index(i,j,k)],p[index(i,j,k)]);
				rho[index(i,j,k)] = getRho(v[0][index(i,j,k)],v[1][index(i,j,k)],v[2][index(i,j,k)],p[index(i,j,k)]);
			}
		}
	}
	// Sedovs blast energy is corrected (4m side)
	e[index(N/2,j,k)]     = blast;
	e[index(N/2 - 1,j,k)] = blast;
	e[index(i,N/2,k)]     = blast;
	e[index(i,N/2 - 1,k)] = blast;
	e[index(i,j,N/2)]     = blast;
	e[index(i,j,N/2 - 1)] = blast;
}


/*
 * Gets energy in terms of other variables (pressure p, velocity v, density rho)
 *
 */
FLOAT stateEq(FLOAT vx, FLOAT vy, FLOAT vz, FLOAT p, FLOAT rho){
	// Which is the state equation?
	return 0;
}

/*
 * Gets rho density in terms of other variables (pressure p, velocity v, density rho)
 *
 */
FLOAT getRho(FLOAT vx, FLOAT vy, FLOAT vz, FLOAT p, FLOAT rho){
	// How to calculate rho?
	return 0;
}

/*
 * Gets sound speed in terms of other variables (pressure p, velocity v, density rho)
 *
 */
FLOAT getSoundSpeed(FLOAT vx, FLOAT vy, FLOAT vz, FLOAT p, FLOAT rho){
	// How to calculate sound speed?
	return 0;
}













