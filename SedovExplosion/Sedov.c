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
FLOAT P = 
























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







