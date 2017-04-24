/*
 * This program simulates a Sedov explosion in a cube with side 256m at 300K and 1atm
 * The blast liberates 10**10J from a cube of 2m (4?) in the center of the cube
 * Spatial resolution is 2m (1?) in each direction
 * Method used is finite volumes to solve Euler equations
 * This code STOPS when the shock wave is of 120m of radius
 * Prints radial density at r= 10,60,120m
 * We will assume rho_ini = 1kg/m3
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <time.h>

//------------------------------------------------------------------------------
//   DEFINITIONS & METHODS
//------------------------------------------------------------------------------
/*
 * FLOAT: data precision
 * atm  : atmospheric pressure in Pa
 * cv   : Constant volume heat capacity of the gas (divided by constant R)
 * gamma: Adiabatic constant of the gas
 */

#define man "./sedov.x "
#define FLOAT float
#define atm 101325
#define cv (3.0/2.0)
#define gamma (5.0/3.0)
void allocateAll();
int indexx(int i, int j, int k);
void __init__();
FLOAT getSoundSpeed(FLOAT vx, FLOAT vy, FLOAT vz, FLOAT p, FLOAT rho);

//------------------------------------------------------------------------------
//   GLOBAL VARIABLES
//------------------------------------------------------------------------------

// Box size in m 
int L = 256;
// Spatial resolution in m
int resol = 2;
// Size of the tridimentional arrays
int N;
// Blast energy liberated in joules/kg (from little square at the center)
FLOAT blast;
// Gas temperature in Kelvin
FLOAT T;
// Gas Pressure in Pascals
FLOAT P;
// Gas initial density
FLOAT rho_ini;
// Velocity v
// Pressure p, internal energy e per unit mass, density rho, total energy E per unit volume
FLOAT e;
FLOAT** v;
FLOAT* p;
FLOAT* E;
FLOAT* rho;
/*
 * AUXILIAR VECTORS THAT FOLLOW WAVE EQ Ut + (d/dxi)Fi = 0
 *      | rho    |           |  rho*vi                   |                                        
 *      | rho*vx |           |  rho*vi*vx + delta_i,0*p  |       
 * U =  | rho*vy |      Fi = |  rho*vi*vy + delta_i,1*p  |       
 *      | rho*vz |           |  rho*vi*vz + delta_i,2*p  |       
 *      | rho*E  |           |  (rho*E+p)*vi             |      
 *
 */             
FLOAT** U;
FLOAT*** F;

//------------------------------------------------------------------------------
//   MAIN
//------------------------------------------------------------------------------


int main(int argc, char** argv){
	// Temperature is 300K
	T = 300.0;
	// Initial gas density is 1kg/m3
	rho_ini = 1.0;
	// Pressure is atmospheric
	P = 1.0*atm;
	// Blast is defined
	blast = pow(10.0,10.0);
	// Internal energy depends on pressure and volume: e = cv*n*T/m = P*(V/m)*(cv/R) = P*cv/rho
	e = P*pow(1.0*L,3.0)*cv/rho_ini;
	// N is defined 
	N = ((FLOAT)L)/((FLOAT)resol);
	// Variables are allocated
	allocateAll();
	// Variables are initialised
	__init__();
	
	return 0;
}






//------------------------------------------------------------------------------
//  AUXILIAR METHODS
//------------------------------------------------------------------------------


/*
 * Allocates cubic nxnxn matrix into unidimentional FLOAT array arr
 *
 */
void allocateAll(){
	v = malloc(3*sizeof(FLOAT*));
	U = malloc(5*sizeof(FLOAT*));
	F = malloc(5*sizeof(FLOAT**));
	int i,j;
	for( i = 0; i < 3; i++ ){
		v[i] = malloc(N*N*N*sizeof(FLOAT));
		F[i] = malloc(5*sizeof(FLOAT*));
		for( j = 0; j < 5; j++ ){
			F[i][j] = malloc(N*N*N*sizeof(FLOAT));
		}
	}
	for( j = 0; j < 5; j++ ){
		U[j] = malloc(N*N*N*sizeof(FLOAT));
	}	
	p = malloc(N*N*N*sizeof(FLOAT));
	rho = malloc(N*N*N*sizeof(FLOAT));
	E = malloc(N*N*N*sizeof(FLOAT));
}

/*
 * Gets unidimentional index given position i,j,k
 * i,j,k < N
 *
 */
int indexx(int i, int j, int k){
	return (int) N*N*i + N*j + k;
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
	for( k = 0; k < N*N*N; k++ ){
		p[k] = P;
		v[0][k] = 0.0;
		v[1][k] = 0.0;
		v[2][k] = 0.0;
		rho[k]  = rho_ini;
		// Initial velocity is 0, so Total energy E is only internal
		// E = rho(e+(1/2)*|v^2|)
		E[k]    = rho_ini*e;
	}

	// Sedovs blast energy is corrected for the blast (4m side)
	E[indexx(N/2  ,N/2  ,N/2  )] = rho_ini*blast;
	E[indexx(N/2-1,N/2  ,N/2  )] = rho_ini*blast;
	E[indexx(N/2  ,N/2-1,N/2  )] = rho_ini*blast;
	E[indexx(N/2  ,N/2  ,N/2-1)] = rho_ini*blast;
	E[indexx(N/2-1,N/2-1,N/2  )] = rho_ini*blast;
	E[indexx(N/2  ,N/2-1,N/2-1)] = rho_ini*blast;
	E[indexx(N/2-1,N/2  ,N/2-1)] = rho_ini*blast;
	E[indexx(N/2-1,N/2-1,N/2-1)] = rho_ini*blast;

	for( k = 0; k < N*N*N; k++ ){
		// Defining F
		for( j = 0; j < 3; j++ ){
			F[j][0][k] = rho[k]*v[j][k];
			for( i = 1; i < 4; i++ ){
				F[j][i][k] = rho[k]*v[j][k]*v[i-1][k] + (j==(i-1))*p[k];//delta_j,i-1 = (j==(i-1))
			}
			F[j][4][k] = (rho[k]*E[k]+p[k])*v[j][k];
		}
		// Defining U
		U[0][k] = rho[k];
		for( i = 1; i < 4; i++ ){
				U[i][k] = rho[k]*v[i-1][k];
		}
		U[4][k] = rho[k]*E[k];
	}
}



/*
 * Gets Vectors U,F in terms of the vectors rho,E,p,v
 *
 */
FLOAT getSoundSpeed(FLOAT vx, FLOAT vy, FLOAT vz, FLOAT p, FLOAT rho){
	// How to calculate sound speed?
	return 0;
}













