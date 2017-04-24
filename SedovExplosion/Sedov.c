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
#define bool int
#define true 1
#define false 0
void allocateAll();
int indexx(int i, int j, int k);
void __init__();
FLOAT getDT();
void allocate3d(FLOAT*** vec);
void evolveU(float dt,float dx);
void updateF();

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
FLOAT**** v;
FLOAT*** p;
FLOAT*** E;
FLOAT*** rho;
/*
 * AUXILIAR VECTORS THAT FOLLOW WAVE EQ Ut + (d/dxi)Fi = 0
 *      | rho    |           |  rho*vi                   |                                        
 *      | rho*vx |           |  rho*vi*vx + delta_i,0*p  |       
 * U =  | rho*vy |      Fi = |  rho*vi*vy + delta_i,1*p  |       
 *      | rho*vz |           |  rho*vi*vz + delta_i,2*p  |       
 *      | rho*E  |           |  (rho*E+p)*vi             |      
 *
 */             
FLOAT**** U;
FLOAT***** F;

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
	   //
	  //  TIME EVOLUTION
	 //
	// Time variable t, deltas in position and time: dx,dt
	FLOAT time = 0;
	FLOAT dx = resol;
	FLOAT dt = getDT();
	bool go_on = true;
	do{	
		// Evolves vector U
		evolveU(dt,dx);
		// Updates Fs in terms of U
		updateF();
		

		time+=dt;
	} while(go_on);
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
	v = malloc(3*sizeof(FLOAT***));
	U = malloc(5*sizeof(FLOAT***));
	F = malloc(5*sizeof(FLOAT****));
	int i,j;
	for( i = 0; i < 3; i++ ){
		allocate3d(v[i]);
		F[i] = malloc(5*sizeof(FLOAT***));
		for( j = 0; j < 5; j++ ){
			allocate3d(F[i][j]);
		}
	}
	for( j = 0; j < 5; j++ ){
		allocate3d(U[j]);
	}	
	allocate3d(p);
	allocate3d(rho);
	allocate3d(E);
}

/*
 * Initializes variables according to the corresponding initial conditions
 *
 */
void __init__(){

	int i;
	int j;
	int k;
	int l;
	int m;
	int n;
	// Fills variables as if it was a normal ideal gas
	for( i = 0; i < N; i++){
		for( j = 0; j < N; j++){
			for( k = 0; k < N; k++){
				p[i][j][k] = P;
				v[0][i][j][k] = 0.0;
				v[1][i][j][k] = 0.0;
				v[2][i][j][k] = 0.0;
				rho[i][j][k]  = rho_ini;
				// Initial velocity is 0, so Total energy E is only internal
				// E = rho(e+(1/2)*|v^2|)
				E[i][j][k]    = rho_ini*e;
			}
		}
	}

	// Sedovs blast energy is corrected for the blast (4m side)
	E[N/2  ][N/2  ][N/2  ] = rho_ini*blast;
	E[N/2-1][N/2  ][N/2  ] = rho_ini*blast;
	E[N/2  ][N/2-1][N/2  ] = rho_ini*blast;
	E[N/2  ][N/2  ][N/2-1] = rho_ini*blast;
	E[N/2  ][N/2-1][N/2-1] = rho_ini*blast;
	E[N/2-1][N/2  ][N/2-1] = rho_ini*blast;
	E[N/2-1][N/2-1][N/2  ] = rho_ini*blast;
	E[N/2-1][N/2-1][N/2-1] = rho_ini*blast;
	for( l = 0; l < N; l++){
		for( m = 0; m < N; m++){
			for( n = 0; n < N; n++){
				// Defining F
				for( j = 0; j < 3; j++ ){
					F[j][0][l][m][n] = rho[l][m][n]*v[j][l][m][n];
					for( i = 1; i < 4; i++ ){
						F[j][i][l][m][n] = rho[l][m][n]*v[j][l][m][n]*v[i-1][l][m][n] + (j==(i-1))*p[l][m][n];
						//delta_j,i-1 = (j==(i-1))
					}
					F[j][4][l][m][n] = (rho[l][m][n]*E[l][m][n]+p[l][m][n])*v[j][l][m][n];
				}
				// Defining U
				U[0][l][m][n] = rho[l][m][n];
				for( i = 1; i < 4; i++ ){
						U[i][l][m][n] = rho[l][m][n]*v[i-1][l][m][n];
				}
				U[4][l][m][n] = rho[l][m][n]*E[l][m][n];
			}
		}
	}
}

/*
 * Evolves vector U given dx and dt according to the finite volume method
 * Finite volume is performed over a volume-centered cubic uniform grid
 * Volume integrals are approximated given variables constant over finite volumes
 * Surface integrals are approximated given variables constant over square surfaces
 * Values of variables on surfaces are calculated as the mean of the adjacent volumes
 * Boundary conditions are taken as free boundary :
 * (directional derivative perpendicular to the boundary surface vanishes) (This vanishes the respective surface integral)
 *
 */
void evolveU(float dt,float dx){
	// Creates new array to avoid over-evolving
	int i,j,k,l;
	FLOAT**** U_new = malloc(5*sizeof(FLOAT***));
	for( i = 0; i < 5; i++){
		allocate3d(U_new[i]);
	}
	for( l = 0; l < 5; l++){
		for( i = 1; i < N-1; i++){
			for( j = 1; j < N-1; j++){
				// CENTRAL VOLUMES		
				for( k = 1; k < N-1; k++){
					U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
					-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);
				}
				// BOUNDARY CONDITIONS (EDGES of the cube)
				//////////////////////////////////////////////// DOWN
				k = 0;			
				U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
				-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]+(F[2][l][i][j][k]));//6
				////////////////////////////////////////////////// UP				
				k = N-1;
				U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]-(F[2][l][i][j][k])
				-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);//3
				////////////////////////////////////////////////// SOUTH (INDICES ARE PERMUTATED)				
				k = 0;
				U_new[l][j][k][i]=U[l][j][k][i]-(0.5*dt/dx)*(F[0][l][j+1][k][i]+F[1][l][j][k+1][i]+F[2][l][j][k][i+1]
				-F[0][l][j-1][k][i]+(F[1][l][j][k][i])-F[2][l][j][k][i-1]);//5
				////////////////////////////////////////////////// NORTH			
				k = N-1;
				U_new[l][j][k][i]=U[l][j][k][i]-(0.5*dt/dx)*(F[0][l][j+1][k][i]-(F[1][l][j][k][i])+F[2][l][j][k][i+1]
				-F[0][l][j-1][k][i]-F[1][l][j][k-1][i]-F[2][l][j][k][i-1]);//2				
				////////////////////////////////////////////////// WEST			
				k = 0;
				U_new[l][k][i][j]=U[l][k][i][j]-(0.5*dt/dx)*(F[0][l][k+1][i][j]+F[1][l][k][i+1][j]+F[2][l][k][i][j+1]
				+(F[0][l][k][i][j])-F[1][l][k][i-1][j]-F[2][l][k][i][j-1]);//4
				////////////////////////////////////////////////// EAST			
				k = N-1;
				U_new[l][k][i][j]=U[l][k][i][j]-(0.5*dt/dx)*(-(F[0][l][k][i][j])+F[1][l][k][i+1][j]+F[2][l][k][i][j+1]
				-F[0][l][k][i][j]-F[1][l][k][i-1][j]-F[2][l][k][i][j-1]);//1
			}
			// BOUNDARY CONDITIONS (EDGES of the cube)
			//////////////////////////////////////////////// DOWN
			k = 0;
			j = 0;//SOUTH	
			U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
			-F[0][l][i-1][j][k]+(F[1][l][i][j][k])+(F[2][l][i][j][k]));//6,5
			j = N-1;//NORTH
			U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]-(F[1][l][i][j][k])+F[2][l][i][j][k+1]
			-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]+(F[2][l][i][j][k]));//6,2
			////////////////////////////////////////////////// UP				
			k = N-1;
			j = 0;// SOUTH
			U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]-(F[2][l][i][j][k])
			-F[0][l][i-1][j][k]+(F[1][l][i][j][k])-F[2][l][i][j][k-1]);//3,5
			j = N-1;// NORTH
			U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]-(F[1][l][i][j][k])-(F[2][l][i][j][k])
			-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);//3,2
			////////////////////////////////////////////////// SOUTH (INDICES ARE PERMUTATED)				
			k = 0;
			j = 0;// WEST
			U_new[l][j][k][i]=U[l][j][k][i]-(0.5*dt/dx)*(F[0][l][j+1][k][i]+F[1][l][j][k+1][i]+F[2][l][j][k][i+1]
			+(F[0][l][j][k][i])+(F[1][l][j][k][i])-F[2][l][j][k][i-1]);//5,4
			j = N-1;// EAST
			U_new[l][j][k][i]=U[l][j][k][i]-(0.5*dt/dx)*(-(F[0][l][j][k][i])+F[1][l][j][k+1][i]+F[2][l][j][k][i+1]
			-F[0][l][j-1][k][i]+(F[1][l][j][k][i])-F[2][l][j][k][i-1]);//5,1
			////////////////////////////////////////////////// NORTH			
			k = N-1;
			j = 0;// WEST
			U_new[l][j][k][i]=U[l][j][k][i]-(0.5*dt/dx)*(F[0][l][j+1][k][i]-(F[1][l][j][k][i])+F[2][l][j][k][i+1]
			+(F[0][l][j][k][i])-F[1][l][j][k-1][i]-F[2][l][j][k][i-1]);//2,4
			j = N-1;// EAST	
			U_new[l][j][k][i]=U[l][j][k][i]-(0.5*dt/dx)*(-(F[0][l][j][k][i])-(F[1][l][j][k][i])+F[2][l][j][k][i+1]
			-F[0][l][j-1][k][i]-F[1][l][j][k-1][i]-F[2][l][j][k][i-1]);//2,1				
			////////////////////////////////////////////////// WEST			
			k = 0;
			j = 0;// DOWN
			U_new[l][k][i][j]=U[l][k][i][j]-(0.5*dt/dx)*(F[0][l][k+1][i][j]+F[1][l][k][i+1][j]+F[2][l][k][i][j+1]
			+(F[0][l][k][i][j])-F[1][l][k][i-1][j]+(F[2][l][k][i][j]));//4,6
			j = N-1;// UP
			U_new[l][k][i][j]=U[l][k][i][j]-(0.5*dt/dx)*(F[0][l][k+1][i][j]+F[1][l][k][i+1][j]-(F[2][l][k][i][j])
			+(F[0][l][k][i][j])-F[1][l][k][i-1][j]-F[2][l][k][i][j-1]);//4,3
			////////////////////////////////////////////////// EAST			
			k = N-1;
			j = 0;// DOWN
			U_new[l][k][i][j]=U[l][k][i][j]-(0.5*dt/dx)*(-(F[0][l][k][i][j])+F[1][l][k][i+1][j]+F[2][l][k][i][j+1]
			-F[0][l][k][i][j]-F[1][l][k][i-1][j]+(F[2][l][k][i][j]));//1,6
			j = N-1;// UP
			U_new[l][k][i][j]=U[l][k][i][j]-(0.5*dt/dx)*(-(F[0][l][k][i][j])+F[1][l][k][i+1][j]-(F[2][l][k][i][j])
			-F[0][l][k][i][j]-F[1][l][k][i-1][j]-F[2][l][k][i][j-1]);//1,3
		}
	// BOUNDARY CONDITIONS (CORNERS of the cube)
	//////////////////////////////////////////////// DOWN
	k = 0;
	j = 0;/////////////////SOUTH
	i = 0;//WEST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
	+(F[0][l][i][j][k])+(F[1][l][i][j][k])+(F[2][l][i][j][k]));//6,5,4
	i = N-1;//EAST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(-(F[0][l][i][j][k])+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
	-F[0][l][i-1][j][k]+(F[1][l][i][j][k])+(F[2][l][i][j][k]));//6,5,1
	j = N-1;///////////////NORTH
	i = 0  ;//WEST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]-(F[1][l][i][j][k])+F[2][l][i][j][k+1]
	+(F[0][l][i][j][k])-F[1][l][i][j-1][k]+(F[2][l][i][j][k]));//6,2,4
	i = N-1;//EAST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(-(F[0][l][i][j][k])-(F[1][l][i][j][k])+F[2][l][i][j][k+1]
	-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]+(F[2][l][i][j][k]));//6,2,1
	////////////////////////////////////////////////// UP
	k = N-1;
	j = 0;//////////////////SOUTH
	i = 0;//WEST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]-(F[2][l][i][j][k])
	+(F[0][l][i][j][k])+(F[1][l][i][j][k])-F[2][l][i][j][k-1]);//3,5,4
	i = N-1;//EAST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(-(F[0][l][i][j][k])+F[1][l][i][j+1][k]-(F[2][l][i][j][k])
	-F[0][l][i-1][j][k]+(F[1][l][i][j][k])-F[2][l][i][j][k-1]);//3,5,1
	j = N-1;////////////////NORTH
	i = 0  ;//WEST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]-(F[1][l][i][j][k])-(F[2][l][i][j][k])
	+(F[0][l][i][j][k])-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);//3,2,4			
	i = N-1;//EAST
	U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(-(F[0][l][i][j][k])-(F[1][l][i][j][k])-(F[2][l][i][j][k])
	-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);//3,2,1				
	}
}


void updateF(){

}


/*
 * Allocates 3d vector
 *
 */
void allocate3d(FLOAT*** vec){
	int i,j;
	vec = malloc(N*sizeof(FLOAT**));
	for( i = 1; i < N-1; i++){
			vec[i] = malloc(N*sizeof(FLOAT*));
			for( j = 1; j < N-1; j++){
				vec[i][j] = malloc(N*sizeof(FLOAT));
		}
	}
}

/*
 * Gets Vectors delta in time according to the sound speed boundary which defines stability
 *
 */
FLOAT getDT(){
	// How to calculate sound speed?
	return 0;
}













