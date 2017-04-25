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
FLOAT getDT(FLOAT dx);
FLOAT*** allocate3d();
FLOAT evolve(float dt,float dx);


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
FLOAT p_ini;
// Gas initial density
FLOAT rho_ini;
// Velocity v
// Pressure p, internal energy e per unit mass, density rho, total energy E per unit volume
FLOAT e_ini;

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
	p_ini = 1.0*atm;
	// Blast is defined J/kg
	blast = pow(10.0,10.0);
	// Internal energy depends on pressure and volume: p = rho*e*(gamma-1)
	e_ini = p_ini/(rho_ini*(gamma-1));
	// N is defined (2 more for phantom) 
	N = ((FLOAT)L)/((FLOAT)resol) + 2;
	// Variables are allocated
	allocateAll();
	printf("aal iis weell (ALLOCATION)N=%d\n",N);
	// Variables are initialised
	__init__();
	printf("aal iis weell (INITIALIZATION)\n");

	   ///////////////////////////
	  //  TIME EVOLUTION       //
	 ///////////////////////////

	// Time variable t, deltas in position and time: dx,dt
	FLOAT time = 0;
	FLOAT dx = resol;
	FLOAT dt = dx/sqrt(rho_ini*gamma*(gamma-1)*blast);
	bool go_on = true;

	do{	
		printf("%f,%f\n",time,dt);
		// Evolves vector U, actualizes dt
		dt = evolve(dt,dx);
		printf("aal iis weell (EVOLVE)\n");

		time+=dt;

	} while(go_on);
	return 0;
}






//------------------------------------------------------------------------------
//  AUXILIAR METHODS
//------------------------------------------------------------------------------


/*
 * Allocates cubic variables U and Fs
 *
 */
void allocateAll(){
	U = malloc(5*sizeof(FLOAT***));
	F = malloc(3*sizeof(FLOAT****));
	int i,j;
	for( i = 0; i < 3; i++ ){
		F[i] = malloc(5*sizeof(FLOAT***));
		for( j = 0; j < 5; j++ ){
			F[i][j] = allocate3d();
		}
	}
	for( j = 0; j < 5; j++ ){
		U[j] = allocate3d();
	}	
}

/*
 * Initializes variables according to the corresponding initial conditions
 *
 */
void __init__(){

	int i,j,l,m,n;
	// Initial velocity and Energy density are defined E = rho(e+0.5v**2)
	FLOAT E_ini = rho_ini*e_ini;
	FLOAT* v = malloc(3*sizeof(FLOAT));
	v[0] = 0;
	v[1] = 0;
	v[2] = 0;

	for( l = 0; l < N; l++){
		for( m = 0; m < N; m++){
			for( n = 0; n < N; n++){

				// Defining F
				for( j = 0; j < 3; j++ ){
					F[j][0][l][m][n] = rho_ini*v[j];
					for( i = 1; i < 4; i++ ){
						F[j][i][l][m][n] = rho_ini*v[j]*v[i-1] + (j==(i-1))*p_ini;
						//delta_j,i-1 = (j==(i-1))
					}
					F[j][4][l][m][n] = (rho_ini*E_ini+p_ini)*v[j];
				}

				// Defining U
				U[0][l][m][n] = rho_ini;
				for( i = 1; i < 4; i++ ){
						U[i][l][m][n] = rho_ini*v[i-1];
				}
				U[4][l][m][n] = rho_ini*E_ini;

			}
		}
	}

	// Corrects for the blast
	E_ini = rho_ini*blast;
	FLOAT p_blast =  (gamma-1)*blast;
	// Defining F
	for( j = 0; j < 3; j++ ){
		for( i = 1; i < 4; i++ ){
			F[j][i][N/2][N/2][N/2] = rho_ini*v[j]*v[i-1] + (j==(i-1))*p_blast;
			//delta_j,i-1 = (j==(i-1))
		}
		F[j][4][N/2][N/2][N/2] = (rho_ini*E_ini+p_ini)*v[j];
	}

	U[4][N/2][N/2][N/2] = rho_ini*E_ini;
	free(v);
}



/*
 * Evolves vector U and actualizes vector F afterwards, given dx and dt according to the finite volume method
 * Finite volume is performed over a volume-centered cubic uniform grid
 * Volume integrals are approximated given variables constant over finite volumes
 * Surface integrals are approximated given variables constant over square surfaces
 * Values of variables on surfaces are calculated as the mean of the adjacent volumes
 * Boundary conditions are taken as free boundary :
 * (directional derivative perpendicular to the boundary surface vanishes) (This vanishes the respective surface integral)
 * !!!RETURN!!! This method returns the time difference dt for the other time step
 * delta in time is calculated according to the sound speed boundary which defines stability
 *
 */
FLOAT evolve(float dt,float dx){
	// Creates new array to avoid over-evolving
	int i,j,k,l,m;
	FLOAT**** U_new = malloc(5*sizeof(FLOAT***));

	for( i = 0; i < 5; i++){
		U_new[i] = allocate3d();
	}

	// Temporal variables to deduce maximum velocity and sound speed
	FLOAT cmax = -1;
	FLOAT vmax = -1;


	  ////////////////////////////
	 //   Evolution of Unew    //
	////////////////////////////
	

	for( l = 0; l < 5; l++){
		for( i = 1; i < N-1; i++){
			for( j = 1; j < N-1; j++){
				for( k = 1; k < N-1; k++){

					U_new[l][i][j][k]=U[l][i][j][k]-(0.5*dt/dx)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
					-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);

				}
			}
		}		
	}


	  //////////////////////////////////////////////////////////////////////////
	 // Passes from U_new to U, Updates F and calculates maximum sound speed //
	//////////////////////////////////////////////////////////////////////////

 
	// Temporal variables
	FLOAT rho,e,E,p;
	FLOAT* v = malloc(3*sizeof(FLOAT));

	for( i = 1; i < N-1; i++){
		for( j = 1; j < N-1; j++){
			for( k = 1; k < N-1; k++){
				for( l = 0; l < 5; l++){

					// U is actualized
					U[l][i][j][k] = U_new[l][i][j][k];	

				}

				// From auxiliar variable U to usual variables rho,p,E,v
				rho  = U[0][i][j][k];
				v[0] = U[1][i][j][k]/rho;
				v[1] = U[2][i][j][k]/rho;
				v[2] = U[3][i][j][k]/rho;
				E    = U[4][i][j][k];
				e    = (E/rho) - 0.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
				p    = (gamma-1)*rho*e;

				// Calculates sound speed and verifies maximum for sound speed and velocities
				//FLOAT c_now = sqrt((gamma+1)*(E[i][j][k]+p[i][j][k]/rho[i][j][k])); DANGER
				FLOAT c_now = sqrt((gamma+1)*(E+p));
				if( c_now > cmax ){
					cmax = c_now;
				}
				if( abs(v[0]) > vmax ){
					vmax = abs(v[0]);
				}
				if( abs(v[1]) > vmax ){
					vmax = abs(v[1]);
				}
				if( abs(v[2]) > vmax ){
					vmax = abs(v[2]);
				}
	

				// From usual variables to auxiliar Fi
				for( l = 0; l < 3; l++ ){

					F[l][0][i][j][k] = rho*v[l];

					for( m = 1; m < 4; m++ ){

						F[l][m][i][j][k] = rho*v[l]*v[m-1] + (l==(m-1))*p;

					}

					F[l][4][i][j][k] = (rho*E+p)*v[l];
				}
			}
		}
	}


	  /////////////////////////////////
	 // Free memory for U_new and v //
	/////////////////////////////////


	for( l = 0; l < 5; l++){	
		for( i = 0; i < N; i++){
			for( j = 0; j < N; j++){
				free(U_new[l][i][j]);
			}
			free(U_new[l][i]);
		}
		free(U_new[l]);
	}
	free(U_new);
	free(v);
	
	// RETURNS delta time
	return dx/(vmax+cmax);
}


/*
 * Allocates 3d vector
 *
 */
FLOAT*** allocate3d(){
	int i,j;
	FLOAT*** vec = malloc(N*sizeof(FLOAT**));
	for( i = 0; i < N; i++){
		vec[i] = malloc(N*sizeof(FLOAT*));
		for( j = 0; j < N; j++){
			vec[i][j] = malloc(N*sizeof(FLOAT));
		}
	}
	return vec;
}














