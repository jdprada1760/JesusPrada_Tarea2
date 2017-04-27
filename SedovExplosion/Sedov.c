/*
 * This program simulates a Sedov explosion in a cube with side 256m at 300K and 1atm
 * The blast liberates 10**10J from a cube of 2m in the center of the cube
 * Spatial resolution is 2m  in each direction
 * Method used is finite volumes to solve Euler equations
 * This code STOPS when the shock wave is of 120m of radius
 * Prints radial density at r= 10,60,120m
 * We will assume gass treated is nitrogen as it composes most of atmospheric gas
 * EULER EQUATION FORMALISM IS OBTAINED FROM ****** http://197.14.51.10:81/pmb/GENIE_DES_PROCEDES/An%20Introduction%20to%20Scientific%20Computing%20Twelve%20Computational%20Projects%20Solved%20with%20MATLAB.pdf  ***** page 215
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
// 1atm in Pa
#define atm 101325
// Boltzmann constant divided by NITROGEN molecular mass (Nitrogen assumed)
#define R 296.803 
// Adiabatic coefficient
#define gamma (5.0/3.0)
// Boolean definitions
#define bool int
#define true 1
#define false 0


void allocateAll();
void __init__();
FLOAT getDT(FLOAT dx);
FLOAT*** allocate3d();
FLOAT evolve(FLOAT inv_vel);
void get_radius();
void save();
// Actual radius of the blast wave
FLOAT radius = -1;

//------------------------------------------------------------------------------
//   GLOBAL VARIABLES
//------------------------------------------------------------------------------

// Box size in m 
int L = 256;
// Spatial resolution in m
int resol = 2;
// Size of the tridimentional arrays and its half
int N, halfN;
// Blast energy liberated in joules/kg (from little square at the center)
FLOAT blast;
// Gas temperature in Kelvin
FLOAT T;
// Gas Pressure in Pascals
FLOAT p_ini;
// Gas initial density
FLOAT rho_ini;

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
// File to write 
FILE* density;
FILE* energy;

//------------------------------------------------------------------------------
//   MAIN
//------------------------------------------------------------------------------


int main(int argc, char** argv){

	// Temperature is 300K
	T = 300.0;
	// Pressure is atmospheric
	p_ini = 1.0*atm;
	// Initial gas density is calculated withstate equations rho*R*T = p (((((R is normal constant over molecular mass)))))
	rho_ini = p_ini/(R*T);
	// Blast is defined J/kg
	blast = pow(10.0,10.0);
	// N is defined (2 more for phantom) 
	N = ((FLOAT)L)/((FLOAT)resol) + 2;
	halfN = N/2;
	// Variables are allocated
	allocateAll();
	printf("aal iis weell (ALLOCATION)N=%d\n",N);
	// Variables are initialised
	__init__();
	printf("aal iis weell (INITIALIZATION)\n");
	// File to write
	density = fopen("Density_sedov.data","w");
	energy = fopen("energy_sedov.data","w");


	   ///////////////////////////
	  //  TIME EVOLUTION       //
	 ///////////////////////////

	// Time variable t, deltas in position and time: dx,dt
	FLOAT time = 0;
	FLOAT dx = resol;
	// Stores inverse of velocity ((((initial velocity is equal to sqrt(gamma*R*T))))
	FLOAT p_blast = rho_ini*((gamma-1)*blast);
	FLOAT inv_vel = 1.0/sqrt(gamma*p_blast/rho_ini);
	bool go_on = true;
	int ji = 0;
	do{	
		// Evolves vector U, actualizes dt
		printf("%f____%f____%f,\n",time,dx*inv_vel,radius);
		inv_vel = evolve(inv_vel);
		printf("aal iis weell (EVOLVE)\n");
		// 1/inv_vel = vel = dx/dt
		time+=dx*inv_vel;
		save();

		if(ji>20){ go_on = false;}
		ji++;

	} while(go_on);
	fclose(density);
	fclose(energy);
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
	FLOAT E_ini = p_ini/(gamma-1);
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
					F[j][4][l][m][n] = (E_ini+p_ini)*v[j];
				}

				// Defining U
				U[0][l][m][n] = rho_ini;
				for( i = 1; i < 4; i++ ){
						U[i][l][m][n] = rho_ini*v[i-1];
				}
				U[4][l][m][n] = E_ini;

			}
		}
	}

	// Corrects for the blast (Assume energy is created by overdensity)
	FLOAT E_blast = rho_ini*blast;
	FLOAT p_blast = rho_ini*((gamma-1)*blast);
	// Defining F
	for( j = 0; j < 3; j++ ){
		for( i = 1; i < 4; i++ ){
			F[j][i][N/2][N/2][N/2] = rho_ini*v[j]*v[i-1] + (j==(i-1))*p_blast;
			//delta_j,i-1 = (j==(i-1))
		}
		F[j][4][N/2][N/2][N/2] = (E_blast+p_blast)*v[j];
	}
	U[4][N/2][N/2][N/2] = E_blast;
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
 * !!!RETURN!!! This method returns the inverse of maximum velocity for the other time step
 *
 */
FLOAT evolve(FLOAT inv_vel){
	// Creates new array to avoid over-evolving
	int i,j,k,l,m;
	FLOAT**** U_new = malloc(5*sizeof(FLOAT***));

	for( i = 0; i < 5; i++){
		U_new[i] = allocate3d();
	}

	// Temporal variables to deduce maximum velocity and sound
	FLOAT cmax = -1;
	FLOAT vmax = -1;


	  ////////////////////////////
	 //   Evolution of Unew    //
	////////////////////////////
	

	for( l = 0; l < 5; l++){
		for( i = 1; i < N-1; i++){
			for( j = 1; j < N-1; j++){
				for( k = 1; k < N-1; k++){

					U_new[l][i][j][k]=U[l][i][j][k]-(0.5*inv_vel)*(F[0][l][i+1][j][k]+F[1][l][i][j+1][k]+F[2][l][i][j][k+1]
					-F[0][l][i-1][j][k]-F[1][l][i][j-1][k]-F[2][l][i][j][k-1]);

				}
			}
		}		
	}


	  //////////////////////////////////////////////////////////////////////////
	 // Passes from U_new to U, Updates F and calculates maximum sound speed //
	//////////////////////////////////////////////////////////////////////////

 
	// Temporal variables
	FLOAT rho,E,p;
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
				p    = (gamma-1)*(E-0.5*rho*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2)));

				// Calculates sound speed and verifies maximum for sound speed and velocities
				FLOAT c_now = sqrt((gamma)*(p/rho));
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
	
	// RETURNS inverse of maximum velocity
	return 1.0/(vmax+cmax);
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


/*
 * Averages over density spherically and saves it
 *
 */
void get_radius(){
	
	int i;
	int maxind;
	int maxrho = -1;
	for ( i = N/2; i < N-2; i++){
		fprintf(density,"%f,",U[0][i][N/2][N/2]);
		fprintf(energy,"%f,",U[4][i][N/2][N/2]);
		if( U[0][i][N/2][N/2] > maxrho ){
			maxrho = U[0][i][N/2][N/2];
			maxind = i;
		}
	}
	fprintf(density,"%f\n",U[0][N-2][N/2][N/2]);
	fprintf(energy,"%f\n",U[4][N-2][N/2][N/2]);
	// abs(maxind - N/2)*dx is the radius
	radius = abs(maxind - N/2)*resol;

}

void save(){
	// indie is the index 
	
	int i,j,k, indie;
	// Auxiliar arrays to calculate rho_prom
	FLOAT* rhoprom = malloc((halfN+1)*sizeof(FLOAT));
	int* count = malloc((halfN+1)*sizeof(int));
	for ( i = 0; i < (halfN+1); i++){
		//printf("JIJIJI___%d___%d\n",i, halfN+1);
		rhoprom[i] = 0;
		count[i] = 0;	
	}
	
	for ( i = 0; i < N; i++){
		for ( j = 0; j < N; j++){
			for ( k = 0; k < N; k++){
				indie = (0.5*(sqrt((pow((2*i-N),2)+pow((2*j-N),2)+pow((2*k-N),2)))));
				if(indie < halfN+1){
					// schock wave does not reach sides of the cube. Not even its corners
					rhoprom[indie] += U[0][i][j][k];
					count[indie] += 1;
					//printf("JIJIJI___%d___%d\n",indie, halfN+1);
				}
			}
		}
	}
	for ( i = 0; i <  halfN; i++){
		//printf("num___%d\n",i);
		if( count[i] != 0 ){
			rhoprom[i]  = rhoprom[i]/count[i];
		}
		else{
			rhoprom[i]  = 0;
		}
		fprintf(density,"%f,",rhoprom[i]);
		//printf("JIJIJI___%d___%d\n",i, halfN);
	}
	if( count[halfN] != 0 ){
		rhoprom[halfN]  = rhoprom[halfN]/count[halfN];
	}
	else{
		rhoprom[halfN]  = 0;
	}
	fprintf(density,"%f\n",rhoprom[halfN]);

	
	
	//free(rhoprom);
	//free(count);
	//fprintf(energy,"%f\n",U[4][N-2][N/2][N/2]);
	// abs(maxind - N/2)*dx is the radius
	//radius = abs(maxind - N/2)*resol;

}












