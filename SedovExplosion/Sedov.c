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
FLOAT get_axis_flux( int i, int j, int k, int axis, int n );
FLOAT flux_limiter( FLOAT r );
// Actual radius of the blast wave
FLOAT radius = -1;
FLOAT time;

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
FILE* times;

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
	blast = pow(10.0,8.0);
	// N is defined (2 more for phantom 
	N = ((FLOAT)L)/((FLOAT)resol) + 4;
	halfN = N/2;
	// Variables are allocated
	allocateAll();
	printf("aal iis weell (ALLOCATION)N=%d\n",N);
	// Variables are initialised
	__init__();
	printf("aal iis weell (INITIALIZATION)\n");
	// File to write
	density = fopen("Density_sedov.data","w");
	times = fopen("Times_sedov.data","w");


	   ///////////////////////////
	  //  TIME EVOLUTION       //
	 ///////////////////////////

	// Time variable t, deltas in position and time: dx,dt
	time = 0;
	FLOAT dx = resol;
	// Stores inverse of velocity ((((initial velocity is equal to sqrt(gamma*R*T))))
	FLOAT p_blast = rho_ini*((gamma-1)*blast);
	FLOAT inv_vel = 1.0/sqrt(gamma*p_blast/rho_ini);
	bool go_on = true;
	int ji = 0;
	
	do{	
		// Evolves vector U, actualizes dt
		printf("%f____%f____%f__________________%f,%f\n",time,dx*inv_vel,radius,U[0][halfN+5][halfN][halfN],U[1][halfN+5][halfN][halfN]);
		inv_vel = evolve(inv_vel);
		printf("aal iis weell (EVOLVE)\n");
		// 1/inv_vel = vel = dx/dt
		time+=dx*inv_vel;
		//save();
		get_radius();

		if(ji>40){ go_on = false;}
		ji++;

	} while(go_on);

	fclose(density);
	fclose(times);
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
 * Surface integrals are calculated according to TVD scheme with Roe approximations and superbee flux limiter
 * Boundary conditions are taken as free boundary :
 * (directional derivative perpendicular to the boundary surface vanishes) -> ghost cells needed
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
	
	for( i = 2; i < N-2; i++){
		for( j = 2; j < N-2; j++){
			for( k = 2; k < N-2; k++){
				for( l = 0; l < 5; l++){
					FLOAT temp = get_axis_flux( i, j, k, 0, l )+get_axis_flux( i, j, k, 1, l )+get_axis_flux( i, j, k, 2, l );
					U_new[l][i][j][k]=U[l][i][j][k]-(inv_vel)*temp;

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

	for( i = 2; i < N-2; i++){
		for( j = 2; j < N-2; j++){
			for( k = 2; k < N-2; k++){
				for( l = 0; l < 5; l++){

					// U is actualized
					U[l][i][j][k] = U_new[l][i][j][k];	

				}

				// From auxiliar variable U to usual variables rho,p,E,v
				rho  = U[0][i][j][k];
				//if( rho < 0 ){ rho = 0; }
				v[0] = U[1][i][j][k]/rho;
				v[1] = U[2][i][j][k]/rho;
				v[2] = U[3][i][j][k]/rho;
				E    = U[4][i][j][k];
				//if( E < 0 ){ E = 0; }
				p    = (gamma-1)*(E-0.5*rho*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2)));
				//if( p < 0 ){ p = 0; }

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
						//printf("%d,%d,%d\n",(l==(m-1)),l,m-1);

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
	return 10.0/(vmax+cmax);
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
	for ( i = 2; i < N-2; i++){
		fprintf(density,"%f,",U[0][i][N/2][N/2]);
		//fprintf(energy,"%f,",U[4][i][N/2][N/2]);
		if( U[0][i][N/2][N/2] > maxrho ){
			maxrho = U[0][i][N/2][N/2];
			maxind = i;
		}
	}
	fprintf(density,"%f\n",U[0][N-2][N/2][N/2]);
	//fprintf(energy,"%f\n",U[4][N-2][N/2][N/2]);
	// abs(maxind - N/2)*dx is the radius
	radius = abs(maxind - N/2)*resol;
	fprintf(times,"%f\n",time);

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
	fprintf(times,"%f\n",time);

	
	
	free(rhoprom);
	free(count);
	//fprintf(energy,"%f\n",U[4][N-2][N/2][N/2]);
	// abs(maxind - N/2)*dx is the radius
	//radius = abs(maxind - N/2)*resol;

}

// gets area flux according to TVD-Roe scheme with superbee flux limiter
// n is component from vectors u and F. Axis is the axis in which flux is going to be calculated
// https://tspace.library.utoronto.ca/bitstream/1807/14332/1/MQ45872.pdf
FLOAT get_axis_flux( int i, int j, int k, int axis, int n ){
	// indices of the points ll(-2), l(-1), c(0), r(+1), rr(+2)
	int** indices = malloc(5*sizeof(int*));
	int l;
	for( l = 0; l < 5; l++ ){
		indices[l] = malloc(3*sizeof(int));
		indices[l][0] = i;
		indices[l][1] = j;
		indices[l][2] = k;
		indices[l][axis] += (l-2);
	}
	// f in i + 1/2
	FLOAT rho_1 = sqrt(U[0][indices[2][0]][indices[2][1]][indices[2][2]]);
	FLOAT rho_2 = sqrt(U[0][indices[3][0]][indices[3][1]][indices[3][2]]);
	FLOAT f_12 = (rho_1*F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]] + rho_2*F[axis][n][indices[3][0]][indices[3][1]][indices[3][2]])/(rho_1+rho_2);
	// f in i + 3/2
	rho_1 = sqrt(U[0][indices[3][0]][indices[3][1]][indices[3][2]]);
	rho_2 = sqrt(U[0][indices[4][0]][indices[4][1]][indices[4][2]]);
	FLOAT f_32 = (rho_1*F[axis][n][indices[3][0]][indices[3][1]][indices[3][2]] + rho_2*F[axis][n][indices[4][0]][indices[4][1]][indices[4][2]])/(rho_1+rho_2);
	// f in i - 1/2
	rho_1 = sqrt(U[0][indices[1][0]][indices[1][1]][indices[1][2]]);
	rho_2 = sqrt(U[0][indices[2][0]][indices[2][1]][indices[2][2]]);
	FLOAT f_m12 = (rho_1*F[axis][n][indices[1][0]][indices[1][1]][indices[1][2]] + rho_2*F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]])/(rho_1+rho_2);
	//  f in i - 3/2
	rho_1 = sqrt(U[0][indices[0][0]][indices[0][1]][indices[0][2]]);
	rho_2 = sqrt(U[0][indices[1][0]][indices[1][1]][indices[1][2]]);
	FLOAT f_m32 = (rho_1*F[axis][n][indices[0][0]][indices[0][1]][indices[0][2]] + rho_2*F[axis][n][indices[1][0]][indices[1][1]][indices[1][2]])/(rho_1+rho_2);

	// Calculates the ratios for flux limiter
	// ratio + in i-1/2
	FLOAT pr_m12 = (F[axis][n][indices[3][0]][indices[3][1]][indices[3][2]]-f_12)/(F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]]-f_m12);
	// ratio + in i-3/2
	FLOAT pr_m32 = (F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]]-f_m12)/(F[axis][n][indices[1][0]][indices[1][1]][indices[1][2]]-f_m32);
	// ratio - in 1+1/2
	FLOAT mr_12 = (F[axis][n][indices[1][0]][indices[1][1]][indices[1][2]]-f_m12)/(F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]]-f_12);
	// ratio + in i+3/2
	FLOAT mr_32 = (F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]]-f_12)/(F[axis][n][indices[3][0]][indices[3][1]][indices[3][2]]-f_32);
	
	// Contribution from r+ and r-
	FLOAT temp1 = (1.0 + 0.5*flux_limiter(pr_m12) -0.5*flux_limiter(pr_m32)/(pr_m32))*(F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]]-f_m12);
	FLOAT temp2 = (1.0 + 0.5*flux_limiter(mr_12) -0.5*flux_limiter(mr_32)/(mr_32))*(f_12-F[axis][n][indices[2][0]][indices[2][1]][indices[2][2]]);
	
	for( l = 0; l < 5; l++ ){
		free(indices[l]);
	}
	free(indices);
	
	if( pr_m32 == 0 || isnan(pr_m32) || isnan(pr_m12)  ){ temp1 = 0;}
	if( mr_32 == 0 || isnan(mr_32) || isnan(mr_12) ){ temp2 = 0;}
	return temp1+temp2;


}


// Flux limiter, superbee scheme max(0,min(1,2r),min(2,r))
FLOAT flux_limiter( FLOAT r ){
	//FLOAT min1, min2;
	if( r >= 2 ){
		//min1 = 1;
		//min2 = 2;
		return 2;
	}
	else if ( r>= 0.5 ){
		//min1 = 1;
		//min2 = r;
		if( 1 >= r ){
			return 1;
		}
		else{
			return r;
		}
	} 
	else if( r>= 0 ){
		//min1 = r;
		//min2 = 2*r;
		return 2*r;
		
	}
	return 0;		

}



