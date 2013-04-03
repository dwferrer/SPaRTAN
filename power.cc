/*
 * power.cc
 *	faster, more accurate power spectrum from input particles
 *  Created on: Sep 27, 2012
 *      Author: dferrer
 */


#include <gsl/gsl_rng.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include "particlestructure.h"



int wrap(int input, int max){
	int result = input % max;
	while(result <0){
		result += max;
	}
	return result;
}

#define squ(x) (x)*(x)
#define arr(a,x,y,z) a[gridN1D*gridN1D *x + gridN1D*y + z]
#define den(x,y,z) arr(density,x,y,z)

extern "C"{

void tsc(particlestructure * positions, float * density, long long int NP, int gridN1D,float boxsize){
	long long int n;
	for(n = 0; n < NP; n++){
		float px = positions[n].position.x/boxsize * gridN1D;
		float py = positions[n].position.y/boxsize * gridN1D;
		float pz = positions[n].position.z/boxsize * gridN1D;

		//round to nearest cell center (we offset the grid .5 so we can use floor instead of round)
		int ix = floor(px+.5);
		int iy = floor(py+.5);
		int iz = floor(pz+.5);

		//calculate distance to cell center
		float dx = ix - px;
		float dy = iy - py;
		float dz = iz - pz;

		//find the tsc weights for each dimension
		float wx = .75 -       squ(dx);
		float wxm1 = .5 * squ(.5 + dx);
		float wxp1 = .5 * squ(.5 - dx);
		float wy = .75 -       squ(dy);
		float wym1 = .5 * squ(.5 + dy);
		float wyp1 = .5 * squ(.5 - dy);
		float wz = .75 -       squ(dz);
		float wzm1 = .5 * squ(.5 + dz);
		float wzp1 = .5 * squ(.5 - dz);

		//find the wrapped x,y,z grid locations of the points we need to change
		int ixm1 =wrap(ix-1,gridN1D);
		int iym1 =wrap(iy-1,gridN1D);
		int izm1 =wrap(iz-1,gridN1D);
		int ixw = wrap(ix,gridN1D);
		int iyw = wrap(iy,gridN1D);
		int izw = wrap(iz,gridN1D);
		int ixp1 =wrap(ix+1,gridN1D);
		int iyp1 =wrap(iy+1,gridN1D);
		int izp1 =wrap(iz+1,gridN1D);

		//change the 27 cells that the cloud touches
		den(ixm1,iym1,izm1) += wxm1*wym1*wzm1;
		den(ixw, iym1,izm1) += wx  *wym1*wzm1;
		den(ixp1,iym1,izm1) += wxp1*wym1*wzm1;

		den(ixm1,iyw ,izm1) += wxm1*wy  *wzm1;
		den(ixw, iyw ,izm1) += wx  *wy  *wzm1;
		den(ixp1,iyw ,izm1) += wxp1*wy  *wzm1;

		den(ixm1,iyp1,izm1) += wxm1*wyp1*wzm1;
		den(ixw, iyp1,izm1) += wx  *wyp1*wzm1;
		den(ixp1,iyp1,izm1) += wxp1*wyp1*wzm1;

		den(ixm1,iym1,izw ) += wxm1*wym1*wz  ;
		den(ixw, iym1,izw ) += wx  *wym1*wz  ;
		den(ixp1,iym1,izw ) += wxp1*wym1*wz  ;

		den(ixm1,iyw ,izw ) += wxm1*wy  *wz  ;
		den(ixw, iyw ,izw ) += wx  *wy  *wz  ;
		den(ixp1,iyw ,izw ) += wxp1*wy  *wz  ;

		den(ixm1,iyp1,izw ) += wxm1*wyp1*wz  ;
		den(ixw, iyp1,izw ) += wx  *wyp1*wz  ;
		den(ixp1,iyp1,izw ) += wxp1*wyp1*wz  ;

		den(ixm1,iym1,izp1) += wxm1*wym1*wzp1;
		den(ixw, iym1,izp1) += wx  *wym1*wzp1;
		den(ixp1,iym1,izp1) += wxp1*wym1*wzp1;

		den(ixm1,iyw ,izp1) += wxm1*wy  *wzp1;
		den(ixw, iyw ,izp1) += wx  *wy  *wzp1;
		den(ixp1,iyw ,izp1) += wxp1*wy  *wzp1;

		den(ixm1,iyp1,izp1) += wxm1*wyp1*wzp1;
		den(ixw, iyp1,izp1) += wx  *wyp1*wzp1;
		den(ixp1,iyp1,izp1) += wxp1*wyp1*wzp1;
	}

}

void fftAndWindow(float * density, long long int NP, int gridN1D,float boxsize){
	fftwf_complex *input = new fftw_complex[gridN1D];

	fftw_plan plan = fftwf_plan_dft_3d(gridN1D,gridN1D,gridN1D,input,input,-1,FFTW_MEASURE);
}

float3 * randIC(int N1D) { //generate uniform random initial conditions
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();

	T = gsl_rng_default; //The mersene twister
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,time(0));

	float3 * pos = new float3[N1D*N1D*N1D];
	long long int n;
	for (n = 0; n < N1D*N1D*N1D; n++){
		pos[n].x = gsl_rng_uniform(r);
		pos[n].z = gsl_rng_uniform(r);
		pos[n].y = gsl_rng_uniform(r);
	}

	return pos;
}

float * zeros(int n1d){
	float * z= new float[n1d*n1d*n1d];
	long long int n;
	for (n = 0; n<n1d*n1d*n1d; n++)
		z[n] = 0;
	return z;
}
}


