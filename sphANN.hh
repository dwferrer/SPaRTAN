/*
 * sphAnn.hh
 *
 *  Created on: Nov 4, 2012
 *      Author: dferrer
 *
 *      basic sph functions using the ANN library for nearest neighbor searching
 */
#ifndef SPHANN_HH
#define SPHANN_HH
#include <vector>
#include <cassert>
#include "threevector.hh"



//#include "Cover_Tree.h"
#include "StopWatch.cc"
#include <math.h>
#include <fenv.h>
#include "ANN/ANN.h"

#define NSMOOTH 60
#include "particlestruct2.hh"

#define pi 3.14159265

#define cube(x) ((x)*(x)*(x))
float kernel(float x, float h){ //cubic spline kernel
	float q = x/h;
	if (q > 2) return 0;
	assert(q >=0);
	if (q <= 1) return 1.0/(4.0*pi*cube(h))* (cube(2.0-q) -4.0*cube(1.0-q));


	return 1.0/(4.0*pi*cube(h)) * cube(2-q);
}

ANNkd_tree* kdtree;
ANNpointArray pts;
void setupCoverTree(particlestructure &p,int N){
	std::cout<<"Making kd tree with size: "<<N<<"...";
	pts = new ANNpoint[N];
	for(int i =0 ; i < N; i++) pts[i] = (ANNcoord *)(p.pos +i);
	kdtree = new ANNkd_tree (pts, N,3);
	std::cout<<"Done.\n";
	std::cout.flush();
}
void destroyCoverTree(){
	delete kdtree;
	delete pts;
	pts = NULL;
	kdtree = NULL;
}


float smoothingLength(particlestructure &p, float3 pos){ //get the smoothing length for a generic position

	ANNidxArray NN = new ANNidx[NSMOOTH];
	ANNdistArray NNd = new ANNdist[NSMOOTH];
	ANNpoint querry = (float *)&pos;
	kdtree->annkSearch(querry,NSMOOTH,NN,NNd,0.0);

	float maxdist = 0;
	for(int i = 0; i < NSMOOTH; i++){
		float dist = NNd[i];
		if (dist > maxdist) maxdist = dist;
	}
	float h = sqrt(maxdist)/2.0;
	assert(h >0);
	delete[] NN;
	delete[] NNd;
	return h;
}

float smoothingLength(particlestructure &p, float3 pos, size_t *NN, float *NNd){ //get the smoothing length for a generic position
	float maxdist = 0;
	for(int i = 0; i < NSMOOTH; i++){
		float dist = NNd[i];
		if (dist > maxdist) maxdist = dist;
	}
	float h = sqrt(maxdist)/2.0;
	assert(h >0);
	return h;
}

float smoothingLength(particlestructure &p, int k){ //get the smoothing length for particle k
	if (p.clean[k]) return p.h[k];

	float maxdist = 0;
	for(int i = 0; i < NSMOOTH; i++){
		float dist = p.getNNdist(k)[i];
		if (dist > maxdist) maxdist = dist;
	}
	float h = sqrt(maxdist)/2.0;
	p.h[k] = h;
	assert(h >0);

	return h;
}

void updateSmoothingLenghts(particlestructure &p, long long int np){
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < np; i++){
		p.h[i] = smoothingLength(p,i);
	}
}

void updateNN(particlestructure &p){
//#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < p.count; i++){
		kdtree->annkSearch(pts[i],NSMOOTH,p.getNN(i),p.getNNdist(i),0.0);
		//for(int j = 0; j < NSMOOTH-1; j++) assert(p.getNN(i)[j] == p.getNN(i)[j+1] ); //Make sure we got discrete points
	}
}

float density(particlestructure &p, float3 position){
	ANNidxArray NN = new ANNidx[NSMOOTH];
	ANNdistArray NNd = new ANNdist[NSMOOTH];
	ANNpoint querry = (float *)&position;
	kdtree->annkSearch(querry,NSMOOTH,NN,NNd,0.0);

	float result = 0;
	for(int i = 0; i < NSMOOTH; i++)
	{
		float x = NNd[i];
		result += p.mass[NN[i]] * kernel(x,smoothingLength(p,NN[i]));
	}
	//result06+= p.mass *kernel(0,p.h);
	//if(createdNN) delete NN;

	assert(result > 0);
	delete[] NN;
	delete[] NNd;
	return result;
}

float density(particlestructure &p, float3 position, size_t * NN, float * NNd){

	float result = 0;
	for(int i = 0; i < NSMOOTH; i++)
	{
		float x = NNd[i];
		result += p.mass[NN[i]] * kernel(x,p.h[NN[i]]);
	}

	assert(result > 0);

	return result;
}

float density(particlestructure &p,int k){
	if(p.clean[k]) return p.density[k];

	float result = 0;
	bool  createdNN = false;




	for(int i = 0; i < NSMOOTH; i++)
	{
		float x = p.getNNdist(k)[i];
		result += p.mass[p.getNN(k)[i]] * kernel(x,smoothingLength(p,p.getNN(k)[i]));
	}
	//result06+= p.mass *kernel(0,p.h);
	//if(createdNN) delete NN;

	assert(result > 0);
	p.density[k] = result;
	p.clean[k] = true;
	return result;
}


#ifndef kbr
#define kbr (366.57) //actually kb/mh in [AU^2/(ka^2 K)
#endif


#ifndef c0
#define c0 40.05 //au/ka
#endif

#ifndef rhoc
#define rhoc (6.73 *pow(10,-2))
#endif

float pressure(particlestructure &p, int k){ //units of [10^-6Msun /(AU ka^2)]
	if(p.T[k] == 0.0) return 0;
	float rho = density(p,k);
	return c0*c0 * rho * ( 1+ pow(rho/rhoc,2/5));
	//return kbr * density(p,NN) * temperature(p) /mu(p);
}

float temperature(particlestructure &p, int k){
	return pressure(p,k)/(kbr*2.0*density(p,k));
}


float sq(float x){
	return x*x;
}

float Fij(float r, float h){ //derivative calculated from mathematica
	float q = r/h;
	if (q >= 2) return 0;
	assert(q >=0);
	if(q == 0) return 0;
	if (q <=1) return -3/(4*pi*sq(cube(h))) * (4 * h - 3 *r);
	return -3/(4*pi*sq(cube(h))) * sq(2*h -r)/r;
}

float3 delWij(const float3 rij, const float h){
	float r = rij.norm();
	return rij *Fij(r,h);
}

float3 delbarWij(particlestructure &p, int i, int j){  //fully symmetric kernel for momentum update equation
	float3 rij = (p.pos[i] - p.pos[j]);
	return .5 * (delWij(rij,smoothingLength(p,i))+delWij(rij,smoothingLength(p,j)));

}
//#define gamma 1.66666667

float cs(particlestructure &p, int k){
	float rho = density(p,k);
	float gamma = 4.0/3 + .06666667 * rho/(rho + rhoc);
	return sqrt(gamma * pressure(p,k)/density(p,k));
}

float mu(particlestructure &p, int k){
	return .1*smoothingLength(p,k) *cs(p,k), density(p,k) ;
}

#ifndef alpha
#define alpha 1.0
#endif

#ifndef beta
#define beta 2.0
#endif

float PIij(particlestructure &p,int i, int j){

	float rhoij = .5*(density(p,i)+density(p,j));
	float hij = .5*(smoothingLength(p,i)+smoothingLength(p,j));
	float cij = .5*(cs(p,i)+cs(p,j));

	float3 rij = p.pos[i] - p.pos[j];
	float3 vij = p.velocity[i] - p.velocity[j];

	float rijvij = rij.dot(vij);

	float delvij = rijvij/(rij.norm2() +.01 *sq(hij));

	float nu = hij/rhoij *(alpha*cij -beta*delvij);
	assert(nu ==nu);
	assert(delvij == delvij);


	return -nu * delvij;
}

#ifndef RHOMIN
#define RHOMIN 1e-18
#endif
float3 accel(particlestructure &p, int i){ //fluid driven acceleration units are [AU /ka^2]

	 int *NN = (p.getNN(i));




	float3 result = float3(0,0,0);
	if(p.T[i] == 0) return result;
	if(density(p,i) < RHOMIN) return result; //we don't really have a fluid any more at this point
	for(int j = 0; j < NSMOOTH; j++){ //from springel+ hernquist, 2002
		if (p.T[NN[j]] ==0) continue;
		if (density(p,NN[j]) < RHOMIN)
		result += p.mass[NN[j]] *(   pressure(p,i)/sq(density(p,i))  +  pressure(p,NN[j])/sq(density(p,NN[j])) + PIij(p,i,NN[j]) )*delbarWij(p,i,NN[j]);
		assert(result == result); //check for nans
	}
	assert(result == result); //check for nans
	assert(!isnanf(result.x));
	assert(!isnanf(result.y));
	assert(!isnanf(result.z));
	return result;
}

#ifndef CP
#define CP 5.0/2.0
#endif

#ifndef kbi
#define kbi (1.0/3.08556 *pow(10,55)) //actually kb in [AU^2/(ka^2 K)
#endif

float cp(particlestructure &p,int i){
	return CP *kbr;
}

#ifndef KAPPA
#define KAPPA 1.0
#endif

float kappa (particlestructure &p, int i){ //thermal conductivity
	return 1.0;
}

float dTCond(particlestructure &p, int i){
	int *NN = (p.getNN(i));
	float *rij = p.getNNdist(i);
	float result = 0;
	for(int j = 0; j < NSMOOTH; j++){// from Monoghan,2005
		result += p.mass[NN[j]] *4 *kappa(p, i)*kappa(p,NN[j])/(kappa(p, i)+kappa(p,NN[j]))  * (p.T[i]-p.T[NN[j]]) *Fij(rij[j],.5*(smoothingLength(p,i)+smoothingLength(p,NN[j])) );
	}
	result /= cp(p, i);
	result *= kbi;
	return result;
}
float dot(float3 a, float3 b){
	return a.x*b.x+a.y*b.y+a.z*b.z;

}


float dTad(particlestructure &p, int i){
	int *NN = (p.getNN(i));
	float *rij = p.getNNdist(i);
		float result = 0;
		if (p.T[i] == 0) return result;
		for(int j = 0; j < NSMOOTH; j++){// from springel,2006
			float3 vij = p.velocity[i]-p.velocity[NN[j]];
			result += .5*p.mass[NN[j]] *(   pressure(p,i)/sq(density(p,i))  +  pressure(p,NN[j])/sq(density(p,NN[j]))  )*dot(vij,delbarWij(p,i,NN[j]));
		}
		result /= cp(p,i);

		return result;
}


#endif










