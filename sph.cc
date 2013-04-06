/*
 * sph.cc
 *
 *  Created on: Nov 4, 2012
 *      Author: dferrer
 *
 *      Simple SPH test program to begin testing for stars project
 */

#include <vector>
#include <cassert>
#include "threevector.hh"
#include "particlestructure.h"

#include "Cover_Tree.h"
#include "StopWatch.cc"
#include <math.h>
#include <fenv.h>


#define pi 3.14159265

#define cube(x) ((x)*(x)*(x))
float kernel(float x, float h){ //cubic spline kernel
	float q = x/h;
	if (q > 2) return 0;
	assert(q >=0);
	if (q <= 1) return 1.0/(4.0*pi*cube(h))* (cube(2.0-q) -4.0*cube(1.0-q));


	return 1.0/(4.0*pi*cube(h)) * cube(2-q);
}

CoverTree<particlestructure> *particlecovertree;

void setupCoverTree(std::vector<particlestructure> &p,float maxsize = 8.0){
	std::cout<<"Making cover tree with size: "<<p.size()<<"...";
	particlecovertree = new CoverTree<particlestructure> (maxsize,p);
	std::cout<<"Done.\n";
	std::cout.flush();
}
void destroyCoverTree(){
	delete particlecovertree;
	particlecovertree = NULL;
}

#define NSMOOTH 60
float smoothingLength(particlestructure &p){
	if (p.clean) return p.h;
	bool  createdNN = false;


	float maxdist = 0;
	for(int i = 0; i < NSMOOTH; i++){
		float dist = (p.NN[i]->position-p.position).norm2();
		if (dist > maxdist) maxdist = dist;
	}
	float h = sqrt(maxdist)/2.0;
	p.h = h;
	assert(h >0);
	//if(createdNN) deleteNN;

	return h;
}

void updateSmoothingLenghts(std::vector<particlestructure> &p, long long int np){
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < np; i++){
		p[i].h = smoothingLength(p[i]);
	}
}

void updateNN(std::vector<particlestructure> &p){
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < p.size(); i++){
		std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(p[i],NSMOOTH);
		//#pragma unroll(NSMOOTH)
		for(int j = 0; j < NSMOOTH; j++)
			p[i].NN[j] =  &(p[NN[j].address]);
	}
}

float density(particlestructure &p){
	if(p.clean) return p.density;

	float result = 0;
	bool  createdNN = false;




	for(int i = 0; i < NSMOOTH; i++)
	{
		float x = (p.position-p.NN.at(i)->position).norm();
		result += p.NN.at(i)->mass * kernel(x,smoothingLength(p));
	}
	//result06+= p.mass *kernel(0,p.h);
	//if(createdNN) delete NN;

	assert(result > 0);
	p.density = result;
	p.clean = true;
	return result;
}


/*float mu(particlestructure &p){
	return 1.0;
}*/

#ifndef kbr
#define kbr (366.57) //actually kb/mh in [AU^2/(ka^2 K)
#endif




#ifndef c0
#define c0 40.05 //au/ka
#endif

#ifndef rhoc
#define rhoc (6.73 *pow(10,-8))
#endif

float pressure(particlestructure &p){ //units of [Msun /(AU ka^2)]
	if(p.T == 0.0) return 0;
	float rho = density(p);
	return c0*c0 * rho * ( 1+ pow(rho/rhoc,2/5));
	//return kbr * density(p,NN) * temperature(p) /mu(p);
}

float temperature(particlestructure &p){
	return pressure(p)/(kbr*2.0*density(p));
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

float3 delWij(float3 rij,float h){
	float r = rij.norm();
	return rij *Fij(r,h);
}

float3 delbarWij(particlestructure &p_i,particlestructure &p_j){  //fully symmetric kernel for momentum update equation
	float3 rij = (p_i.position - p_j.position);
	return .5 * (delWij(rij,smoothingLength(p_i))+delWij(rij,smoothingLength(p_j)));

}
#define gamma 1.66666667

float cs(particlestructure &p){
	return sqrt(gamma * pressure(p)/density(p));
}

float mu(particlestructure &p){
	return .1*smoothingLength(p) *cs(p), density(p) ;
}

#ifndef alpha
#define alpha 1.0
#endif

#ifndef beta
#define beta 2.0
#endif

float PIij(particlestructure &p_i,particlestructure &p_j){

	float rhoij = .5*(density(p_i)+density(p_j));
	float hij = .5*(smoothingLength(p_i)+smoothingLength(p_j));
	float cij = .5*(cs(p_i)+cs(p_j));

	float3 rij = p_i.position - p_j.position;
	float3 vij = p_i.velocity - p_j.velocity;

	float rijvij = rij.dot(vij);

	float delvij = rijvij/(rij.norm2() +.01 *sq(hij));

	float nu = hij/rhoij *(alpha*cij -beta*delvij);
	assert(nu ==nu);
	assert(delvij == delvij);


	return -nu * delvij;
}



float3 accel(particlestructure &p){ //fluid driven acceleration units are [AU /ka^2]

	std::vector<particlestructure *> &NN = p.NN;




	float3 result = float3(0,0,0);
	if(p.T == 0) return result;
	for(int j = 0; j < NSMOOTH; j++){ //from springel+ hernquist, 2002
		if (NN[j]->T ==0) continue;
		result += -NN[j]->mass *(   pressure(p)/sq(density(p))  +  pressure(*NN[j])/sq(density(*NN[j])) + PIij(p,*NN[j]) )*delbarWij(p,*NN[j]);
		assert(result == result); //check for nans
	}
	//result += -p.mass *(   pressure(p,&NN)/sq(density(p,&NN))  +  pressure(p,&NN)/sq(density(p,&NN))  )*delbarWij(p,p);
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

float cp(particlestructure &p){
	return CP *kbr;
}

#ifndef KAPPA
#define KAPPA 1.0
#endif

float kappa (particlestructure &p){ //thermal conductivity
	return 1.0;
}

float dTCond(particlestructure &p){
	std::vector<particlestructure*> &NN = p.NN;
	float result = 0;
	for(int j = 0; j < NSMOOTH; j++){// from Monoghan,2005
		float rij = (p.position-NN[j]->position).norm();
		result += NN[j]->mass *4 *kappa(p)*kappa(*NN[j])/(kappa(p)+kappa(*NN[j]))  * (p.T-NN[j]->T) *Fij(rij,.5*(smoothingLength(p)+smoothingLength(*NN[j])) );
	}
	result /= cp(p);
	result *= kbi;
	return result;
}
float dot(float3 a, float3 b){
	return a.x*b.x+a.y*b.y+a.z*b.z;

}


float dTad(particlestructure &p){
		std::vector<particlestructure *> &NN = p.NN;
		float result = 0;
		for(int j = 0; j < NSMOOTH; j++){// from springel,2006
			float rij = (p.position-NN[j]->position).norm();
			float3 vij = p.velocity-NN[j]->velocity;

			result += .5*NN[j]->mass *(   pressure(p)/sq(density(p))  +  pressure(*NN[j])/sq(density(*NN[j]))  )*dot(vij,delbarWij(p,*NN[j]));
		}
		result /= cp(p);

		return result;
}


#include "gravity.hh"
#include "ic.hh"
#include "sink.hh"
#include "timestep.hh"
#include "micro.hh"
#include "reversible.hh"
#include "pyinterface.hh"









