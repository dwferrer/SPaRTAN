/*
 * mcrt.hh
 * Monte-carlo radiative transfer
 *  Created on: Apr 5, 2013
 *      Author: dferrer
 */

#ifndef MCRT_HH_
#define MCRT_HH_

#define FLANN_USE_CUDA
#include <flann/flann.hpp>
#include <gsl/gsl_rng.h>

#include "sphANN.hh"
#include <math.h>

#define Clight 63239726.3451

gsl_rng * mcrng;

gsl_rng * rngsetup(){
	gsl_rng * r;
	const gsl_rng_type * T;
	gsl_rng_env_setup();

	T = gsl_rng_default; //The mersene twister
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,1234567890);
	return r;
}

struct lightpacket{
	float3 *position;
	size_t *NN;
	float  *NNd;
	float3  direction;
	float   tau;
	float   wavelength; //in microns
	float   lastdl;
};

void swaplp(lightpacket * l1, lightpacket * l2){
	lightpacket temp = *l1;
	memcpy(l2,l1,sizeof(lightpacket));
	memcpy(l1,&temp,sizeof(lightpacket));
}

float temp(float density){
	return 10;
}

//units are AU^2 /(1e-6 solar mass)
float opacity(float density, float T, float wavelength){ //we use the hildebrand dust opacity for now. We retain the temperature and density dependence in case a more complex fit is desired
	float result = 0;
	if (wavelength <= 250) result  = 2.5e3 * 1/.1125484 * wavelength;
	else result  = 6.25e5 * 1/.1125484 * sq(wavelength);
	return result;
}


float scattering(float density, float T, float wavelength){
	return .1; //fixme: place holder!!!!
}

float opticaldepth(float pathlength, float opacity, float density){
	return opacity * density * pathlength;
}

float radiativedl(float h, float rho){
	return std::min(h/2,rad/10);
}


void propogateLightPacket(particlestructure &p,lightpacket &l){

	*(l.position) += .5 *l.direction * l.lastdl;

	float localdensity = density(p,*l.position,l.NN,l.NNd);
	float localh = smoothingLength(p, *l.position,l.NN,l.NNd);

	float dl = radiativedl(localh, localdensity);

	float dtau = opticaldepth(.5 *(dl + l.lastdl),opacity(localdensity,temp(localdensity),l.wavelength),l.wavelength);


	l.tau -=dtau;

	*(l.position) += .5 *l.direction *dl;
	l.lastdl = dl;
}

#define MAXX (1.5*rad)


void populateLightPacket_plane(lightpacket &l,int d, int tb){ //put a bunch of light packets in the top/bottom(tb) plane perpendicular to d. basic rt test function

		float * pos = (float *) l.position;
		float *dir = (float *) &l.direction;
		for(int i = 0; i < 3; i++){
			pos[i] = 2*(gsl_rng_uniform(mcrng)-.5) *MAXX;
			dir[i] = 0;
		}
		pos[d] = tb *MAXX;
		dir[d] = 1;

		l.wavelength = 1000;
		l.tau = 0;
}

void populateLightPackets_plane(particlestructure &p, lightpacket * l,flann::Matrix<float> &pos,flann::Matrix<size_t>   NN,flann::Matrix<float> &NNd,long long int start,  long long int count, long long int total, long long int np){ //for now we just use the plane generator
	for(int i = 0; i < count; i++ ){
		l[i].position = (float3 *)(pos[start +i]);

		int c = (start +i)/np;
		int d = c/3;
		int tb = (c%2) *2 -1;
		populateLightPacket_plane(l[i],d,tb);
	}
}


// macro for accessing the array of boundary points
#define boundary(arr,coord,tb,x,y) (arr[n1d*n1d *(2*coord +tb) +n1d*x + y])

int wrap(float x, int n1d){

	int result = remainder(x,float(n1d));
	while (result >= n1d) result -= n1d;
	while (result < 0) result += n1d;
	return result;

}

//cull points if they are outside the boundary and update the boundary image
bool cull(lightpacket &l,float * result, int *counts, int n1d){
	float * pos = (float* ) l.position;

	for(int d = 0; d < 3; d++){
		if (pos[d] > -MAXX && pos[d] < MAXX) continue;
		else{
			int tb = -1;
			if(pos[d] > -MAXX) tb = 1;
			int x = wrap(pos[wrap(d+1,3)]/(2 *MAXX) *n1d,n1d);
			int y = wrap(pos[wrap(d+2,3)]/(2 *MAXX) *n1d,n1d);
			boundary(result,d,tb,x,y) += l.tau;
			boundary(counts,d,tb,x,y) += 1;
			return true;
		}
	}


	return false;

}

bool scatterabsorb(particlestructure &p, lightpacket &l){
	return false;
}


void doRT(particlestructure &p, long long int npackets, float * result,int * counts, int n1d, long long int blocksize = 1<<20){

	long long int np = 6 *npackets;
	assert(blocksize < npackets);
	mcrng = rngsetup();
	flann::Matrix<float> pos((float *)(p.pos),p.count,3);

	flann::KDTreeCuda3dIndexParams params;
	flann::SearchParams sp;
	flann::Index<flann::L2<float> > flannindex( pos, params );
	flannindex.buildIndex();

	flann::Matrix<float> querry(new float[3*blocksize],blocksize,3);
	flann::Matrix<size_t>   NN(new size_t[NSMOOTH *blocksize],blocksize,NSMOOTH);
	flann::Matrix<float> NNd(new float[NSMOOTH *blocksize],blocksize,NSMOOTH);
	long long int done = 0;

	lightpacket *l = new lightpacket[blocksize];

	populateLightPackets_plane(p,l,querry,NN,NNd,0,blocksize,np,npackets);
	for(int i =0; i < blocksize; i++){
		l[i].NN  = (size_t *)NN[i];
		l[i].NNd = (float *)NNd[i];
	}

	long long int lcount = blocksize;
	while (lcount >0){
		std::printf("Propogating light packets. %lld are done and %lld remain.",done,np-done);
		flannindex.knnSearch(querry,NN,NNd,NSMOOTH,sp); //get the nearest neighbors and distances for each light packet

		#pragma omp parallel for schedule(dynamic,1)
		for(int i =0; i < lcount; i++){
			propogateLightPacket(p,l[i]);
		}

		//cull, scatter, and absorb light packets
		std::vector<lightpacket *> toswap;
		for (int i = 0; i < lcount; i++){
			if(cull(l[i],result,counts,n1d) || scatterabsorb(p,l[i])){
				toswap.push_back(&(l[i]));
				done ++;
			}
		}

		//swap out the removed lightpackets for new ones
		for(int i = 0; i < toswap.size(); i++) swaplp(toswap[i],&(l[lcount -i -1]));

		int tofill = toswap.size();
		if (done + tofill > np){
			tofill = np-done;
			lcount -=toswap.size() - tofill;
		}
		populateLightPackets_plane(p,l,pos,NN,NNd,done,tofill,np,npackets);
	}

	delete [] querry.ptr();
	delete [] NN.ptr();
	delete [] NNd.ptr();
}



#endif /* MCRT_HH_ */
