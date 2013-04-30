/*
 * particlestruct2.hh
 *
 *  Created on: Apr 7, 2013
 *      Author: dferrer
 */

#ifndef PARTICLESTRUCT2_HH_
#define PARTICLESTRUCT2_HH_
#include <cassert>
struct gravdata{
	float3 position;
	float h;
};

class particlestructure{
public:
	long long int count;
	float3 * pos;
	float *h;
	float3 *velocity;
	long long int *id;
	int *address;

	unsigned short *microlevel;

	float *density;
	bool *clean; //if the smoothing length and density need to be re-calculated

	float *mass;
	float *T; //temperature

	ANNidx * NN;
	ANNdist * NNdist;

	particlestructure(int NP){
		pos = new float3[NP];
		h = new float[NP];
		velocity = new float3[NP];
		id = new long long int[NP];
		address = new int[NP];
		microlevel = new unsigned short[NP];
		density = new float[NP];
		clean = new bool[NP];
		mass = new float[NP];
		T = new float[NP];
		NN = new ANNidx *[NP * NSMOOTH];
		NNdist = new ANNdist *[NP *NSMOOTH];
	}

	~particlestructure(){
		delete[] pos;
		delete[] h;
		delete[] velocity;
		delete[] id;
		delete[] address;
		delete[] microlevel;
		delete[] density;
		delete[] clean;
		delete[] mass;
		delete[] T;
		delete[] NN;
		delete[] NNdist;

	}

	ANNidxArray  getNN(int i){
		return (ANNidxArray )(NN + NSMOOTH*i);
	}

	float * getNNdist(int i){
		return (float *) (NNdist + NSMOOTH*i);
	}

	void swap(long long int i, long long int j){
		assert(i < count);
		assert(j < count);

		float3 tpos;
		float th;
		float3 tvelocity;
		long long int tid;
		int taddress;

		unsigned short tmicrolevel;

		float tdensity;
		bool tclean;

		float tmass;
		float tT;
		float tNN;

		tpos = pos[i];
		th = h[i];
		tvelocity = velocity[i];
		tid = id[i];
		taddress = taddress[i];
		tmicrolevel = tmicrolevel[i];
		tdensity = density[i];
		tclean = clean[i];
		tmass = mass[i];
		tT = T[i];
		tNN = NN[i];

		pos[i] = pos[j];
		velocity[i] = velocity[j];
		id[i] = id[j];
		address[i] = address[j];
		microlevel[i] = microlevel[j];
		density[i] = density[j];
		clean[i] = clean[j];
		mass[i] = mass[j];
		tT[i] = T[j];
		NN[i] = NN[j];

		pos[i] = tpos;
		velocity[i] = tvelocity;
		id[i] = tid;
		address[i] = taddress;
		microlevel[i] = tmicrolevel;
		density[i] = tdensity;
		clean[i] = tclean;
		mass[i] = tmass;
		tT[i] = T;
		NN[i] = tNN;

	}

};

#endif /* PARTICLESTRUCT2_HH_ */
