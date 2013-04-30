/*
 * pyinterface.hh
 *
 *  Created on: Nov 15, 2012
 *      Author: dferrer
 */





#ifndef PYINTERFACE_HH_
#define PYINTERFACE_HH_

#include "sphANN.hh"
#include "gravity.hh"
#include "ic.hh"
#include "sink.hh"
#include "timestep.hh"
#include "micro.hh"
#include "reversible.hh"
#include "mcrt.hh"

extern "C"{

particlestructure* getPS(int n1d){
	particlestructure * res = new particlestructure(n1d);

	int count = makeHomogeneousSphere(*res,n1d,rad);
	std::cout<<"Sphere Count: "<<count<<"\n";

	return res;
}

void freeps(particlestructure *ps){
	delete ps;
}

int getsize(particlestructure *ps){

	return ps->count;
}

void getPos(particlestructure *ps, float * result){
	float3 * res = (float3 *) result;
	int size = ps->count;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = ps->pos[i];
	}
}

void getVel(particlestructure *ps, float * result){
	float3 * res = (float3 *) result;
	int size = ps->count;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = ps->velocity[i];
	}
}

void setPos(particlestructure *ps, float * input){

	float3 * in = (float3 *) input;
	int size = ps->count;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		ps->velocity[i] = in[i];
	}


}

void setVel(particlestructure *ps, float * input){

	float3 * in = (float3 *) input;
	int size = ps->count;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		ps->velocity[i] = in[i];
	}

}

void setT(particlestructure * ps, float * input){
	int size = ps->count;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		ps->T[i] = input[i];
	}
}

particlestructure* injectIC(long long int n, float3 * pos, float3 * vel, float * T, float * m){
	particlestructure* p = new particlestructure(n);
	std::cout<<"N: "<<n<<"\n";
	std::cout<<"Size: "<<n*sizeof(particlestructure) *n<<"\n";

	for(int i = 0; i< n; i++){
		p->mass[i] = m[i];
		p->pos[i] = pos[i];
		p->velocity[i] = pos[i];
		p->T[i] = T[i];
		p->clean[i] = false;
		p->density[i] = 0;
		p->h[i] = 0;
		p->microlevel[i] = 0;
		p->id[i] = i;
		p->address[i] = i;
	}
	p->count = n;
	return p;
}

void getAcc(particlestructure *ps, float * result){
	float3 * res = (float3 *) result;
	int size = ps->count;
	setupCoverTree(*ps,size);
	updateNN(*ps);
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = accel(*ps,i);
	}
}

void getg(particlestructure *ps, float * result){

	float3 * res = (float3 *) result;
	int size = ps->count;
	std::vector<float3> g;
	g.resize(size);
	for (int i = 0; i < g.size(); i++){
		g[i].x = 0;
		g[i].y = 0;
		g[i].z = 0;
	}

	directGrav(*ps,g,FS);

#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = g[i];
	}

}

void getT(particlestructure *ps, float * result){
	int size = ps->count;
	setupCoverTree(*ps,size);
	updateNN(*ps);
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
			result[i] = temperature(*ps,i);
		}
	destroyCoverTree();
}

void getm(particlestructure *ps, float * result){
	int size = ps->count;
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
			result[i] = ps->mass[i];
		}
}

void getP(particlestructure *ps, float * result){
	int size = ps->count;
	setupCoverTree(*ps,size);
	updateNN(*ps);
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
			result[i] = pressure(*ps,i);

		}
#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++) ps->clean[i] = false;
	destroyCoverTree();
}

void getdT(particlestructure *ps, float * result){
	int size = ps->count;
		setupCoverTree(*ps,size);
		#pragma omp parallel for schedule(dynamic,1)
			for(int i = 0; i < size; i++){
				result[i] = dTad(*ps,i);
			}
}

void getDens(particlestructure *ps, float * result){
	int size = ps->count;
	setupCoverTree(*ps,size);
	updateNN(*ps);
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
//			std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(ps->at(i),100);
			result[i] = density(*ps,i);
		}
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++) ps->clean[i] = false;
	destroyCoverTree();
}

void geth(particlestructure *ps, float * result){
	int size = ps->count;
	setupCoverTree(*ps,size);
	updateNN(*ps);
	std::cout<<"N: "<<size<<"\n";
	std::cout.flush();
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
//			std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(ps->at(i),100);
			result[i] = smoothingLength(*ps,i);
		}
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++) ps->clean[i] = false;
	destroyCoverTree();
}

void getMicro(particlestructure *ps, unsigned short * result){
	int size = ps->count;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		result[i] = ps->microlevel[i];
	}
}

float stepforward(particlestructure *ps, bool first = false, bool last = false, float dt = 0){

	return timestep(*ps,first,last,dt);
}

float domicro(std::vector<particlestructure> *ps,float *dt){
	return 0;//microstep(*ps,*dt);
}

float rts(particlestructure *ps,bool first, bool last, float3 * lastacc, float * lfs ){//reversible time step
	ctsetup = new StopWatch("ctsetup");
	nnupdate = new StopWatch("nnupdate");
	grav = new StopWatch("cudagrav");
	sph = new StopWatch("sph");
	kdk = new StopWatch("kdk");
	getrdt = new StopWatch("rdt");
	std::vector<float3> la;
	la.resize(ps->count);
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i<ps->count; i++){
		la[i] = lastacc[i];
	}

	float result = reversibletimestep(*ps,la, *lfs, first,last);
	delete ctsetup;
	delete nnupdate;
	delete grav;
	delete sph;
	delete kdk;
	delete getrdt;
	return result;
}
}


#endif /* PYINTERFACE_HH_ */
