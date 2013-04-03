/*
 * pyinterface.hh
 *
 *  Created on: Nov 15, 2012
 *      Author: dferrer
 */

#ifndef PYINTERFACE_HH_
#define PYINTERFACE_HH_


extern "C"{

std::vector<particlestructure> * getPS(int n1d){
	std::vector<particlestructure> * res = new std::vector<particlestructure>;
	res->reserve(cube(n1d));

	int count = makeHomogeneousSphere(*res,n1d,rad);
	std::cout<<"Sphere Count: "<<count<<"\n";

	return res;
}

void freeps(std::vector<particlestructure> *ps){
	delete ps;
}

int getsize(std::vector<particlestructure> *ps){

	return ps->size();
}

void getPos(std::vector<particlestructure> *ps, float * result){
	float3 * res = (float3 *) result;
	int size = ps->size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = ps->at(i).position;
	}
}

void getVel(std::vector<particlestructure> *ps, float * result){
	float3 * res = (float3 *) result;
	int size = ps->size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = ps->at(i).velocity;
	}
}

void setPos(std::vector<particlestructure> *ps, float * input){

	float3 * in = (float3 *) input;
	int size = ps->size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		ps->at(i).position = in[i];
	}


}

void setVel(std::vector<particlestructure> *ps, float * input){

	float3 * in = (float3 *) input;
	int size = ps->size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		ps->at(i).velocity = in[i];
	}

}

void setT(std::vector<particlestructure> *ps, float * input){
	int size = ps->size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		ps->at(i).T = input[i];
	}
}

std::vector<particlestructure> * injectIC(long long int n, float3 * pos, float3 * vel, float * T, float * m){
	std::vector<particlestructure> * ps = new std::vector<particlestructure>;
	std::cout<<"N: "<<n<<"\n";
	std::cout<<"Size: "<<n*sizeof(particlestructure)<<"\n";
	ps->reserve(n);

	for(int i = 0; i< n; i++){
		particlestructure p;
		p.mass = m[i];
		p.position = pos[i];
		p.velocity = pos[i];
		p.T = T[i];
		p.NN.resize(NSMOOTH);
		p.clean = false;
		p.density = 0;
		p.h = 0;
		p.microlevel = 0;
		p.id = i;
		p.address = i;
		ps->push_back(p);
	}
	return ps;
}

void getAcc(std::vector<particlestructure> *ps, float * result){
	float3 * res = (float3 *) result;
	int size = ps->size();
	setupCoverTree(*ps,sq(4*rad));
	updateNN(*ps);
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		res[i] = accel(ps->at(i));
	}
}

void getg(std::vector<particlestructure> *ps, float * result){

	float3 * res = (float3 *) result;
	int size = ps->size();
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

void getT(std::vector<particlestructure> *ps, float * result){
	int size = ps->size();
	setupCoverTree(*ps,sq(40*rad));
	updateNN(*ps);
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
			result[i] = temperature(ps->at(i));
		}
}

void getm(std::vector<particlestructure> *ps, float * result){
	int size = ps->size();
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
			result[i] = ps->at(i).mass;
		}
}

void getP(std::vector<particlestructure> *ps, float * result){
	int size = ps->size();
	setupCoverTree(*ps,sq(4*rad));
	updateNN(*ps);
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
			result[i] = pressure(ps->at(i));

		}
#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++) ps->at(i).clean = false;
	destroyCoverTree();
}

void getdT(std::vector<particlestructure> *ps, float * result){
	int size = ps->size();
		setupCoverTree(*ps,sq(4*rad));
		#pragma omp parallel for schedule(dynamic,1)
			for(int i = 0; i < size; i++){
				result[i] = dTad(ps->at(i));
			}
}

void getDens(std::vector<particlestructure> *ps, float * result){
	int size = ps->size();
	setupCoverTree(*ps,sq(40*rad));
	updateNN(*ps);
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
//			std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(ps->at(i),100);
			result[i] = density(ps->at(i));
		}
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++) ps->at(i).clean = false;
	destroyCoverTree();
}

void geth(std::vector<particlestructure> *ps, float * result){
	int size = ps->size();
	setupCoverTree(*ps,sq(40*rad));
	updateNN(*ps);
	std::cout<<"N: "<<size<<"\n";
	std::cout.flush();
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++){
//			std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(ps->at(i),100);
			result[i] = smoothingLength(ps->at(i));
		}
	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < size; i++) ps->at(i).clean = false;
	destroyCoverTree();
}

void getMicro(std::vector<particlestructure> *ps, unsigned short * result){
	int size = ps->size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		result[i] = ps->at(i).microlevel;
	}
}

float stepforward(std::vector<particlestructure> *ps, bool first = false, bool last = false, float dt = 0){

	return timestep(*ps,first,last,dt);
}

float domicro(std::vector<particlestructure> *ps,float *dt){
	return microstep(*ps,*dt);
}

float rts(std::vector<particlestructure> *ps,bool first, bool last, float3 * lastacc, float * lfs ){//reversible time step
	std::vector<float3> la;
	la.resize(ps->size());
	#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i<ps->size(); i++){
		la[i] = lastacc[i];
	}
	return reversibletimestep(*ps,la, *lfs, first,last);
}
}


#endif /* PYINTERFACE_HH_ */
