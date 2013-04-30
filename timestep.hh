/*
 * timestep.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef TIMESTEP_HH_
#define TIMESTEP_HH_


float rdt(particlestructure &p,std::vector<float3> &a){
        float dt = 1000.0;
#pragma omp parallel for schedule(dynamic,1)
        for (int i = 0; i < p.count; i++){
        		float temp = density(p,i);
        		float hi = smoothingLength(p,i);
                float dti = .3 * hi /(cs(p,i) + sqrt(a[i].norm()*hi));
                if (dt >= dti) dt = dti;
        }
        return dt;
}

void kick (particlestructure &p,std::vector<float3> &acc,float dt){
	assert(p.count == acc.size());
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.count; i++){
		assert(p.velocity[i] == p.velocity[i]); //check for nans;
		assert(acc[i] == acc[i]);
		p.velocity[i] += dt* acc[i];
	}
}

void updateT(particlestructure &p,std::vector<float> &dT,float dt){
	assert(p.count == dT.size());
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.count; i++){
		p.T[i] += dt* dT[i];
		if(p.T[i] <0) p.T[i] = 0;
	}
}

void drift(particlestructure &p, float dt){
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.count; i++){
		assert(p.velocity[i] == p.velocity[i]);
		p.position[i] += dt* p.velocity[i];
		p.clean[i] = false;
	}
}

float getdt(particlestructure &p,int i, float3 &acc){

	float result = .3 * smoothingLength(p,i) /(cs(p,i) + sqrt(acc.norm()*smoothingLength(p,i)));
	return result;

}

float timestep(particlestructure &p, bool first = false,bool last = false, float lastdt = 0.0){
	int np = p.count;
	std::vector<float3> acc;


//	std::vector<float> dT;
	acc.resize(np);
//	dT.resize(np);

	setupCoverTree(p,p.count);
	updateNN(p);
	float maxa = 0;
	float minh = 4*rad;

#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < np; i++){
		acc[i] = accel(p,i);
	}

	for(int i = 0; i < np; i++){
		if (p.h[i] < minh) minh = p.h[i];
	}

	float fs = minh/(NSMOOTH *5);

	directGrav(p,acc,fs);


	for(int i = 0; i < np; i++){
		if(acc[i].norm() > maxa) maxa = acc[i].norm();
	}

//#pragma omp parallel for schedule(dynamic,1)
//	for(int i = 0; i < np; i++){
//		dT[i] = dTad(p[i]);
//
//	}

	float dt = rdt(p,acc);


//	float dtT = dt; //maximum temperature timestep

//	for (int i = 0; i < np; i++){
//		if ( std::abs(.5 *p[i].T/dT[i]) < dtT) dtT = std::abs(.5 *p[i].T/dT[i]);
//	}


//	if (dtT < dt) dt = dtT;
//	updateT(p,dT,dt);


	kick(p,acc,dt/2+lastdt/2);
	drift(p,dt);

	if (last) kick(p,acc,dt/2);
//	std::cout<<"dT 0:" <<dT[0]<<"\n";
//	std::cout<<"dT 0:" <<dT[1]<<"\n";
	destroyCoverTree();
	return dt;
}

#endif /* TIMESTEP_HH_ */
