/*
 * reversible.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef REVERSIBLE_HH_
#define REVERSIBLE_HH_

#include "StopWatch.cc"
StopWatch *ctsetup;
StopWatch *nnupdate;
StopWatch *grav;
StopWatch *sph;
StopWatch *kdk;
StopWatch *getrdt;
/*float rdt(std::vector<particlestructure> &p,std::vector<float3> &a){
	float dt = 1000.0;
	for (int i = 0; i < p.size(); i++){
		float dti = .3 * smoothingLength(p[i]) /(cs(p[i]) + sqrt(a[i].norm()*smoothingLength(p[i])));
		if (dt >= dti) dt = dti;
	}
	return dt;
}*/

//A more time reversible variable time step integrator based on hut, 1995. TODO:has no temperature integrator as of yet
float reversibletimestep(std::vector<particlestructure> &p,std::vector<float3> &lastacc, float &lastfs, bool first = false,bool last = false){
	int np = p.size();
	std::vector<float3> acc;


	std::vector<float> dT;
	acc.resize(np);



	float maxa = 0;
	float minh = 4*rad;
	if (first){ //we need to get the initial acceleration for the first timestep
		setupCoverTree(p,sq(40*rad));
		updateNN(p);

	#pragma omp parallel for schedule(dynamic,1)
		for(int i = 0; i < np; i++){
			lastacc[i] = accel(p[i]);
		}


		for(int i = 0; i < np; i++){
			if (p[i].h < minh) minh = p[i].h;
		}

		float fs = minh/(NSMOOTH *5);
		lastfs = fs;

		cudaGrav(p,lastacc,fs);

		destroyCoverTree();
	}
	for(int i = 0; i < np; i++){
		if(lastacc[i].norm() > maxa) maxa = lastacc[i].norm();
	}

	getrdt->Start();
	float lastdt = rdt(p,lastacc);
	getrdt->Stop();

	kdk->Start();
	kick(p,lastacc,lastdt/2);
	drift(p,lastdt);
	kdk->Stop();

	ctsetup->Start();
	setupCoverTree(p,sq(40*rad));
	ctsetup->Stop();

	nnupdate->Start();
	updateNN(p);
	nnupdate->Stop();

	//if(doaccretion(p)){
//		std::cout<<"\nAccreted\n\n";
//		destroyCoverTree();
//		setupCoverTree(p,sq(40*rad));
//		updateNN(p);
//	}
	sph->Start();
	for(int i = 0; i < np; i++){
		acc[i] = accel(p[i]);
	}
	sph->Stop();

	for(int i = 0; i < np; i++){
		if (p[i].h < minh) minh = p[i].h;
	}

	float fs = minh/(NSMOOTH *5);

	grav->Start();
	cudaGrav(p,acc,fs);
	grav->Stop();

	for(int i = 0; i < np; i++){
		if(acc[i].norm() > maxa) maxa = acc[i].norm();
	}
	getrdt->Start();
	float dt1 = rdt(p,acc);
	getrdt->Stop();

	float dt = .5*(lastdt+dt1);

	float delt = dt1-lastdt;
	kdk->Start();
	kick(p,lastacc,dt/2);
	drift(p,delt/2);

	kick(p,lastacc,-lastdt/2);
	kick(p,acc,dt/2);
	kdk->Stop();
	lastacc = acc;
	lastfs = fs;
	destroyCoverTree();
	return dt;

}

#endif /* REVERSIBLE_HH_ */
