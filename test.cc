/*
 * test.cc
 *
 *  Created on: Nov 5, 2012
 *      Author: dferrer
 */
#include "pyinterface.cc"
#include "StopWatch.cc"
#include <cassert>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <fenv.h>


void testCoverTree(){
	std::cout<<"...Starting Cover Tree Test...\n";
	int size = 113078;
	particlestructure particles(size);
	particles.count = size;
	for(int i = 0; i <size; i++ ){
		float3 x (i/(1.0*((float)(size))),0,0);
		particles.pos[i] = x;
	}
	for (int i = 1; i < size; i++){
		assert( particles.pos[i-1].x <  particles.pos[i].x);
		}
	StopWatch setuptime("Cover Tree Setup Time");
	setuptime.Start();
	setupCoverTree(particles, size);
	setuptime.Stop();
	std::cout<<"Cover Tree Setup Time: "<<setuptime.Elapsed()<<" s \n";

	StopWatch runtime("32 Nearest Neighbors Recovery and Iteration Time");
	runtime.Start();
	updateNN(particles);
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < size; i++){
		for(int j = 0; j < NSMOOTH; j++){
			float x = particles.pos[particles.getNN(i)[j]].x;
			float highbound = particles.pos[i].x + (1.01*NSMOOTH)/((float)(size));
			float lowbound = particles.pos[i].x - (1.01*NSMOOTH)/((float)(size));
			if (x >= highbound || x <= lowbound ){
				std::cout<<"Problem! Assertion  "<<lowbound<<" < " <<x <<" < "<<highbound<<" is going to fail.\n";
			}
			assert(x < highbound);
			assert(x > lowbound);
		}
	}
	runtime.Stop();
	std::cout<<"\nRun time: "<<runtime.Elapsed()<<" s\n";
	std::cout<<"Rate: "<< size/runtime.Elapsed()<<" particles/s\n\n";
	std::cout<<"...............Passed CoverTree Test................\n";
}



void testDirect(){
	std::cout<<"\n...Starting Direct Gravity Test...\n";
	long int size = 30000;//113078;
	particlestructure particles(size);
	std::cout<<"Size: "<<size<<"\n";

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();

	T = gsl_rng_default; //The mersene twister
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,time(0));

	std::vector<float3> acc;
        acc.resize(size);

	std::vector<float3> acc2;
        acc2.resize(size);

	for(int i = 0; i < size; i++){
		float3 pos(gsl_rng_uniform(r),gsl_rng_uniform(r),gsl_rng_uniform(r));
		particles.pos[i] = pos;
		particles.mass[i] = 1.0e-8f/size;
		acc[i].x = 0;
		acc[i].y = 0;
		acc[i].z = 0;
		acc2[i].x = 0;
		acc2[i].y = 0;
		acc2[i].z = 0;
	}
	StopWatch runtime("Direct Runtime");
	runtime.Start();
	cudaGrav(particles,acc,.0001);
	runtime.Stop();
	std::cout<<"\nCuda Run time: "<<runtime.Elapsed()<<" s\n";
	std::cout<<"Rate: "<< (double)(size*size)/(pow(10,9)*runtime.Elapsed())<<" Gdirect/s \n\n";
	runtime.Clear();
	runtime.Start();
        directGrav(particles,acc2,.0001);
        runtime.Stop();
        std::cout<<"\nCPU Run time: "<<runtime.Elapsed()<<" s\n";
        std::cout<<"Rate: "<< (double)(size*size)/(pow(10,9)*runtime.Elapsed())<<" Gdirect/s \n\n";

	for(int i = 0; i <size; i++){
		if( (acc[i]-acc2[i]).norm() >= 1e-5){
			printf("Error! Particle %d has cudaAcc: (%f,%f,%f) and cpuAcc: (%f,%f,%f) which differ by: %f. \n",i,acc[i].x,acc[i].y,acc[i].z,acc2[i].x,acc2[i].y,acc2[i].z, (acc[i]-acc2[i]).norm());
			assert((acc[i]-acc2[i]).norm() < 1e-5);
		}
	}
	std::cout<<"...............Passed Direct Test................\n";


}


void testTimeStep(){
	std::cout<<"\n...Starting Single Timestep Test...\n";
	long int size = cube(40);
	particlestructure particles(size);
	std::cout<<"Size: "<<size<<"\n";

	int count = makeHomogeneousSphere(particles,40,rad);
	std::cout<<"Sphere count: "<<count<<"\n";
	StopWatch runtime("Single Timestep");
	runtime.Start();
	float dt  = timestep(particles,true);
	runtime.Stop();
	std::cout<<"dt: "<<dt<<"\n";
	std::cout<<"\nRun time: "<<runtime.Elapsed()<<" s\n";
	std::cout<<"Rate: "<< (double)(count)/(runtime.Elapsed())<<" particle/s \n\n";
	std::cout<<"...............Passed Single Time Step Test................\n";

}

void testrTimeStep(){
	ctsetup = new StopWatch("ctsetup");
	nnupdate = new StopWatch("nnupdate");
	grav = new StopWatch("cudagrav");
	sph = new StopWatch("sph");
	kdk = new StopWatch("kdk");
	getrdt = new StopWatch("rdt");
	int steps = 10;
	std::cout<<"\n...Starting Single Timestep Test...\n";
	long int size = cube(60);
	particlestructure particles(size);

	std::cout<<"Size: "<<size<<"\n";

	int count = makeHomogeneousSphere(particles,60,rad);
	std::cout<<"Sphere count: "<<count<<"\n";
	StopWatch runtime("Single Timestep");
	std::vector<float3> a(count,(float3)(0,0,0));
	float fs =0;
	runtime.Start();
	float dt  = reversibletimestep(particles,a,fs,true,false);
	for(int i = 0; i < steps -1; i++) {
		std::cout<<"Step: "<<1+i<<"\n";
		reversibletimestep(particles,a,fs,false,false);
	}
	runtime.Stop();
	std::cout<<"dt: "<<dt<<"\n";
	std::cout<<"\nRun time: "<<runtime.Elapsed()<<" s\n";
	std::cout<<"Rate: "<< steps *(double)(count)/(runtime.Elapsed())<<" particle/s \n\n";
	printf("Timing Breakdown: \n\tCoverTree Setup: %f\n\tNN Update: %f\n\tGravity: %f\n\tSPH: %f\n\tKDK: %f\n\tChoose TimeStep:%f\n",
			ctsetup->Elapsed(),nnupdate->Elapsed(),grav->Elapsed(),sph->Elapsed(),kdk->Elapsed(),getrdt->Elapsed());


	std::cout<<"...............Passed Single Time Step Test................\n";
	delete ctsetup;
	delete nnupdate;
	delete grav;
	delete sph;
	delete kdk;
	delete getrdt;
}



void testMicro(){
	std::cout<<"\n...Starting Microstep Test...\n";
	int n1d = 60;
	long int size = cube(n1d);
	particlestructure particles(size);
	std::cout<<"Size: "<<size<<"\n";

	int count = makeHomogeneousSphere(particles,n1d,rad);
	std::cout<<"Sphere count: "<<count<<"\n";
	StopWatch runtime("Single Timestep");
	runtime.Start();

	float dt  = timestep(particles,1,0,0);
	float smalldt = dt;
	//for (int i = 0; i < 100; i++) dt  += microstep(particles,smalldt);
	runtime.Stop();
	std::cout<<"Total dt: "<<dt<<"\n";
	std::cout<<"Last Small dt: "<<smalldt<<"\n";
	std::cout<<"\nRun time: "<<runtime.Elapsed()<<" s\n";
	std::cout<<"Rate: "<< 101*(double)(count)/(runtime.Elapsed())<<" particle/s \n\n";
	std::cout<<"...............Passed Single Time Step Test................\n";
}

void pythonTest(){
	int n1d = 20;
	int size = cube(20);
	particlestructure * ps = getPS(20);


	int s = ps->count;
	float * g = new float[3*s];
	for(int i = 0; i < 3*s; i++) g[i] = 0;
	getg(ps,g);
	delete[] g;


	float dt = stepforward(ps,1,0,0.0);
	float totaldt  = dt;
	float smalldt = dt;
	for(int i = 0; i < 30; i++){
		dt = 1;//domicro(ps,&smalldt);
		totaldt += dt;
		std::cout<<"i: "<<"\n";
		std::cout<<"dt: "<<dt<<"\n";
		std::cout<<"microdt: "<<smalldt<<"\n";
		std::cout<<"total time: "<<totaldt<<" ka\n";
	}
	dt = stepforward(ps,0,1);
	totaldt += dt;
	std::cout<<"total time: "<<totaldt<<" ka\n";
	delete ps;

}

int main(){
	feenableexcept(FE_INVALID | FE_DIVBYZERO);
	//testCoverTree();
	testDirect();
	//testrTimeStep();
	//pythonTest();

	//testMicro();

	return 0;

}
