/*
 * gravity.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef GRAVITY_HH_
#define GRAVITY_HH_

#ifndef G
#define G (39.42f) //G *10(-6)Msun*(1 ka)^2/(1 au)^3. change for other units
#endif
#ifndef FS
#define FS 1.0
#endif

//get the gravitational force on each particle via the direct double sum. Set cfs to be non zero to force a global softening length.
//level is the maximum micro-stepping  level that the calculation should be done for. This is currently always less than 10, so omitting this
//argument will do the calculation for every particle in p
//TODO: Write this using trees
void directGrav(particlestructure &p, std::vector<float3> &acc,float cfs = 0){
	assert(p.count == acc.size());
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < p.count; i++) {
		{
			float fs;
			if (cfs == 0) fs = sq(p.h[i]);
			else fs = cfs;
			for (int j = 0; j < p.count ; j++){
				float3 xij = (p.pos[i]-p.pos[j]);
				float rij = xij.norm2();
				float ex = (float)(-3.0/2.0);
				acc[i] += -G *p.mass[j] * (std::pow(rij+fs,ex)) * xij;
			}
			assert(acc[i]==acc[i]); //check for NaNs
		}
	}
}

void gpugravity(float * pos, float *accel, long long int N); //provided by gravity.cu

void cudaGrav(particlestructure &p, std::vector<float3> &acc,float cfs = 0){

	long long  int np = p.count;
	int leftover =1024 -( np%1024); //FIXME:This should actually be the number of threads
	//first we need to copy the particles into good cuda order
	float * pos = new float[(np+leftover)*4];
	float *a = new float[(np+leftover)*4];

	#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < np; i++){
		float fs;
		if (cfs == 0) fs = sq(p.h[i]);
		else fs = cfs;
		pos[4*i + 0] = p.pos[i].x;
		pos[4*i + 1] = p.pos[i].y;
		pos[4*i + 2] = p.pos[i].z;
		pos[4*i + 3] = p.mass[i];
		a[4*i +0] = acc[i].x/G;
		a[4*i +1] = acc[i].y/G;
		a[4*i +2] = acc[i].z/G;
		a[4*i +3] = fs;
	}
	//we need to pad the remaining data for the cuda code so it doesn't give spurious answers
	for(int i = np; i <np + leftover; i++){
		pos[4*i + 0] = 0;
                pos[4*i + 1] = 0;
                pos[4*i + 2] = 0;
                pos[4*i + 3] = 0;
                a[4*i +0] = 0;
                a[4*i +1] = 0;
                a[4*i +2] = 0;
                a[4*i +3] = 1; //avoids an overflow in divide
	}
	gpugravity(pos,a,np+leftover);

#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < np; i++){
		acc[i].x = G*a[4*i +0];
		acc[i].y = G*a[4*i +1];
		acc[i].z = G*a[4*i +2];
	}




}
#endif /* GRAVITY_HH_ */
