/*
 * gravity.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef GRAVITY_HH_
#define GRAVITY_HH_

#ifndef G
#define G 3.942 *pow(10,7) //G *Msun*(1 ka)^2/(1 au)^3. change for other units
#endif
#ifndef FS
#define FS 1.0
#endif

//get the gravitational force on each particle via the direct double sum. Set cfs to be non zero to force a global softening length.
//level is the maximum micro-stepping  level that the calculation should be done for. This is currently always less than 10, so omitting this
//argument will do the calculation for every particle in p
//TODO: Write this using trees or move to the GPU
void directGrav(std::vector<particlestructure> &p, std::vector<float3> &acc,float cfs = 0){
	assert(p.size() == acc.size());
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < p.size(); i++) {
		{
			float fs;
			if (cfs == 0) fs = sq(p[i].h /(10 *NSMOOTH));
			else fs = cfs;
			for (int j = 0; j < p.size() ; j++){
				float3 xij = (p[i].position-p[j].position);
				float rij = xij.norm2();
				float ex = (float)(-3.0/2.0);
				acc[i] += -G *p[j].mass * (std::pow(rij+sq(fs),ex)) * xij;
				/*if (acc[i].x >10000000){
					std::cout<<"h:"<<fs<<"\n";
					std::cout<<"j: "<<j<<"\n";
					std::cout<<"pos i: "<<p[i].position<<"\n";
					std::cout<<"pos j: "<<p[j].position<<"\n";
					std::cout<<"xij: "<<xij<<"\n";
					std::cout<<"id i"<<p[i].id<<"\n";
					std::cout<<"id j"<<p[j].id<<"\n";
					assert(1==0);
				}*/
			}
			assert(acc[i]==acc[i]); //check for NaNs
		}
	}
}


#endif /* GRAVITY_HH_ */
