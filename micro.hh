/*
 * micro.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef MICRO_HH_
#define MICRO_HH_

#define MINMICROLEVEL 3
#define MAXMICROLEVEL (10+ MINMICROLEVEL)

void microupdateNN(std::vector<particlestructure *> &p,std::vector<particlestructure > &pf){
#pragma omp parallel for schedule(dynamic,1)
        for(int i = 0; i < p.size(); i++){
                std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(*p[i],NSMOOTH);
                for(int j = 0; j < NSMOOTH; j++)
                        p[i]->NN[j] =  &(pf[NN[j].address]);
        }


}

void microGrav(std::vector<particlestructure *>  &p1, std::vector<particlestructure> &p2, std::vector<float3> &acc,float cfs = 0){
	assert(cfs ==cfs);
	int p1size = p1.size();
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < p1size; i++) {

			float fs;
			if (cfs == 0) fs = sq(p1[i]->h /(10 *NSMOOTH));
			else fs = cfs;

			if(!(acc[p1[i]->address] == acc[p1[i]->address])){
				std::cout<<"address: "<<p1[i]->address<<"\n";
				std::cout<<"p[i]"<<p1[i]->mass<<"\n";
			}

			assert(acc[p1[i]->address] == acc[p1[i]->address]);
			for (int j = 0; j < p2.size(); j++){
				float3 xij = (p1[i]->position-p2[j].position);
				float rij = xij.norm();
				assert(rij ==rij);

				acc[p1[i]->address] += -G *p2[j].mass * (pow(sq(rij)+sq(fs),-3.0/2.0)) * xij;
				if(!(acc[p1[i]->address] == acc[p1[i]->address])){
					std::cout<<"address: "<<p1[i]->address<<"\n";
					std::cout<<"rij: "<<rij<<"\n";
					std::cout<<"xij: "<<xij<<"\n";
					std::cout<<"fs :"<<fs<<"\n";
					std::cout<<"p2 mass"<<p2[j].mass<<"\n";
					std::cout<<"result: "<< -G *p2[j].mass * (pow(sq(rij)+sq(fs),-3.0/2.0)) * xij<<"\n";
				}
				assert(acc[p1[i]->address] == acc[p1[i]->address]); //check for NaN's
			}

	}
}


//Setup the microstepping levels and lists levellists should be an array of arrays of pointers to particle structure elements, each of length np
//levelcounts is an array holding the number of particles in each microstepping level
unsigned short setupMicro(std::vector<particlestructure> &p, std::vector<float3> &acc,/*std::vector<float> &dT,*/float fs, float globaldt, std::vector<std::vector<particlestructure *> > &levellists){
	long int np = p.size();

	for(int i = 0; i < np; i++){
		float a = acc[i].norm();
		float dt = getdt(p[i],acc[i]);

//		float dtT = std::abs(.3 *p[i].T/dT[i]);
//		if(p[i].T >0 && dtT <dt) dt = dtT;
		float ratio = dt/globaldt;
		assert(ratio >=1.0);

		unsigned short level = floor(log2(ratio))+MINMICROLEVEL;
		if (level >MAXMICROLEVEL) level = MAXMICROLEVEL;

		p[i].microlevel = level;
		levellists[level].push_back(&(p[i]));
	}

	unsigned short maxlevel;
	long int totalcount = 0;
	for(maxlevel = 0; maxlevel < MAXMICROLEVEL; maxlevel++){
		totalcount+= levellists[maxlevel].size();
		if(((float)(totalcount))/np > .9) break;

	}
	//put the particles on the un-needed levels on to the max level list

	for(int l = maxlevel+1; l <= MAXMICROLEVEL; l++){
		for(int i = 0; i < levellists[l].size(); i++)
			levellists[maxlevel].push_back(levellists[l][i]);
		levellists[l].clear();
	}
	return maxlevel;
}

//********************* New KDK operators for microstepping******************************
void microkick(std::vector<particlestructure *>  &p, std::vector<float3> &acc, float dt){
	//#pragma omp parallel for schedule(dynamic, 1)
		for ( int i = 0; i < p.size(); i++){
			assert(p[i]->velocity == p[i]->velocity); //check for nans;

			int idx = p[i]->address;
			assert(acc[idx] == acc[idx]);
			p[i]->velocity += dt* acc[idx];
		}
}

void microdrift(std::vector<particlestructure *> &p, float dt){
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.size(); i++){
		if(!(p[i]->velocity == p[i]->velocity)) std::cout<<p[i]->address;
		assert(p[i]->velocity == p[i]->velocity);
		if(p[i]->position.norm() > 2*rad && p[i]->position.dot(p[i]->velocity) >0) p[i]->velocity *=-1;
		p[i]->position += dt* p[i]->velocity;
		p[i]->clean = false;
	}
}

void microupdateT(std::vector<particlestructure *> &p,std::vector<float> &dT,float dt){
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.size(); i++){
		int idx = p[i]->address;
		p[i]->T += dt* dT[idx];
		//assert(p[i]->T >0);
		if(p[i]->T <0) p[i]->T = 0;
	}
}

void fullkick(std::vector<particlestructure>  &p, std::vector<float3> &acc, float dt){
	#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < acc.size(); i++){
		assert(p[i].velocity == p[i].velocity); //check for nans;
		assert(acc[i] == acc[i]);
		float mydt = dt * (pow(2,p[i].microlevel-MINMICROLEVEL));
		p[i].velocity += mydt* acc[i];
	}
}

void fulldrift(std::vector<particlestructure> &p, float dt){
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.size(); i++){
		assert(p[i].velocity == p[i].velocity);
		float mydt = dt * (pow(2,p[i].microlevel-MINMICROLEVEL));;
		if(p[i].position.norm() > 2*rad) p[i].velocity *=-1;
		p[i].position += mydt* p[i].velocity;
		p[i].clean = false;
	}
}

void fullupdateT(std::vector<particlestructure> &p,std::vector<float> &dT,float dt){
	assert(p.size() == dT.size());
#pragma omp parallel for schedule(dynamic, 1)
	for ( int i = 0; i < p.size(); i++){
		float mydt = dt * (pow(2,p[i].microlevel-MINMICROLEVEL));;
		p[i].T += mydt* dT[i];
		if(p[i].T <0) p[i].T = 0;
	}
}

//***************************************************************************************

void recalcMicrolevel(int level, std::vector<std::vector<particlestructure *> > &levellists,std::vector<float3> &acc,/*std::vector<float> &dT,*/float fs, float globaldt, unsigned short maxlevel){
	std::vector<int> toErase;
	for(int i = 0; i <levellists[level].size(); i++){
		particlestructure * p = levellists[level][i];

		float a = acc[p->address].norm();
		float dt = getdt(*p,acc[p->address]);
//		float dtT = std::abs(.3 *p->T/dT[p->address]);

//		if(p->T >0 && dtT <dt) dt = dtT;
		float ratio = dt/globaldt;

		double lev = floor(log2(ratio))+MINMICROLEVEL;
		if(!(lev >= 0)){
			std::cerr<<"Assertion is going to fail! \n";
			std::cerr<<"lev: "<<lev <<" < 0\n";
//			std::cerr<<"ratio: "<<ratio<<"\n";
			std::cerr<<"dt: "<<dt<<"\n";
			std::cerr<<"global dt: "<<globaldt<<"\n";
//			std::cerr<<"Temp dt: "<<dtT<<"\n";
			std::cerr<<"T: "<<p->T<<"\n";
//			std::cerr<<"dT: "<<dT[p->address]<<"\n\n";
			std::cerr<<"address: "<<p->address<<"\n";
			std::cerr<<"level: "<<level<<"\n";
			std::cerr<<"i: "<<i<<"\n";
			std::cerr<<"a: "<<a<<"\n";
			std::cerr<<"size: "<< levellists[level].size()<<"\n";

		}
		if (lev < 0){
			std::cout<<"A particle has fallen below the lowest micro level!\n";
			lev = 0;
		}
		unsigned short l = lev;


		if(l < level){
			p->microlevel = l;
			toErase.push_back(i);
			if (l > maxlevel) l = maxlevel;
			levellists[l].push_back(p);
		}
	}

	//now we must remove the departed particles from the current level list
	//We do this by swapping the particles that must be removed to the end of the level list,
	//then calling erase on the back of the list
	int nerase = toErase.size();
	int swapstart = levellists[level].size() - nerase;
	//do the swap
	for(int i = 0; i < nerase; i++){
		int idx = toErase[i];
		particlestructure * temp = levellists[level][idx];
		levellists[level][idx] = levellists[level][swapstart +i];
		levellists[level][swapstart +i] = temp;
	}
	//erase the back of the list
	levellists[level].erase(levellists[level].begin()+swapstart,levellists[level].end());
}


//do a timestep with microstepping. returns the full timestep time, and puts the microstep time in lastdt
float microstep(std::vector<particlestructure> &p,float &lastdt){

	int np = p.size();
	std::vector<float3> acc;


//	std::vector<float> dT;
	acc.resize(np);
//	dT.resize(np);

	setupCoverTree(p,sq(4*rad));
	updateNN(p);
	float maxa = 0;
	float minh = 4*rad;
#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < np; i++){
		acc[i] = accel(p[i]);
	}

	for(int i = 0; i < np; i++){
		assert(p[i].h >0);
		if (p[i].h < minh) minh = p[i].h;
	}
//	if(minh < FS) minh = FS;
	float fs = 2*minh/(NSMOOTH);
	assert(fs == fs);

	directGrav(p,acc,fs);


	for(int i = 0; i < np; i++){
		if(acc[i].norm() > maxa) maxa = acc[i].norm();
	}
	float dt = rdt(p,acc);

	//calculate the temperature change
//	#pragma omp parallel for schedule(dynamic,1)
//	for(int i = 0; i < np; i++){
//		dT[i] = dTad(p[i]);

//	}
	//calculate the maximum temperature time step
//	float dtT = dt;
//	for (int i = 0; i < np; i++){
//		if (p[i].T > 0 && std::abs(.3 *p[i].T/dT[i]) < dtT) dtT = std::abs(.3 *p[i].T/dT[i]);
//	}
	float dta = dt;
//	if (dtT< dt) dt = dtT;
//	std::cout<<"dtT: "<<dtT<<"\n";
	std::cout<<"dta: "<<dta<<"\n";
	//set up the microstepping level lists
	std::vector<std::vector<particlestructure *> > levellists;
	levellists.resize(MAXMICROLEVEL+1);
	for(int i = 0; i <= MAXMICROLEVEL; i++) levellists[i].reserve(np);


	unsigned short maxmicrolevel = setupMicro(p,acc,fs, dt,levellists);


	//do the first timestep for every particle
	kick(p,acc,lastdt/2);
	fullkick(p,acc,dt/2);
	fulldrift(p,dt);
	//fullupdateT(p,dT,dt);
	destroyCoverTree();
	setupCoverTree(p,sq(40*rad));

	int timestepcount = 1<<maxmicrolevel;
	std::cout<<"Beginning microstepping with "<<maxmicrolevel-MINMICROLEVEL<<" levels at dt: "<<dt<<"...\n";
	for (int i = 0; i < levellists.size(); i++) std::cout<<"Level "<<i<<": "<<levellists[i].size()<<"\n";

	//std::cout.flush();
	std::cout<<"LL: "<<levellists.size()<<"\n";


	for(int i = 1; i < timestepcount; i++){
		//feenableexcept(FE_DIVBYZERO|FE_INVALID);
		std::cout<<"Step: "<<i<<"\n";
		std::cout.flush();
		for(int l = maxmicrolevel; l >= 0; l-=1){

			if(i % (1<<l) ==0){
				int npl = levellists[l].size();
				if (npl != 0) {
					std::cout<<"Doing micro for level: "<<l<<"\n";
//					std::vector<float3> lacc;
//					std::vector<float> ldT;
//					lacc.resize(npl);
//					ldT.resize(npl);
					std::cout<<"Getting accelerations\n";
					#pragma omp parallel for schedule(dynamic,1)
					for( int j = 0; j <levellists[l].size(); j++){
						acc[levellists[l][j]->address] = accel(*levellists[l][j]) ;
						assert(acc[levellists[l][j]->address] == acc[levellists[l][j]->address]);
//						dT[levellists[l][j]->address] = dTad(*levellists[l][j]);
					}
					std::cout<<"Doing micrograv\n";
					microGrav(levellists[l],p,acc,fs);
					float ldt = dt* pow(2,l-MINMICROLEVEL);
					assert(ldt == ldt);




					recalcMicrolevel(l, levellists,acc,fs, dt,maxmicrolevel);
					for( int j = 0; j <levellists[l].size(); j++) {
						assert(acc[levellists[l][j]->address] == acc[levellists[l][j]->address]);
						assert(levellists[l][j]->velocity == levellists[l][j]->velocity);
					}
					std::cout<<"New level occupancy:\n";
					for (int ix = 0; ix < levellists.size(); ix++) std::cout<<"Level "<<ix<<": "<<levellists[ix].size()<<"\n";
					//we need to remove the particles from the covertree since they are about to change
					std::cout<<"Changing covertree\n";
					for( int j = 0; j <levellists[l].size(); j++) particlecovertree->remove(*levellists[l][j]);


					for( int j = 0; j <levellists[l].size(); j++) {
											assert(acc[levellists[l][j]->address] == acc[levellists[l][j]->address]);
											assert(levellists[l][j]->velocity == levellists[l][j]->velocity);
											assert(levellists[l][j]->position == levellists[l][j]->position);
					}
					std::cout<<"Doing kdk\n";
					microkick(levellists[l],acc,ldt);
					microdrift(levellists[l],ldt);
//					microupdateT(levellists[l],dT,ldt);
					std::cout<<"unlceaning\n";
					#pragma omp parallel for schedule(dynamic,1)
					for(int j = 0; j <p.size(); j++) p[i].clean = false;
					std::cout<<"re-adding to covertree\n";
					//re-add the particles to the covertree
					for( int j = 0; j <levellists[l].size(); j++) particlecovertree->insert(*levellists[l][j]);
					microupdateNN(levellists[l],p);
//					destroyCoverTree();
//					setupCoverTree(p,sq(40*rad));

				}
			}

		}
	}
	//std::cout<<"Done Microstepping\n";
	//std::cout.flush();

	//we need to sync the time of every particle now by doing a fullkick again, then a constant timestep kick back dt/2
	fullkick(p,acc,dt/2);
	kick(p,acc,-dt/2);


	destroyCoverTree();
	lastdt = dt;
	return dt *(1<<(maxmicrolevel-MINMICROLEVEL));


}


#endif /* MICRO_HH_ */
