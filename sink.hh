/*
 * sink.hh
 *
 *  Created on: Nov 28, 2012
 *      Author: dferrer
 */

#ifndef SINK_HH_
#define SINK_HH_

bool isBound(particlestructure &p1, particlestructure &p2){ //check if p2 is gravitationally bound to p1
	float ke = .5 * p2.mass *(p2.velocity-p1.velocity).norm2();
	float pe = G * p1.mass*p2.mass/(p2.position-p1.position).norm();
	return ke <=pe;
}

bool momentumcheck(particlestructure &p1, particlestructure &p2){ //check if the specific angular momentum of p2 is low enough to be accreted
	float h = ((p2.position-p1.position).cross(p2.velocity-p1.velocity)).norm();
	float h0 = sqrt(1*G*(p1.mass+p2.mass)); //specific angular momentum of circular orbit at 1 AU
	return h <= h0;
}

#define rhoa .000067 //the accretion density

bool doaccretion(std::vector<particlestructure> &p){
	bool accreted = false;
	for (int i = 0; i < p.size(); i++){
		if(density(p[i])> rhoa) p[i].T = 0; //turn the particle into a sink particle
		if(p[i].T == 0){
			std::vector<long int> toAccrete;
			toAccrete.reserve(500);
			std::vector<particlestructure> NN = particlecovertree->kNearestNeighbors(p[i],500);
			for(int j = 0; j < NN.size(); j++){ //gather particles to be accreted
				if(NN[j].id != p[i].id && isBound(p[i],NN[j]) && momentumcheck(p[i],NN[j])){
					accreted = true;
					toAccrete.push_back(NN[j].address);
				}
			}

			//add the accreted particle to the sink particle;
			float3 cmpos = float3(0,0,0);
			float3 cmvel = float3(0,0,0);
			double totalmass = p[i].mass;
			cmpos+= p[i].mass * p[i].position;
			cmvel+= p[i].mass * p[i].velocity;
			for(int j = 0; j < toAccrete.size(); j++){
				cmpos+=p[toAccrete[j]].mass* p[toAccrete[j]].position;
				cmvel+= p[toAccrete[j]].mass * p[toAccrete[j]].velocity;
				totalmass+= p[toAccrete[j]].mass;
			}
			cmpos/= totalmass;
			cmvel /= totalmass;
			p[i].mass = totalmass;
			p[i].velocity = cmvel;
			p[i].position = cmpos;

			//perform swap
			int nerase = toAccrete.size();
			int swapstart = p.size() - nerase;
			//do the swap
			for(int i = 0; i < nerase; i++){
				int idx = toAccrete[i];
				particlestructure temp = p[idx];
				p[idx] = p[swapstart +i];
				p[idx].address = idx;
				p[swapstart +i] = temp;
				p[swapstart + i].address = swapstart+i;
			}

			p.erase(p.begin()+swapstart,p.end());




		}
	}
	return accreted;
}



#endif /* SINK_HH_ */
