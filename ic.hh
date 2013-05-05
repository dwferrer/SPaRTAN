/*
 * ic.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef IC_HH_
#define IC_HH_

const float rad = 3870.0; //AU


#ifndef T0
#define T0 10 //K
#endif

#ifndef M0
#define M0 (50*1e6)

int makeHomogeneousSphere(particlestructure &p, int n1d,float radius){
	int particlesadded = 0;
	float r2 = radius*radius;
	std::vector<float> phi;
	phi.reserve(cube(n1d));
	float3 x,y,z;
	x.x = 1.0;
	x.y = 0;
	x.z = 0;
	y.x = 0;
	y.y = 1;
	y.z = 0;
	z.x = 0;
	z.y = 0;
	z.z  =1;
	float3 omega = .001*z;

	for(int i = -n1d/2; i< n1d/2; i++){
		for(int j = -n1d/2; j < n1d/2; j++){
			for(int k = -n1d/2; k < n1d/2; k++){
				float rijk2 = i*i+j*j+k*k;
				rijk2 *= sq(2*radius/(n1d));
				if (rijk2 <= r2) {
					//std::cout<<"(i,j,k): ("<<i<<" , "<<j<<" , "<<k<<")\n";

					p.pos[particlesadded] = (2*(float)(i)/((float)(n1d)))*radius *x +(2*(float)(j)/((float)(n1d)))*radius*y+(2*(float)(k)/((float)(n1d)))*radius*z ;
					//std::cout<<"p: "<<p.pos<<"\n";
					p.velocity[particlesadded] = omega.cross(p.pos[particlesadded]);
					p.id[particlesadded] = particlesadded;
					p.address[particlesadded] = particlesadded;
					p.microlevel[particlesadded] = 0;
					p.h[particlesadded] = 0;
					p.density[particlesadded] = 0;



					p.clean[particlesadded] = false;
					p.T[particlesadded] = T0;



					float myphi;
					if (p.pos[particlesadded].norm() == 0) myphi = 0;
					else myphi =2*sq(p.pos[particlesadded].z/p.pos[particlesadded].norm())-1;

					phi.push_back(myphi);
					particlesadded++;
				}
			}

		}
	}
	p.count = particlesadded;

#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < particlesadded; i ++){
		p.mass[i] = M0/((float)(particlesadded));//*(1+ .01*phi[i]);
		//if (p[i].pos.norm2() ==0) p[i].mass = 1;
	}


	return particlesadded;
}



#endif /* IC_HH_ */
