/*
 * ic.hh
 *
 *  Created on: Nov 19, 2012
 *      Author: dferrer
 */

#ifndef IC_HH_
#define IC_HH_

float rad = 3870.0; //AU
int n1d = 20;


#ifndef T0
#define T0 10 //K
#endif

int makeHomogeneousSphere(std::vector<particlestructure> &p, int n1d,float radius){
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
	p.resize(n1d*n1d*n1d);

	for(int i = -n1d/2; i< n1d/2; i++){
		for(int j = -n1d/2; j < n1d/2; j++){
			for(int k = -n1d/2; k < n1d/2; k++){
				float rijk2 = i*i+j*j+k*k;
				rijk2 *= sq(2*radius/(n1d));
				if (rijk2 <= r2) {
					//std::cout<<"(i,j,k): ("<<i<<" , "<<j<<" , "<<k<<")\n";
					particlestructure pijk;
					pijk.position = (2*(float)(i)/((float)(n1d)))*radius *x +(2*(float)(j)/((float)(n1d)))*radius*y+(2*(float)(k)/((float)(n1d)))*radius*z ;
					//std::cout<<"pijk: "<<pijk.position<<"\n";
					pijk.velocity = omega.cross(pijk.position);
					pijk.id = particlesadded;
					pijk.address = particlesadded;
					pijk.microlevel = 0;
					pijk.h = 0;
					pijk.density = 0;



					pijk.clean = false;
					pijk.T = T0;


					p[particlesadded] = pijk;
					p[particlesadded].NN.resize(NSMOOTH);
					float myphi =2*sq(pijk.position.z/pijk.position.norm())-1;
					if (pijk.position.norm() == 0) myphi = 0;
					phi.push_back(myphi);
					particlesadded++;
				}
			}

		}
	}
	p.resize(particlesadded);

#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < particlesadded; i ++){
		p[i].mass = 50.0/((float)(particlesadded)) *(1+ .01*phi[i]);
		if (p[i].position.norm2() ==0) p[i].mass = 1;
	}


	return particlesadded;
}



#endif /* IC_HH_ */
