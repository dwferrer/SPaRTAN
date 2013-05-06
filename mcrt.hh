/*
 * mcrt.hh
 * Monte-carlo radiative transfer
 *  Created on: Apr 5, 2013
 *      Author: dferrer
 */

#ifndef MCRT_HH_
#define MCRT_HH_

#define Clight 63239726.3451

struct lightpacket{
	float3 position;
	float3 direction;
	float tau;
	float wavelength; //in micron
};

float temp(float density){
	return 10;
}

float opacity(float density, float T, float wavelength){
	return .1; //fixme: place holder!!!!
}


float scattering(float density, float T, float wavelength){
	return .1; //fixme: place holder!!!!
}

float opticaldepth(float pathlength, float opacity, float density){
	return opacity * density * pathlength;
}

float radiativedl(float h, float rho){
	return h/2;
}

void propogateLightPacket(particlestructure &p,lightpacket &l, float &lastdl){

	l.position += .5 *l.direction * lastdl;

	float localdensity = density(p,l.position);
	float localh = smoothingLength(p, l.position);

	float dl = radiativedl(localh, localdensity);

	float dtau = opticaldepth(.5 *(dl + lastdl),opacity(localdensity,temp(localdensity),l.wavelength),l.wavelength);

	l.tau -=dtau;

	l.position += .5 *l.direction *dl;

}


#endif /* MCRT_HH_ */
