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

void propogateLightPacket(particlestructure &p,lightpacket l, float lastdl){

	float3 currentpos = l.position + l.direction *Clight * lastdl;

	float localdensity = density(p,currentpos);
	float localh = smoothingLength(p, currentpos);

	float dl = radiativedl(localh, localdensity);

	float dtau = opticaldepth(dl,opacity(localdensity,temp(localdensity),l.wavelength),l.wavelength);






}


#endif /* MCRT_HH_ */
