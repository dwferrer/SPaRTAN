/*
 * particlestructure.h
 *
 *  Created on: Nov 4, 2012
 *      Author: dferrer
 */

#ifndef PARTICLESTRUCTURE_H_
#define PARTICLESTRUCTURE_H_

#include "threevector.hh"
class particlestructure{
public:
	float3 position;
	float3 velocity;
	long int id;
	long int address;

	unsigned short microlevel;

	float h;
	float density;
	bool clean; //if the smoothing length and density need to be re-calculated

	float mass;
	float T; //temperature

	 std::vector <particlestructure *> NN;


	float distance (const particlestructure &p) const;
	bool operator==(const particlestructure &p) const;
	void print();
};

 float particlestructure::distance(const particlestructure &p) const{
	return (position-p.position).norm2();
}

bool particlestructure::operator==(const particlestructure &p) const{
	return position==p.position;
}


#endif /* PARTICLESTRUCTURE_H_ */
