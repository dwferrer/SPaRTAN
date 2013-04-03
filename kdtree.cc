/*
 * kdtree.cc
 *
 *  Created on: Nov 4, 2012
 *      Author: dferrer
 */

#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cfloat>

#include "particlestructure.h"

#include <threevector.hh>

struct KDNode{
	KDNode * parent;
	KDNode * smallchild;
	KDNode * bigchild;
	particlestructure p;
	int depth;
	long long int nchildren;
};

class SortOperator {
public:
	int axis;

    bool operator() (KDNode &a, KDNode &b) ;
};


class KDTree{  //Container for underlying k-d tree. Constructed from a list of positions.
public:
	KDNode * root;

	KDTree(KDNode * poslist, int poslistlen);
	KDNode * makeKD(KDNode * poslist,int poslistlen, KDNode * parent,int depth, SortOperator *sorters);


};

bool SortOperator::operator() (KDNode &a, KDNode &b) {
    bool result = false;
	if (axis ==0){
    	result = a.p.position.x < b.p.position.x;
    }
	else if (axis==1){
		result = a.p.position.y < b.p.position.y;
	}
	else{
		result = a.p.position.z < b.p.position.z;
	}
    return result;
}

KDTree::KDTree(KDNode *poslist,int poslistlen){
	SortOperator *sorter;
	sorter = new SortOperator[3];
	sorter[0].axis = 0;
	sorter[1].axis = 1;
	sorter[2].axis = 2;
	root = makeKD(poslist,poslistlen,0,0,sorter);

	delete[] sorter;
}

KDNode * KDTree::makeKD(KDNode * poslist,int poslistlen, KDNode * parent,int depth, SortOperator *sorters){

	//std::cout<<"Depth: "<<depth <<" len: "<<poslistlen<<"\n";
	if(poslistlen <1 ) return 0; //empty branch

	int axis = depth % 3; //We want to cycle through each axis as we decrease in depth

	std::sort(poslist,poslist+poslistlen,sorters[axis]); //We sort the list to get the median. This makes tree construction O(N log^2 N).
	int midpoint = poslistlen/2;
	KDNode *current; //begin constructing a new KDNode to return
	long long int nchild = 0;
	current = new KDNode;
	current->depth = depth;
	current->parent = parent;
	current->p.position = poslist[midpoint];
	//std::cout<<poslist[midpoint].x<<" , "<<poslist[midpoint].y<<" , "<<poslist[midpoint].z<<"\n";
	current->smallchild = makeKD(poslist, midpoint,current,depth+1,sorters);
	if(current->smallchild != 0) nchild += current->smallchild->nchildren +1;
	current->bigchild = makeKD(poslist + midpoint+1, poslistlen -midpoint -1,current,depth+1,sorters);
	if(current->bigchild != 0) nchild += current->bigchild->nchildren +1;
	return current;

}

float getNNLength(int NN, KDNode *x){ //get the length to the N nearest neighbors
	int foundNN = x->nchildren;
	KDNode *prev = x;
	//first we move up the tree until
	for(KDNode * current = x->parent; x != 0 && foundNN <NN; current = current->parent) {
		prev = current;
		foundNN = current->nchildren;
	}

	found

	return 0;
}







