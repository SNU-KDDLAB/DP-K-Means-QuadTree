#ifndef DPQT_H
#define DPQT_H

#include <iostream>
#include <vector>
#include "util_gen.h"

using std::vector;

int partId(int dims, double* pivot, double* p) {
	int pid = 0;
	for (int d = 0; d < dims; ++d)
		if (p[d] >= pivot[d])
			pid += 1 << d;
	return pid;
}

void newBound(int dims, double** bound, double** nBound, int partId) {
	for (int d = 0; d < dims; ++d) {
		if ((partId >> d) % 2 == 0) {
			nBound[0][d] = bound[0][d];
			nBound[1][d] = (bound[0][d] + bound[1][d]) / 2;
		}
		else {
			nBound[0][d] = (bound[0][d] + bound[1][d]) / 2;
			nBound[1][d] = bound[1][d];
		}
	}
}

bool contain(int dims, double** bound, double** bound2) { // bound contains bound2
	for (int d = 0; d < dims; ++d) {
		if (bound[0][d] > bound2[0][d])
			return false;
		if (bound[1][d] < bound2[1][d])
			return false;
	}
	return true;
}

double overlap(int dims, double** bound, double** bound2) {
	double overlap = 1;
	double len;
	for (int d = 0; d < dims; ++d) {
		len = MIN(bound[1][d], bound2[1][d]) - MAX(bound[0][d], bound2[0][d]);
		if (len <= 0)
			return 0;
		else
			overlap *= len;
	}
	return overlap;
}

double area(int dims, double** bound) {
	double area = 1;
	for (int d = 0; d < dims; ++d) {
		area *= bound[1][d] - bound[0][d];
	}
	return area;
}

class DPQTNode;

class DPQTTree {
public:
	DPQTNode* root;
	int maxH;
	Noise noiseType;
	double eps;
	int seed;
	double T;
	int dims;
	int n;
	double** points;
	DPQTTree(int n, double** points, double** bound, double T, int dims, int maxH, Noise noiseType, double eps, int seed);
	~DPQTTree();
};

class DPQTNode {

public:
	double remain_budget;
	double** bound;
	int count;
	double noisyCount;
	// non-leaf
	double* pivot;
	int nChildren;
	int h;
	DPQTNode** children;
	DPQTTree* tree;

	void set(double** bound, int h) {
		int dims = tree->dims;
		this->h = h;
		this->bound = new double*[2];
		this->bound[0] = new double[dims];
		this->bound[1] = new double[dims];
		for (int d = 0; d < dims; ++d) {
			this->bound[0][d] = bound[0][d];
			this->bound[1][d] = bound[1][d];
		}
	}

	double rangeQuery(double** bound) {
		if (contain(tree->dims, bound, this->bound))
			return noisyCount;
		if (0 == nChildren)
			return noisyCount * overlap(tree->dims, this->bound, bound) / area(tree->dims, this->bound);
		double q = 0;
		for (int j = 0; j < nChildren; ++j) {
			if (NULL != children[j])
				q += children[j]->rangeQuery(bound);
		}
		return q;
	}
	void getLeaves(vector<DPQTNode*>* leaves){
		if (0 == nChildren) 
			leaves->push_back(this);
		else
			for (int j = 0; j < nChildren; ++j)
				if (NULL != children[j])
					children[j]->getLeaves(leaves);
	}
private:

	bool stopCondition() {
		double density = noisyCount / area(tree->dims, this->bound);
		//std::cout<<"split density = "<<density<<" T condition = "<<tree->T<<std::endl;
		//return h >= tree->maxH || density > tree->T || noisyCount < 0; // density based
		return h >= tree->maxH || noisyCount <= tree->T;
	}
public:
	DPQTNode(double** domain, DPQTTree* tree) { //root
		this->tree = tree;
		// root constructor
		set(domain, 0);

		double budget = findBudget(tree->noiseType, tree->eps, h, tree->maxH, tree->dims);
		double noise = laplace(0, 1 / budget, tree->seed);
		count = tree->n;
		noisyCount = count + noise;
		if (stopCondition()) {
			nChildren = 0;
			return;
		}
		nChildren = 1 << tree->dims;
		vector<int>** splits = new vector<int>*[nChildren];
		computePivot(splits);
		for (int i = 0; i < tree->n; ++i) {
			int pid = partId(tree->dims, pivot, tree->points[i]);
			if (NULL == splits[pid])
				splits[pid] = new vector<int>();
			splits[pid]->push_back(i);
		}
		makeChildren(splits, 1, tree->eps - budget);
	}
private:
	DPQTNode(vector<int>* ids, double** bound, int h, double Rbudget, DPQTTree* tree) {
		this->tree = tree;
		// non-root constructor
		set(bound, h);
		double budget = findBudget(tree->noiseType, tree->eps, h, tree->maxH, tree->dims);
		double noise = laplace(0, 1 / budget, tree->seed);

		count = ids->size();
		noisyCount = count + noise;
		if (stopCondition()) {
			nChildren = 0;
			return;
		}
		nChildren = 1 << tree->dims;
		vector<int>** splits = new vector<int>*[nChildren];
		computePivot(splits);
		for (int idx = 0; idx < count; ++idx) {
			int i = (*ids)[idx];
			int pid = partId(tree->dims, pivot, tree->points[i]);
			if (NULL == splits[pid])
				splits[pid] = new vector<int>();
			splits[pid]->push_back(i);
		}
		makeChildren(splits, h + 1, Rbudget - budget);
	}

	void computePivot(vector<int>** splits) {
		int dims = tree->dims;
		this->pivot = new double[dims];
		for (int d = 0; d < dims; ++d)
			this->pivot[d] = (bound[0][d] + bound[1][d]) / 2;
		for (int j = 0; j < nChildren; ++j)
			splits[j] = NULL;
	}

	void makeChildren( vector<int>** splits, int h, double budget) {
		double** nBound = new double*[2];
		nBound[0] = new double[tree->dims];
		nBound[1] = new double[tree->dims];
		children = new DPQTNode*[nChildren];
		if (Noise::rel == tree->noiseType)
			noisyCount = 0;
		for (int j = 0; j < nChildren; ++j) {
			if (NULL == splits[j])
				children[j] = NULL;
			else{
				newBound(tree->dims, bound, nBound, j);
				children[j] = new DPQTNode(splits[j], nBound, h, budget, tree);
				delete splits[j];
				if (Noise::rel == tree->noiseType)
					noisyCount += children[j]->noisyCount;
			}
		}
		delete[] nBound[0];
		delete[] nBound[1];
		delete[] nBound;
		delete[] splits;
	}
public:
	~DPQTNode(){
		delete[] bound[0];
		delete[] bound[1];
		delete[] bound;
		
		for (int i = 0; i < nChildren; ++i){
			if (children[i] != NULL)
				delete children[i];
		}
		if (nChildren != 0) {
			delete[] pivot;
			delete[] children;
		}
	}

};

DPQTTree::DPQTTree(int n, double** points, double** bound, double T, int dims, int maxH, Noise noiseType, double eps, int seed){
	this->n = n;
	this->points = points;
	this->T = T;
	this->dims = dims;
	this->maxH = maxH;
	this->noiseType = noiseType;
	this->eps = eps;
	this->seed = seed;
	this->root = new DPQTNode(bound, this);
}

DPQTTree::~DPQTTree(){
	delete this->root;
}

#endif