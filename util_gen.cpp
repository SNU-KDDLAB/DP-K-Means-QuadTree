#include "util_gen.h"
#include "basic.h"
#include <random>
#include <cmath>

#define relBudget 1.4

/*double geoNoise(double eps, int depth, int tot_depth, int seed){
	double epsi = geoBudget(eps, depth, tot_depth);
	return laplace(0, 1/epsi, seed);
}*/

double uniformNoise(double eps, int tot_depth, int seed){
	return laplace(0, double(tot_depth + 1)/eps, seed);
}

double relativeNoise(double eps, int depth, int tot_depth, int dims, int seed) {
	double epsi = relativeBudget(eps, depth, tot_depth, dims);
	return laplace(0, 1 / epsi, seed);
}

double geoBudget(double eps, int depth, int tot_depth){
	return eps* pow(2, (tot_depth - depth) / 3.0) * (pow(2, 1 / 3.0) - 1) / (pow(2, (tot_depth + 1) / 3.0) - 1);
}

double uniformBudget(double eps, int tot_depth){
	return eps/(tot_depth + 1);
}

double relativeBudget(double eps, int depth, int tot_depth, int dims) {
	return eps * pow(relBudget, depth * dims) / (pow(relBudget, (tot_depth + 1) * dims) - 1) * (pow(relBudget, dims) - 1);
}

double findBudget(Noise noiseType, double eps, int depth, int tot_depth, int dims) {
	if(noiseType == Noise::uni)
		return uniformBudget(eps, tot_depth);
	if(noiseType == Noise::geo)
		return geoBudget(eps, depth, tot_depth);
	if(noiseType == Noise::rel)
		return relativeBudget(eps, depth, tot_depth, dims);
	return -1;
}